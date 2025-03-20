// SPDX-FileCopyrightText: Copyright (c) Kitware Inc.
// SPDX-License-Identifier: BSD-3-Clause

#include "vtkPVGeometryFilter.h"

#include "vtkAMRInformation.h"
#include "vtkAffineArray.h"
#include "vtkAlgorithmOutput.h"
#include "vtkBoundingBox.h"
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkCellGrid.h"
#include "vtkCommand.h"
#include "vtkCompositeDataSet.h"
#include "vtkConstantArray.h"
#include "vtkConvertToPartitionedDataSetCollection.h"
#include "vtkDataObjectMeshCache.h"
#include "vtkDataObjectTreeRange.h"
#include "vtkExplicitStructuredGrid.h"
#include "vtkExplicitStructuredGridSurfaceFilter.h"
#include "vtkFeatureEdges.h"
#include "vtkFloatArray.h"
#include "vtkGarbageCollector.h"
#include "vtkGenericDataSet.h"
#include "vtkGenericGeometryFilter.h"
#include "vtkGeometryFilter.h"
#include "vtkHyperTreeGrid.h"
#include "vtkHyperTreeGridFeatureEdges.h"
#include "vtkHyperTreeGridGeometry.h"
#include "vtkImageData.h"
#include "vtkInformation.h"
#include "vtkInformationIntegerVectorKey.h"
#include "vtkInformationKey.h"
#include "vtkInformationVector.h"
#include "vtkMath.h"
#include "vtkMultiProcessController.h"
#include "vtkNew.h"
#include "vtkObjectFactory.h"
#include "vtkOutlineSource.h"
#include "vtkOverlappingAMR.h"
#include "vtkPVTrivialProducer.h"
#include "vtkPartitionedDataSet.h"
#include "vtkPartitionedDataSetCollection.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkPolyDataNormals.h"
#include "vtkPolygon.h"
#include "vtkRange.h"
#include "vtkRecoverGeometryWireframe.h"
#include "vtkRectilinearGrid.h"
#include "vtkRectilinearGridOutlineFilter.h"
#include "vtkSmartPointer.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkStructuredGrid.h"
#include "vtkStructuredGridOutlineFilter.h"
#include "vtkTimerLog.h"
#include "vtkTriangleFilter.h"
#include "vtkUniformGrid.h"
#include "vtkUnsignedIntArray.h"
#include "vtkUnstructuredGrid.h"
#include "vtkUnstructuredGridGeometryFilter.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <memory>
#include <numeric>
#include <string>
#include <vector>

namespace details
{
static constexpr const char* ORIGINAL_FACE_IDS = "RecoverWireframeOriginalFaceIds";
static constexpr const char* TEMP_ORIGINAL_IDS = "__original_ids__";

//----------------------------------------------------------------------------
void AddOriginalIds(vtkDataSetAttributes* attributes, vtkIdType size)
{
  vtkNew<vtkAffineArray<vtkIdType>> ids;
  ids->SetBackend(std::make_shared<vtkAffineImplicitBackend<vtkIdType>>(1, 0));
  ids->SetNumberOfTuples(size);
  ids->SetName(details::TEMP_ORIGINAL_IDS);
  attributes->AddArray(ids);
}

//----------------------------------------------------------------------------
/**
 * Add an ids array on underlying PointData and CellData.
 * This is used by the vtkDataObjectMeshCache to forward
 * attributes data from a new input to the cached mesh.
 */
void AddTemporaryOriginalIdsArrays(vtkDataObject* object)
{
  auto dataSet = vtkDataSet::SafeDownCast(object);
  auto dataTree = vtkDataObjectTree::SafeDownCast(object);
  if (dataTree)
  {
    auto options = vtk::DataObjectTreeOptions::TraverseSubTree |
      vtk::DataObjectTreeOptions::SkipEmptyNodes | vtk::DataObjectTreeOptions::VisitOnlyLeaves;
    for (vtkDataObject* dataLeaf : vtk::Range(dataTree, options))
    {
      auto leafDataSet = vtkDataSet::SafeDownCast(dataLeaf);
      if (leafDataSet)
      {
        details::AddOriginalIds(leafDataSet->GetPointData(), leafDataSet->GetNumberOfPoints());
        details::AddOriginalIds(leafDataSet->GetCellData(), leafDataSet->GetNumberOfCells());
      }
    }
  }
  else if (dataSet)
  {
    details::AddOriginalIds(dataSet->GetPointData(), dataSet->GetNumberOfPoints());
    details::AddOriginalIds(dataSet->GetCellData(), dataSet->GetNumberOfCells());
  }
}

//----------------------------------------------------------------------------
/**
 * Cleanup the temporary array, as we do not want it to exists outside
 * of this filter.
 */
void CleanupTemporaryOriginalIds(vtkDataObject* object)
{
  auto dataSet = vtkDataSet::SafeDownCast(object);
  auto dataTree = vtkDataObjectTree::SafeDownCast(object);

  if (dataTree)
  {
    auto options = vtk::DataObjectTreeOptions::TraverseSubTree |
      vtk::DataObjectTreeOptions::SkipEmptyNodes | vtk::DataObjectTreeOptions::VisitOnlyLeaves;
    for (vtkDataObject* dataLeaf : vtk::Range(dataTree, options))
    {
      auto leafDataSet = vtkDataSet::SafeDownCast(dataLeaf);
      if (leafDataSet)
      {
        leafDataSet->GetPointData()->RemoveArray(details::TEMP_ORIGINAL_IDS);
        leafDataSet->GetCellData()->RemoveArray(details::TEMP_ORIGINAL_IDS);
      }
    }
  }
  else if (dataSet)
  {
    dataSet->GetPointData()->RemoveArray(details::TEMP_ORIGINAL_IDS);
    dataSet->GetCellData()->RemoveArray(details::TEMP_ORIGINAL_IDS);
  }
}

};

template <typename T>
void GetValidWholeExtent(T* ds, const int wholeExt[6], int validWholeExt[6])
{
  if (wholeExt[0] <= wholeExt[1] && wholeExt[2] <= wholeExt[3] && wholeExt[4] <= wholeExt[5])
  {
    std::copy(wholeExt, wholeExt + 6, validWholeExt);
  }
  else
  {
    ds->GetExtent(validWholeExt);
  }
}

//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkPVGeometryFilter);
//----------------------------------------------------------------------------
vtkCxxSetObjectMacro(vtkPVGeometryFilter, Controller, vtkMultiProcessController);

//----------------------------------------------------------------------------
vtkInformationKeyMacro(vtkPVGeometryFilter, POINT_OFFSETS, IntegerVector);
vtkInformationKeyMacro(vtkPVGeometryFilter, VERTS_OFFSETS, IntegerVector);
vtkInformationKeyMacro(vtkPVGeometryFilter, LINES_OFFSETS, IntegerVector);
vtkInformationKeyMacro(vtkPVGeometryFilter, POLYS_OFFSETS, IntegerVector);
vtkInformationKeyMacro(vtkPVGeometryFilter, STRIPS_OFFSETS, IntegerVector);

//----------------------------------------------------------------------------
vtkPVGeometryFilter::vtkPVGeometryFilter()
{
  this->OutlineFlag = 0;
  this->UseOutline = 1;
  this->GenerateFeatureEdges = false;
  this->BlockColorsDistinctValues = 7;
  // generating cell normals by default really slows down paraview
  // it is especially noticeable with the OpenGL2 backend.  Leaving
  // it on for the old backend as some tests rely on the cell normals
  // to be there as they use them for other purposes/etc.
  this->GenerateCellNormals = false;
  this->GeneratePointNormals = false;
  this->Splitting = 1;
  this->FeatureAngle = 30.0;
  this->Triangulate = false;
  this->NonlinearSubdivisionLevel = 1;

  this->GeometryFilter = vtkSmartPointer<vtkGeometryFilter>::New();
  // we're prepping geometry for rendering
  // fast mode might generate wrong results because of certain assumptions that don't always hold
  // for that reason, the default is to have fast mode off
  this->GeometryFilter->SetFastMode(false);
  // Always preserve ghost interfaces. This is necessary to ensure that
  // if ghost cells are present, they can be used to generate point normals.
  this->GeometryFilter->SetRemoveGhostInterfaces(false);
  this->GenericGeometryFilter = vtkSmartPointer<vtkGenericGeometryFilter>::New();
  this->UnstructuredGridGeometryFilter = vtkSmartPointer<vtkUnstructuredGridGeometryFilter>::New();
  this->RecoverWireframeFilter = vtkSmartPointer<vtkRecoverGeometryWireframe>::New();
  this->FeatureEdgesFilter = vtkSmartPointer<vtkFeatureEdges>::New();
  this->PolyDataNormals = vtkSmartPointer<vtkPolyDataNormals>::New();
  this->PolyDataNormals->SetComputePointNormals(this->GeneratePointNormals);
  this->PolyDataNormals->SetComputeCellNormals(this->GenerateCellNormals);
  this->PolyDataNormals->SetSplitting(this->Splitting);
  this->PolyDataNormals->SetFeatureAngle(this->FeatureAngle);
  this->PolyDataNormals->AutoOrientNormalsOff();
  this->PolyDataNormals->ConsistencyOff();
  this->PolyDataNormals->NonManifoldTraversalOff();
  this->PolyDataNormals->FlipNormalsOff();
  this->PolyDataNormals->SetOutputPointsPrecision(vtkAlgorithm::DEFAULT_PRECISION);

  // Setup a callback for the internal readers to report progress.
  this->GeometryFilter->AddObserver(
    vtkCommand::ProgressEvent, this, &vtkPVGeometryFilter::HandleGeometryFilterProgress);
  this->GenericGeometryFilter->AddObserver(
    vtkCommand::ProgressEvent, this, &vtkPVGeometryFilter::HandleGeometryFilterProgress);
  this->UnstructuredGridGeometryFilter->AddObserver(
    vtkCommand::ProgressEvent, this, &vtkPVGeometryFilter::HandleGeometryFilterProgress);
  this->RecoverWireframeFilter->AddObserver(
    vtkCommand::ProgressEvent, this, &vtkPVGeometryFilter::HandleGeometryFilterProgress);
  this->PolyDataNormals->AddObserver(
    vtkCommand::ProgressEvent, this, &vtkPVGeometryFilter::HandleGeometryFilterProgress);

  this->Controller = nullptr;
  this->SetController(vtkMultiProcessController::GetGlobalController());
  this->GenerateProcessIds = (this->Controller && this->Controller->GetNumberOfProcesses() > 1);

  this->OutlineSource = vtkSmartPointer<vtkOutlineSource>::New();

  this->PassThroughCellIds = 1;
  this->PassThroughPointIds = 1;

  this->HideInternalAMRFaces = true;
  this->UseNonOverlappingAMRMetaDataForOutlines = true;

  this->MeshCache->SetConsumer(this);
  this->MeshCache->AddOriginalIds(vtkDataObject::POINT, details::TEMP_ORIGINAL_IDS);
  this->MeshCache->AddOriginalIds(vtkDataObject::CELL, details::TEMP_ORIGINAL_IDS);
}

//----------------------------------------------------------------------------
vtkPVGeometryFilter::~vtkPVGeometryFilter()
{
  this->SetController(nullptr);
}

//----------------------------------------------------------------------------
int vtkPVGeometryFilter::RequestDataObject(vtkInformation* vtkNotUsed(request),
  vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
  auto input = vtkDataObject::GetData(inputVector[0], 0);
  int outputType = -1;

  if (!input)
  {
    vtkErrorMacro("Missing input data");
    return 0;
  }

  if (input->IsA("vtkDataSet") || input->IsA("vtkGenericDataSet") || input->IsA("vtkCellGrid") ||
    input->IsA("vtkHyperTreeGrid"))
  {
    outputType = VTK_POLY_DATA;
  }
  else if (input->IsA("vtkMultiBlockDataSet"))
  {
    // Some developers have sub-classed vtkMultiBlockDataSet, in which case,
    // we try to preserve the type.
    outputType = input->GetDataObjectType();
  }
  else if (input->IsA("vtkCompositeDataSet"))
  {
    outputType = VTK_PARTITIONED_DATA_SET_COLLECTION;
  }

  return vtkDataObjectAlgorithm::SetOutputDataObject(
           outputType, outputVector->GetInformationObject(0), /*exact*/ true)
    ? 1
    : 0;
}

//----------------------------------------------------------------------------
void vtkPVGeometryFilter::HandleGeometryFilterProgress(vtkObject* caller, unsigned long, void*)
{
  vtkAlgorithm* algorithm = vtkAlgorithm::SafeDownCast(caller);
  // This limits progress for only the GeometryFilter.
  if (algorithm)
  {
    const double progress = algorithm->GetProgress();
    if (progress > 0.0 && progress < 1.0)
    {
      this->UpdateProgress(progress);
    }
    if (this->AbortExecute)
    {
      algorithm->SetAbortExecute(1);
    }
  }
}

//----------------------------------------------------------------------------
int vtkPVGeometryFilter::CheckAttributes(vtkDataObject* input)
{
  if (auto ds = vtkDataSet::SafeDownCast(input))
  {
    return ds->CheckAttributes();
  }
  else if (auto cds = vtkCompositeDataSet::SafeDownCast(input))
  {
    auto iter = vtk::TakeSmartPointer(cds->NewIterator());
    for (iter->GoToFirstItem(); !iter->IsDoneWithTraversal(); iter->GoToNextItem())
    {
      vtkDataObject* curDs = iter->GetCurrentDataObject();
      if (curDs)
      {
        return this->CheckAttributes(curDs);
      }
    }
  }
  return 0;
}

//----------------------------------------------------------------------------
void vtkPVGeometryFilter::ExecuteAMRBlockOutline(
  const double bounds[6], vtkPolyData* output, const bool extractface[6])
{
  // we generate outline faces, so that front face/back face culling works if
  // needed BUG #0011065. We only do this for AMR datasets for now, but we can
  // extend to all types of datasets, if needed.

  auto points = vtkSmartPointer<vtkPoints>::New();
  points->SetNumberOfPoints(8);

  points->SetPoint(0, bounds[0], bounds[2], bounds[4]);
  points->SetPoint(1, bounds[1], bounds[2], bounds[4]);
  points->SetPoint(2, bounds[0], bounds[3], bounds[4]);
  points->SetPoint(3, bounds[1], bounds[3], bounds[4]);
  points->SetPoint(4, bounds[0], bounds[2], bounds[5]);
  points->SetPoint(5, bounds[1], bounds[2], bounds[5]);
  points->SetPoint(6, bounds[0], bounds[3], bounds[5]);
  points->SetPoint(7, bounds[1], bounds[3], bounds[5]);

  auto lines = vtkSmartPointer<vtkCellArray>::New();
  lines->Allocate(lines->EstimateSize(12, 2));

  // xmin face
  if (extractface[0])
  {
    vtkIdType pts[4] = { 0, 4, 6, 2 };
    lines->InsertNextCell(4, pts);
  }
  // xmax face
  if (extractface[1])
  {
    vtkIdType pts[4] = { 1, 3, 7, 5 };
    lines->InsertNextCell(4, pts);
  }
  // ymin face
  if (extractface[2])
  {
    vtkIdType pts[4] = { 0, 1, 5, 4 };
    lines->InsertNextCell(4, pts);
  }
  // ymax face
  if (extractface[3])
  {
    vtkIdType pts[4] = { 2, 6, 7, 3 };
    lines->InsertNextCell(4, pts);
  }
  // zmin face
  if (extractface[4])
  {
    vtkIdType pts[4] = { 0, 2, 3, 1 };
    lines->InsertNextCell(4, pts);
  }
  // zmax face
  if (extractface[5])
  {
    vtkIdType pts[4] = { 4, 5, 7, 6 };
    lines->InsertNextCell(4, pts);
  }
  lines->Squeeze();

  output->SetPoints(points);
  output->SetPolys(lines);

  this->OutlineFlag = 1;
}

//----------------------------------------------------------------------------
void vtkPVGeometryFilter::ExecuteAMRBlock(
  vtkUniformGrid* input, vtkPolyData* output, const bool extractface[6])
{
  assert(input != nullptr && output != nullptr && this->UseOutline == 0);
  if (input->GetNumberOfCells() > 0)
  {
    int extent[6];
    input->GetExtent(extent);
    this->GeometryFilter->StructuredExecute(input, output, extent, const_cast<bool*>(extractface));
  }
  this->OutlineFlag = 0;
}

//----------------------------------------------------------------------------
void vtkPVGeometryFilter::ExecuteBlock(vtkDataObject* input, vtkPolyData* output, int doCommunicate,
  int updatePiece, int updateNumPieces, int updateGhosts, const int* wholeExtent)
{
  // Copy field data from the input block to the output block
  output->GetFieldData()->PassData(input->GetFieldData());

  if (auto imageData = vtkImageData::SafeDownCast(input))
  {
    this->ImageDataExecute(imageData, output, doCommunicate, updatePiece, wholeExtent);
  }
  else if (auto structuredGrid = vtkStructuredGrid::SafeDownCast(input))
  {
    this->StructuredGridExecute(
      structuredGrid, output, updatePiece, updateNumPieces, updateGhosts, wholeExtent);
  }
  else if (auto rectilinearGrid = vtkRectilinearGrid::SafeDownCast(input))
  {
    this->RectilinearGridExecute(
      rectilinearGrid, output, updatePiece, updateNumPieces, updateGhosts, wholeExtent);
  }
  else if (auto unstructuredGridBase = vtkUnstructuredGridBase::SafeDownCast(input))
  {
    this->UnstructuredGridExecute(unstructuredGridBase, output, doCommunicate);
  }
  else if (auto polyData = vtkPolyData::SafeDownCast(input))
  {
    this->PolyDataExecute(polyData, output, doCommunicate);
  }
  else if (auto hyperTreeGrid = vtkHyperTreeGrid::SafeDownCast(input))
  {
    // If feature edges are requested, vtkHyperTreeGridFeatureEdges will generate them,
    // so no need to compute the HTG geometry.
    if (!this->GenerateFeatureEdges)
    {
      this->HyperTreeGridExecute(hyperTreeGrid, output, doCommunicate);
    }
  }
  else if (auto explicitStructuredGrid = vtkExplicitStructuredGrid::SafeDownCast(input))
  {
    this->ExplicitStructuredGridExecute(explicitStructuredGrid, output, doCommunicate, wholeExtent);
  }
  else if (auto dataset = vtkDataSet::SafeDownCast(input))
  {
    this->DataSetExecute(dataset, output, doCommunicate);
  }
  else if (auto genericDataSet = vtkGenericDataSet::SafeDownCast(input))
  {
    this->GenericDataSetExecute(genericDataSet, output, doCommunicate);
  }
  else if (auto cellGrid = vtkCellGrid::SafeDownCast(input))
  {
    this->CellGridExecute(cellGrid, output, doCommunicate);
  }
}

//----------------------------------------------------------------------------
void vtkPVGeometryFilter::UpdateCache(vtkDataObject* output)
{
  this->MeshCache->UpdateCache(output);
  details::CleanupTemporaryOriginalIds(output);
}

//----------------------------------------------------------------------------
bool vtkPVGeometryFilter::UseCacheIfPossible(vtkDataObject* input, vtkDataObject* output)
{
  details::AddTemporaryOriginalIdsArrays(input);
  this->MeshCache->SetOriginalDataObject(input);

  vtkDataObjectMeshCache::Status status = this->MeshCache->GetStatus();
  if (status.enabled())
  {
    this->MeshCache->CopyCacheToDataObject(output);
    details::CleanupTemporaryOriginalIds(output);
    return true;
  }

  return false;
}

//----------------------------------------------------------------------------
int vtkPVGeometryFilter::RequestData(
  vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
  auto input = vtkDataObject::GetData(inputVector[0], 0);
  auto dataObjectOutput = vtkDataObject::GetData(outputVector, 0);

  vtkSmartPointer<vtkDataObject> modifiedInput = input;
  const bool isCachingSupported = this->MeshCache->IsSupportedData(input);
  if (isCachingSupported)
  {
    // create a copy as we add some temporary array.
    modifiedInput.TakeReference(input->NewInstance());
    modifiedInput->ShallowCopy(input);
    if (this->UseCacheIfPossible(modifiedInput, dataObjectOutput))
    {
      return 1;
    }
  }

  if (input->IsA("vtkCompositeDataSet"))
  {
    vtkTimerLog::MarkStartEvent("vtkPVGeometryFilter::RequestData");
    vtkGarbageCollector::DeferredCollectionPush();
    if (input->IsA("vtkUniformGridAMR"))
    {
      this->RequestAMRData(request, inputVector, outputVector);
    }
    else
    {
      this->RequestDataObjectTree(request, inputVector, outputVector);
    }
    vtkTimerLog::MarkStartEvent("vtkPVGeometryFilter::GarbageCollect");
    vtkGarbageCollector::DeferredCollectionPop();
    vtkTimerLog::MarkEndEvent("vtkPVGeometryFilter::GarbageCollect");
    vtkTimerLog::MarkEndEvent("vtkPVGeometryFilter::RequestData");
  }
  else
  {
    vtkPolyData* output = vtkPolyData::GetData(outputVector, 0);
    assert(output != nullptr);

    int procid = 0;
    int numProcs = 1;
    if (this->Controller)
    {
      procid = this->Controller->GetLocalProcessId();
      numProcs = this->Controller->GetNumberOfProcesses();
    }
    int* wholeExtent =
      vtkStreamingDemandDrivenPipeline::GetWholeExtent(inputVector[0]->GetInformationObject(0));

    auto inputHTG = vtkHyperTreeGrid::SafeDownCast(input);
    if (this->GenerateFeatureEdges && inputHTG)
    {
      this->GenerateFeatureEdgesHTG(inputHTG, output);
    }
    else
    {
      this->ExecuteBlock(modifiedInput, output, 1, procid, numProcs, 0, wholeExtent);
      this->CleanupOutputData(output);
    }
  }

  if (isCachingSupported)
  {
    this->UpdateCache(dataObjectOutput);
  }
  return 1;
}

//----------------------------------------------------------------------------
void vtkPVGeometryFilter::GenerateFeatureEdgesHTG(vtkHyperTreeGrid* input, vtkPolyData* output)
{
  vtkNew<vtkHyperTreeGridFeatureEdges> featureEdgesFilter;
  vtkNew<vtkHyperTreeGrid> htgCopy;
  htgCopy->ShallowCopy(input);
  featureEdgesFilter->SetInputData(htgCopy);
  featureEdgesFilter->Update();
  output->ShallowCopy(featureEdgesFilter->GetOutput());
  output->RemoveGhostCells();
  if (this->GenerateProcessIds)
  {
    this->GenerateProcessIdsArrays(output);
  }
}

//----------------------------------------------------------------------------
void vtkPVGeometryFilter::GenerateProcessIdsArrays(vtkPolyData* output)
{
  const unsigned int procId =
    this->Controller ? static_cast<unsigned int>(this->Controller->GetLocalProcessId()) : 0;

  const vtkIdType numPoints = output->GetNumberOfPoints();
  if (numPoints > 0)
  {
    vtkNew<vtkConstantArray<unsigned int>> pointsProcArray;
    pointsProcArray->SetNumberOfValues(numPoints);
    pointsProcArray->ConstructBackend(procId);
    pointsProcArray->SetName("vtkProcessId");
    output->GetPointData()->AddArray(pointsProcArray);
  }

  const vtkIdType numCells = output->GetNumberOfCells();
  if (numCells > 0)
  {
    vtkNew<vtkConstantArray<unsigned int>> cellsProcArray;
    cellsProcArray->SetNumberOfValues(numCells);
    cellsProcArray->ConstructBackend(procId);
    cellsProcArray->SetName("vtkProcessId");
    output->GetCellData()->AddArray(cellsProcArray);
  }
}

//----------------------------------------------------------------------------
void vtkPVGeometryFilter::CleanupOutputData(vtkPolyData* output)
{
  if (this->GenerateFeatureEdges)
  {
    this->FeatureEdgesFilter->SetInputData(output);
    this->FeatureEdgesFilter->Update();
    output->ShallowCopy(this->FeatureEdgesFilter->GetOutput());
  }
  if (this->GenerateCellNormals || this->GeneratePointNormals)
  {
    this->ExecuteNormalsComputation(output);
  }
  output->RemoveGhostCells();
  if (this->GenerateProcessIds && output)
  {
    this->GenerateProcessIdsArrays(output);
  }
}

//----------------------------------------------------------------------------
void vtkPVGeometryFilter::AddCompositeIndex(vtkPolyData* pd, unsigned int index)
{
  vtkNew<vtkConstantArray<unsigned int>> pointsCIArray;
  pointsCIArray->SetNumberOfValues(pd->GetNumberOfPoints());
  pointsCIArray->ConstructBackend(index);
  pointsCIArray->SetName("vtkCompositeIndex");
  pd->GetPointData()->AddArray(pointsCIArray);

  vtkNew<vtkConstantArray<unsigned int>> cellsCIArray;
  cellsCIArray->SetNumberOfValues(pd->GetNumberOfCells());
  cellsCIArray->ConstructBackend(index);
  cellsCIArray->SetName("vtkCompositeIndex");
  pd->GetCellData()->AddArray(cellsCIArray);
}

//----------------------------------------------------------------------------
void vtkPVGeometryFilter::AddBlockColors(vtkDataObject* pd, unsigned int index)
{
  vtkNew<vtkUnsignedIntArray> blockColorsArray;
  blockColorsArray->SetNumberOfTuples(1);
  blockColorsArray->SetValue(0, index % this->BlockColorsDistinctValues);
  blockColorsArray->SetName("vtkBlockColors");
  pd->GetFieldData()->AddArray(blockColorsArray);
}

//----------------------------------------------------------------------------
void vtkPVGeometryFilter::AddHierarchicalIndex(
  vtkPolyData* pd, unsigned int level, unsigned int index)
{
  vtkNew<vtkConstantArray<unsigned int>> dslevel;
  dslevel->SetNumberOfValues(pd->GetNumberOfCells());
  dslevel->ConstructBackend(level);
  dslevel->SetName("vtkAMRLevel");
  pd->GetCellData()->AddArray(dslevel);

  vtkNew<vtkConstantArray<unsigned int>> dsindex;
  dsindex->SetNumberOfValues(pd->GetNumberOfCells());
  dsindex->ConstructBackend(index);
  dsindex->SetName("vtkAMRIndex");
  pd->GetCellData()->AddArray(dsindex);
}

//----------------------------------------------------------------------------
int vtkPVGeometryFilter::RequestAMRData(
  vtkInformation*, vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
  vtkTimerLog::MarkStartEvent("vtkPVGeometryFilter::RequestAMRData");

  // STEP 0: Acquire input & output object
  auto output = vtkPartitionedDataSetCollection::GetData(outputVector, 0);
  if (!output)
  {
    vtkErrorMacro("Output vtkPartitionedDataSetCollection is nullptr.");
    return 0;
  }

  vtkUniformGridAMR* amr = vtkUniformGridAMR::GetData(inputVector[0], 0);
  if (!amr)
  {
    vtkErrorMacro("Input vtkUniformGridAMR is nullptr.");
    return 0;
  }
  output->CopyStructure(amr);

  // STEP 1: Check Attributes
  vtkTimerLog::MarkStartEvent("vtkPVGeometryFilter::CheckAttributes");
  if (this->CheckAttributes(amr))
  {
    vtkErrorMacro("CheckAttributes() failed!");
    return 0;
  }
  vtkTimerLog::MarkEndEvent("vtkPVGeometryFilter::CheckAttributes");

  // STEP 2: Loop through data, determine if they are visible and call
  // execute block to get the polydata to render.
  // We use different implementations for overlapping and non-overlapping amrs.
  vtkOverlappingAMR* overlappingAMR = vtkOverlappingAMR::SafeDownCast(amr);

  double bounds[6];
  amr->GetBounds(bounds);
  const vtkBoundingBox amrBBox(bounds);
  if (this->Controller && this->Controller->GetNumberOfProcesses() > 1)
  {
    // Since bounds are not necessarily synced up, especially for non-overlapping
    // AMR datasets, we sync them up across all processes.
    vtkBoundingBox recvAmrBBox;
    this->Controller->AllReduce(amrBBox, recvAmrBBox);
    recvAmrBBox.GetBounds(bounds);
  }

  for (unsigned int level = 0; level < amr->GetNumberOfLevels(); ++level)
  {
    const unsigned int num_datasets = amr->GetNumberOfDataSets(level);
    for (unsigned int partitionIdx = 0; partitionIdx < num_datasets; ++partitionIdx)
    {
      vtkUniformGrid* ug = amr->GetDataSet(level, partitionIdx);
      if (!ug && ((this->UseOutline == 0) || (!overlappingAMR)))
      {
        // if this->UseOutline == 0, we need uniform grid to be present.

        // if this->UseOutline ==1, we need ug only for non-overlapping AMR. For
        // overlapping AMR, we can generate outline using the meta-data
        // available.
        continue;
      }

      if (overlappingAMR && !this->UseNonOverlappingAMRMetaDataForOutlines && !ug)
      {
        // for non-overlapping AMR, if we were told to not use meta-data, don't.
        continue;
      }

      double data_bounds[6];
      double error_margin = 0.01;

      // we have different mechanisms for determining if any of the faces of the
      // block are visible and what faces are visible based on the type of amr.
      if (overlappingAMR)
      {
        // for overlappingAMR, we use the meta-data to determine AMR bounds.
        overlappingAMR->GetAMRInfo()->GetBounds(level, partitionIdx, data_bounds);
        double data_spacing[3];
        overlappingAMR->GetAMRInfo()->GetSpacing(level, data_spacing);
        error_margin = vtkMath::Norm(data_spacing);
      }
      else if (ug)
      {
        // for non-overlapping AMR, we use the bounds from the heavy-data itself.
        ug->GetBounds(data_bounds);

        double data_spacing[3];
        ug->GetSpacing(data_spacing);
        error_margin = vtkMath::Norm(data_spacing);
      }
      else
      {
        continue; // skip block.
      }

      bool extractface[6] = { true, true, true, true, true, true };
      for (int cc = 0; this->HideInternalAMRFaces && cc < 6; cc++)
      {
        const double delta = std::abs(data_bounds[cc] - bounds[cc]);
        extractface[cc] = (delta < error_margin);
      }

      if (!(extractface[0] || extractface[1] || extractface[2] || extractface[3] ||
            extractface[4] || extractface[5]))
      {
        // we are not extracting a single face. nothing to do here.
        continue;
      }

      const vtkNew<vtkPolyData> outputPartition;
      if (this->UseOutline)
      {
        this->ExecuteAMRBlockOutline(data_bounds, outputPartition, extractface);
        // don't process attribute arrays when generating outlines.
      }
      else
      {
        this->ExecuteAMRBlock(ug, outputPartition, extractface);
        // add atttribute arrays when not generating outlines
        this->CleanupOutputData(outputPartition);
        this->AddCompositeIndex(outputPartition, amr->GetCompositeIndex(level, partitionIdx));
        this->AddHierarchicalIndex(outputPartition, level, partitionIdx);
        // we don't call this->AddBlockColors() for AMR dataset since it doesn't
        // make sense, nor can be supported since all datasets merged into a
        // single polydata for rendering.
      }
      output->SetPartition(level, partitionIdx, outputPartition);
    }
  }

  vtkTimerLog::MarkEndEvent("vtkPVGeometryFilter::RequestAMRData");
  return 1;
}

//----------------------------------------------------------------------------
int vtkPVGeometryFilter::RequestDataObjectTree(
  vtkInformation*, vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
  vtkTimerLog::MarkStartEvent("vtkPVGeometryFilter::RequestDataObjectTree");

  auto output = vtkDataObjectTree::GetData(outputVector, 0);
  if (!output)
  {
    vtkErrorMacro("Output vtkPartitionedDataSetCollection is nullptr.");
    return 0;
  }

  auto realInput = vtkDataObjectTree::GetData(inputVector[0], 0);
  if (!realInput)
  {
    vtkErrorMacro("Input vtkDataObjectTree is nullptr.");
    return 0;
  }
  vtkSmartPointer<vtkDataObjectTree> tempInput;
  if (realInput->IsA("vtkPartitionedDataSetCollection") || realInput->IsA("vtkMultiBlockDataSet"))
  {
    tempInput = realInput;
  }
  else
  {
    vtkNew<vtkConvertToPartitionedDataSetCollection> converter;
    converter->SetInputDataObject(realInput);
    converter->SetContainerAlgorithm(this);
    converter->Update();
    tempInput = converter->GetOutput();
  }
  output->CopyStructure(tempInput);

  vtkTimerLog::MarkStartEvent("vtkPVGeometryFilter::CheckAttributes");
  if (this->CheckAttributes(tempInput))
  {
    return 0;
  }
  vtkTimerLog::MarkEndEvent("vtkPVGeometryFilter::CheckAttributes");

  // create a copy as we add some temporary array.
  vtkSmartPointer<vtkDataObjectTree> input;
  input.TakeReference(realInput->NewInstance());
  input->ShallowCopy(realInput);

  vtkTimerLog::MarkStartEvent("vtkPVGeometryFilter::ExecuteCompositeDataSet");

  auto inIter = vtk::TakeSmartPointer(tempInput->NewTreeIterator());
  inIter->VisitOnlyLeavesOn();
  inIter->SkipEmptyNodesOn();

  // get a block count for progress scaling.
  unsigned int totalNumberOfBlocks = 0;
  for (inIter->InitTraversal(); !inIter->IsDoneWithTraversal(); inIter->GoToNextItem())
  {
    ++totalNumberOfBlocks;
  }

  int* wholeExtent =
    vtkStreamingDemandDrivenPipeline::GetWholeExtent(inputVector[0]->GetInformationObject(0));
  int numInputs = 0;
  for (inIter->InitTraversal(); !inIter->IsDoneWithTraversal(); inIter->GoToNextItem())
  {
    vtkDataObject* block = inIter->GetCurrentDataObject();
    if (!block)
    {
      continue;
    }

    vtkNew<vtkPolyData> tmpOut;
    auto blockHTG = vtkHyperTreeGrid::SafeDownCast(block);
    if (this->GenerateFeatureEdges && blockHTG)
    {
      this->GenerateFeatureEdgesHTG(blockHTG, tmpOut);
    }
    else
    {
      this->ExecuteBlock(block, tmpOut, 0, 0, 1, 0, wholeExtent);
      this->CleanupOutputData(tmpOut);
    }
    // skip empty nodes.
    if (tmpOut->GetNumberOfPoints() > 0)
    {
      output->SetDataSet(inIter, tmpOut);
      this->AddCompositeIndex(tmpOut, inIter->GetCurrentFlatIndex());
    }
    this->UpdateProgress(static_cast<float>(++numInputs) / totalNumberOfBlocks);
  }
  vtkTimerLog::MarkEndEvent("vtkPVGeometryFilter::ExecuteCompositeDataSet");

  auto outIter = vtk::TakeSmartPointer(output->NewTreeIterator());
  if (this->Controller && this->Controller->GetNumberOfProcesses() > 1)
  {
    // When running in parallel, processes may have nullptr-leaf nodes at
    // different locations. To make our life easier in subsequent filtering such as
    // vtkAllToNRedistributeCompositePolyData or vtkKdTreeManager we ensure that
    // all nullptr-leafs match up across processes i.e. if any leaf is non-nullptr on
    // any process, then all other processes add empty polydatas for that leaf.
    std::vector<unsigned char> non_null_leaves;
    non_null_leaves.reserve(totalNumberOfBlocks); // just an estimate.
    outIter->VisitOnlyLeavesOn();
    outIter->SkipEmptyNodesOn();
    for (outIter->InitTraversal(); !outIter->IsDoneWithTraversal(); outIter->GoToNextItem())
    {
      non_null_leaves.resize(outIter->GetCurrentFlatIndex() + 1, 0);
      non_null_leaves[outIter->GetCurrentFlatIndex()] = static_cast<unsigned char>(1);
    }

    const int count = static_cast<int>(non_null_leaves.size());
    int reduced_size;
    this->Controller->AllReduce(&count, &reduced_size, 1, vtkCommunicator::MAX_OP);
    assert(reduced_size >= static_cast<int>(non_null_leaves.size()));
    non_null_leaves.resize(reduced_size, 0);
    // if reduced_size ==0, then all processes have no non-nullptr-leaves, so
    // nothing special to do here.
    if (reduced_size != 0)
    {
      std::vector<unsigned char> reduced_non_null_leaves;
      reduced_non_null_leaves.resize(reduced_size, 0);
      this->Controller->AllReduce(non_null_leaves.data(), reduced_non_null_leaves.data(),
        reduced_size, vtkCommunicator::MAX_OP);

      outIter->SkipEmptyNodesOff();
      outIter->VisitOnlyLeavesOn();
      for (outIter->InitTraversal(); !outIter->IsDoneWithTraversal(); outIter->GoToNextItem())
      {
        const unsigned int index = outIter->GetCurrentFlatIndex();
        if (outIter->GetCurrentDataObject() == nullptr &&
          index < static_cast<unsigned int>(reduced_non_null_leaves.size()) &&
          reduced_non_null_leaves[index] != 0)
        {
          vtkNew<vtkPolyData> trivalInput;
          this->AddCompositeIndex(trivalInput, index);
          output->SetDataSet(outIter, trivalInput);
        }
      }
    }
  }

  unsigned int block_id = 0;
  if (auto outputPDC = vtkPartitionedDataSetCollection::SafeDownCast(output))
  {
    // To avoid paraview/paraview#20908, we use a 2-level approach.
    for (block_id = 0; block_id < outputPDC->GetNumberOfPartitionedDataSets(); ++block_id)
    {
      auto datasets = vtkCompositeDataSet::GetDataSets(outputPDC->GetPartitionedDataSet(block_id));
      for (auto dataset : datasets)
      {
        this->AddBlockColors(dataset, block_id);
      }
    }
  }
  else // vtkMultiBlockDataSet
  {
    // At this point, all ranks have consistent tree structure with leaf nodes non-nullptr
    // at exactly same locations. This is a good point to assign block colors.
    outIter->SkipEmptyNodesOff();
    outIter->VisitOnlyLeavesOff();
    for (outIter->InitTraversal(); !outIter->IsDoneWithTraversal(); outIter->GoToNextItem())
    {
      auto block = outIter->GetCurrentDataObject();
      if (!block)
      {
        ++block_id;
      }
      else if (auto pds = vtkPartitionedDataSet::SafeDownCast(block))
      {
        for (unsigned int i = 0; i < pds->GetNumberOfPartitions(); ++i)
        {
          if (auto ds = pds->GetPartition(i))
          {
            this->AddBlockColors(ds, block_id);
          }
          outIter->GoToNextItem();
        }
        // MB increments the iterator by the number of partitions, while PDC increments it by 1
        block_id += pds->GetNumberOfPartitions();
      }
      else if (!block->IsA("vtkMultiBlockDataSet"))
      {
        this->AddBlockColors(block, block_id);
        ++block_id;
      }
    }
  }

  if (block_id > 0)
  {
    // Add block colors to root-node's field data to keep it from being flagged as partial.
    this->AddBlockColors(output, 0);
  }

  vtkTimerLog::MarkEndEvent("vtkPVGeometryFilter::RequestDataObjectTree");
  return 1;
}

//----------------------------------------------------------------------------
// We need to change the mapper.  Now it always flat shades when cell normals
// are available.
void vtkPVGeometryFilter::ExecuteCellNormals(vtkPolyData* output, int doCommunicate)
{
  // Do not generate cell normals if any of the processes
  // have lines, verts or strips.
  vtkCellArray* aPrim;
  int skip = 0;
  aPrim = output->GetVerts();
  if (aPrim && aPrim->GetNumberOfCells())
  {
    skip = 1;
  }
  aPrim = output->GetLines();
  if (aPrim && aPrim->GetNumberOfCells())
  {
    skip = 1;
  }
  aPrim = output->GetStrips();
  if (aPrim && aPrim->GetNumberOfCells())
  {
    skip = 1;
  }
  if (this->Controller && doCommunicate)
  {
    int reduced_skip = 0;
    if (!this->Controller->AllReduce(&skip, &reduced_skip, 1, vtkCommunicator::MAX_OP))
    {
      vtkErrorMacro("Failed to reduce correctly.");
      skip = 1;
    }
    else
    {
      skip = reduced_skip;
    }
  }
  if (skip)
  {
    return;
  }

  aPrim = output->GetPolys();
  const vtkIdType numPolys = aPrim ? aPrim->GetNumberOfCells() : 0;
  if (numPolys != output->GetNumberOfCells())
  {
    vtkErrorMacro("Number of numPolys does not match output.");
    return;
  }

  vtkNew<vtkFloatArray> cellNormals;
  cellNormals->SetName("cellNormals");
  cellNormals->SetNumberOfComponents(3);
  cellNormals->SetNumberOfTuples(numPolys);

  if (aPrim)
  {
    vtkPoints* p = output->GetPoints();
    vtkNew<vtkIdList> tempPtIds;
    vtkIdType npts;
    const vtkIdType* pts;
    double polyNorm[3];
    for (vtkIdType cellId = 0; cellId < numPolys; cellId++)
    {
      aPrim->GetCellAtId(cellId, npts, pts, tempPtIds);
      vtkPolygon::ComputeNormal(p, static_cast<int>(npts), pts, polyNorm);
      cellNormals->SetTuple(cellId, polyNorm);
    }
  }

  output->GetCellData()->AddArray(cellNormals);
  output->GetCellData()->SetActiveNormals(cellNormals->GetName());
}

//----------------------------------------------------------------------------
void vtkPVGeometryFilter::ExecuteNormalsComputation(vtkPolyData* output)
{
  const auto cellNormals = vtkFloatArray::FastDownCast(output->GetCellData()->GetNormals());
  const bool hasCellNormals = this->GenerateCellNormals ? cellNormals != nullptr : true;
  const auto pointNormals = vtkFloatArray::FastDownCast(output->GetPointData()->GetNormals());
  const bool hasPointNormals = this->GeneratePointNormals ? pointNormals != nullptr : true;
  if (hasPointNormals && hasCellNormals)
  {
    return;
  }
  this->PolyDataNormals->SetInputData(output);
  this->PolyDataNormals->Update();
  this->PolyDataNormals->SetInputData(nullptr);
  output->ShallowCopy(this->PolyDataNormals->GetOutput());
}

//----------------------------------------------------------------------------
void vtkPVGeometryFilter::DataSetExecute(vtkDataSet* input, vtkPolyData* output, int doCommunicate)
{
  double bds[6];
  int procid = 0;

  if (!doCommunicate && input->GetNumberOfPoints() == 0)
  {
    return;
  }

  if (this->Controller)
  {
    procid = this->Controller->GetLocalProcessId();
  }

  input->GetBounds(bds);
  const vtkBoundingBox dataSetBBox(bds);
  if (procid && doCommunicate)
  {
    // Satellite node
    vtkBoundingBox recvBbox;
    this->Controller->Reduce(dataSetBBox, recvBbox, 0);
  }
  else
  {
    if (this->Controller && doCommunicate)
    {
      vtkBoundingBox recvBBox;
      this->Controller->Reduce(dataSetBBox, recvBBox, 0);
      recvBBox.GetBounds(bds);
    }

    if (bds[1] >= bds[0] && bds[3] >= bds[2] && bds[5] >= bds[4])
    {
      // only output in process 0.
      this->OutlineSource->SetBounds(bds);
      this->OutlineSource->Update();

      output->SetPoints(this->OutlineSource->GetOutput()->GetPoints());
      output->SetLines(this->OutlineSource->GetOutput()->GetLines());
    }
  }
}

//----------------------------------------------------------------------------
void vtkPVGeometryFilter::GenericDataSetExecute(
  vtkGenericDataSet* input, vtkPolyData* output, int doCommunicate)
{
  double bds[6];
  int procid = 0;

  if (!this->UseOutline)
  {
    this->OutlineFlag = 0;

    // Geometry filter
    this->GenericGeometryFilter->SetInputData(input);
    this->GenericGeometryFilter->Update();
    output->ShallowCopy(this->GenericGeometryFilter->GetOutput());

    return;
  }

  // Just outline
  this->OutlineFlag = 1;

  if (!doCommunicate && input->GetNumberOfPoints() == 0)
  {
    return;
  }

  if (this->Controller)
  {
    procid = this->Controller->GetLocalProcessId();
  }

  input->GetBounds(bds);
  const vtkBoundingBox dataSetBBox(bds);
  if (procid && doCommunicate)
  {
    // Satellite node
    vtkBoundingBox recvBBox;
    this->Controller->Reduce(dataSetBBox, recvBBox, 0);
  }
  else
  {
    if (doCommunicate)
    {
      vtkBoundingBox recvBBox;
      this->Controller->Reduce(dataSetBBox, recvBBox, 0);
      recvBBox.GetBounds(bds);
    }

    // only output in process 0.
    this->OutlineSource->SetBounds(bds);
    this->OutlineSource->Update();

    output->SetPoints(this->OutlineSource->GetOutput()->GetPoints());
    output->SetLines(this->OutlineSource->GetOutput()->GetLines());
  }
}

//----------------------------------------------------------------------------
void vtkPVGeometryFilter::CellGridExecute(
  vtkCellGrid* input, vtkPolyData* output, int doCommunicate)
{
  double bounds[6];
  int procid = 0;

  if (!doCommunicate && input->GetNumberOfCells() == 0)
  {
    return;
  }

  if (this->Controller)
  {
    procid = this->Controller->GetLocalProcessId();
  }

  input->GetBounds(bounds);
  const vtkBoundingBox dataSetBBox(bounds);
  if (procid && doCommunicate)
  {
    // Satellite node
    vtkBoundingBox recvBbox;
    this->Controller->Reduce(dataSetBBox, recvBbox, 0);
  }
  else
  {
    if (this->Controller && doCommunicate)
    {
      vtkBoundingBox recvBBox;
      this->Controller->Reduce(dataSetBBox, recvBBox, 0);
      recvBBox.GetBounds(bounds);
    }

    if (bounds[1] >= bounds[0] && bounds[3] >= bounds[2] && bounds[5] >= bounds[4])
    {
      // only output in process 0.
      this->OutlineSource->SetBounds(bounds);
      this->OutlineSource->Update();

      output->SetPoints(this->OutlineSource->GetOutput()->GetPoints());
      output->SetLines(this->OutlineSource->GetOutput()->GetLines());
    }
  }
}

//----------------------------------------------------------------------------
void vtkPVGeometryFilter::ImageDataExecute(
  vtkImageData* input, vtkPolyData* output, int doCommunicate, int updatePiece, const int* ext)
{
  // If doCommunicate is false, use extent because the block is
  // entirely contained in this process.
  if (!doCommunicate)
  {
    ext = input->GetExtent();
  }

  if (!this->UseOutline)
  {
    if (input->GetNumberOfCells() > 0)
    {
      this->GeometryFilter->StructuredExecute(
        input, output, const_cast<int*>(ext), nullptr, nullptr);
    }
    this->OutlineFlag = 0;
    return;
  }
  this->OutlineFlag = 1;

  // Otherwise, let OutlineSource do all the work
  if (ext[1] >= ext[0] && ext[3] >= ext[2] && ext[5] >= ext[4] &&
    (updatePiece == 0 || !doCommunicate))
  {
    double* spacing = input->GetSpacing();
    double* origin = input->GetOrigin();
    double bounds[6];
    bounds[0] = spacing[0] * ((float)ext[0]) + origin[0];
    bounds[1] = spacing[0] * ((float)ext[1]) + origin[0];
    bounds[2] = spacing[1] * ((float)ext[2]) + origin[1];
    bounds[3] = spacing[1] * ((float)ext[3]) + origin[1];
    bounds[4] = spacing[2] * ((float)ext[4]) + origin[2];
    bounds[5] = spacing[2] * ((float)ext[5]) + origin[2];

    this->OutlineSource->SetBounds(bounds);
    this->OutlineSource->Update();
    auto outline = this->OutlineSource->GetOutput();

    output->SetPoints(outline->GetPoints());
    output->SetLines(outline->GetLines());
    output->SetPolys(outline->GetPolys());
  }
  else
  {
    vtkNew<vtkPoints> pts;
    output->SetPoints(pts);
  }
}

//----------------------------------------------------------------------------
void vtkPVGeometryFilter::StructuredGridExecute(vtkStructuredGrid* input, vtkPolyData* output,
  int updatePiece, int updateNumPieces, int updateGhosts, const int* wholeExtentArg)
{
  int wholeExtent[6];
  ::GetValidWholeExtent(input, wholeExtentArg, wholeExtent);

  if (!this->UseOutline)
  {
    if (input->GetNumberOfCells() > 0)
    {
      this->GeometryFilter->StructuredExecute(input, output, wholeExtent, nullptr, nullptr);
    }
    this->OutlineFlag = 0;
    return;
  }
  this->OutlineFlag = 1;

  vtkNew<vtkPVTrivialProducer> producer;
  producer->SetOutput(input);
  producer->SetWholeExtent(
    wholeExtent[0], wholeExtent[1], wholeExtent[2], wholeExtent[3], wholeExtent[4], wholeExtent[5]);

  vtkNew<vtkStructuredGridOutlineFilter> outline;
  outline->SetInputConnection(producer->GetOutputPort());
  outline->UpdatePiece(updatePiece, updateNumPieces, updateGhosts);
  output->CopyStructure(outline->GetOutput());
}

//----------------------------------------------------------------------------
void vtkPVGeometryFilter::RectilinearGridExecute(vtkRectilinearGrid* input, vtkPolyData* output,
  int vtkNotUsed(updatePiece), int vtkNotUsed(updateNumPieces), int vtkNotUsed(updateGhosts),
  const int* wholeExtent)
{
  if (!this->UseOutline)
  {
    if (input->GetNumberOfCells() > 0)
    {
      this->GeometryFilter->StructuredExecute(
        input, output, const_cast<int*>(wholeExtent), nullptr, nullptr);
    }
    this->OutlineFlag = 0;
    return;
  }
  this->OutlineFlag = 1;

  vtkNew<vtkPVTrivialProducer> producer;
  producer->SetOutput(input);
  producer->SetWholeExtent(
    wholeExtent[0], wholeExtent[1], wholeExtent[2], wholeExtent[3], wholeExtent[4], wholeExtent[5]);

  vtkNew<vtkRectilinearGridOutlineFilter> outline;
  outline->SetInputConnection(producer->GetOutputPort());
  outline->Update();
  output->CopyStructure(outline->GetOutput());
}

//----------------------------------------------------------------------------
void vtkPVGeometryFilter::UnstructuredGridExecute(
  vtkUnstructuredGridBase* input, vtkPolyData* output, int doCommunicate)
{
  if (!this->UseOutline)
  {
    this->OutlineFlag = 0;

    bool handleSubdivision = (this->Triangulate != 0) && (input->GetNumberOfCells() > 0);
    if (!handleSubdivision && (this->NonlinearSubdivisionLevel > 0))
    {
      // Check to see if the data actually has nonlinear cells.  Handling
      // nonlinear cells adds unnecessary work if we only have linear cells.
      if (input->GetNumberOfCells() > 0)
      {
        auto helper = vtkGeometryFilterHelper::CharacterizeUnstructuredGrid(input);
        if (!helper->IsLinear)
        {
          handleSubdivision = true;
        }
        delete helper;
      }
    }

    vtkSmartPointer<vtkIdTypeArray> facePtIds2OriginalPtIds;

    auto inputClone = vtk::TakeSmartPointer(input->NewInstance());
    inputClone->ShallowCopy(input);
    input = inputClone;

    if (handleSubdivision)
    {
      // Use the vtkUnstructuredGridGeometryFilter to extract 2D surface cells
      // from the geometry.  This is important to extract an appropriate
      // wireframe in vtkRecoverGeometryWireframe.  Also, at the time of this
      // writing vtkGeometryFilter only properly subdivides 2D cells past level 1.
      this->UnstructuredGridGeometryFilter->SetInputData(input);
      // Let the vtkUnstructuredGridGeometryFilter record from which point and
      // cell each face comes from in the standard vtkOriginalCellIds array.
      this->UnstructuredGridGeometryFilter->SetPassThroughCellIds(this->PassThroughCellIds);
      this->UnstructuredGridGeometryFilter->SetPassThroughPointIds(this->PassThroughPointIds);
      this->UnstructuredGridGeometryFilter->SetMatchBoundariesIgnoringCellOrder(
        this->MatchBoundariesIgnoringCellOrder);
      // Turn off ghost cell clipping. This ensures that ghost cells are retained
      // and handed to the GeometryFilter to ensure only valid faces are
      // generated. If this weren't here, then the GeometryFilter would
      // generate boundary faces between normal cells and where the ghost cells
      // used to be, which is not correct.
      this->UnstructuredGridGeometryFilter->DuplicateGhostCellClippingOff();
      // Disable point merging as it may prevent the correct visualization
      // of non-continuous attributes.
      this->UnstructuredGridGeometryFilter->MergingOff();
      this->UnstructuredGridGeometryFilter->Update();

      this->UnstructuredGridGeometryFilter->SetInputData(nullptr);
      // Feed the extracted surface as the input to the rest of the processing.
      input->ShallowCopy(this->UnstructuredGridGeometryFilter->GetOutput());

      // Keep a handle to the vtkOriginalPointIds array.  We might need it.
      facePtIds2OriginalPtIds =
        vtkIdTypeArray::SafeDownCast(input->GetPointData()->GetArray("vtkOriginalPointIds"));

      // Flag the data set surface filter to record original cell ids, but do it
      // in a specially named array that vtkRecoverGeometryWireframe will later
      // use.  Note that because the data set comes from
      // UnstructuredGridGeometryFilter, the ids will represent the faces rather
      // than the original cells, which is important.
      this->GeometryFilter->PassThroughCellIdsOn();
      this->GeometryFilter->SetOriginalCellIdsName(details::ORIGINAL_FACE_IDS);

      if (this->PassThroughPointIds)
      {
        // vtkGeometryFilter is going to strip the vtkOriginalPointIds
        // created by the vtkUnstructuredGridGeometryFilter because it
        // cannot interpolate the ids.  Make the vtkGeometryFilter make
        // its own original ids array.  We will resolve them later.
        this->GeometryFilter->PassThroughPointIdsOn();
      }
    }

    if (input->GetNumberOfCells() > 0)
    {
      this->GeometryFilter->UnstructuredGridExecute(input, output);
    }

    if (this->Triangulate && (output->GetNumberOfPolys() > 0))
    {
      // Triangulate the polygonal mesh if requested to avoid rendering
      // issues of non-convex polygons.
      vtkNew<vtkTriangleFilter> triangleFilter;
      triangleFilter->SetInputData(output);
      triangleFilter->Update();
      output->ShallowCopy(triangleFilter->GetOutput());
    }

    if (handleSubdivision && !this->GenerateFeatureEdges)
    {
      // Restore state of GeometryFilter.
      this->GeometryFilter->SetPassThroughCellIds(this->PassThroughCellIds);
      this->GeometryFilter->SetOriginalCellIdsName(nullptr);
      this->GeometryFilter->SetPassThroughPointIds(this->PassThroughPointIds);

      this->GeometryFilter->SetMatchBoundariesIgnoringCellOrder(
        this->MatchBoundariesIgnoringCellOrder);

      // Now use vtkRecoverGeometryWireframe to create an edge flag attribute
      // that will cause the wireframe to be rendered correctly.
      vtkNew<vtkPolyData> nextStageInput;
      nextStageInput->ShallowCopy(output); // Yes output is correct.
      this->RecoverWireframeFilter->SetInputData(nextStageInput.Get());
      this->RecoverWireframeFilter->SetCellIdsAttribute(details::ORIGINAL_FACE_IDS);
      this->RecoverWireframeFilter->Update();
      this->RecoverWireframeFilter->SetInputData(nullptr);

      // Get what should be the final output.
      output->ShallowCopy(this->RecoverWireframeFilter->GetOutput());

      if (this->PassThroughPointIds)
      {
        // The output currently has a vtkOriginalPointIds array that maps points
        // to the data containing only the faces.  Correct this to point to the
        // original data set.
        auto polyPtIds2FacePtIds =
          vtkIdTypeArray::SafeDownCast(output->GetPointData()->GetArray("vtkOriginalPointIds"));
        if (!facePtIds2OriginalPtIds || !polyPtIds2FacePtIds)
        {
          vtkErrorMacro(<< "Missing original point id arrays.");
          return;
        }
        const vtkIdType numPts = polyPtIds2FacePtIds->GetNumberOfTuples();
        vtkNew<vtkIdTypeArray> polyPtIds2OriginalPtIds;
        polyPtIds2OriginalPtIds->SetName("vtkOriginalPointIds");
        polyPtIds2OriginalPtIds->SetNumberOfComponents(1);
        polyPtIds2OriginalPtIds->SetNumberOfTuples(numPts);
        for (vtkIdType polyPtId = 0; polyPtId < numPts; polyPtId++)
        {
          const vtkIdType facePtId = polyPtIds2FacePtIds->GetValue(polyPtId);
          vtkIdType originalPtId = -1;
          if (facePtId >= 0)
          {
            originalPtId = facePtIds2OriginalPtIds->GetValue(facePtId);
          }
          polyPtIds2OriginalPtIds->SetValue(polyPtId, originalPtId);
        }
        output->GetPointData()->AddArray(polyPtIds2OriginalPtIds.Get());
      }
    }

    output->GetCellData()->RemoveArray(details::ORIGINAL_FACE_IDS);
    return;
  }

  this->OutlineFlag = 1;

  this->DataSetExecute(input, output, doCommunicate);
}

//----------------------------------------------------------------------------
void vtkPVGeometryFilter::PolyDataExecute(
  vtkPolyData* input, vtkPolyData* output, int doCommunicate)
{
  if (!this->UseOutline)
  {
    this->OutlineFlag = 0;
    output->ShallowCopy(input);
    if (this->PassThroughCellIds)
    {
      vtkNew<vtkIdTypeArray> originalCellIds;
      originalCellIds->SetName("vtkOriginalCellIds");
      originalCellIds->SetNumberOfComponents(1);
      vtkNew<vtkIdTypeArray> originalFaceIds;
      originalFaceIds->SetName(details::ORIGINAL_FACE_IDS);
      originalFaceIds->SetNumberOfComponents(1);
      vtkCellData* outputCD = output->GetCellData();
      outputCD->AddArray(originalCellIds.Get());
      if (this->Triangulate)
      {
        outputCD->AddArray(originalFaceIds.Get());
      }
      const vtkIdType numberOfCells = output->GetNumberOfCells();
      originalCellIds->SetNumberOfValues(numberOfCells);
      originalFaceIds->SetNumberOfValues(numberOfCells);
      std::iota(originalCellIds->GetPointer(0), originalCellIds->GetPointer(numberOfCells), 0);
      std::iota(originalFaceIds->GetPointer(0), originalFaceIds->GetPointer(numberOfCells), 0);
    }
    if (this->PassThroughPointIds)
    {
      vtkNew<vtkIdTypeArray> originalPointIds;
      originalPointIds->SetName("vtkOriginalPointIds");
      originalPointIds->SetNumberOfComponents(1);
      vtkPointData* outputPD = output->GetPointData();
      outputPD->AddArray(originalPointIds.Get());
      const vtkIdType numberOfPoints = output->GetNumberOfPoints();
      originalPointIds->SetNumberOfValues(numberOfPoints);
      std::iota(originalPointIds->GetPointer(0), originalPointIds->GetPointer(numberOfPoints), 0);
    }

    if (this->Triangulate)
    {
      // Triangulate the polygonal mesh.
      vtkNew<vtkTriangleFilter> triangleFilter;
      triangleFilter->SetInputData(output);
      triangleFilter->Update();

      // Now use vtkRecoverGeometryWireframe to create an edge flag attribute
      // that will cause the wireframe to be rendered correctly.
      this->RecoverWireframeFilter->SetInputData(triangleFilter->GetOutput());
      this->RecoverWireframeFilter->Update();
      this->RecoverWireframeFilter->SetInputData(nullptr);

      // Get what should be the final output.
      output->ShallowCopy(this->RecoverWireframeFilter->GetOutput());

      output->GetCellData()->RemoveArray(details::ORIGINAL_FACE_IDS);
    }
    return;
  }

  this->OutlineFlag = 1;
  this->DataSetExecute(input, output, doCommunicate);
}

//----------------------------------------------------------------------------
void vtkPVGeometryFilter::HyperTreeGridExecute(
  vtkHyperTreeGrid* input, vtkPolyData* output, int doCommunicate)
{
  if (!this->UseOutline)
  {
    this->OutlineFlag = 0;

    vtkNew<vtkHyperTreeGridGeometry> internalFilter;
    vtkNew<vtkHyperTreeGrid> htgCopy;
    htgCopy->ShallowCopy(input);
    internalFilter->SetInputData(htgCopy);
    internalFilter->SetPassThroughCellIds(this->PassThroughCellIds);
    internalFilter->SetOriginalCellIdArrayName("vtkOriginalCellIds");
    internalFilter->Update();
    output->ShallowCopy(internalFilter->GetOutput());
    return;
  }

  this->OutlineFlag = 1;
  double bds[6];
  int procid = 0;
  if (!doCommunicate && input->GetNumberOfCells() == 0)
  {
    return;
  }

  if (this->Controller)
  {
    procid = this->Controller->GetLocalProcessId();
  }

  input->GetBounds(bds);
  const vtkBoundingBox dataSetBBox(bds);
  if (procid && doCommunicate)
  {
    // Satellite node
    vtkBoundingBox recvBBox;
    this->Controller->Reduce(dataSetBBox, recvBBox, 0);
  }
  else
  {
    if (this->Controller && doCommunicate)
    {
      vtkBoundingBox recvBBox;
      this->Controller->Reduce(dataSetBBox, recvBBox, 0);
      recvBBox.GetBounds(bds);
    }

    if (bds[1] >= bds[0] && bds[3] >= bds[2] && bds[5] >= bds[4])
    {
      // only output in process 0.
      this->OutlineSource->SetBounds(bds);
      this->OutlineSource->Update();

      output->SetPoints(this->OutlineSource->GetOutput()->GetPoints());
      output->SetLines(this->OutlineSource->GetOutput()->GetLines());
    }
  }
}

//----------------------------------------------------------------------------
void vtkPVGeometryFilter::ExplicitStructuredGridExecute(
  vtkExplicitStructuredGrid* input, vtkPolyData* out, int doCommunicate, const int* wholeExtent)
{
  vtkNew<vtkPVTrivialProducer> producer;
  producer->SetOutput(input);
  producer->SetWholeExtent(
    wholeExtent[0], wholeExtent[1], wholeExtent[2], wholeExtent[3], wholeExtent[4], wholeExtent[5]);
  producer->Update();

  if (!this->UseOutline)
  {
    this->OutlineFlag = 0;

    vtkNew<vtkExplicitStructuredGridSurfaceFilter> internalFilter;
    internalFilter->SetPassThroughPointIds(this->PassThroughPointIds);
    internalFilter->SetPassThroughCellIds(this->PassThroughCellIds);
    internalFilter->SetInputConnection(producer->GetOutputPort());
    internalFilter->Update();
    out->ShallowCopy(internalFilter->GetOutput());
    return;
  }
  auto in = vtkExplicitStructuredGrid::SafeDownCast(producer->GetOutputDataObject(0));

  this->OutlineFlag = 1;
  this->DataSetExecute(in, out, doCommunicate);
}

//----------------------------------------------------------------------------
int vtkPVGeometryFilter::FillInputPortInformation(int port, vtkInformation* info)
{
  if (!this->Superclass::FillInputPortInformation(port, info))
  {
    return 0;
  }

  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
  info->Append(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkGenericDataSet");
  info->Append(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkCompositeDataSet");
  info->Append(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkHyperTreeGrid");
  info->Append(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkCellGrid");
  return 1;
}

//-----------------------------------------------------------------------------
void vtkPVGeometryFilter::ReportReferences(vtkGarbageCollector* collector)
{
  this->Superclass::ReportReferences(collector);
  vtkGarbageCollectorReport(collector, this->GeometryFilter, "GeometryFilter");
  vtkGarbageCollectorReport(collector, this->GenericGeometryFilter, "GenericGeometryFilter");
  vtkGarbageCollectorReport(
    collector, this->UnstructuredGridGeometryFilter, "UnstructuredGridGeometryFilter");
  vtkGarbageCollectorReport(collector, this->RecoverWireframeFilter, "RecoverWireframeFilter");
}

//----------------------------------------------------------------------------
void vtkPVGeometryFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);

  os << indent << "OutlineFlag: " << (this->OutlineFlag ? "on" : "off") << endl;
  os << indent << "UseOutline: " << (this->UseOutline ? "on" : "off") << endl;
  os << indent << "GenerateFeatureEdges: " << (this->GenerateFeatureEdges ? "on" : "off") << endl;
  os << indent << "BlockColorsDistinctValues: " << this->BlockColorsDistinctValues << endl;
  os << indent << "GenerateCellNormals: " << (this->GenerateCellNormals ? "on" : "off") << endl;
  os << indent << "GeneratePointNormals: " << (this->GeneratePointNormals ? "on" : "off") << endl;
  os << indent << "Splitting: " << (this->Splitting ? "on" : "off") << endl;
  os << indent << "FeatureAngle: " << this->FeatureAngle << endl;
  os << indent << "Triangulate: " << (this->Triangulate ? "on" : "off") << endl;
  os << indent << "NonlinearSubdivisionLevel: " << this->NonlinearSubdivisionLevel << endl;
  os << indent << "MatchBoundariesIgnoringCellOrder: " << this->MatchBoundariesIgnoringCellOrder
     << endl;
  os << indent << "Controller: " << this->Controller << endl;
  os << indent << "PassThroughCellIds: " << (this->PassThroughCellIds ? "on" : "off") << endl;
  os << indent << "PassThroughPointIds: " << (this->PassThroughPointIds ? "on" : "off") << endl;
  os << indent << "GenerateProcessIds: " << (this->GenerateProcessIds ? "on" : "off") << endl;
  os << indent << "HideInternalAMRFaces: " << (this->HideInternalAMRFaces ? "on" : "off") << endl;
  os << indent << "UseNonOverlappingAMRMetaDataForOutlines: "
     << (this->UseNonOverlappingAMRMetaDataForOutlines ? "on" : "off") << endl;
}

//----------------------------------------------------------------------------
void vtkPVGeometryFilter::SetGenerateFeatureEdges(bool val)
{
  if (this->GenerateFeatureEdges != val)
  {
    this->GenerateFeatureEdges = val;
    this->Modified();
  }
}

//----------------------------------------------------------------------------
void vtkPVGeometryFilter::SetGenerateCellNormals(int val)
{
  if (this->GenerateCellNormals != static_cast<bool>(val))
  {
    this->GenerateCellNormals = val;
    if (this->PolyDataNormals)
    {
      this->PolyDataNormals->SetComputeCellNormals(this->GenerateCellNormals);
    }
    this->Modified();
  }
}

//----------------------------------------------------------------------------
void vtkPVGeometryFilter::SetGeneratePointNormals(bool val)
{
  if (this->GeneratePointNormals != val)
  {
    this->GeneratePointNormals = val;
    if (this->PolyDataNormals)
    {
      this->PolyDataNormals->SetComputePointNormals(this->GeneratePointNormals);
    }
    this->Modified();
  }
}

//----------------------------------------------------------------------------
void vtkPVGeometryFilter::SetSplitting(bool val)
{
  if (this->Splitting != val)
  {
    this->Splitting = val;
    if (this->PolyDataNormals)
    {
      this->PolyDataNormals->SetSplitting(this->Splitting);
    }
    this->Modified();
  }
}

//----------------------------------------------------------------------------
void vtkPVGeometryFilter::SetFeatureAngle(double val)
{
  if (this->FeatureAngle != val)
  {
    this->FeatureAngle = val;
    if (this->PolyDataNormals)
    {
      this->PolyDataNormals->SetFeatureAngle(this->FeatureAngle);
    }
    this->Modified();
  }
}

//----------------------------------------------------------------------------
void vtkPVGeometryFilter::SetPassThroughCellIds(int newvalue)
{
  this->PassThroughCellIds = newvalue;
  if (this->GeometryFilter)
  {
    this->GeometryFilter->SetPassThroughCellIds(this->PassThroughCellIds);
  }
  if (this->GenericGeometryFilter)
  {
    this->GenericGeometryFilter->SetPassThroughCellIds(this->PassThroughCellIds);
  }
}

//----------------------------------------------------------------------------
void vtkPVGeometryFilter::SetPassThroughPointIds(int newvalue)
{
  this->PassThroughPointIds = newvalue;
  if (this->GeometryFilter)
  {
    this->GeometryFilter->SetPassThroughPointIds(this->PassThroughPointIds);
  }
}

//-----------------------------------------------------------------------------
void vtkPVGeometryFilter::SetNonlinearSubdivisionLevel(int newvalue)
{
  if (this->NonlinearSubdivisionLevel != newvalue)
  {
    this->NonlinearSubdivisionLevel = newvalue;
    if (this->GeometryFilter)
    {
      this->GeometryFilter->SetNonlinearSubdivisionLevel(this->NonlinearSubdivisionLevel);
    }
    this->Modified();
  }
}

//-----------------------------------------------------------------------------
void vtkPVGeometryFilter::SetMatchBoundariesIgnoringCellOrder(int newvalue)
{
  if (this->MatchBoundariesIgnoringCellOrder != newvalue)
  {
    this->MatchBoundariesIgnoringCellOrder = newvalue;
    if (this->GeometryFilter)
    {
      this->GeometryFilter->SetMatchBoundariesIgnoringCellOrder(
        this->MatchBoundariesIgnoringCellOrder);
    }
    this->Modified();
  }
}
