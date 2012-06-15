/*=========================================================================

  Program:   ParaView
  Module:    vtkPVMetaClipDataSet.cxx

  Copyright (c) Kitware, Inc.
  All rights reserved.
  See Copyright.txt or http://www.paraview.org/HTML/Copyright.html for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkPVMetaClipDataSet.h"

#include "vtkAlgorithm.h"
#include "vtkPVClipDataSet.h"
#include "vtkDataObject.h"
#include "vtkExtractGeometry.h"
#include "vtkImplicitFunction.h"
#include "vtkInformation.h"
#include "vtkNew.h"
#include "vtkObjectFactory.h"
#include "vtkSmartPointer.h"

class vtkPVMetaClipDataSet::vtkInternals
{
public:
  vtkNew<vtkPVClipDataSet> Clip;
  vtkNew<vtkExtractGeometry> ExtractCells;

  vtkInternals()
  {
    this->ExtractCells->SetExtractInside(1);
    this->ExtractCells->SetExtractOnlyBoundaryCells(0);
    this->ExtractCells->SetExtractBoundaryCells(1);
  }

};
//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkPVMetaClipDataSet);
//----------------------------------------------------------------------------
vtkPVMetaClipDataSet::vtkPVMetaClipDataSet()
{
  // Setup default configuration
  this->SetOutputType(VTK_UNSTRUCTURED_GRID);

  this->Internal = new vtkInternals();

  this->RegisterFilter(this->Internal->Clip.GetPointer());
  this->RegisterFilter(this->Internal->ExtractCells.GetPointer());

  this->Superclass::SetActiveFilter(0);
}

//----------------------------------------------------------------------------
vtkPVMetaClipDataSet::~vtkPVMetaClipDataSet()
{
  delete this->Internal; this->Internal = NULL;
}

//----------------------------------------------------------------------------
void vtkPVMetaClipDataSet::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}

//----------------------------------------------------------------------------
void vtkPVMetaClipDataSet::SetImplicitFunction(vtkImplicitFunction* func)
{
  this->Internal->Clip->SetClipFunction(func);
  this->Internal->ExtractCells->SetImplicitFunction(func);
  this->Modified();
}

//----------------------------------------------------------------------------
void vtkPVMetaClipDataSet::SetValue( double value)
{
  this->Internal->Clip->SetValue(value);
  this->Modified();
}
//----------------------------------------------------------------------------
void vtkPVMetaClipDataSet::SetUseValueAsOffset(int value)
{
  this->Internal->Clip->SetUseValueAsOffset(value);
  this->Modified();
}
//----------------------------------------------------------------------------
void vtkPVMetaClipDataSet::PreserveInputCells(int keepCellAsIs)
{
  this->SetActiveFilter(keepCellAsIs);
}
//----------------------------------------------------------------------------
void vtkPVMetaClipDataSet::SetInputArrayToProcess(int idx, int port, int connection,
                            int fieldAssociation,
                            const char *name)
{
  this->Internal->Clip->SetInputArrayToProcess(idx, port, connection, fieldAssociation, name);
  this->Modified();
}
//----------------------------------------------------------------------------

void vtkPVMetaClipDataSet::SetInputArrayToProcess(int idx, int port, int connection,
                            int fieldAssociation,
                            int fieldAttributeType)
{
  this->Internal->Clip->SetInputArrayToProcess(idx, port, connection, fieldAssociation, fieldAttributeType);
  this->Modified();
}
//----------------------------------------------------------------------------

void vtkPVMetaClipDataSet::SetInputArrayToProcess(int idx, vtkInformation *info)
{
  this->Internal->Clip->SetInputArrayToProcess(idx, info);
  this->Modified();
}//----------------------------------------------------------------------------

void vtkPVMetaClipDataSet::SetInsideOut(int insideOut)
{
  this->Internal->Clip->SetInsideOut(insideOut);
  this->Internal->ExtractCells->SetExtractInside(insideOut);
  this->Modified();
}
