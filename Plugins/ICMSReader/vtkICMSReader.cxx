#include "vtkICMSReader.h"
 
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkInformationVector.h"
#include "vtkInformation.h"
#include "vtkDataObject.h"
#include "vtkSmartPointer.h"
#include "vtkVertexGlyphFilter.h"
 
#include <sstream>

vtkStandardNewMacro(vtkICMSReader);
 
vtkICMSReader::vtkICMSReader()
{
  this->FileName = NULL;
  this->SetNumberOfInputPorts(0);
  this->SetNumberOfOutputPorts(1);
}
 
int vtkICMSReader::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *outputVector)
{
 
  // get the info object
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
 
  // get the ouptut
   vtkPolyData *output = vtkPolyData::SafeDownCast(
            outInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();

  // Open file
  std::ifstream in(FileName);

  std::string title;
  std::getline(in,title);
  // Check for empty file
  if(in.eof())
    return 0;

  int np, nf, numVar;
  in >> np >> nf >> numVar;

  std::string line;
  std::getline(in,line); // read end of line
  std::getline(in,line); // read variable list

  // Check for unstructured ICMS
  if(np > 0)
    return 0;

  np *= -1;

  polydata->Allocate(nf);

  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  points->SetNumberOfPoints(np);

  float xyz[3];
  for(int i=0; i<np; i++)
  {
    in >> xyz[0]  >> xyz[1]  >> xyz[2];
    getline( in, line );  // read end of line
    points->InsertPoint(i, xyz[0], xyz[1], xyz[2]);
  }

  polydata->SetPoints(points);
 
  std::vector<int> faceNodes;
  //faceNodes.reserve( 10 );
  vtkIdType cellList[100];

  for(int i=0; i<nf; i++)
  {
    //faceNodes.clear();
    getline(in, line);  // read end of line
    //util::IoUtils::removeTrailingWhiteSpace(line);
    std::stringstream istr(  line );

    int numPts = 0;

    while(!istr.eof()){
      int node = 0;
      istr >> node;
      if(node == 0)
        break;

      //faceNodes.push_back(node-1);
      cellList[numPts] = node-1;
      numPts++;

    }
    int vtkCellType = VTK_EMPTY_CELL;
    if (numPts == 3)
      vtkCellType = VTK_TRIANGLE;
    else if (numPts == 4)
      vtkCellType = VTK_QUAD;
    else
      vtkCellType = VTK_POLYGON;

    polydata->InsertNextCell(vtkCellType, numPts, cellList);

  }

  output->ShallowCopy(polydata);
 
  return 1;
}
 
void vtkICMSReader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
 
  os << indent << "File Name: "
      << (this->FileName ? this->FileName : "(none)") << "\n";
}

