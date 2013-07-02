#pragma once
 
#include "vtkPolyDataAlgorithm.h"
#include "vtkMultiBlockDataSetAlgorithm.h"

//class vtkICMSReader : public vtkPolyDataAlgorithm
class vtkICMSReader : public vtkMultiBlockDataSetAlgorithm
{
public:
  //vtkTypeMacro(vtkICMSReader,vtkPolyDataAlgorithm);
  vtkTypeMacro(vtkICMSReader,vtkMultiBlockDataSetAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);
 
  static vtkICMSReader *New();
 
  // Description:
  // Specify file name of the .abc file.
  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);
 
protected:
  vtkICMSReader();
  ~vtkICMSReader(){}
 
  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
 
private:
  vtkICMSReader(const vtkICMSReader&);  // Not implemented.
  void operator=(const vtkICMSReader&);  // Not implemented.
 
  char* FileName;
};
 

