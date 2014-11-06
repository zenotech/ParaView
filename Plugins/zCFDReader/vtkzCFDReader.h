/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkPOpenFOAMReader.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkPOpenFOAMReader - reads a decomposed dataset in OpenFOAM format
// .SECTION Description
// vtkPOpenFOAMReader creates a multiblock dataset. It reads
// parallel-decomposed mesh information and time dependent data.  The
// polyMesh folders contain mesh information. The time folders contain
// transient data for the cells. Each folder can contain any number of
// data files.

// .SECTION Thanks
// This class was developed by Takuya Oshima at Niigata University,
// Japan (oshima@eng.niigata-u.ac.jp).

#ifndef __vtkzCFDReader_h
#define __vtkzCFDReader_h

#include <map>
#include <set>

//#include "vtkOpenFOAMReader.h"
#include "vtkPolyDataAlgorithm.h"
#include "vtkMultiBlockDataSetAlgorithm.h"
#include "vtkStringArray.h"

class vtkDataArraySelection;
class vtkMultiProcessController;

class vtkzCFDReader : public vtkMultiBlockDataSetAlgorithm
{
public:

  //ETX
  static vtkzCFDReader *New();
  vtkTypeMacro(vtkzCFDReader, vtkMultiBlockDataSetAlgorithm);

  void PrintSelf(ostream &os, vtkIndent indent);

  // Description:
  // Set and get case type. 0 = decomposed case, 1 = reconstructed case.
  //void SetCaseType(const int t);
  //vtkGetMacro(CaseType, caseType);
  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);

  vtkSetStringMacro(CaseFileName);
  vtkGetStringMacro(CaseFileName);


  void SetRefresh() { this->Refresh = true; this->Modified(); }

  // Description:
  // Get the number of Patches (including Internal Mesh) available in the input.
  int GetNumberOfZoneArrays(void)
  { return this->GetNumberOfSelectionArrays(this->ZoneDataArraySelection); }

  // Description:
  // Get/Set whether the Patch with the given name is to
  // be read.
  int GetZoneArrayStatus(const char *name)
  { return this->GetSelectionArrayStatus(this->ZoneDataArraySelection, name); }
  void SetZoneArrayStatus(const char *name, int status)
  { this->SetSelectionArrayStatus(this->ZoneDataArraySelection, name,
    status); }

  // Description:
  // Get the name of the Patch with the given index in
  // the input.
  const char *GetZoneArrayName(int index)
  { return this->GetSelectionArrayName(this->ZoneDataArraySelection, index); }


  int GetNumberOfCellArrays(void)
  { return this->GetNumberOfSelectionArrays(this->CellDataArraySelection); }

  int GetCellArrayStatus(const char *name)
  { return this->GetSelectionArrayStatus(this->CellDataArraySelection, name); }
  void SetCellArrayStatus(const char *name, int status)
  { this->SetSelectionArrayStatus(this->CellDataArraySelection, name,
    status); }

  const char *GetCellArrayName(int index)
   { return this->GetSelectionArrayName(this->CellDataArraySelection, index); }



  // Description:
  // Set and get the controller.
  virtual void SetController(vtkMultiProcessController *);
  vtkGetObjectMacro(Controller, vtkMultiProcessController);

protected:
  vtkzCFDReader();
  ~vtkzCFDReader();

  void ReadPython(std::map<int,std::string> &zoneToBc);

  int RequestInformation(vtkInformation *, vtkInformationVector **,
    vtkInformationVector *);
  int RequestData(vtkInformation *, vtkInformationVector **,
    vtkInformationVector *);

  int GetNumberOfSelectionArrays(vtkDataArraySelection *s);

  int GetSelectionArrayStatus(vtkDataArraySelection *s,
      const char *name);

  void SetSelectionArrayStatus(vtkDataArraySelection *s,
      const char* name, int status);

  const char *GetSelectionArrayName(vtkDataArraySelection *s,
      int index);

private:
  vtkMultiProcessController *Controller;
  //caseType CaseType;
  unsigned long MTimeOld;
  int NumProcesses;
  int ProcessId;

  char *FileName;
  char *CaseFileName;
  vtkStdString *FileNameOld;

  vtkStdString *ProblemName;
  vtkStdString *CaseName;

  // refresh flag
  bool Refresh;

  vtkDataArraySelection *ZoneDataArraySelection;
  std::map<int,std::set<int> > selectionToZone;

  vtkDataArraySelection *CellDataArraySelection;
  std::map<int,std::set<int> > selectionToData;

  vtkzCFDReader(const vtkzCFDReader &); // Not implemented.
  void operator=(const vtkzCFDReader &); // Not implemented.

  void GatherMetaData();
  void BroadcastStatus(int &);
  void Broadcast(vtkStringArray *);
  void AllGather(vtkStringArray *);
  void AllGather(vtkDataArraySelection *);

  bool ReadStatusFile(const char *fileName);

};

#endif
