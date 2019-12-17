
#include <set>
#include <map>

#include "vtkzCFDReader.h"

#include "vtkAppendCompositeDataLeaves.h"
#include "vtkCharArray.h"
#include "vtkCollection.h"
#include "vtkDataArraySelection.h"
#include "vtkDirectory.h"
#include "vtkDoubleArray.h"
#include "vtkFieldData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkIntArray.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkUnstructuredGrid.h"
#include "vtkMultiProcessController.h"
#include "vtkObjectFactory.h"
#include "vtkSortDataArray.h"
#include "vtkStdString.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkStringArray.h"
#include <vtkCellType.h>
#include <vtkCellLocator.h>
#include <vtkIntArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>

#include <vtksys/SystemTools.hxx>

#include <vtk_jsoncpp.h>

#include <hdf5/hdffile.hpp>
#include <hdf5/hdfdataset.hpp>

#include <Python.h>
#include <boost/python.hpp>



namespace detail
{

  struct BoundaryCondition
  {

  enum {
    NONE=0x0,INTERIOR=0x2,WALL=0x3,INFLOW=0x4,OUTFLOW=0x5,SYMMETRY=0x7,FARFIELD=0x9,PERIODIC=0x12
    };

    static std::string toString(int bc);

  };

  std::string BoundaryCondition::toString(int bc)
  {
    std::string name = "UNKNOWN";

    switch(bc)
    {
    case NONE:
    case INTERIOR:
      name = "interior";
      break;
    case WALL:
      name = "wall";
      break;
    case SYMMETRY:
      name = "symmetry";
      break;
    case FARFIELD:
      name = "farfield";
      break;
    case INFLOW:
      name = "inflow";
      break;
    case OUTFLOW:
      name = "outflow";
      break;
    case PERIODIC:
      name = "periodic";
      break;
    }
    return name;
  }

  /** Cell type encoder
   *  The first five bits are used to encode the fundamental cell type
   *  The next 11bits encode the number of faces (max 2047 faces)
   *  The remaining bits encode the order
   */
  struct CellType
  {
    enum{
      TETRA=0x1,PRISM=0x2,PYRA=0x4,HEX=0x8,POLY=0x10,
      NUM_TYPE=0x5,TYPE_MASK=0x1F,NUM_FACE_MASK=0x7FF
    };

    static int getType(int cellType)
    {
      return cellType & TYPE_MASK;
    }

    static int getNumFaces(int cellType)
    {
      return ( cellType >> NUM_TYPE ) & NUM_FACE_MASK;
    }

    static int encode(int cellType, int numFaces)
    {
      return (numFaces << 5) + cellType;
    }

  };

}

vtkStandardNewMacro(vtkzCFDReader);
vtkCxxSetObjectMacro(vtkzCFDReader, Controller, vtkMultiProcessController);

//-----------------------------------------------------------------------------
vtkzCFDReader::vtkzCFDReader()
{

  this->ZoneDataArraySelection = vtkDataArraySelection::New();
  this->CellDataArraySelection = vtkDataArraySelection::New();

  this->FileName = NULL;
  this->CaseFileName = NULL;
  this->FileNameOld = new vtkStdString;
  this->ProblemName = new vtkStdString;
  this->CaseName = new vtkStdString;

  this->Refresh = false;

  this->NumProcesses = 1;
  this->ProcessId = 0;

//#ifdef PARAVIEW_USE_MPI
  this->Controller = NULL;
  this->SetController(vtkMultiProcessController::GetGlobalController());

  if (this->Controller == NULL)
    {
    this->NumProcesses = 1;
    this->ProcessId = 0;
    }
  else
    {
    this->NumProcesses = this->Controller->GetNumberOfProcesses();
    this->ProcessId = this->Controller->GetLocalProcessId();
    }
//#endif
  //this->CaseType = RECONSTRUCTED_CASE;
  this->MTimeOld = 0;

  this->SetNumberOfInputPorts(0);
  this->SetNumberOfOutputPorts(1);

}

//-----------------------------------------------------------------------------
vtkzCFDReader::~vtkzCFDReader()
{
  delete this->FileNameOld;

  delete this->ProblemName;
  delete this->CaseName;

  this->ZoneDataArraySelection->Delete();

  this->SetController(NULL);
}

//-----------------------------------------------------------------------------
void vtkzCFDReader::PrintSelf(ostream &os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
  //os << indent << "Case Type: " << this->CaseType << endl;
  os << indent << "MTimeOld: " << this->MTimeOld << endl;
  os << indent << "Number of Processes: " << this->NumProcesses << endl;
  os << indent << "Process Id: " << this->ProcessId << endl;
  os << indent << "Controller: " << this->Controller << endl;
}

//-----------------------------------------------------------------------------
/*
void vtkzCFDReader::SetCaseType(const int t)
{
  if (this->CaseType != t)
    {
    this->CaseType = static_cast<caseType>(t);
    this->Refresh = true;
    this->Modified();
    }
}
*/
//-----------------------------------------------------------------------------

void vtkzCFDReader::ReadPython(std::map<int,std::string> &zoneToBc)
{

  using namespace boost::python;

  const char *env = vtksys::SystemTools::GetEnv("ZCFD_HOME");

  if(env)
  {
    PyObject *obj = PySys_GetObject("path");
    PyObject *pPath = PyUnicode_FromString(env);
    PyList_Append(obj,pPath);
    Py_DECREF(pPath);

    PySys_SetObject("path",obj);

    //PySys_SetPath(const_cast<char*>(env));
  }

  std::cout << "Reading Case File " <<  *CaseName << " " << std::string(env) << std::endl;

  //Initialize python
  Py_Initialize();

  PyRun_SimpleString("import sys\nprint sys.path");

  //Get the main module
  object main_module = import("__main__");
  object main_namespace = main_module.attr("__dict__");

  try
  {
    object perror = exec_file( (const char *)(*CaseName), main_namespace, main_namespace);
  }
  catch(const boost::python::error_already_set&)
  {
    PyObject *e, *v, *t;
    PyErr_Fetch(&e, &v, &t);

    //Get error message
    Py_ssize_t size;
    /*const char *pStrErrorMessage = PyUnicode_AsUTF8AndSize(v, &size);

    std::cout << "Python Exception: " << pStrErrorMessage << std::endl;*/

    // A NULL e means that there is not available Python
    // exception
    if (!e) return;
/*
    // See if the exception was an AttributeError. If so,
    // throw a C++ version of that exception
    if (PyErr_GivenExceptionMatches(perror.ptr(), e))
    {
      // We construct objects now since we plan to keep
      // ownership of the references.
      object e_obj(handle<>(allow_null(e));
      object v_obj(handle<>(allow_null(v));
      object t_obj(handle<>(allow_null(t));

      throw AttributeException(e_obj, v_obj, t_obj);
    }
*/
    // We didn't do anything with the Python exception,
    // and we never took ownership of the refs, so it's
    // safe to simply pass them back to PyErr_Restore
    PyErr_Restore(e, v, t);

    // Rethrow the exception (or whatever...this
    // is just an example.)
    throw;
  }

  extract<dict> cppdict_ext(main_namespace["parameters"]);

  assert(cppdict_ext.check());

  dict param = cppdict_ext();

  std::cout << "Read Case File" << std::endl;

  //object r = param.get("reference");
  //std::string s = extract<std::string>(r);
  //std::cout << s;

  bool found = true;
  int count = 1;
  while(found)
  {
    std::ostringstream key;
    key << "BC_" << count;
    if(param.has_key(key.str()))
    {
      object ob = param.get(key.str());
      dict bc = extract<dict>(ob);

      ob = bc.get("type");
      std::string bctype = extract<std::string>(ob);

      if(bc.has_key("zone"))
      {
        object ob = bc.get("zone");
        list l = extract<list>(ob);
        for(int i = 0; i < len(l); ++i)
        {
          object ob = l[i];
          int z = extract<int>(ob);

          zoneToBc[z] = bctype;
        }
      }
    }
    else
      found = false;
    count++;
  }

}

bool vtkzCFDReader::ReadStatusFile(const char *fileName)
{
  std::string status;

  std::string settingsFileName(fileName);
  std::ifstream settingsFile(settingsFileName.c_str(), ios::in | ios::binary | ios::ate);
  vtkDebugMacro("Attempting to load settings file '" << fileName << "'");
  if (settingsFile.is_open())
    {
    std::streampos size = settingsFile.tellg();
    settingsFile.seekg(0, ios::beg);
    int stringSize = size;
    char * settingsString = new char[stringSize+1];
    settingsFile.read(settingsString, stringSize);
    settingsString[stringSize] = '\0';
    settingsFile.close();

    status = settingsString;

    delete[] settingsString;
    }
  std::cout << status << std::endl;

  // Parse the user settings
  Json::Reader reader;
  Json::Value val;
  bool success = reader.parse(status, val, false);

  if(val.isMember("num processor"))
  {
    std::cout << val.get("num processor",Json::Value(0)) << std::endl;
  }
  bool problemNameSet=false;
  if(val.isMember("problem"))
  {

    std::cout << val.get("problem",Json::Value("UNKNOWN")) << std::endl;

    *this->ProblemName = vtksys::SystemTools::GetFilenamePath(fileName) + "/";

    *this->ProblemName += vtkStdString(val.get("problem",Json::Value("UNKNOWN")).asString());


    std::cout << "Problem Name: " << *this->ProblemName << std::endl;

    if(!vtksys::SystemTools::FileExists(*this->ProblemName))
    {
      *this->ProblemName += ".h5";
      problemNameSet = vtksys::SystemTools::FileExists(*this->ProblemName);
    }
    else
      problemNameSet=true;
  }

  if(val.isMember("case"))
  {
    *this->CaseName = vtksys::SystemTools::GetFilenamePath(fileName) + "/";

    *this->CaseName += vtkStdString(val.get("case",Json::Value("UNKNOWN")).asString());


    std::cout << "Case Name: " << *this->CaseName << std::endl;

    if(!vtksys::SystemTools::FileExists(*this->CaseName))
    {
      *this->CaseName += ".py";
      if(!vtksys::SystemTools::FileExists(*this->CaseName))
      {
        *this->CaseName = "";
      }
    }

  }

/*
  for(Json::Value::iterator itr = val.begin(); itr != val.end(); ++itr)
  {
    std::cout << *itr << std::endl;
  }
*/
  return problemNameSet;
}

int vtkzCFDReader::RequestInformation(vtkInformation *vtkNotUsed(request),
    vtkInformationVector **vtkNotUsed(inputVector), vtkInformationVector *outputVector)
{

  const char *env = vtksys::SystemTools::GetEnv("ZCFD_HOME");
  if(env == 0)
  {
    vtkErrorMacro("ZCFD_HOME env has to be specified!");
    return 0;
  }
  std::string zcfdhome(env);
  std::cout << "ZCFD_HOME: " << zcfdhome << std::endl;

  if (!this->FileName || strlen(this->FileName) == 0)
  {
    vtkErrorMacro("FileName has to be specified!");
    return 0;
  }

  if (*FileNameOld != this->FileName || this->Refresh)
  {
    // retain selection status when just refreshing a case
    //if (*this->FileNameOld != "" && *this->FileNameOld != this->FileName)
    {
      // clear selections
      this->ZoneDataArraySelection->RemoveAllArrays();
    }

    std::cout << "RequestInformation: Loading Array selection" << std::endl;

    *this->FileNameOld = vtkStdString(this->FileName);

    // Check if h5 or status file
    std::string extension = vtksys::SystemTools::GetFilenameLastExtension(this->FileName);

    std::cout << "Extension: " << extension << std::endl;

    if(extension == ".txt") // _status.txt
    {
      // Read status file
      bool problemNameSet = ReadStatusFile(this->FileName);
      if(!problemNameSet)
        return 0;
    }
    else if(extension == ".h5") // Mesh file
    {
      *this->ProblemName = vtkStdString(this->FileName);
    }

    vtksys::SystemTools::FileExists(this->FileName);
    vtksys::SystemTools::GetFilenameLastExtension(this->FileName);
    vtksys::SystemTools::GetFilenamePath(FileName);
    vtksys::SystemTools::GetFilenameName(FileName);
    vtksys::SystemTools::GetFilenameWithoutExtension(FileName);

    if(CaseFileName != NULL)
    {
      *this->CaseName = vtkStdString(this->CaseFileName);
      std::cout << "RequestInformation: Case Name " << std::endl;
    }

    //this->InternalMeshSelectionStatus
    //    = this->Parent->GetZonehArrayStatus(internalMeshName.c_str());
    std::map<int,std::string> zoneToBc;
    if(*CaseName != "")
    {
      ReadPython(zoneToBc);
    }

    // Read zone indexes from h5 file
    hdf::HDFFile<> file(*this->ProblemName, hdf::HDFFile<>::readonly);
    boost::shared_ptr<hdf::HDFGroup<> > meshg = file.openGroup("mesh",true);
    int totalNumFaces, totalNumCells;

    meshg->openAttribute("numFaces")->readData(totalNumFaces);
    meshg->openAttribute("numCells")->readData(totalNumCells);

    int numCells = totalNumCells;
    int numAllFaces = totalNumFaces;


    std::vector<int> zoneInfo;
    std::vector<int> faceBC;

    faceBC.resize(numAllFaces);
    zoneInfo.resize(numAllFaces*2);

    meshg->openDataset("faceBC")->readData(faceBC);
    meshg->openDataset("faceInfo")->readData(zoneInfo);

    std::map<std::string,std::set<int> > bc;
    std::set<int> zone;
    std::map<int,std::vector<int> > bcToZone;

    assert(faceBC.size() == zoneInfo.size()/2);

    // Count boundary faces
    int numBoundaryFaces = 0;
    for (int i = 0; i < faceBC.size(); i++)
    {
      int z = zoneInfo[i*2];
      if (faceBC[i] != detail::BoundaryCondition::NONE)
      {
        numBoundaryFaces++;
        zone.insert(z);
      }

      if(zoneToBc.find(z) != zoneToBc.end())
      {
        bc[zoneToBc[z]].insert(z);
      }
      else
      {
        std::string str = detail::BoundaryCondition::toString(faceBC[i]);
        if(str != "interior")
          bc[str].insert(z);
      }
    }

    std::set<int> cellZoneSet;
    {
      std::vector<int> cellZones(numCells,0);
      try{
        // Read cells zones

        meshg->openDataset("cellZone")->readData(cellZones);

      }
      catch(...){
      }
      for(int i = 0; i < cellZones.size(); ++i){
        cellZoneSet.insert(cellZones[i]);
      }
    }
    for(std::set<int>::iterator itr = cellZoneSet.begin(); itr != cellZoneSet.end(); ++itr){
      char buf[1024];
      std::sprintf(buf, "%d", *itr);

      bc["interior-"+std::string(buf)].insert(*itr);
    }

    int count = 0;

    selectionToZone.clear();

    std::cout << "RequestInformation: Populating Selection " << GetNumberOfSelectionArrays(ZoneDataArraySelection) << std::endl;

    // For each boundary condition
    for(std::map<std::string,std::set<int> >::iterator it = bc.begin(); it != bc.end(); ++it)
    {
      this->ZoneDataArraySelection->AddArray(it->first.c_str());
      selectionToZone[count] = it->second;
      std::cout << it->first << " " << *(it->second.begin()) << std::endl;
      count++;
    }

    // For each zone create a set
    int numZones = zone.size();
    for(std::set<int>::iterator it = zone.begin(); it != zone.end(); ++it)
    {
      char buf[1024];
      std::sprintf(buf, "%d", *it);
      //this->ZoneDataArraySelection->AddArray(("Zone "+std::to_string(*it)).c_str());
      this->ZoneDataArraySelection->AddArray(("Zone "+std::string(buf)).c_str());
      selectionToZone[count].insert(*it);
      count++;
    }

    double timeRange[2];
    timeRange[0] = timeRange[1] = 0.0;
    outputVector->GetInformationObject(0)->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(), timeRange, 0);
    outputVector->GetInformationObject(0)->Set(vtkStreamingDemandDrivenPipeline::TIME_RANGE(), timeRange, 2);

    this->Refresh = false;

  }

  /*

  if (this->CaseType == RECONSTRUCTED_CASE)
    {
    int ret = 1;
    if (this->ProcessId == 0)
      {
      ret = this->Superclass::RequestInformation(request, inputVector,
          outputVector);
      }
    if (this->NumProcesses > 1)
      {
      // if there was an error in process 0 abort all processes
      this->BroadcastStatus(ret);
      if (ret == 0)
        {
        vtkErrorMacro(<< "The master process returned an error.");
        return 0;
        }

      vtkDoubleArray *timeValues;
      if (this->ProcessId == 0)
        {
        timeValues = this->Superclass::GetTimeValues();
        }
      else
        {
        timeValues = vtkDoubleArray::New();
        }
      this->Controller->Broadcast(timeValues, 0);
      if (this->ProcessId != 0)
        {
        this->Superclass::SetTimeInformation(outputVector, timeValues);
        timeValues->Delete();
        this->Superclass::Refresh = false;
        }
      this->GatherMetaData(); // pvserver deadlocks without this
      }

    return ret;
    }

  if (!this->Superclass::FileName || strlen(this->Superclass::FileName) == 0)
    {
    vtkErrorMacro("FileName has to be specified!");
    return 0;
    }

  if (*this->Superclass::FileNameOld != this->Superclass::FileName
      || this->Superclass::ListTimeStepsByControlDict
          != this->Superclass::ListTimeStepsByControlDictOld
      || this->Superclass::Refresh)
    {
    // retain selection status when just refreshing a case
    if (*this->Superclass::FileNameOld != "" && *this->Superclass::FileNameOld != this->Superclass::FileName)
      {
      // clear selections
      this->Superclass::CellDataArraySelection->RemoveAllArrays();
      this->Superclass::PointDataArraySelection->RemoveAllArrays();
      this->Superclass::LagrangianDataArraySelection->RemoveAllArrays();
      this->Superclass::PatchDataArraySelection->RemoveAllArrays();
      }

    *this->Superclass::FileNameOld = vtkStdString(this->FileName);
    this->Superclass::Readers->RemoveAllItems();
    this->Superclass::NumberOfReaders = 0;

    vtkStringArray *procNames = vtkStringArray::New();
    vtkDoubleArray *timeValues;

    // recreate case information
    vtkStdString masterCasePath, controlDictPath;
    this->Superclass::CreateCasePath(masterCasePath, controlDictPath);

    this->Superclass::CreateCharArrayFromString(this->Superclass::CasePath,
        "CasePath", masterCasePath);

    int ret = 1;
    if (this->ProcessId == 0)
      {
      // search and list processor subdirectories
      vtkDirectory *dir = vtkDirectory::New();
      if (!dir->Open(masterCasePath.c_str()))
        {
        vtkErrorMacro(<< "Can't open " << masterCasePath.c_str());
        dir->Delete();
        this->BroadcastStatus(ret = 0);
        return 0;
        }
      vtkIntArray *procNos = vtkIntArray::New();
      for (int fileI = 0; fileI < dir->GetNumberOfFiles(); fileI++)
        {
        const vtkStdString subDir(dir->GetFile(fileI));
        if (subDir.substr(0, 9) == "processor")
          {
          const vtkStdString procNoStr(subDir.substr(9, vtkStdString::npos));
          char *conversionEnd;
          const int procNo = strtol(procNoStr.c_str(), &conversionEnd, 10);
          if (procNoStr.c_str() + procNoStr.length() == conversionEnd && procNo
              >= 0)
            {
            procNos->InsertNextValue(procNo);
            procNames->InsertNextValue(subDir);
            }
          }
        }
      procNos->Squeeze();
      procNames->Squeeze();
      dir->Delete();

      // sort processor subdirectories by processor numbers
      vtkSortDataArray::Sort(procNos, procNames);
      procNos->Delete();

      // get time directories from the first processor subdirectory
      if (procNames->GetNumberOfTuples() > 0)
        {
        vtkOpenFOAMReader *masterReader = vtkOpenFOAMReader::New();
        masterReader->SetFileName(this->FileName);
        masterReader->SetParent(this);
        if (!masterReader->MakeInformationVector(outputVector, procNames
        ->GetValue(0)) || !masterReader->MakeMetaDataAtTimeStep(true))
          {
          procNames->Delete();
          masterReader->Delete();
          this->BroadcastStatus(ret = 0);
          return 0;
          }
        this->Superclass::Readers->AddItem(masterReader);
        timeValues = masterReader->GetTimeValues();
        masterReader->Delete();
        }
      else
        {
        timeValues = vtkDoubleArray::New();
        this->Superclass::SetTimeInformation(outputVector, timeValues);
        }
      }
    else
      {
      timeValues = vtkDoubleArray::New();
      }

    if (this->NumProcesses > 1)
      {
      // if there was an error in process 0 abort all processes
      this->BroadcastStatus(ret);
      if (ret == 0)
        {
        vtkErrorMacro(<< "The master process returned an error.");
        timeValues->Delete(); // don't have to care about process 0
        return 0;
        }

      this->Broadcast(procNames);
      this->Controller->Broadcast(timeValues, 0);
      if (this->ProcessId != 0)
        {
        this->Superclass::SetTimeInformation(outputVector, timeValues);
        timeValues->Delete();
        }
      }

    if (this->ProcessId == 0 && procNames->GetNumberOfTuples() == 0)
      {
      timeValues->Delete();
      }

    // create reader instances for other processor subdirectories
    // skip processor0 since it's already created
    for (int procI = (this->ProcessId ? this->ProcessId : this->NumProcesses); procI
        < procNames->GetNumberOfTuples(); procI += this->NumProcesses)
      {
      vtkOpenFOAMReader *subReader = vtkOpenFOAMReader::New();
      subReader->SetFileName(this->FileName);
      subReader->SetParent(this);
      // if getting metadata failed simply delete the reader instance
      if (subReader->MakeInformationVector(NULL, procNames->GetValue(procI))
          && subReader->MakeMetaDataAtTimeStep(true))
        {
        this->Superclass::Readers->AddItem(subReader);
        }
      else
        {
        vtkWarningMacro(<<"Removing reader for processor subdirectory "
            << procNames->GetValue(procI).c_str());
        }
      subReader->Delete();
      }

    procNames->Delete();

    this->GatherMetaData();
    this->Superclass::Refresh = false;
    }

  outputVector->GetInformationObject(0)->Set(
    CAN_HANDLE_PIECE_REQUEST(),
    1);
*/
  return 1;
}

//-----------------------------------------------------------------------------
int vtkzCFDReader::RequestData(vtkInformation *vtkNotUsed(request),
    vtkInformationVector **vtkNotUsed(inputVector), vtkInformationVector *outputVector)
{

  vtkInformation* outInfo = outputVector->GetInformationObject(0);
  vtkMultiBlockDataSet
      *output =
          vtkMultiBlockDataSet::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  int nSteps = 0;
  double requestedTimeValue(0);
  if (outInfo->Has(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP()))
    {
    requestedTimeValue
        = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());
    nSteps = outInfo->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
    }

  // Read request zones
  hdf::HDFFile<> file(*this->ProblemName, hdf::HDFFile<>::readonly);
  boost::shared_ptr<hdf::HDFGroup<> > meshg = file.openGroup("mesh",false);
  int totalNumFaces, totalNumCells;

  meshg->openAttribute("numFaces")->readData(totalNumFaces);
  meshg->openAttribute("numCells")->readData(totalNumCells);

  int numCells = totalNumCells;
  int numAllFaces = totalNumFaces;


  std::vector<int> zoneInfo;
  //std::vector<int> faceBC;

  //faceBC.resize(numAllFaces);
  zoneInfo.resize(numAllFaces*2);

  //meshg->openDataset("faceBC")->readData(faceBC);
  meshg->openDataset("faceInfo")->readData(zoneInfo);

  std::vector<int> faceTypeGobalPtr;
  faceTypeGobalPtr.reserve(numAllFaces+1);

  meshg->openDataset("faceType")->readData(faceTypeGobalPtr);
  faceTypeGobalPtr.push_back(0);
  for(int i=faceTypeGobalPtr.size()-1;i>0;--i)
    faceTypeGobalPtr[i] = faceTypeGobalPtr[i-1];
  faceTypeGobalPtr[0] = 0;
  for(int i=1;i<faceTypeGobalPtr.size();++i)
    faceTypeGobalPtr[i] += faceTypeGobalPtr[i-1];

  std::vector<int> faceNodes;
  faceNodes.reserve(faceTypeGobalPtr[numAllFaces]);
  meshg->openDataset("faceNodes")->readData(faceNodes);

  std::vector<double> vertex;
  meshg->openDataset("nodeVertex")->readData(vertex);

  std::vector<int> faceCell(numAllFaces*2);
  meshg->openDataset("faceCell")->readData(faceCell);

  // Read cell zones
  std::vector<int> cellZones(numCells,0);
  try{
    // Read cells zones

    meshg->openDataset("cellZone")->readData(cellZones);

  }
  catch(...){
  }

  std::vector<int> cellType(totalNumCells);
  // Check for cellFace
  try
  {
    meshg->openDataset("cellType")->readData(cellType);
  }
  catch(...)
  {
    std::cout << "zCFDReader: cellType missing. Treating all cells as polyhedral " << std::endl;
    //std::vector<hsize_t> dimsf(1);
    //dimsf[0] = totalNumCells;
    //hdf::Slab<1> filespace(dimsf);
    //boost::shared_ptr<hdf::HDFDataSet<> > dataset
    //  = meshg->createDataset<int>("cellType", filespace);

    //Read cell face
    std::vector<int> cellFaceCount(totalNumCells,0);
    for(int i=0;i<numAllFaces;++i)
    {
      int left = faceCell[i*2];
      int right = faceCell[i*2 + 1];
      if(left < totalNumCells)
        cellFaceCount[left]++;
      if(right < totalNumCells)
        cellFaceCount[right]++;
    }
    for(int i=0;i<totalNumCells;++i)
    {
      cellType[i] = detail::CellType::encode(detail::CellType::POLY,cellFaceCount[i]);
    }
    //dataset->writeData(cellType);
  }
  std::vector<int> cellFacePtr(totalNumCells+1,0);
  for(int i=0;i<totalNumCells;++i)
  {
    cellFacePtr[i+1] = cellFacePtr[i] + detail::CellType::getNumFaces(cellType[i]);
  }
  std::vector<int> cellFace(cellFacePtr[totalNumCells]);
  try
  {
    meshg->openDataset("cellFace")->readData(cellFace);
  }
  catch(...)
  {
    // Need to create cellFaces
    std::cout << "zCFDReader: Creating cell face information " << std::endl;
    //std::vector<hsize_t> dimsf(1);
    //dimsf[0] = cellFace.size();
    //hdf::Slab<1> filespace(dimsf);
    //boost::shared_ptr<hdf::HDFDataSet<> > dataset
    //  = meshg->createDataset<int>("cellFace", filespace);

    std::vector<int> count(totalNumCells,0);
    for(int i = 0; i < numAllFaces; ++i)
    {
      for(int j = 0; j < 2; ++j)
      {
        int c = faceCell[i*2 + j];
        if(c < totalNumCells)
        {
          cellFace[cellFacePtr[c]+count[c]] = i;
          count[c]++;
        }
      }
    }
    //dataset->writeData(cellFace);
  }

  for(int i=0;i<GetNumberOfSelectionArrays(ZoneDataArraySelection);++i)
  {
    std::string name = GetSelectionArrayName(ZoneDataArraySelection,i);

    if(name.find("interior") != 0 && GetSelectionArrayStatus(ZoneDataArraySelection,name.c_str()))
    {
      // User has selected zone for extraction

      // Get list of faces
      std::set<int> &zoneSet = selectionToZone[i];
      std::vector<int> faceList;
      for(int j=0;j<zoneInfo.size()/2;++j)
      {
        int z = zoneInfo[j*2];
        if(zoneSet.find(z) != zoneSet.end())
        {
          faceList.push_back(j);
          //std::cout << "Zones: " << z << std::endl;
        }
      }
      int numFaces = faceList.size();

      std::cout << "Faces: " << numFaces << std::endl;

      vtkSmartPointer<vtkPolyData> polyData =
          vtkSmartPointer<vtkPolyData>::New();
      polyData->Allocate(numFaces); // set number of faces

      vtkSmartPointer<vtkFloatArray> var =
          vtkSmartPointer<vtkFloatArray>::New();
      var->SetName("zone");
      var->SetNumberOfComponents(1);
      var->SetNumberOfValues(numFaces);

      vtkIdType cellList[100];
      int count = 0;
      int faceCount = 0;
      std::vector<int> nodeMap(vertex.size()/3, -1);
      std::vector<int> reverseMap;

      // Extract faces
      for(int f=0;f<faceList.size();++f)
      {
        int face=faceList[f];
        int z = zoneInfo[face*2];
        int start = faceTypeGobalPtr[face];
        int end   = faceTypeGobalPtr[face+1];
        int numPts = end - start;
        assert(numPts < 100);

        for (int j = start; j < end; ++j)
        {
          int n = faceNodes[j];
          int mapIndex = -1;
          if (nodeMap[n] == -1)
          {
            mapIndex = count;
            nodeMap[n] = count;
            reverseMap.push_back(n);
            count++;
          }
          else
          {
            mapIndex = nodeMap[n];
          }
          cellList[j - start] = mapIndex;
        }

        int vtkCellType = VTK_EMPTY_CELL;
        if (numPts == 3)
          vtkCellType = VTK_TRIANGLE;
        else if (numPts == 4)
          vtkCellType = VTK_QUAD;
        else
          vtkCellType = VTK_POLYGON;

        polyData->InsertNextCell(vtkCellType, numPts, cellList);

        float va = z;
        var->InsertTuple(faceCount,&va);
        faceCount++;
      }

      polyData->GetCellData()->SetActiveScalars("zone");
      polyData->GetCellData()->SetScalars(var);

      vtkSmartPointer<vtkPoints> newPts = vtkSmartPointer<vtkPoints>::New();
      newPts->SetNumberOfPoints(count);

      float xyz[3];
      for (int v = 0; v < 3; ++v)
        xyz[v] = 0.0;

      for (int i = 0; i < reverseMap.size(); ++i)
      {
        int n = reverseMap[i];

        for (int v = 0; v < 3; ++v)
          xyz[v] = vertex[3*n+v];

        newPts->InsertPoint(i, xyz[0], xyz[1], xyz[2]);
      }
      polyData->SetPoints(newPts);

      // Extract nodes
      //vtkUnstructuredGrid *ugrid = vtkUnstructuredGrid::New();

      //vtkMultiBlockDataSet *subOutput = vtkMultiBlockDataSet::New();

      const int blockI = output->GetNumberOfBlocks();
      output->SetBlock(blockI, polyData);
      output->GetMetaData(blockI)->Set(vtkCompositeDataSet::NAME(), name.c_str());
      //polyData->Delete();

    }
    else if(name.find("interior") == 0 && GetSelectionArrayStatus(ZoneDataArraySelection,name.c_str()))
    {

      // Get list of cells
      std::set<int> &zoneSet = selectionToZone[i];

      std::vector<int> cellList;
      cellList.reserve(totalNumCells);
      for(std::set<int>::iterator it = zoneSet.begin(); it != zoneSet.end(); ++it){
        for(int i = 0; i < totalNumCells; ++i){
          if(cellZones[i] == *it){
            cellList.push_back(i);
          }
        }
      }

      // Add interior
      vtkSmartPointer<vtkUnstructuredGrid> uGrid = vtkSmartPointer<
          vtkUnstructuredGrid>::New();
      uGrid->Allocate(cellList.size());

      for(int j = 0; j < cellList.size(); ++j)
      {
        int i = cellList[j];

        std::set<vtkIdType> cellNodeList;
        std::vector<vtkIdType> faceStream;

        int start = cellFacePtr[i];
        int end  = cellFacePtr[i+1];

        for(int f = start; f < end; ++f)
        {
          int face = cellFace[f];
          int z = zoneInfo[face*2];

          int start = faceTypeGobalPtr[face];
          int end   = faceTypeGobalPtr[face+1];
          int numPts = end - start;
          assert(numPts < 100);
          faceStream.push_back(numPts);

          for (int j = start; j < end; ++j)
          {
            int n = faceNodes[j];
            cellNodeList.insert(n);
            faceStream.push_back(n);
          }
        }
        std::vector<vtkIdType> cNodeList(cellNodeList.begin(),
            cellNodeList.end());
        vtkIdType nFaces = end - start;
        vtkIdType nNodes = cNodeList.size();

        uGrid->InsertNextCell(VTK_POLYHEDRON, nNodes, &cNodeList[0], nFaces,
            &faceStream[0]);
      }

      vtkSmartPointer<vtkPoints> newPts = vtkSmartPointer<vtkPoints>::New();
      newPts->SetNumberOfPoints(vertex.size());

      float xyz[3];
      for (int v = 0; v < 3; ++v)
        xyz[v] = 0.0;

      for (int i = 0; i < vertex.size()/3; ++i)
      {
        for (int v = 0; v < 3; ++v)
          xyz[v] = vertex[i*3+v];

        newPts->InsertPoint(i, xyz[0], xyz[1], xyz[2]);
      }
      uGrid->SetPoints(newPts);

      const int blockI = output->GetNumberOfBlocks();
      output->SetBlock(blockI, uGrid);
      output->GetMetaData(blockI)->Set(vtkCompositeDataSet::NAME(), name.c_str());
    }
  }


  /*
  if (this->CaseType == RECONSTRUCTED_CASE)
    {
    int ret = 1;
    if (this->ProcessId == 0)
      {
      ret = this->Superclass::RequestData(request, inputVector, outputVector);
      }
    this->BroadcastStatus(ret);
    this->GatherMetaData();
    return ret;
    }

  vtkInformation* outInfo = outputVector->GetInformationObject(0);
  vtkMultiBlockDataSet
      *output =
          vtkMultiBlockDataSet::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  int ret = 1;
  if (this->Superclass::Readers->GetNumberOfItems() > 0)
    {
    int nSteps = 0;
    double requestedTimeValue(0);
    if (outInfo->Has(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP()))
      {
      requestedTimeValue
          = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());
      nSteps = outInfo->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
      if (nSteps > 0)
        {
        outInfo->Set(vtkDataObject::DATA_TIME_STEP(), requestedTimeValue);
        }
      }

    vtkAppendCompositeDataLeaves *append = vtkAppendCompositeDataLeaves::New();
    // append->AppendFieldDataOn();

    vtkOpenFOAMReader *reader;
    this->Superclass::CurrentReaderIndex = 0;
    this->Superclass::Readers->InitTraversal();
    while ((reader
        = vtkOpenFOAMReader::SafeDownCast(this->Superclass::Readers->GetNextItemAsObject()))
        != NULL)
      {
      // even if the child readers themselves are not modified, mark
      // them as modified if "this" has been modified, since they
      // refer to the property of "this"
      if ((nSteps > 0 && reader->SetTimeValue(requestedTimeValue))
          || this->MTimeOld != this->GetMTime())
        {
        reader->Modified();
        }
      if (reader->MakeMetaDataAtTimeStep(false))
        {
        append->AddInputConnection(reader->GetOutputPort());
        }
      }

    this->GatherMetaData();

    if (append->GetNumberOfInputConnections(0) == 0)
      {
      output->Initialize();
      ret = 0;
      }
    else
      {
      // reader->RequestInformation() and RequestData() are called
      // for all reader instances without setting UPDATE_TIME_STEPS
      append->Update();
      output->ShallowCopy(append->GetOutput());
      }
    append->Delete();

    // known issue: output for process without sub-reader will not have CasePath
    output->GetFieldData()->AddArray(this->Superclass::CasePath);
    }
  else
    {
    this->GatherMetaData();
    // page 322 of The ParaView Guide says the output must be initialized
    output->Initialize();
    }

  this->Superclass::UpdateStatus();
  this->MTimeOld = this->GetMTime();
*/
  return 1;
}

int vtkzCFDReader::GetNumberOfSelectionArrays(vtkDataArraySelection *s)
{
  return s->GetNumberOfArrays();
}
int vtkzCFDReader::GetSelectionArrayStatus(vtkDataArraySelection *s,
    const char *name)
{
  return s->ArrayIsEnabled(name);
}
void vtkzCFDReader::SetSelectionArrayStatus(vtkDataArraySelection *s,
    const char* name, int status)
{
  unsigned long int mTime = s->GetMTime();
  if (status)
    {
    s->EnableArray(name);
    }
  else
    {
    s->DisableArray(name);
    }
  if (mTime != s->GetMTime()) // indicate that the pipeline needs to be updated
    {
    this->Modified();
    }
}
const char *vtkzCFDReader::GetSelectionArrayName(vtkDataArraySelection *s,
    int index)
{
  return s->GetArrayName(index);
}
//-----------------------------------------------------------------------------
void vtkzCFDReader::BroadcastStatus(int &status)
{
  if (this->NumProcesses > 1)
    {
    this->Controller->Broadcast(&status, 1, 0);
    }
}

//-----------------------------------------------------------------------------
void vtkzCFDReader::GatherMetaData()
{
  if (this->NumProcesses > 1)
    {
    /*
    this->AllGather(this->Superclass::PatchDataArraySelection);
    this->AllGather(this->Superclass::CellDataArraySelection);
    this->AllGather(this->Superclass::PointDataArraySelection);
    this->AllGather(this->Superclass::LagrangianDataArraySelection);
    // omit removing duplicated entries of LagrangianPaths as well
    // when the number of processes is 1 assuming there's no duplicate
    // entry within a process
    this->AllGather(this->Superclass::LagrangianPaths);
    */
    }
}

//-----------------------------------------------------------------------------
// Broadcast a vtkStringArray in process 0 to all processes
void vtkzCFDReader::Broadcast(vtkStringArray *sa)
{
  vtkIdType lengths[2];
  if (this->ProcessId == 0)
    {
    lengths[0] = sa->GetNumberOfTuples();
    lengths[1] = 0;
    for (int strI = 0; strI < sa->GetNumberOfTuples(); strI++)
      {
      lengths[1] += static_cast<vtkIdType>(sa->GetValue(strI).length()) + 1;
      }
    }
  this->Controller->Broadcast(lengths, 2, 0);
  char *contents = new char [lengths[1]];
  if (this->ProcessId == 0)
    {
    for (int strI = 0, idx = 0; strI < sa->GetNumberOfTuples(); strI++)
      {
      const int len = static_cast<int>(sa->GetValue(strI).length()) + 1;
      memmove(contents + idx, sa->GetValue(strI).c_str(), len);
      idx += len;
      }
    }
  this->Controller->Broadcast(contents, lengths[1], 0);
  if (this->ProcessId != 0)
    {
    sa->Initialize();
    sa->SetNumberOfTuples(lengths[0]);
    for (int strI = 0, idx = 0; strI < lengths[0]; strI++)
      {
      sa->SetValue(strI, contents + idx);
      idx += static_cast<int>(sa->GetValue(strI).length()) + 1;
      }
    }
  delete [] contents;
}

//-----------------------------------------------------------------------------
// AllGather vtkStringArray from and to all processes
void vtkzCFDReader::AllGather(vtkStringArray *s)
{
  vtkIdType length = 0;
  for (int strI = 0; strI < s->GetNumberOfTuples(); strI++)
    {
    length += static_cast<vtkIdType>(s->GetValue(strI).length()) + 1;
    }
  vtkIdType *lengths = new vtkIdType [this->NumProcesses];
  this->Controller->AllGather(&length, lengths, 1);
  vtkIdType totalLength = 0;
  vtkIdType *offsets = new vtkIdType [this->NumProcesses];
  for (int procI = 0; procI < this->NumProcesses; procI++)
    {
    offsets[procI] = totalLength;
    totalLength += lengths[procI];
    }
  char *allContents = new char [totalLength], *contents = new char [length];
  for (int strI = 0, idx = 0; strI < s->GetNumberOfTuples(); strI++)
    {
    const int len = static_cast<int>(s->GetValue(strI).length()) + 1;
    memmove(contents + idx, s->GetValue(strI).c_str(), len);
    idx += len;
    }
  this->Controller->AllGatherV(contents, allContents, length, lengths, offsets);
  delete [] contents;
  delete [] lengths;
  delete [] offsets;
  s->Initialize();
  for (int idx = 0; idx < totalLength; idx += static_cast<int>(strlen(allContents + idx)) + 1)
    {
    const char *str = allContents + idx;
    // insert only when the same string is not found
    if (s->LookupValue(str) == -1)
      {
      s->InsertNextValue(str);
      }
    }
  s->Squeeze();
  delete [] allContents;
}

//-----------------------------------------------------------------------------
// AllGather vtkDataArraySelections from and to all processes
void vtkzCFDReader::AllGather(vtkDataArraySelection *s)
{
  vtkIdType length = 0;
  for (int strI = 0; strI < s->GetNumberOfArrays(); strI++)
    {
    length += static_cast<vtkIdType>(strlen(s->GetArrayName(strI))) + 2;
    }
  vtkIdType *lengths = new vtkIdType [this->NumProcesses];
  this->Controller->AllGather(&length, lengths, 1);
  vtkIdType totalLength = 0;
  vtkIdType *offsets = new vtkIdType [this->NumProcesses];
  for (int procI = 0; procI < this->NumProcesses; procI++)
    {
    offsets[procI] = totalLength;
    totalLength += lengths[procI];
    }
  char *allContents = new char [totalLength], *contents = new char [length];
  for (int strI = 0, idx = 0; strI < s->GetNumberOfArrays(); strI++)
    {
    const char *arrayName = s->GetArrayName(strI);
    contents[idx] = s->ArrayIsEnabled(arrayName);
    const int len = static_cast<int>(strlen(arrayName)) + 1;
    memmove(contents + idx + 1, arrayName, len);
    idx += len + 1;
    }
  this->Controller->AllGatherV(contents, allContents, length, lengths, offsets);
  delete [] contents;
  delete [] lengths;
  delete [] offsets;
  // do not RemoveAllArray so that the previous arrays are preserved
  // s->RemoveAllArrays();
  for (int idx = 0; idx < totalLength; idx += static_cast<int>(strlen(allContents + idx + 1)) + 2)
    {
    const char *arrayName = allContents + idx + 1;
    s->AddArray(arrayName);
    if (allContents[idx] == 0)
      {
      s->DisableArray(arrayName);
      }
    else
      {
      s->EnableArray(arrayName);
      }
    }
  delete [] allContents;
}
