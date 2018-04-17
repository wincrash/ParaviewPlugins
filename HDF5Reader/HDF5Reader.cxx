#include "HDF5Reader.h"
 
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkInformationVector.h"
#include "vtkInformation.h"
#include "vtkDataObject.h"
#include "vtkSmartPointer.h"
#include "vtkVertexGlyphFilter.h"
 
#include <algorithm>
#include "hdf_wrapper.hpp"
#include "vtkUnstructuredGrid.h"
#include <iostream>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkFieldData.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkTriangle.h>
#include <vtkUnstructuredGrid.h>
#include <vtkVertex.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkPolyData.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkVariantArray.h>
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
vtkIdType HDF5Reader::findClosestSolutionIndex(double Time)
{
    double dt=1e24;
    vtkIdType ti=0;//time index
    for(vtkIdType i=0;i<times->GetNumberOfTuples();i++)
    {
        if(dt > std::abs(Time - times->GetValue(i) ) )
        {
            ti=i;
            dt=std::abs(Time-times->GetValue(i));
        }
    }

  return ti;
}



std::vector<double> GetData(h5cpp::Group &gr, std::string NAME,int kiekis)
{
    std::vector<double> data;
    if (gr.exists(NAME)) {
        auto dt = gr.open_dataset(NAME);
        hsize_t dim[2];
        dt.get_dataspace().get_dims(dim);
        data.resize(dim[0] * kiekis);
        dt.read(&data[0]);
    }
    return data;
}
void AddToVTKScalar(auto *location,std::vector<double>&array,std::string name)
{
    vtkDoubleArray*data = vtkDoubleArray::New();
    data->SetNumberOfTuples(array.size());
    data->SetNumberOfComponents(1);
    data->SetName(name.c_str());
        data->SetNumberOfValues(array.size());
    for (int ii = 0; ii < array.size(); ii++) {
        data->SetTuple1(ii, array[ii]);
    }
    location->AddArray(data);
    data->Delete();
}

void AddToVTKVector(auto *location,std::vector<double>&array,std::string name)
{
    vtkDoubleArray*data = vtkDoubleArray::New();
    data->SetNumberOfTuples(array.size()/4);
    data->SetNumberOfComponents(3);
    data->SetNumberOfValues(array.size()/4.0*3);
    data->SetName(name.c_str());
    for (int ii = 0; ii < array.size()/4; ii++) {
        data->SetTuple3(ii, array[ii*4],array[ii*4+1],array[ii*4+2]);
    }
    location->AddArray(data);
    data->Delete();
}








vtkStandardNewMacro(HDF5Reader);
 


HDF5Reader::HDF5Reader()
{
  this->FileName = NULL;
    this->DirectoryName = NULL;
  this->SetNumberOfInputPorts(0);
  this->SetNumberOfOutputPorts(1);
    this->times = vtkSmartPointer<vtkDoubleArray>::New();
    this->timesNames = vtkSmartPointer<vtkVariantArray>::New();
}
 
int HDF5Reader::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *outputVector)
{
 
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
    double requestedTime = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());
    std::cout<<"  "<<findClosestSolutionIndex(requestedTime)<<"\n";


 
  // get the ouptut
   vtkPolyData *output = vtkPolyData::SafeDownCast(
            outInfo->Get(vtkDataObject::DATA_OBJECT()));
 
  // Here is where you would read the data from the file. In this example,
  // we simply create a point.
 std::cout<<"iddd "<<findClosestSolutionIndex(requestedTime)<<"\n";
std::string tempFilename=filenames[findClosestSolutionIndex(requestedTime)];
std::cout<<"filename "<<tempFilename<<"\n";

  //vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
  h5cpp::File file(tempFilename, "a");

  auto gr = file.root();
    h5cpp::File file1(filenames[findClosestSolutionIndex(0)], "a");

  auto gr0 = file1.root();

  //auto gr = root.require_group(root.get_link_name(findClosestSolutionIndex(requestedTime)));
  //auto gr0=root.require_group(root.get_link_name(0));
  std::vector<double>oldPosition=GetData(gr0,"POSITIONS",4);
  file1.close();

  std::vector<double> STEP;
  STEP.push_back(gr.attrs().get<int>("STEP"));
  std::vector<double> TIME;
  TIME.push_back(gr.attrs().get<double>("TIME"));
  //vtkPolyData*output = vtkPolyData::New();
  vtkPoints*points = vtkPoints::New();
  std::vector<double> POSITIONS=GetData(gr,"POSITIONS",4);
  points->SetNumberOfPoints(POSITIONS.size()/4);
  std::vector<double> DEFORM;
  DEFORM.resize(POSITIONS.size());
  for(int i=0;i<POSITIONS.size()/4;i++)
  {
      points->SetPoint(i,POSITIONS[i*4],POSITIONS[i*4+1],POSITIONS[i*4+2]);
      DEFORM[i*4]=-oldPosition[i*4]+POSITIONS[i*4];
      DEFORM[i*4+1]=-oldPosition[i*4+1]+POSITIONS[i*4+1];
      DEFORM[i*4+2]=-oldPosition[i*4+2]+POSITIONS[i*4+2];
      DEFORM[i*4+3]=0;
  }

  output->SetPoints(points);
  vtkPointData* pdata = output->GetPointData();
  vtkCellData* cdata = output->GetCellData();
 vtkFieldData* fdata = output->GetFieldData();

     AddToVTKVector(pdata,DEFORM,"DISPLACEMENT");


  AddToVTKScalar(fdata,STEP,"STEP");
  AddToVTKScalar(fdata,TIME,"TIME");
  std::vector<double> data;
  data=GetData(gr,"UNIQUE_RADIUS",1);
  AddToVTKScalar(fdata,data,"UNIQUE_RADIUS");

  data=GetData(gr,"FORCE",4);
  AddToVTKVector(pdata,data,"FORCE");

  data=GetData(gr,"VELOCITY",4);
  AddToVTKVector(pdata,data,"VELOCITY");
 data=GetData(gr,"ACCELERATION",4);
  AddToVTKVector(output->GetPointData(),data,"ACCELERATION");
  data=GetData(gr,"PARTICLE_MATERIAL",1);
  AddToVTKScalar(output->GetPointData(),data,"PARTICLE_MATERIAL");
  data=GetData(gr,"PARTICLE_TYPE",1);
  AddToVTKScalar(output->GetPointData(),data,"PARTICLE_TYPE");
  data=GetData(gr,"RADIUS",1);
  AddToVTKScalar(output->GetPointData(),data,"RADIUS");
  data=GetData(gr,"PARTICLE_FIX",1);
  AddToVTKScalar(output->GetPointData(),data,"PARTICLE_FIX");
  data=GetData(gr,"PARTICLE_TYPE",1);
  AddToVTKScalar(output->GetPointData(),data,"PARTICLE_TYPE");

  if(gr.exists("BOND_PARTICLE_1") &&gr.exists("BOND_PARTICLE_2")){
      data=GetData(gr,"BOND_STATE",1);
      AddToVTKScalar(cdata,data,"BOND_STATE");

      data=GetData(gr,"BOND_F_LIMIT_N",1);
      AddToVTKScalar(cdata,data,"BOND_F_LIMIT_N");
      data=GetData(gr,"BOND_F_LIMIT_T",1);
      AddToVTKScalar(cdata,data,"BOND_F_LIMIT_T");
      data=GetData(gr,"BOND_F_N",4);
      AddToVTKVector(cdata,data,"BOND_F_N");
      std::vector<double> particle1=GetData(gr,"BOND_PARTICLE_1",1);
      std::vector<double> particle2=GetData(gr,"BOND_PARTICLE_2",1);
      std::vector<double> distance;
      std::vector<double> enlongation;
      distance.resize(particle1.size(),0);
      enlongation.resize(particle1.size(),0);
        vtkCellArray*lines=vtkCellArray::New();
      for(int i=0;i<particle1.size();i++)
      {
          lines->InsertNextCell(2);
          int id1=particle1[i];
          int id2=particle2[i];

          lines->InsertCellPoint(id1);
          lines->InsertCellPoint(id2);
          double p1[3];
          double p2[3];
          double oldp1[3];
          double oldp2[3];
          p1[0]=POSITIONS[id1*4];
          p1[1]=POSITIONS[id1*4+1];
          p1[2]=POSITIONS[id1*4+2];

          p2[0]=POSITIONS[id2*4];
          p2[1]=POSITIONS[id2*4+1];
          p2[2]=POSITIONS[id2*4+2];

          oldp1[0]=oldPosition[id1*4];
          oldp1[1]=oldPosition[id1*4+1];
          oldp1[2]=oldPosition[id1*4+2];

          oldp2[0]=oldPosition[id2*4];
          oldp2[1]=oldPosition[id2*4+1];
          oldp2[2]=oldPosition[id2*4+2];


          double d1=sqrt((p1[0]-p2[0])*(p1[0]-p2[0])+(p1[1]-p2[1])*(p1[1]-p2[1])+(p1[2]-p2[2])*(p1[2]-p2[2]));
          double oldd1=sqrt((oldp1[0]-oldp2[0])*(oldp1[0]-oldp2[0])+(oldp1[1]-oldp2[1])*(oldp1[1]-oldp2[1])+(oldp1[2]-oldp2[2])*(oldp1[2]-oldp2[2]));
          distance[i]=d1;
          enlongation[i]=oldd1-d1;
      }
      AddToVTKScalar(cdata,distance,"DISTANCE");
      AddToVTKScalar(cdata,enlongation,"ENLONGATION");
      output->SetLines(lines);
      lines->Delete();

  }else
  {
      vtkCellArray*verts=vtkCellArray::New();
      for(int i=0;i<POSITIONS.size()/4;i++)
      {
          verts->InsertNextCell(1);
          verts->InsertCellPoint(i);

      }
      output->SetVerts(verts);
      verts->Delete();

  }
  points->Delete();



file.close();
  return 1;
}
 
void HDF5Reader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
 
  os << indent << "File Name: "
      << (this->FileName ? this->FileName : "(none)") << "\n";
}

std::vector<std::string> split(const std::string& str, const std::string& delim)
{
    std::vector<std::string> tokens;
    size_t prev = 0, pos = 0;
    do
    {
        pos = str.find(delim, prev);
        if (pos == std::string::npos) pos = str.length();
        std::string token = str.substr(prev, pos-prev);
        if (!token.empty()) tokens.push_back(token);
        prev = pos + delim.length();
    }
    while (pos < str.length() && prev < str.length());
    return tokens;
}

int HDF5Reader::RequestInformation(vtkInformation *vtkNotUsed(request), vtkInformationVector **vtkNotUsed(inputVector), vtkInformationVector *outputVector)
{
    if (!this->FileName || strlen(this->FileName) == 0)
      {
      std::cout<<"FileName has to be specified!";
      //this->SetErrorCode(vtkErrorCode::NoFileNameError);
      return 0;
      }
    filenames.resize(0,"");
    std::cout<<fs::path(FileName).parent_path();

    for(auto& p: fs::directory_iterator(fs::path(FileName).parent_path()))
    {
        if (p.path().string().find(".h5") != std::string::npos) {

            filenames.push_back(p.path().string());
        }
    }
    std::sort(filenames.begin(),filenames.end());

    for (int i = 0; i < filenames.size(); i++) {
            h5cpp::File file(FileName, "a");
            auto root = file.root();
            double time=root.attrs().get<double>("TIME");
            std::cout<<time<<"\n";
            this->timesNames->InsertNextValue(root.attrs().get<double>("TIME"));
           // this->times->InsertNextValue(root.attrs().get<double>("TIME"));
            this->times->InsertNextValue(i);
            file.close();
    }


  /*  h5cpp::File file(FileName, "a");
    auto root = file.root();

    for (int i = 0; i < root.size(); i++) {
        std::cout << i<<"\t" << root.get_link_name(i) << std::endl;
        auto gr = root.require_group(root.get_link_name(i));
         this->timesNames->InsertNextValue(gr.attrs().get<double>("TIME"));
        this->times->InsertNextValue(gr.attrs().get<double>("TIME"));
    }

    file.close();
*/
    vtkIdType nsteps = times->GetNumberOfTuples();
    double * trange = new double[2];
    trange[0]=times->GetValue(0);
    trange[1]=times->GetValue(nsteps-1);
    vtkInformation * outInfo = outputVector->GetInformationObject(0);
    outInfo->Set(
                vtkStreamingDemandDrivenPipeline::TIME_STEPS(),
                times->GetPointer(0), nsteps
                );
    outInfo->Set(
                vtkStreamingDemandDrivenPipeline::TIME_RANGE(),
                trange,2);

    //2do create a list of fields (in current time only?)
    //2do create a list of blocks i.e. patches,zones regions etc


    return 1;
}
