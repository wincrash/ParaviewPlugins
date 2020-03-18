#include "HDF5ReaderConvex.h"

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
#include <vtkUnstructuredGrid.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkVariantArray.h>
#include <experimental/filesystem>
#include <vtkMultiBlockDataSet.h>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/lexical_cast.hpp>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridAlgorithm.h>
#include <algorithm>
#include <functional>
#include <iostream>
#include <iterator>
#include <numeric>
#include <vector>
namespace fs = std::experimental::filesystem;
vtkIdType HDF5ReaderConvex::findClosestSolutionIndex(double Time)
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
        int rank=dt.get_dataspace().get_rank();
        hsize_t dim[2];
        dt.get_dataspace().get_dims(dim);
        if(rank==1)
        data.resize(dim[0] );
        else
        data.resize(dim[0] *dim[1] );
        dt.read(&data[0]);
    }
    return data;
}

std::vector<int> GetDataINT(h5cpp::Group &gr, std::string NAME,int kiekis)
{
    std::vector<int> data;
    if (gr.exists(NAME)) {
        auto dt = gr.open_dataset(NAME);
        hsize_t dim[2];
        dt.get_dataspace().get_dims(dim);
        data.resize(dim[0] * kiekis);
        dt.read(&data[0]);
    }
    return data;
}


void AddToVTKScalar(vtkFieldData *location,std::vector<double>&array,std::string name)
{
    if(array.size()>0)
    {
        //vtkSmartPointer<vtkDoubleArray>data = vtkSmartPointer<vtkDoubleArray>::New();
        vtkDoubleArray*data=vtkDoubleArray::New();
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
}



void AddToVTKScalar(vtkFieldData *location,std::vector<int>&array,std::string name)
{
    if(array.size()>0)
    {
        //vtkSmartPointer<vtkDoubleArray>data = vtkSmartPointer<vtkDoubleArray>::New();
        vtkDoubleArray*data=vtkDoubleArray::New();
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
}

void AddToVTKVector(vtkFieldData *location,std::vector<double>&array,std::string name)
{
    if(array.size()>0)
    {
        //vtkSmartPointer<vtkDoubleArray>data = vtkSmartPointer<vtkDoubleArray>::New();
        vtkIntArray*data=vtkIntArray::New();
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
}








vtkStandardNewMacro(HDF5ReaderConvex);



HDF5ReaderConvex::HDF5ReaderConvex()
{
    this->FileName = NULL;
    this->DirectoryName = NULL;
    this->SetNumberOfInputPorts(0);
    this->SetNumberOfOutputPorts(1);
    this->times = vtkSmartPointer<vtkDoubleArray>::New();
    this->timesNames = vtkSmartPointer<vtkVariantArray>::New();
}

int HDF5ReaderConvex::RequestData(
        vtkInformation *vtkNotUsed(request),
        vtkInformationVector **vtkNotUsed(inputVector),
        vtkInformationVector *outputVector)
{

    vtkInformation *outInfo = outputVector->GetInformationObject(0);


    double requestedTime = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());
    std::cout<<"  "<<findClosestSolutionIndex(requestedTime)<<"\n";


    vtkUnstructuredGrid *output = vtkUnstructuredGrid::SafeDownCast(
                outInfo->Get(vtkDataObject::DATA_OBJECT()));

    // Here is where you would read the data from the file. In this example,
    // we simply create a point.
    std::cout<<"iddd "<<findClosestSolutionIndex(requestedTime)<<"\n";
    std::string tempFilename=filenames[findClosestSolutionIndex(requestedTime)];
    std::cout<<"filename "<<tempFilename<<"\n";

    //For displacement
    h5cpp::File file1(filenames[findClosestSolutionIndex(0)], "r");
    auto gr0 = file1.root();
    for(size_t i=0;i<gr0.size();i++)
        std::cout<<"Input faile esantys datasetai : "<<i<<" name = "<<gr0.get_link_name(i)<<"\n";
    std::vector<double>oldCENTROID=GetData(gr0,"CENTROID",1);

    file1.close();


    h5cpp::File file(tempFilename, "r");
    auto gr = file.root();

    std::vector<double> STEP;
    STEP.push_back(gr.attrs().get<int>("STEP"));
    std::vector<double> TIME;
    TIME.push_back(gr.attrs().get<double>("TIME"));


    std::vector<double> CENTROID=GetData(gr,"CENTROID",1);
    std::vector<double> VERTICES=GetData(gr,"VERTICES",1);
 //   std::vector<int> VERTICES_COUNT=GetDataINT(gr,"VERTICES_COUNT",1);
    std::vector<int> FACE_IDS=GetDataINT(gr,"FACE_IDS",1);

  //  std::vector<double> VELOCITY=GetData(gr,"VELOCITY",1);
   // std::vector<double> ANGULAR_VELOCITY=GetData(gr,"ANGULAR_VELOCITY",1);
   // std::vector<double> FORCE=GetData(gr,"FORCE",1);
   // std::vector<double> MOMENTUM=GetData(gr,"MOMENTUM",1);
   // std::vector<int> MATERIAL=GetDataINT(gr,"MATERIAL",1);
   // std::vector<int> FIX=GetDataINT(gr,"FIX",1);
    std::vector<double> BOUNDING_RADIUS=GetData(gr,"BOUNDING_RADIUS",1);


    std::vector<double> DEFORM;
    DEFORM.resize(CENTROID.size(),0);
    std::transform(CENTROID.begin(), CENTROID.end(), oldCENTROID.begin(), DEFORM.begin(), [&](double l, double r)
    {
        return (l - r);
    });


    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    output->Allocate();

    int k=0;
    for(int i=0;i<VERTICES_COUNT.size();i++)
    {
        double x=CENTROID[i*4+0];
        double y=CENTROID[i*4+1];
        double z=CENTROID[i*4+2];

        for(int a=0;a<VERTICES_COUNT[i];a++)
        {
            points->InsertNextPoint(VERTICES[k+0]+x,VERTICES[k+1]+y,VERTICES[k+2]+z);
            k=k+4;
        }
    }

    std::vector<int> VERTICES_START;
    VERTICES_START.resize(VERTICES_COUNT.size(),0);
    std::exclusive_scan(VERTICES_COUNT.begin(), VERTICES_COUNT.end(),VERTICES_START.begin(),0);

    int i=0;
    for(int particleID=0;particleID<CENTROID.size()/4;particleID++)
    {
        int nfaces=FACE_IDS[i];i++;
        vtkSmartPointer<vtkCellArray> dodechedronFaces = vtkSmartPointer<vtkCellArray>::New();
        for(int k=0;k<nfaces;k++)
        {
            int Vcount=FACE_IDS[i];i++;
            dodechedronFaces->InsertNextCell(Vcount);
            for(int v=0;v<Vcount;v++)
            {
                dodechedronFaces->InsertCellPoint(FACE_IDS[i]+VERTICES_START[particleID]);i++;
            }
        }
        int VertexCount=VERTICES_COUNT[particleID];
        vtkIdType dodechedronPointsIds[VertexCount];
        for(int b=0;b<VertexCount;b++)
        {
            dodechedronPointsIds[b]=VERTICES_START[particleID]+b;
        }
        output->InsertNextCell(VTK_POLYHEDRON,
                               VertexCount, dodechedronPointsIds,
                               nfaces, dodechedronFaces->GetPointer());
    }



    output->SetPoints(points);
    vtkCellData* cdata = output->GetCellData();
    vtkFieldData* fdata = output->GetFieldData();

    AddToVTKVector(cdata,DEFORM,"DISPLACEMENT");
    AddToVTKVector(cdata,CENTROID,"CENTROID");
    AddToVTKScalar(cdata,BOUNDING_RADIUS,"BOUNDING_RADIUS");
    AddToVTKScalar(fdata,STEP,"STEP");
    AddToVTKScalar(fdata,TIME,"TIME");


    //AddToVTKVector(cdata,VELOCITY,"VELOCITY");
    //AddToVTKVector(cdata,ANGULAR_VELOCITY,"ANGULAR_VELOCITY");
    //AddToVTKVector(cdata,FORCE,"FORCE");
   // AddToVTKVector(cdata,MOMENTUM,"MOMENTUM");

   //AddToVTKScalar(cdata,VERTICES_COUNT,"VERTICES_COUNT");
  // AddToVTKScalar(cdata,MATERIAL,"MATERIAL");
  // AddToVTKScalar(cdata,FIX,"FIX");


    file.close();
    return 1;
}

void HDF5ReaderConvex::PrintSelf(ostream& os, vtkIndent indent)
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

int HDF5ReaderConvex::RequestInformation(vtkInformation *vtkNotUsed(request), vtkInformationVector **vtkNotUsed(inputVector), vtkInformationVector *outputVector)
{
    if (!this->FileName || strlen(this->FileName) == 0)
    {
        std::cout<<"FileName has to be specified!";
        return 0;
    }
    filenames.resize(0,"");
    std::string fileStringStart=fs::path(FileName).filename().string().substr(0,fs::path(FileName).filename().string().length()-12);
    std::cout<<fs::path(FileName).parent_path()<<"   "<<FileName<<"  "<<fileStringStart<<"\n";

    for(auto& p: fs::directory_iterator(fs::path(FileName).parent_path()))
    {
        if (p.path().string().find(".ch5") != std::string::npos)
            if (p.path().string().find(fileStringStart) != std::string::npos) {

                filenames.push_back(p.path().string());
            }
    }
    std::sort(filenames.begin(),filenames.end());

    for (int i = 0; i < filenames.size(); i++) {
        this->times->InsertNextValue(i);
    }


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




    return 1;
}
