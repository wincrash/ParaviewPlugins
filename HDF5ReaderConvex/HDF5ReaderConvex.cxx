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

#include <boost/array.hpp>
typedef boost::array<double,4> REAL4;
typedef boost::array<int,4> INT4;
typedef boost::array<int,2> INT2;
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
        vtkIntArray*data=vtkIntArray::New();
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
        vtkDoubleArray*data=vtkDoubleArray::New();
        data->SetNumberOfTuples(array.size());
        data->SetNumberOfComponents(3);
        data->SetNumberOfValues(array.size());
        data->SetName(name.c_str());
        for (int ii = 0; ii < array.size(); ii++) {
            data->SetTuple1(ii, array[ii]);
        }
        location->AddArray(data);
        data->Delete();
    }
}

void AddToVTKVector(vtkFieldData *location,std::vector<REAL4>&array,std::string name)
{
    if(array.size()>0)
    {
        //vtkSmartPointer<vtkDoubleArray>data = vtkSmartPointer<vtkDoubleArray>::New();
        std::cout<<"bandoma ideti "<<name<<" Su dydziu "<<array.size()<<"\n";
        vtkDoubleArray*data=vtkDoubleArray::New();
        // data->SetNumberOfTuples(array.size());
        data->SetNumberOfComponents(3);
        data->SetName(name.c_str());
        for (int ii = 0; ii < array.size(); ii++) {
            REAL4 t=array[ii];
            data->InsertNextTuple3( t[0],t[1],t[2]);
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

    auto readDataSetREAL4=[](auto &gr,auto &array,std::string name,size_t N)
    {
        std::cout<<"Reading from file REAL4 size="<<N<<" array name="<<name<<"\n";
        for(size_t ii=0;ii<gr.size();ii++)
            if(gr.get_link_name(ii).compare(name)==0){
                auto dt = gr.open_dataset(name.c_str());
                array.resize(N,REAL4({1,2,3,4}));
                std::vector<double> tempArray(N*4);
                dt.read(&tempArray[0]);
                for(size_t i=0;i<N;i++)
                {
                    array[i]=REAL4({tempArray[i*4],tempArray[i*4+1],tempArray[i*4+2],tempArray[i*4+3]});
                }
            }
    };

    auto readDataSetREAL=[](auto &gr,auto &array,std::string name,size_t N)
    {
        std::cout<<"Reading from file REAL size="<<N<<" array name="<<name<<"\n";
        for(size_t ii=0;ii<gr.size();ii++)
            if(gr.get_link_name(ii).compare(name)==0){
                auto dt = gr.open_dataset(name.c_str());
                array.resize(N);
                std::vector<double> tempArray(N);
                dt.read(&tempArray[0]);
                for(size_t i=0;i<N;i++)
                {
                    array[i]=tempArray[i];
                }
            }
    };

    auto readDataSetINT2=[](auto &gr,auto &array,std::string name,size_t N)
    {
        std::cout<<"Reading from file INT2 size="<<N<<" array name="<<name<<"\n";
        for(size_t ii=0;ii<gr.size();ii++)
            if(gr.get_link_name(ii).compare(name)==0){
                auto dt = gr.open_dataset(name.c_str());
                std::cout<<"Reading from file INT2 size="<<N<<" array name="<<name<<"\n";
                array.resize(N);
                std::cout<<"Reading from file INT2 size="<<N<<" array name="<<name<<"\n";
                std::vector<int> tempArray(N*2);
                std::cout<<"Reading from file INT2 size="<<N<<" array name="<<name<<"\n";
                dt.read(&tempArray[0]);
                std::cout<<"Reading from file INT2 size="<<N<<" array name="<<name<<"\n";
                for(size_t i=0;i<N;i++)
                {
                    array[i]=INT2({tempArray[i*2],tempArray[i*2+1]});
                }
                std::cout<<"Reading from file INT2 size="<<N<<" array name="<<name<<"\n";
            }
    };
    auto readDataSetINT4=[](auto &gr,auto &array,std::string name,size_t N)
    {
        std::cout<<"Reading from file INT4 size="<<N<<" array name="<<name<<"\n";
        for(size_t ii=0;ii<gr.size();ii++)
            if(gr.get_link_name(ii).compare(name)==0){
                auto dt = gr.open_dataset(name.c_str());
                array.resize(N);
                std::vector<int> tempArray(N*4);
                dt.read(&tempArray[0]);
                for(size_t i=0;i<N;i++)
                {
                    array[i]=INT4({tempArray[i*4],tempArray[i*4+1],tempArray[i*4+2],tempArray[i*4+3]});
                }
            }
    };
    auto readDataSetINT=[](auto &gr,auto &array,std::string name,size_t N)
    {
        std::cout<<"Reading from file INT size="<<N<<" array name="<<name<<"\n";
        for(size_t ii=0;ii<gr.size();ii++)
            if(gr.get_link_name(ii).compare(name)==0){
                auto dt = gr.open_dataset(name.c_str());

                std::vector<int> tempArray(N);
                array.resize(N);
                dt.read(&tempArray[0]);
                for(size_t i=0;i<N;i++)
                {
                    array[i]=tempArray[i];
                }
            }
    };

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

    h5cpp::File file(tempFilename, "r");
    auto gr = file.root();


    for(size_t i=0;i<gr.size();i++)
        std::cout<<"Input faile esantys datasetai : "<<i<<" name = "<<gr.get_link_name(i)<<"\n";

    size_t PARTICLE_COUNT=gr.attrs().get<int>("PARTICLE_COUNT");
    size_t MAX_FACE_COUNT=gr.attrs().get<int>("MAX_FACE_COUNT");
    size_t MAX_VERTICES_COUNT=gr.attrs().get<int>("MAX_VERTICES_COUNT");
    size_t MAX_FACE_INDEX_COUNT=gr.attrs().get<int>("MAX_FACE_INDEX_COUNT");
    size_t FACE_INDEX_OFFSET=gr.attrs().get<int>("FACE_INDEX_OFFSET");


    std::vector<double> STEP;
    STEP.push_back(gr.attrs().get<int>("STEP"));
    std::vector<double> TIME;
    TIME.push_back(gr.attrs().get<double>("TIME"));


    std::vector<REAL4> CENTROID;
    std::vector<REAL4> VERTICES;
    std::vector<INT2> PARTICLE_GEOMETRY_INFO;
    std::vector<int> FIX;
    readDataSetINT(gr,FIX,"FIX",PARTICLE_COUNT);
    std::vector<int> MATERIAL;
    readDataSetINT(gr,MATERIAL,"MATERIAL",PARTICLE_COUNT);
    std::vector<REAL4> FORCE;
    readDataSetREAL4(gr,FORCE,"FORCE",PARTICLE_COUNT);
    std::vector<REAL4> TORQUE;
    readDataSetREAL4(gr,TORQUE,"TORQUE",PARTICLE_COUNT);
    std::vector<REAL4> VELOCITY;
    readDataSetREAL4(gr,VELOCITY,"VELOCITY",PARTICLE_COUNT);

    std::vector<int> NN_COUNT;
    readDataSetINT(gr,NN_COUNT,"NN_COUNT",PARTICLE_COUNT);

    std::vector<double> BOUNDING_RADIUS;
    readDataSetREAL(gr,BOUNDING_RADIUS,"BOUNDING_RADIUS",PARTICLE_COUNT);

    std::vector<REAL4> ANGULAR_VELOCITY;
    readDataSetREAL4(gr,ANGULAR_VELOCITY,"ANGULAR_VELOCITY",PARTICLE_COUNT);

    readDataSetREAL4(gr,CENTROID,"CENTROID",PARTICLE_COUNT);
    readDataSetREAL4(gr,VERTICES,"VERTICES",PARTICLE_COUNT*MAX_VERTICES_COUNT);
    readDataSetINT2(gr,PARTICLE_GEOMETRY_INFO,"PARTICLE_GEOMETRY_INFO",PARTICLE_COUNT);
    std::vector<int> FACE_IDS=GetDataINT(gr,"FACE_IDS",1);



    vtkPoints *points = vtkPoints::New();
    output->Allocate();

    std::vector<int> VERTICES_COUNT(PARTICLE_COUNT,0);
    for(size_t i=0;i<PARTICLE_COUNT;i++)
    {
        REAL4 center=CENTROID[i];
        size_t verticeCount=PARTICLE_GEOMETRY_INFO[i][0];
        VERTICES_COUNT[i]=verticeCount;
        for(size_t a=0;a<verticeCount;a++)
        {
            REAL4 v=VERTICES[i*MAX_VERTICES_COUNT+a];
            points->InsertNextPoint(v[0]+center[0],v[1]+center[1],v[2]+center[2]);
        }
    }

    std::vector<int> VERTICES_START;
    VERTICES_START.resize(PARTICLE_COUNT,0);
    std::exclusive_scan(VERTICES_COUNT.begin(), VERTICES_COUNT.end(),VERTICES_START.begin(),0);


    for(int particleID=0;particleID<PARTICLE_COUNT;particleID++)
    {
        size_t nfaces=PARTICLE_GEOMETRY_INFO[particleID][1];
        int i=particleID*FACE_INDEX_OFFSET;
        vtkCellArray* dodechedronFaces = vtkCellArray::New();
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


    std::cout<<"output number of cells "<<output->GetNumberOfCells()<<"\n";
    output->SetPoints(points);
    vtkCellData* cdata =output->GetCellData();
    vtkFieldData* fdata = output->GetFieldData();

    //AddToVTKVector(cdata,DEFORM,"DISPLACEMENT");
    AddToVTKVector(cdata,CENTROID,"CENTROID");
    //AddToVTKScalar(cdata,BOUNDING_RADIUS,"BOUNDING_RADIUS");
    AddToVTKScalar(fdata,STEP,"STEP");
    AddToVTKScalar(fdata,TIME,"TIME");



    AddToVTKVector(cdata,VELOCITY,"VELOCITY");
    AddToVTKVector(cdata,ANGULAR_VELOCITY,"ANGULAR_VELOCITY");
    AddToVTKVector(cdata,FORCE,"FORCE");
    AddToVTKVector(cdata,TORQUE,"TORQUE");

    AddToVTKScalar(cdata,VERTICES_COUNT,"VERTICES_COUNT");
    AddToVTKScalar(cdata,MATERIAL,"MATERIAL");
    AddToVTKScalar(cdata,FIX,"FIX");


    AddToVTKScalar(cdata,NN_COUNT,"NN_COUNT");
    AddToVTKScalar(cdata,BOUNDING_RADIUS,"BOUNDING_RADIUS");

    output->Modified();
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
