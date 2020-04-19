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
#include <vtkTetra.h>

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

    std::vector<double> STEP;
    STEP.push_back(gr.attrs().get<int>("STEP"));
    std::vector<double> TIME;
    TIME.push_back(gr.attrs().get<double>("TIME"));




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

    std::vector<int> PP_COUNT;
    readDataSetINT(gr,PP_COUNT,"PP_COUNT",PARTICLE_COUNT);

    std::vector<REAL4> ANGULAR_VELOCITY;
    readDataSetREAL4(gr,ANGULAR_VELOCITY,"ANGULAR_VELOCITY",PARTICLE_COUNT);

    std::vector<double> MASS;
    readDataSetREAL(gr,MASS,"MASS",PARTICLE_COUNT);
    std::vector<double> VOLUME;
    readDataSetREAL(gr,VOLUME,"VOLUME",PARTICLE_COUNT);




    /*
    group.attrs().create<int>("FACE_IDS_COUNT", FACE_IDS_COUNT);
    group.attrs().create<int>("POINT_COUNT", POINT_COUNT);
    WriteToHdfReal4(group, POINTS, "POINTS");

    WriteToHdfInt(group, FACE_COUNT, "FACE_COUNT");
    WriteToHdfInt(group, FACES_START, "FACES_START");
    WriteToHdfInt(group, POINTS_START, "POINTS_START");
    */

int FACE_IDS_COUNT=gr.attrs().get<int>("FACE_IDS_COUNT");
int POINT_COUNT=gr.attrs().get<int>("POINT_COUNT");
    std::vector<REAL4> POINTS;
    readDataSetREAL4(gr,POINTS,"POINTS",POINT_COUNT);
    std::vector<int> FACE_IDS;
    readDataSetINT(gr,FACE_IDS,"FACE_IDS",FACE_IDS_COUNT);
    std::vector<int> FACE_COUNT;
    readDataSetINT(gr,FACE_COUNT,"FACE_COUNT",PARTICLE_COUNT);
    std::vector<int> FACES_START;
    readDataSetINT(gr,FACES_START,"FACES_START",PARTICLE_COUNT);
    std::vector<int> POINTS_START;
    readDataSetINT(gr,POINTS_START,"POINTS_START",PARTICLE_COUNT);
    std::vector<int> POINTS_COUNT;
    readDataSetINT(gr,POINTS_COUNT,"POINTS_COUNT",PARTICLE_COUNT);

    vtkPoints *points = vtkPoints::New();
    output->Allocate(1000,1000);

    int index=0;
    int findex=0;
    for(int i=0;i<PARTICLE_COUNT;i++)
    {
        int point_startas=POINTS_START[i];
        vtkIdType dodechedronPointsIds[POINTS_COUNT[i]];
        for(int z=0;z<POINTS_COUNT[i];z++)
        {
            dodechedronPointsIds[z]=index;
            REAL4 P=POINTS[index++];
            points->InsertNextPoint(P[0],P[1],P[2]);
        }

        vtkSmartPointer<vtkCellArray> dodechedronFaces =vtkSmartPointer<vtkCellArray>::New();
        int face_count=FACE_COUNT[i];
        for (size_t ii = 0; ii < face_count; ii++)
        {
            int fsize=FACE_IDS[findex++];
            dodechedronFaces->InsertNextCell(fsize);
            for(size_t k=0;k<fsize;k++)
                dodechedronFaces->InsertCellPoint(FACE_IDS[findex++]+point_startas);
        }

        output->InsertNextCell(VTK_POLYHEDRON,
                               POINTS_COUNT[i], dodechedronPointsIds,
                               face_count, dodechedronFaces->GetPointer());
    }













    std::cout<<"output number of cells "<<output->GetNumberOfCells()<<"\n";
    output->SetPoints(points);
    vtkCellData* cdata =output->GetCellData();
    vtkFieldData* fdata = output->GetFieldData();

    AddToVTKScalar(fdata,STEP,"STEP");
    AddToVTKScalar(fdata,TIME,"TIME");



    AddToVTKVector(cdata,VELOCITY,"VELOCITY");
    AddToVTKVector(cdata,ANGULAR_VELOCITY,"ANGULAR_VELOCITY");
    AddToVTKVector(cdata,FORCE,"FORCE");
    AddToVTKVector(cdata,TORQUE,"TORQUE");

    AddToVTKScalar(cdata,MATERIAL,"MATERIAL");
    AddToVTKScalar(cdata,FIX,"FIX");
        AddToVTKScalar(cdata,POINTS_COUNT,"POINTS_COUNT");
        AddToVTKScalar(cdata,FACE_COUNT,"FACE_COUNT");


    AddToVTKScalar(cdata,NN_COUNT,"NN_COUNT");
    AddToVTKScalar(cdata,PP_COUNT,"PP_COUNT");
    AddToVTKScalar(cdata,MASS,"MASS");
    AddToVTKScalar(cdata,VOLUME,"VOLUME");

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
