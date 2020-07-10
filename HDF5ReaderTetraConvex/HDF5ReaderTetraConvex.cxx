#include "HDF5ReaderTetraConvex.h"

#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkInformationVector.h"
#include "vtkInformation.h"
#include "vtkDataObject.h"
#include "vtkSmartPointer.h"
#include "vtkVertexGlyphFilter.h"

#include <algorithm>
//#include "hdf_wrapper.hpp"
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
#define GetParticleFix(x) (x >> 24) & 0xFF
#define GetParticleWall(x) (x >> 16) & 0xFF
#define GetParticleMaterial(x) (x >> 8) & 0xFF
#define GetParticleLaisvas(x) x  & 0xFF

#include <boost/array.hpp>
typedef boost::array<double,4> REAL4;
typedef boost::array<int,4> INT4;
typedef boost::array<int,2> INT2;
namespace fs = std::experimental::filesystem;
vtkIdType HDF5ReaderTetraConvex::findClosestSolutionIndex(double Time)
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

#include <hdf5.h>
#include "H5Cpp.h"

typedef struct
{
    unsigned char FIX;
    unsigned char MATERIAL;
    unsigned char WALL;
    unsigned char LAISVAS;
} PARTICLE_INFO_TYPE;
typedef std::vector<PARTICLE_INFO_TYPE> PARTICLE_INFO_ARRAY;

auto WriteParticleInfo=[](H5::H5File* file,std::vector<PARTICLE_INFO_TYPE> &data,size_t COUNT)->void
{
    hsize_t dim[] = {COUNT};
    H5::DataSpace space( 1, dim );
    H5::CompType mtype1( sizeof(PARTICLE_INFO_TYPE) );
    mtype1.insertMember( "FIX", HOFFSET(PARTICLE_INFO_TYPE, FIX), H5::PredType::NATIVE_UCHAR);
    mtype1.insertMember( "MATERIAL", HOFFSET(PARTICLE_INFO_TYPE, MATERIAL), H5::PredType::NATIVE_UCHAR);
    mtype1.insertMember( "WALL", HOFFSET(PARTICLE_INFO_TYPE, WALL), H5::PredType::NATIVE_UCHAR);
    mtype1.insertMember( "LAISVAS", HOFFSET(PARTICLE_INFO_TYPE, LAISVAS), H5::PredType::NATIVE_UCHAR);
    H5::DataSet* dataset = new H5::DataSet(file->createDataSet("PARTICLE_INFO", mtype1, space));
    dataset->write( &data[0], mtype1 );
    delete dataset;
};


auto ReadParticleInfo=[](H5::H5File* file,std::vector<PARTICLE_INFO_TYPE> &data,size_t COUNT)->void
{
    data.resize(COUNT);
    H5::DataSet* dataset = new H5::DataSet (file->openDataSet( "PARTICLE_INFO" ));
    H5::CompType mtype1( sizeof(PARTICLE_INFO_TYPE) );
    mtype1.insertMember( "FIX", HOFFSET(PARTICLE_INFO_TYPE, FIX), H5::PredType::NATIVE_UCHAR);
    mtype1.insertMember( "MATERIAL", HOFFSET(PARTICLE_INFO_TYPE, MATERIAL), H5::PredType::NATIVE_UCHAR);
    mtype1.insertMember( "WALL", HOFFSET(PARTICLE_INFO_TYPE, WALL), H5::PredType::NATIVE_UCHAR);
    mtype1.insertMember( "LAISVAS", HOFFSET(PARTICLE_INFO_TYPE, LAISVAS), H5::PredType::NATIVE_UCHAR);
    dataset->read( &data[0], mtype1 );
    delete dataset;
};

auto ReadIntArray=[](H5::H5File* file,std::vector<int> &data,size_t COUNT,std::string name)->void
{
    data.resize(COUNT);
    H5::DataSet* dataset = new H5::DataSet (file->openDataSet( name ));
    dataset->read( &data[0], H5::PredType::NATIVE_INT );
    delete dataset;
};


auto ReadRealArray=[](H5::H5File* file,std::vector<double> &data,size_t COUNT,std::string name)->void
{
    data.resize(COUNT);
    H5::DataSet* dataset = new H5::DataSet (file->openDataSet( name ));
    dataset->read( &data[0], H5::PredType::NATIVE_DOUBLE );
    delete dataset;
};


auto ReadReal4Array=[](H5::H5File* file,std::vector<REAL4> &data,size_t COUNT,std::string name)->void
{
    data.resize(COUNT);
    H5::DataSet* dataset = new H5::DataSet (file->openDataSet( name ));
    dataset->read( &data[0], H5::PredType::NATIVE_DOUBLE );
    delete dataset;
};


auto ReadIntAttribute=[](H5::H5File*file,std::string attribute_name)->int
{
    H5::DataSpace att_space(H5S_SCALAR);
    auto att=file->openAttribute(attribute_name);
    int val=0;
    att.read(H5::PredType::NATIVE_INT,&val);
    return val;
};
auto ReadDoubleAttribute=[](H5::H5File*file,std::string attribute_name)->double
{
    H5::DataSpace att_space(H5S_SCALAR);
    auto att=file->openAttribute(attribute_name);
    double val=0;
    att.read(H5::PredType::NATIVE_DOUBLE,&val);
    return val;
};

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



vtkStandardNewMacro(HDF5ReaderTetraConvex);


HDF5ReaderTetraConvex::HDF5ReaderTetraConvex()
{
    this->FileName = NULL;
    this->DirectoryName = NULL;
    this->SetNumberOfInputPorts(0);
    this->SetNumberOfOutputPorts(1);
    this->times = vtkSmartPointer<vtkDoubleArray>::New();
    this->timesNames = vtkSmartPointer<vtkVariantArray>::New();
}

int HDF5ReaderTetraConvex::RequestData(
        vtkInformation *vtkNotUsed(request),
        vtkInformationVector **vtkNotUsed(inputVector),
        vtkInformationVector *outputVector)
{

    vtkInformation *outInfo = outputVector->GetInformationObject(0);
    double requestedTime = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());
    vtkUnstructuredGrid *output = vtkUnstructuredGrid::SafeDownCast(
                outInfo->Get(vtkDataObject::DATA_OBJECT()));


    std::string tempFilename=filenames[findClosestSolutionIndex(requestedTime)];

    std::cout<<"cia as "<<tempFilename<<"\n";


    std::vector<REAL4> HOST_POINTS_INIT;
    {
        H5::H5File* fileInit= new H5::H5File( filenames[0], H5F_ACC_RDONLY );
        size_t PARTICLE_COUNT=ReadIntAttribute(fileInit,"PARTICLE_COUNT");
        ReadReal4Array(fileInit,HOST_POINTS_INIT,PARTICLE_COUNT*4,"POINTS");
        delete fileInit;

    }



    H5::H5File* file= new H5::H5File( tempFilename, H5F_ACC_RDONLY );
    size_t PARTICLE_COUNT=ReadIntAttribute(file,"PARTICLE_COUNT");
    size_t TIME_STEP=ReadIntAttribute(file,"STEP");
    double TIME=ReadDoubleAttribute(file,"TIME");

    std::vector<REAL4> HOST_POINTS;
    std::vector<PARTICLE_INFO_TYPE> HOST_INFO;

    ReadReal4Array(file,HOST_POINTS,PARTICLE_COUNT*4,"POINTS");
    ReadParticleInfo(file,HOST_INFO,PARTICLE_COUNT);

    std::vector<int> FIX(PARTICLE_COUNT,0);
    std::vector<int> MATERIAL(PARTICLE_COUNT,0);
    std::vector<int> WALL(PARTICLE_COUNT,0);
    std::vector<int> LAISVAS(PARTICLE_COUNT,0);
     vtkPoints *points = vtkPoints::New();
    points->SetNumberOfPoints(0);
     points->SetDataTypeToDouble();
     output->Allocate(1000,1000);
     output->SetPoints(points);
    std::vector<REAL4> DISPLACEMENT;
    DISPLACEMENT.resize(PARTICLE_COUNT);
    auto CalcCentroid=[](REAL4 P1,REAL4 P2,REAL4 P3,REAL4 P4)->REAL4
    {
        REAL4 center=(REAL4){0,0,0,0};
        for(int i=0;i<3;i++)
        {
            center[i]=(P1[i]+P2[i]+P3[i]+P4[i])/4.0;
        }

        return center;

    };

    for(size_t i=0;i<PARTICLE_COUNT;i++)
    {
        REAL4 P1=HOST_POINTS[i*4+0];
        REAL4 P2=HOST_POINTS[i*4+1];
        REAL4 P3=HOST_POINTS[i*4+2];
        REAL4 P4=HOST_POINTS[i*4+3];


        REAL4 IP1=HOST_POINTS_INIT[i*4+0];
        REAL4 IP2=HOST_POINTS_INIT[i*4+1];
        REAL4 IP3=HOST_POINTS_INIT[i*4+2];
        REAL4 IP4=HOST_POINTS_INIT[i*4+3];
        PARTICLE_INFO_TYPE info=HOST_INFO[i];


        FIX[i]=info.FIX;
        MATERIAL[i]=info.MATERIAL;
        WALL[i]=info.WALL;
        LAISVAS[i]=info.LAISVAS;

        REAL4 INIT_CENTER=CalcCentroid(IP1,IP2,IP3,IP4);
        REAL4 CENTER=CalcCentroid(P1,P2,P3,P4);
        CENTER[0]=CENTER[0]-INIT_CENTER[0];
        CENTER[1]=CENTER[1]-INIT_CENTER[1];
        CENTER[2]=CENTER[2]-INIT_CENTER[2];
        DISPLACEMENT[i]=CENTER;


        int offsetas=points->GetNumberOfPoints();
        points->InsertNextPoint(P1[0],P1[1],P1[2]);
        points->InsertNextPoint(P2[0],P2[1],P2[2]);
        points->InsertNextPoint(P3[0],P3[1],P3[2]);
        points->InsertNextPoint(P4[0],P4[1],P4[2]);
        vtkIdType ptIds[] = {offsetas, offsetas+1, offsetas+2, offsetas+3};
        output->InsertNextCell( VTK_TETRA, 4, ptIds );
    }




    std::cout<<"output number of cells "<<output->GetNumberOfCells()<<"\n";
    output->SetPoints(points);
    vtkCellData* cdata =output->GetCellData();
    vtkFieldData* fdata = output->GetFieldData();
    std::vector<int> STEPAS;
    STEPAS.push_back(TIME_STEP);
    std::vector<double> TIMES;
    TIMES.push_back(TIME);
    AddToVTKScalar(fdata,STEPAS,"STEP");
    AddToVTKScalar(fdata,TIMES,"TIME");

    AddToVTKScalar(cdata,MATERIAL,"MATERIAL");
    AddToVTKScalar(cdata,FIX,"FIX");
    AddToVTKScalar(cdata,WALL,"WALL");
    AddToVTKScalar(cdata,LAISVAS,"LAISVAS");

    AddToVTKVector(cdata,DISPLACEMENT,"DISPLACEMENT");
    if(file->exists("VELOCITY"))
    {
        std::vector<REAL4> HOST_VELOCITY;
        ReadReal4Array(file,HOST_VELOCITY,PARTICLE_COUNT,"VELOCITY");
        AddToVTKVector(cdata,HOST_VELOCITY,"VELOCITY");
    }
    if(file->exists("ANGULAR_VELOCITY"))
    {
        std::vector<REAL4> HOST_ANGULAR_VELOCITY;
        ReadReal4Array(file,HOST_ANGULAR_VELOCITY,PARTICLE_COUNT,"ANGULAR_VELOCITY");
        AddToVTKVector(cdata,HOST_ANGULAR_VELOCITY,"ANGULAR_VELOCITY");
    }

    if(file->exists("FORCE"))
    {
        std::vector<REAL4> HOST_FORCE;
        ReadReal4Array(file,HOST_FORCE,PARTICLE_COUNT,"FORCE");
        AddToVTKVector(cdata,HOST_FORCE,"FORCE");
    }
    if(file->exists("TORQUE"))
    {
        std::vector<REAL4> HOST_TORQUE;
        ReadReal4Array(file,HOST_TORQUE,PARTICLE_COUNT,"TORQUE");
        AddToVTKVector(cdata,HOST_TORQUE,"TORQUE");
    }

    if(file->exists("NN_COUNT"))
    {
        std::vector<int> HOST_NN_COUNT;
        ReadIntArray(file,HOST_NN_COUNT,PARTICLE_COUNT,"NN_COUNT");
        AddToVTKScalar(cdata,HOST_NN_COUNT,"NN_COUNT");
    }
    if(file->exists("PP_COUNT"))
    {
        std::vector<int> HOST_PP_COUNT;
        ReadIntArray(file,HOST_PP_COUNT,PARTICLE_COUNT,"PP_COUNT");
        AddToVTKScalar(cdata,HOST_PP_COUNT,"PP_COUNT");
    }
    if(file->exists("MASS"))
    {
        std::vector<double> HOST_MASS;
        ReadRealArray(file,HOST_MASS,PARTICLE_COUNT,"MASS");
        AddToVTKScalar(cdata,HOST_MASS,"MASS");
    }

    if(file->exists("VOLUME"))
    {
        std::vector<double> HOST_VOLUME;
        ReadRealArray(file,HOST_VOLUME,PARTICLE_COUNT,"VOLUME");
        AddToVTKScalar(cdata,HOST_VOLUME,"VOLUME");
    }
    cdata->SetActiveVectors("VELOCITY");
    cdata->SetActiveScalars("FIX");






    delete file;


    output->Modified();
    return 1;
}

void HDF5ReaderTetraConvex::PrintSelf(ostream& os, vtkIndent indent)
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

int HDF5ReaderTetraConvex::RequestInformation(vtkInformation *vtkNotUsed(request), vtkInformationVector **vtkNotUsed(inputVector), vtkInformationVector *outputVector)
{
    if (!this->FileName || strlen(this->FileName) == 0)
    {
        std::cout<<"FileName has to be specified!";
        return 0;
    }
    filenames.resize(0,"");
    std::string fileStringStart=fs::path(FileName).filename().string().substr(0,fs::path(FileName).filename().string().length()-12);
    std::cout<<"HDF5ReaderTetraConvex "<<fs::path(FileName).parent_path()<<"   "<<FileName<<"  "<<fileStringStart<<"\n";

    for(auto& p: fs::directory_iterator(fs::path(FileName).parent_path()))
    {
        if (p.path().string().find(".th5") != std::string::npos)
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
