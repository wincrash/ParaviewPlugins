#include "HDF5ReaderTetraBonds.h"

#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkInformationVector.h"
#include "vtkInformation.h"
#include "vtkDataObject.h"
#include "vtkSmartPointer.h"
#include "vtkVertexGlyphFilter.h"

#include <algorithm>
//#include "hdf_wrapper.hpp"
#include "vtkPolyData.h"
#include <iostream>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkFieldData.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkTriangle.h>
#include <vtkPolyData.h>
#include <vtkVertex.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkPolyData.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkVariantArray.h>
#include <experimental/filesystem>
#include <vtkMultiBlockDataSet.h>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/lexical_cast.hpp>
#include <vtkPolyData.h>
#include <vtkPolyDataAlgorithm.h>
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
vtkIdType HDF5ReaderTetraBonds::findClosestSolutionIndex(double Time)
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


typedef struct
{
    unsigned char FIX;
    unsigned char MATERIAL;
    unsigned char WALL;
    unsigned char LAISVAS;
} PARTICLE_INFO_TYPE;
typedef std::vector<PARTICLE_INFO_TYPE> PARTICLE_INFO_ARRAY;



typedef struct
{
    int ID1;
    int ID2;
    unsigned char STATE;
    unsigned char TYPE;
    unsigned char ID_A1;
    unsigned char ID_A2;
    unsigned char ID_A3;
    unsigned char ID_B1;
    unsigned char ID_B2;
    unsigned char ID_B3;
}BONDS_TYPE;
typedef std::vector<BONDS_TYPE> BONDS;
#include <hdf5.h>
#include "H5Cpp.h"
auto ReadBondsInfo=[](H5::H5File* file,std::vector<BONDS_TYPE> &data,size_t COUNT)->void
{
    data.resize(COUNT);
    H5::DataSet* dataset = new H5::DataSet (file->openDataSet( "BONDS" ));
    H5::CompType mtype1( sizeof(BONDS_TYPE));
    mtype1.insertMember( "ID1", HOFFSET(BONDS_TYPE, ID1), H5::PredType::NATIVE_INT);
    mtype1.insertMember( "ID2", HOFFSET(BONDS_TYPE, ID2), H5::PredType::NATIVE_INT);
    mtype1.insertMember( "STATE", HOFFSET(BONDS_TYPE, STATE), H5::PredType::NATIVE_UCHAR);
    mtype1.insertMember( "TYPE", HOFFSET(BONDS_TYPE, TYPE), H5::PredType::NATIVE_UCHAR);
    mtype1.insertMember( "ID_A1", HOFFSET(BONDS_TYPE, ID_A1), H5::PredType::NATIVE_UCHAR);
    mtype1.insertMember( "ID_A2", HOFFSET(BONDS_TYPE, ID_A2), H5::PredType::NATIVE_UCHAR);
    mtype1.insertMember( "ID_A3", HOFFSET(BONDS_TYPE, ID_A3), H5::PredType::NATIVE_UCHAR);
    mtype1.insertMember( "ID_B1", HOFFSET(BONDS_TYPE, ID_B1), H5::PredType::NATIVE_UCHAR);
    mtype1.insertMember( "ID_B2", HOFFSET(BONDS_TYPE, ID_B2), H5::PredType::NATIVE_UCHAR);
    mtype1.insertMember( "ID_B3", HOFFSET(BONDS_TYPE, ID_B3), H5::PredType::NATIVE_UCHAR);
    dataset->read( &data[0], mtype1 );
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



vtkStandardNewMacro(HDF5ReaderTetraBonds);


HDF5ReaderTetraBonds::HDF5ReaderTetraBonds()
{
    this->FileName = NULL;
    this->DirectoryName = NULL;
    this->SetNumberOfInputPorts(0);
    this->SetNumberOfOutputPorts(1);
    this->times = vtkSmartPointer<vtkDoubleArray>::New();
    this->timesNames = vtkSmartPointer<vtkVariantArray>::New();
}

int HDF5ReaderTetraBonds::RequestData(
        vtkInformation *vtkNotUsed(request),
        vtkInformationVector **vtkNotUsed(inputVector),
        vtkInformationVector *outputVector)
{

    vtkInformation *outInfo = outputVector->GetInformationObject(0);
    double requestedTime = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());
    vtkPolyData *output = vtkPolyData::SafeDownCast(
                outInfo->Get(vtkDataObject::DATA_OBJECT()));


    std::string tempFilename=filenames[findClosestSolutionIndex(requestedTime)];

    std::cout<<tempFilename<<"\n";


    std::vector<REAL4> HOST_POINTS_INIT;
    {
        H5::H5File* fileInit= new H5::H5File( filenames[0], H5F_ACC_RDONLY );
        size_t PARTICLE_COUNT=ReadIntAttribute(fileInit,"PARTICLE_COUNT");
        ReadReal4Array(fileInit,HOST_POINTS_INIT,PARTICLE_COUNT*4,"POINTS");
        delete fileInit;

    }

    std::vector<REAL4> HOST_POINTS;
    std::vector<PARTICLE_INFO_TYPE> HOST_PARTICLE_INFO;

    H5::H5File* file= new H5::H5File( tempFilename, H5F_ACC_RDONLY );
    size_t PARTICLE_COUNT=ReadIntAttribute(file,"PARTICLE_COUNT");
    size_t BOND_COUNT=ReadIntAttribute(file,"BOND_COUNT");
    size_t TIME_STEP=ReadIntAttribute(file,"STEP");
    double TIME=ReadDoubleAttribute(file,"TIME");

    ReadReal4Array(file,HOST_POINTS,PARTICLE_COUNT*4,"POINTS");



    if(BOND_COUNT>0)
    {
        std::vector<BONDS_TYPE> HOST_BONDS;
        ReadBondsInfo(file,HOST_BONDS,BOND_COUNT);
        vtkPoints *pointsB = vtkPoints::New();
        pointsB->SetNumberOfPoints(0);
        pointsB->SetDataTypeToDouble();
        vtkCellArray*cells=vtkCellArray::New();
        std::vector<int> BOND_STATE(BOND_COUNT,0);
        std::vector<int> BOND_TYPE(BOND_COUNT,0);
        std::vector<REAL4> displacement(BOND_COUNT);
        for(size_t i=0;i<BOND_COUNT;i++)
        {
            auto bond=HOST_BONDS[i];
            BOND_STATE[i]=bond.STATE;
            BOND_TYPE[i]=bond.TYPE;

            auto A1=HOST_POINTS[bond.ID1*4+bond.ID_A1];
            auto A2=HOST_POINTS[bond.ID1*4+bond.ID_A2];
            auto A3=HOST_POINTS[bond.ID1*4+bond.ID_A3];

            auto B1=HOST_POINTS[bond.ID2*4+bond.ID_B1];
            auto B2=HOST_POINTS[bond.ID2*4+bond.ID_B2];
            auto B3=HOST_POINTS[bond.ID2*4+bond.ID_B3];


            auto IA1=HOST_POINTS_INIT[bond.ID1*4+bond.ID_A1];
            auto IA2=HOST_POINTS_INIT[bond.ID1*4+bond.ID_A2];
            auto IA3=HOST_POINTS_INIT[bond.ID1*4+bond.ID_A3];

            auto IB1=HOST_POINTS_INIT[bond.ID2*4+bond.ID_B1];
            auto IB2=HOST_POINTS_INIT[bond.ID2*4+bond.ID_B2];
            auto IB3=HOST_POINTS_INIT[bond.ID2*4+bond.ID_B3];

            REAL4 dis=REAL4({0,0,0,0});
            REAL4 AC=REAL4({0,0,0,0});
            REAL4 BC=REAL4({0,0,0,0});

            REAL4 IAC=REAL4({0,0,0,0});
            REAL4 IBC=REAL4({0,0,0,0});

            for(int k=0;k<3;k++)
            {
                AC[k]=(A1[k]+A2[k]+A3[k])/3.0;
                BC[k]=(B1[k]+B2[k]+B3[k])/3.0;
                IAC[k]=(IA1[k]+IA2[k]+IA3[k])/3.0;
                IBC[k]=(IB1[k]+IB2[k]+IB3[k])/3.0;
                dis[k]=(AC[k]+BC[k])*0.5-(IAC[k]+IBC[k])*0.5;
            }

            displacement[i]=dis;


            cells->InsertNextCell(3);
            cells->InsertCellPoint(pointsB->GetNumberOfPoints());
            cells->InsertCellPoint(pointsB->GetNumberOfPoints()+1);
            cells->InsertCellPoint(pointsB->GetNumberOfPoints()+2);
            pointsB->InsertNextPoint(
                        (A1[0]+B1[0])*0.5,
                        (A1[1]+B1[1])*0.5,
                        (A1[2]+B1[2])*0.5);
            pointsB->InsertNextPoint(
                        (A2[0]+B2[0])*0.5,
                        (A2[1]+B2[1])*0.5,
                        (A2[2]+B2[2])*0.5);
            pointsB->InsertNextPoint(
                        (A3[0]+B3[0])*0.5,
                        (A3[1]+B3[1])*0.5,
                        (A3[2]+B3[2])*0.5);

        }
        output->SetPoints(pointsB);
        output->SetPolys(cells);
        vtkCellData* cbdata =output->GetCellData();
        AddToVTKScalar(cbdata,BOND_STATE,"STATE");
        cbdata->SetActiveScalars("STATE");
        AddToVTKScalar(cbdata,BOND_TYPE,"TYPE");
        AddToVTKVector(cbdata,displacement,"DISPLACEMENT");
        if(file->exists("BONDS_FN"))
        {
            std::vector<REAL4> HOST_BONDS_FN;
            ReadReal4Array(file,HOST_BONDS_FN,BOND_COUNT,"BONDS_FN");
            AddToVTKVector(cbdata,HOST_BONDS_FN,"BONDS_FN");
        }
        if(file->exists("BONDS_FT"))
        {
            std::vector<REAL4> HOST_BONDS_FT;
            ReadReal4Array(file,HOST_BONDS_FT,BOND_COUNT,"BONDS_FT");
            AddToVTKVector(cbdata,HOST_BONDS_FT,"BONDS_FT");
        }


        if(file->exists("BONDS_SIGMA_N"))
        {
            std::vector<double> HOST_BONDS_SIGMA_N;
            ReadRealArray(file,HOST_BONDS_SIGMA_N,BOND_COUNT,"BONDS_SIGMA_N");
            AddToVTKScalar(cbdata,HOST_BONDS_SIGMA_N,"BONDS_SIGMA_N");
        }
        if(file->exists("BONDS_SIGMA_T"))
        {
            std::vector<double> HOST_BONDS_SIGMA_T;
            ReadRealArray(file,HOST_BONDS_SIGMA_T,BOND_COUNT,"BONDS_SIGMA_T");
            AddToVTKScalar(cbdata,HOST_BONDS_SIGMA_T,"BONDS_SIGMA_T");
        }



    }

    delete file;


    output->Modified();
    return 1;
}

void HDF5ReaderTetraBonds::PrintSelf(ostream& os, vtkIndent indent)
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

int HDF5ReaderTetraBonds::RequestInformation(vtkInformation *vtkNotUsed(request), vtkInformationVector **vtkNotUsed(inputVector), vtkInformationVector *outputVector)
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
