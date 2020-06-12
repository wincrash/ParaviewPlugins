#include "HDF5ReaderConvex.h"

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


#define MAX_FACE_COUNT 6
#define MAX_VERTEX_COUNT_IN_FACE 8
#define MAX_VERTICES 12
typedef struct
{
    REAL4 POINTS[MAX_VERTICES];
} CONVEX_POINTS_TYPE;

typedef struct
{
    unsigned char FACE_VERTEX_COUNT[MAX_FACE_COUNT];
    unsigned char IDS[MAX_FACE_COUNT][MAX_VERTEX_COUNT_IN_FACE];
} FACE_INFO_TYPE;


typedef struct
{
    unsigned char VERTEX_COUNT;
    unsigned char FACE_COUNT;
    unsigned char FIX;
    unsigned char MATERIAL;
} PARTICLE_INFO_TYPE;


#define MAX_BOND_FACE_VERTICES 10
typedef struct
{
    int ID1;
    int ID2;
    unsigned char STATE;
    unsigned char TYPE;
    unsigned char COUNT;
    unsigned char PADDING;
    unsigned char POINTS_A[MAX_BOND_FACE_VERTICES];
    unsigned char POINTS_B[MAX_BOND_FACE_VERTICES];
}BONDS_TYPE;

#include <hdf5.h>
#include "H5Cpp.h"


auto ReadBondsInfo=[](H5::H5File* file,std::vector<BONDS_TYPE> &data,size_t COUNT)->void
{
    data.resize(COUNT);
    H5::DataSet* dataset = new H5::DataSet (file->openDataSet( "BONDS" ));
    H5::CompType mtype1( sizeof(BONDS_TYPE));
    hsize_t dim2[] = {MAX_BOND_FACE_VERTICES};
    auto array2_tid = H5Tarray_create(H5T_NATIVE_UCHAR, 1,dim2);
    mtype1.insertMember( "ID1", HOFFSET(BONDS_TYPE, ID1), H5::PredType::NATIVE_INT);
    mtype1.insertMember( "ID2", HOFFSET(BONDS_TYPE, ID2), H5::PredType::NATIVE_INT);
    mtype1.insertMember( "STATE", HOFFSET(BONDS_TYPE, STATE), H5::PredType::NATIVE_UCHAR);
    mtype1.insertMember( "TYPE", HOFFSET(BONDS_TYPE, TYPE), H5::PredType::NATIVE_UCHAR);
    mtype1.insertMember( "COUNT", HOFFSET(BONDS_TYPE, COUNT), H5::PredType::NATIVE_UCHAR);
    mtype1.insertMember( "PADDING", HOFFSET(BONDS_TYPE, PADDING), H5::PredType::NATIVE_UCHAR);
    mtype1.insertMember( "POINTS_A", HOFFSET(BONDS_TYPE, POINTS_A), array2_tid);
    mtype1.insertMember( "POINTS_B", HOFFSET(BONDS_TYPE, POINTS_B), array2_tid);
    dataset->read( &data[0], mtype1 );
    delete dataset;
};

auto ReadConvexPointsInfo=[](H5::H5File* file,std::vector<CONVEX_POINTS_TYPE> &data,size_t COUNT)->void
{
    data.resize(COUNT);
    H5::DataSet* dataset = new H5::DataSet (file->openDataSet( "CONVEX_POINTS" ));
    H5::CompType mtype1( sizeof(CONVEX_POINTS_TYPE));
    hsize_t dim2[] = {MAX_VERTICES,4};
    auto array2_tid = H5Tarray_create(H5T_NATIVE_DOUBLE, 2,dim2);
    mtype1.insertMember( "POINTS", HOFFSET(CONVEX_POINTS_TYPE, POINTS), array2_tid);
    dataset->read( &data[0], mtype1 );
    delete dataset;
};



auto ReadFaceInfo=[](H5::H5File* file,std::vector<FACE_INFO_TYPE> &data,size_t COUNT)->void
{
    data.resize(COUNT);
    H5::DataSet* dataset = new H5::DataSet (file->openDataSet( "FACE_INFO" ));
    H5::CompType mtype1( sizeof(FACE_INFO_TYPE) );
    hsize_t dim1[] = {MAX_FACE_COUNT};
    auto array1_tid = H5Tarray_create(H5T_NATIVE_UCHAR, 1,dim1);
    hsize_t dim2[] = {MAX_FACE_COUNT,MAX_VERTEX_COUNT_IN_FACE};
    auto array2_tid = H5Tarray_create(H5T_NATIVE_UCHAR, 2,dim2);
    mtype1.insertMember( "FACE_VERTEX_COUNT", HOFFSET(FACE_INFO_TYPE, FACE_VERTEX_COUNT), array1_tid);
    mtype1.insertMember( "IDS", HOFFSET(FACE_INFO_TYPE, IDS), array2_tid);
    dataset->read( &data[0], mtype1 );
    delete dataset;
};



auto ReadParticleInfo=[](H5::H5File* file,std::vector<PARTICLE_INFO_TYPE> &data,size_t COUNT)->void
{
    data.resize(COUNT);
    H5::DataSet* dataset = new H5::DataSet (file->openDataSet( "PARTICLE_INFO" ));
    H5::CompType mtype1( sizeof(PARTICLE_INFO_TYPE) );
    mtype1.insertMember( "VERTEX_COUNT", HOFFSET(PARTICLE_INFO_TYPE, VERTEX_COUNT), H5::PredType::NATIVE_UCHAR);
    mtype1.insertMember( "FACE_COUNT", HOFFSET(PARTICLE_INFO_TYPE, FACE_COUNT), H5::PredType::NATIVE_UCHAR);
    mtype1.insertMember( "FIX", HOFFSET(PARTICLE_INFO_TYPE, FIX), H5::PredType::NATIVE_UCHAR);
    mtype1.insertMember( "MATERIAL", HOFFSET(PARTICLE_INFO_TYPE, MATERIAL), H5::PredType::NATIVE_UCHAR);
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



vtkStandardNewMacro(HDF5ReaderConvex);


HDF5ReaderConvex::HDF5ReaderConvex()
{
    this->FileName = NULL;
    this->DirectoryName = NULL;
    this->SetNumberOfInputPorts(0);
    this->SetNumberOfOutputPorts(2);
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
    vtkUnstructuredGrid *output = vtkUnstructuredGrid::SafeDownCast(
                outInfo->Get(vtkDataObject::DATA_OBJECT()));


    vtkInformation *outInfo1 = outputVector->GetInformationObject(1);
    double requestedTime1 = outInfo1->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());
    vtkUnstructuredGrid *BONDAI = vtkUnstructuredGrid::SafeDownCast(
                outInfo1->Get(vtkDataObject::DATA_OBJECT()));

    std::string tempFilename=filenames[findClosestSolutionIndex(requestedTime)];

    std::cout<<tempFilename<<"\n";

    std::vector<CONVEX_POINTS_TYPE> HOST_CONVEX_POINTS;
    std::vector<FACE_INFO_TYPE> HOST_FACE_INFO;
    std::vector<PARTICLE_INFO_TYPE> HOST_PARTICLE_INFO;

    H5::H5File* file= new H5::H5File( tempFilename, H5F_ACC_RDONLY );
    size_t PARTICLE_COUNT=ReadIntAttribute(file,"PARTICLE_COUNT");
    size_t BOND_COUNT=ReadIntAttribute(file,"BOND_COUNT");
    size_t TIME_STEP=ReadIntAttribute(file,"STEP");
    double TIME=ReadDoubleAttribute(file,"TIME");

    ReadConvexPointsInfo(file,HOST_CONVEX_POINTS,PARTICLE_COUNT);
    ReadParticleInfo(file,HOST_PARTICLE_INFO,PARTICLE_COUNT);
    ReadFaceInfo(file,HOST_FACE_INFO,PARTICLE_COUNT);

    std::vector<int> FIX(PARTICLE_COUNT,0);
    std::vector<int> MATERIAL(PARTICLE_COUNT,0);
    vtkPoints *points = vtkPoints::New();
    points->SetNumberOfPoints(0);
    points->SetDataTypeToDouble();
    output->Allocate(1000,1000);

    for(size_t i=0;i<PARTICLE_COUNT;i++)
    {
        CONVEX_POINTS_TYPE CPOINTS=HOST_CONVEX_POINTS[i];
        FACE_INFO_TYPE CFACE=HOST_FACE_INFO[i];
        PARTICLE_INFO_TYPE CINFO=HOST_PARTICLE_INFO[i];
        FIX[i]=CINFO.FIX;
        MATERIAL[i]=CINFO.MATERIAL;
        int pointStart=points->GetNumberOfPoints();
        vtkIdType dodechedronPointsIds[CINFO.VERTEX_COUNT];
        for(size_t z=0;z<CINFO.VERTEX_COUNT;z++)
        {
            dodechedronPointsIds[z]=points->GetNumberOfPoints();
            points->InsertNextPoint(CPOINTS.POINTS[z][0],CPOINTS.POINTS[z][1],CPOINTS.POINTS[z][2]);
        }

        vtkSmartPointer<vtkCellArray> dodechedronFaces =vtkSmartPointer<vtkCellArray>::New();

        for (size_t z = 0; z < CINFO.FACE_COUNT; z++)
        {
            dodechedronFaces->InsertNextCell(CFACE.FACE_VERTEX_COUNT[z]);
            for(size_t k=0;k<CFACE.FACE_VERTEX_COUNT[z];k++)
                dodechedronFaces->InsertCellPoint(CFACE.IDS[z][k]+pointStart);
        }

        output->InsertNextCell(VTK_POLYHEDRON,
                               CINFO.VERTEX_COUNT, dodechedronPointsIds,
                               CINFO.FACE_COUNT, dodechedronFaces->GetPointer());
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
        std::vector<int> HOST_MASS;
        ReadIntArray(file,HOST_MASS,PARTICLE_COUNT,"MASS");
        AddToVTKScalar(cdata,HOST_MASS,"MASS");
    }
    if(file->exists("VOLUME"))
    {
        std::vector<int> HOST_VOLUME;
        ReadIntArray(file,HOST_VOLUME,PARTICLE_COUNT,"VOLUME");
        AddToVTKScalar(cdata,HOST_VOLUME,"VOLUME");
    }
    cdata->SetActiveVectors("VELOCITY");
    cdata->SetActiveScalars("FIX");






    if(BOND_COUNT>0)
    {
        std::vector<BONDS_TYPE> HOST_BONDS;
        ReadBondsInfo(file,HOST_BONDS,BOND_COUNT);
        vtkPoints *pointsB = vtkPoints::New();
        pointsB->SetNumberOfPoints(0);
        pointsB->SetDataTypeToDouble();
        BONDAI->Allocate(1000,1000);
        std::vector<int> BOND_STATE(BOND_COUNT,0);
        std::vector<int> BOND_TYPE(BOND_COUNT,0);
        for(size_t i=0;i<BOND_COUNT;i++)
        {
            auto bond=HOST_BONDS[i];
            BOND_STATE[i]=bond.STATE;
            BOND_TYPE[i]=bond.TYPE;
            auto P1=HOST_CONVEX_POINTS[bond.ID1];
            auto P2=HOST_CONVEX_POINTS[bond.ID2];
            vtkIdType PP[bond.COUNT];

            for(size_t k=0;k<bond.COUNT;k++)
            {
                PP[k]=pointsB->GetNumberOfPoints();
                pointsB->InsertNextPoint(
                            (P1.POINTS[bond.POINTS_A[k]][0]+P2.POINTS[bond.POINTS_B[k]][0])*0.5,
                            (P1.POINTS[bond.POINTS_A[k]][1]+P2.POINTS[bond.POINTS_B[k]][1])*0.5,
                            (P1.POINTS[bond.POINTS_A[k]][2]+P2.POINTS[bond.POINTS_B[k]][2])*0.5);

            }
            BONDAI->InsertNextCell(VTK_POLYGON,bond.COUNT,PP);


        }
         BONDAI->SetPoints(pointsB);
         vtkCellData* cbdata =BONDAI->GetCellData();
         AddToVTKScalar(cbdata,BOND_STATE,"STATE");
         cbdata->SetActiveScalars("STATE");
         AddToVTKScalar(cbdata,BOND_TYPE,"TYPE");
         if(file->exists("BONDS_FN"))
         {
             std::vector<REAL4> HOST_BONDS_FN;
             ReadReal4Array(file,HOST_BONDS_FN,BOND_COUNT,"BONDS_FN");
             for(int zz=0;zz<HOST_BONDS_FN.size();zz++)
             {
                 std::cout<<"zz "<<zz<<"  "<<HOST_BONDS_FN[zz][0]<<"  "<<HOST_BONDS_FN[zz][
                            1]<<"  "<<HOST_BONDS_FN[zz][2]<<"  "<<HOST_BONDS_FN[zz][3]<<"\n";
             }
             AddToVTKVector(cbdata,HOST_BONDS_FN,"BONDS_FN");
         }
         if(file->exists("BONDS_FT"))
         {
             std::vector<REAL4> HOST_BONDS_FT;
             ReadReal4Array(file,HOST_BONDS_FT,BOND_COUNT,"BONDS_FT");
             AddToVTKVector(cbdata,HOST_BONDS_FT,"BONDS_FT");
         }



    }



    delete file;


    BONDAI->Modified();
    output->Modified();
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
