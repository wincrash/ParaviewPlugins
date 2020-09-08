#include "McdemReader.h"

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


typedef double REAL;
typedef boost::array<double,4> REAL4;
typedef boost::array<int,4> INT4;
typedef boost::array<int,2> INT2;
typedef  unsigned char UCHAR ;

typedef boost::array<UCHAR,4> UCHAR4;
typedef boost::array<UCHAR,2> UCHAR2;
namespace fs = std::experimental::filesystem;






auto dot=[](REAL4 vect_A, REAL4 vect_B)->double
{
    return vect_A[0]*vect_B[0]+
            vect_A[1]*vect_B[1]+
            vect_A[2]*vect_B[2]+
            vect_A[3]*vect_B[3];
};
auto cross=[](REAL4 vect_A, REAL4 vect_B)->REAL4
{
    REAL4 cross_P=REAL4({0,0,0,0});
    cross_P[0] = vect_A[1] * vect_B[2] - vect_A[2] * vect_B[1];
    cross_P[1] = vect_A[2] * vect_B[0] - vect_A[0] * vect_B[2];
    cross_P[2] = vect_A[0] * vect_B[1] - vect_A[1] * vect_B[0];
    return cross_P;
};
#define CLIP_EPSILON 1E-6
auto CLIP=[](REAL4 *POINTS,REAL4 *cv,REAL4 cuttingPlane,int COUNT)->int
{
    REAL4 v1=POINTS[0];
    REAL dist2;
    REAL dist1=dot(cuttingPlane,v1)-cuttingPlane[3];
    bool inside=(dist1<=CLIP_EPSILON);
    int curv=0;
    REAL4 v2;
    REAL skirtumas;
    REAL d;
    REAL4 t;

    for (int i=1; i<=COUNT; i++)
    {
        v2=POINTS[i % COUNT];
        dist2=dot(cuttingPlane,v2)-cuttingPlane[3];
        skirtumas=dist1-dist2;
        d=dist1/skirtumas;
        for(int f=0;f<4;f++)// ant gpu sito nereikia
            t[f]=v1[f]+(v2[f]-v1[f])*d;

        // Sutherland-hodgeman clipping
        if (inside && (dist2<=CLIP_EPSILON)) cv[curv++]=v2;	// Both in
        else if ((!inside) && (dist2<=CLIP_EPSILON))		// Coming in
        {
            inside=true;
            if(fabs(skirtumas)<CLIP_EPSILON) continue;
            cv[curv++]=t;
            cv[curv++]=v2;
        } else if (inside && (dist2>CLIP_EPSILON))		// Going out
        {
            inside=false;
            if(fabs(skirtumas)<CLIP_EPSILON) continue;
            cv[curv++]=t;
        }		// Both out
        v1=v2;
        dist1=dist2;
    }
    return curv;
};

auto normalize=[](REAL4 A)->REAL4
{
    REAL ilgis=std::sqrt(dot(A,A));
    A[0]=A[0]/ilgis;
    A[1]=A[1]/ilgis;
    A[2]=A[2]/ilgis;
    A[3]=A[3]/ilgis;
    return A;
};

auto CreatePlaneCUT=[](const REAL4 FacePoint,const REAL4 CENTROID)->REAL4
{
    //REAL4 PLANE=normalize(REAL4({FacePoint[0]-CENTROID[0],FacePoint[1]-CENTROID[1],FacePoint[2]-CENTROID[2],0}));
    REAL4 PLANE=normalize(REAL4({FacePoint[0],FacePoint[1],FacePoint[2],0}));
    //PLANE[0]*=-1;
    //PLANE[1]*=-1;
    //PLANE[2]*=-1;
    PLANE[3]=0;
REAL4 FacePoint1=FacePoint;
    FacePoint1[0]+=CENTROID[0];
    FacePoint1[1]+=CENTROID[1];
    FacePoint1[2]+=CENTROID[2];
    PLANE[3]=dot(FacePoint1,PLANE);
    //PrintREAL4("CENTROID",CENTROID);
   // PrintREAL4("FacePoint",FacePoint);
    //PrintREAL4("CreatePlane ",PLANE);
    return PLANE;
};



auto RotatePointQuarterion=[](REAL4 u,REAL4 v)->REAL4
{
    REAL s = u[3];
    u[3]=0;
    REAL dot1=dot(u, v);
    REAL dot2=dot(u, u);
    REAL4 cross1=cross(u, v);
    REAL4 rez=REAL4({0,0,0,0});
    for(int i=0;i<3;i++)
    {
        rez[i]= 2.0*dot1*u[i]+(s*s-dot2)*v[i]+2.0*s*cross1[i];
    }

    return rez;
};





#define MAX_VERTICES_IN_FACE 32
#define SURFACE_MAX_SIZE 1000
auto CreateFace=[](const REAL4 CENTER,const std::vector<REAL4> FACE_POINTS,int faceCount,int faceID,int &COUNT,REAL4 QUART)->std::vector<REAL4>
{
    REAL4 NORMAL=normalize(FACE_POINTS[faceID]);
    REAL scale=-dot(FACE_POINTS[faceID],NORMAL);
    REAL4 POINTS[MAX_VERTICES_IN_FACE];
    REAL4 POINTS_CV[MAX_VERTICES_IN_FACE];
    COUNT=4;
    //create initial surface
    POINTS[0]=REAL4({scale,-SURFACE_MAX_SIZE,-SURFACE_MAX_SIZE,0});
    POINTS[1]=REAL4({scale,SURFACE_MAX_SIZE,-SURFACE_MAX_SIZE,0});
    POINTS[2]=REAL4({scale,SURFACE_MAX_SIZE,SURFACE_MAX_SIZE,0});
    POINTS[3]=REAL4({scale,-SURFACE_MAX_SIZE,SURFACE_MAX_SIZE,0});
    REAL4 XAXIS=REAL4({1,0,0,0});

    //tikrinama kad nebutu ta pati kryptis ir 180 laipsniu jei
    //taip atsitiktu tiesiog pakeiciame vietoj X i Y krypti
    if(NORMAL[0]<-0.99 ||NORMAL[0]>0.99)
    {
        XAXIS=REAL4({0,1,0,0});
        POINTS[0]=REAL4({-SURFACE_MAX_SIZE,scale,-SURFACE_MAX_SIZE,0});
        POINTS[1]=REAL4({SURFACE_MAX_SIZE,scale,-SURFACE_MAX_SIZE,0});
        POINTS[2]=REAL4({SURFACE_MAX_SIZE,scale,SURFACE_MAX_SIZE,0});
        POINTS[3]=REAL4({-SURFACE_MAX_SIZE,scale,SURFACE_MAX_SIZE,0});
    }


    REAL4 q=REAL4({0,0,0,0});
    std::swap(NORMAL,XAXIS);
    q=cross(NORMAL,XAXIS);
    q[3]=1 + dot(NORMAL, XAXIS);
    q=normalize(q);



    POINTS[0]=RotatePointQuarterion(q,POINTS[0]);
    POINTS[1]=RotatePointQuarterion(q,POINTS[1]);
    POINTS[2]=RotatePointQuarterion(q,POINTS[2]);
    POINTS[3]=RotatePointQuarterion(q,POINTS[3]);


  //  POINTS[0]=RotatePointQuarterion(QUART,POINTS[0]);
   // POINTS[1]=RotatePointQuarterion(QUART,POINTS[1]);
    //POINTS[2]=RotatePointQuarterion(QUART,POINTS[2]);
//    POINTS[3]=RotatePointQuarterion(QUART,POINTS[3]);


    for(int i=0;i<3;i++)
    {
        POINTS[0][i]=POINTS[0][i]+CENTER[i];
        POINTS[1][i]=POINTS[1][i]+CENTER[i];
        POINTS[2][i]=POINTS[2][i]+CENTER[i];
        POINTS[3][i]=POINTS[3][i]+CENTER[i];
    }

    bool state=true;
    for(int i=0;i<faceCount;i++)
    {
        if(i==faceID) continue;
        REAL4 cutting_plane=CreatePlaneCUT(FACE_POINTS[i],CENTER);
        if(state)
        {
            COUNT=CLIP(POINTS,POINTS_CV,cutting_plane,COUNT);
            state=false;
        }else
        {
            COUNT=CLIP(POINTS_CV,POINTS,cutting_plane,COUNT);
            state=true;
        }
        if(COUNT<3) break;
    }


    std::vector<REAL4> RET;

    if(state)
        for(int i=0;i<COUNT;i++)
            RET.push_back(POINTS[i]);
    else
        for(int i=0;i<COUNT;i++)
            RET.push_back(POINTS_CV[i]);
    return RET;
};






#include "HDF5Helper.h"
vtkIdType McdemReader::findClosestSolutionIndex(double Time)
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
        for (size_t ii = 0; ii < array.size(); ii++) {
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
        for (size_t ii = 0; ii < array.size(); ii++) {
            data->SetTuple1(ii, array[ii]);
        }
        location->AddArray(data);
        data->Delete();
    }
}


void AddToVTKScalar(vtkFieldData *location,std::vector<UCHAR>&array,std::string name)
{
    if(array.size()>0)
    {
        //vtkSmartPointer<vtkDoubleArray>data = vtkSmartPointer<vtkDoubleArray>::New();
        vtkIntArray*data=vtkIntArray::New();
        data->SetNumberOfTuples(array.size());
        data->SetNumberOfComponents(1);
        data->SetName(name.c_str());
        data->SetNumberOfValues(array.size());
        for (size_t ii = 0; ii < array.size(); ii++) {
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
        for (size_t ii = 0; ii < array.size(); ii++) {
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
        for (size_t ii = 0; ii < array.size(); ii++) {
            REAL4 t=array[ii];
            data->InsertNextTuple3( t[0],t[1],t[2]);
        }
        location->AddArray(data);
        data->Delete();
    }
}



vtkStandardNewMacro(McdemReader);


McdemReader::McdemReader()
{
    this->FileName = nullptr;
    this->DirectoryName = nullptr;
    this->SetNumberOfInputPorts(0);
    this->SetNumberOfOutputPorts(1);
    this->times = vtkSmartPointer<vtkDoubleArray>::New();
    this->timesNames = vtkSmartPointer<vtkVariantArray>::New();
}

int McdemReader::RequestData(
        vtkInformation *vtkNotUsed(request),
        vtkInformationVector **vtkNotUsed(inputVector),
        vtkInformationVector *outputVector)
{

    vtkInformation *outInfo = outputVector->GetInformationObject(0);
    double requestedTime = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());
    vtkUnstructuredGrid *output = vtkUnstructuredGrid::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
    output->Allocate(1000,1000);

    std::string tempFilename=filenames[findClosestSolutionIndex(requestedTime)];

    std::cout<<"McdemReader "<<tempFilename<<"\n";



    //For displacement
    h5cpp::File file1(filenames[findClosestSolutionIndex(0)], "r");
    auto gr0 = file1.root();
    size_t PARTICLE_COUNT=gr0.attrs().get<int>("PARTICLE_COUNT");
    std::vector<REAL4> OLD_POSITION(PARTICLE_COUNT,REAL4({0,0,0,0}));
    readDataSetREAL4(gr0,OLD_POSITION,"POSITION",PARTICLE_COUNT);
    file1.close();


    h5cpp::File file(tempFilename, "r");
    auto gr=file.root();
    std::vector<double> STEP;
    STEP.push_back(gr.attrs().get<int>("STEP"));
    std::vector<double> TIME;
    TIME.push_back(gr.attrs().get<double>("TIME"));

    std::vector<double> MAX_FACE_COUNT;
    MAX_FACE_COUNT.push_back(gr.attrs().get<int>("MAX_FACE_COUNT"));
    size_t MAX_POSSIBLE_FACES=gr.attrs().get<int>("MAX_POSSIBLE_FACES");
    size_t BOND_COUNT=gr.attrs().get<int>("BOND_COUNT");

    std::vector<UCHAR> FACE_COUNT(PARTICLE_COUNT,0);//UCHAR
    std::vector<REAL4> FACES(PARTICLE_COUNT*MAX_POSSIBLE_FACES,REAL4({0,0,0,0}));//REAL4
    std::vector<REAL4> POSITION(PARTICLE_COUNT,REAL4({0,0,0,0}));//REAL4
    std::vector<REAL4> VELOCITY(PARTICLE_COUNT,REAL4({0,0,0,0}));//REAL4
    std::vector<REAL4> ANGULAR_VELOCITY(PARTICLE_COUNT,REAL4({0,0,0,0}));//REAL4
    std::vector<REAL4> QUARTERION(PARTICLE_COUNT,REAL4({0,0,0,0}));//REAL4
    std::vector<REAL>  INERTIA_TENSOR(PARTICLE_COUNT*9,0.0);//REAL*9
    std::vector<UCHAR> FIX(PARTICLE_COUNT,0);//UCHAR
    std::vector<UCHAR> MATERIAL(PARTICLE_COUNT,0);//UCHAR
    std::vector<UCHAR> WALL(PARTICLE_COUNT,0);//UCHAR
    std::vector<REAL>  RADIUS(PARTICLE_COUNT,0);//REAL


    readDataSetUCHAR(gr,FACE_COUNT,"FACE_COUNT",PARTICLE_COUNT);
    readDataSetREAL4(gr,FACES,"FACES",PARTICLE_COUNT*MAX_POSSIBLE_FACES);
    readDataSetREAL4(gr,POSITION,"POSITION",PARTICLE_COUNT);
    readDataSetREAL4(gr,VELOCITY,"VELOCITY",PARTICLE_COUNT);
    readDataSetREAL4(gr,ANGULAR_VELOCITY,"ANGULAR_VELOCITY",PARTICLE_COUNT);
    readDataSetREAL4(gr,QUARTERION,"QUARTERION",PARTICLE_COUNT);
    readDataSetUCHAR(gr,FIX,"FIX",PARTICLE_COUNT);
    readDataSetUCHAR(gr,MATERIAL,"MATERIAL",PARTICLE_COUNT);
    readDataSetUCHAR(gr,WALL,"WALL",PARTICLE_COUNT);
    readDataSetREAL(gr,RADIUS,"RADIUS",PARTICLE_COUNT);




    vtkPoints*points=vtkPoints::New();
    points->SetDataTypeToDouble();





    for(size_t idx=0;idx<PARTICLE_COUNT;idx++)
    {
        auto QUART=QUARTERION[idx];
       // QUART=REAL4({0.1,0,0, 0.5});
        QUART=normalize(QUART);
        //std::cout<<"QUART "<<QUART[0]<<" "<<QUART[1]<<" "<<QUART[2]<<" "<<QUART[3]<<"\n";
        auto face_count=FACE_COUNT[idx];
        std::vector<REAL4> F;
        F.resize(face_count);
        for(size_t i=0;i<face_count;i++)
        {
            F[i]=RotatePointQuarterion(QUART,FACES[idx*MAX_POSSIBLE_FACES+i]);
        }
        REAL4 CENTER=POSITION[idx];


        std::vector<vtkIdType> pointsIDS;
        vtkSmartPointer<vtkCellArray> dodechedronFaces =vtkSmartPointer<vtkCellArray>::New();

        for(int i=0;i<face_count;i++)
        {
            int COUNT=0;
            std::vector<REAL4> TASKAI=CreateFace(CENTER,F,F.size(),i,COUNT,QUART);
            if(COUNT<3)continue;
            dodechedronFaces->InsertNextCell(COUNT);
            for(size_t z=0;z<COUNT;z++)
            {
                pointsIDS.push_back(points->GetNumberOfPoints());
                dodechedronFaces->InsertCellPoint(points->GetNumberOfPoints());
                points->InsertNextPoint(TASKAI[z][0],TASKAI[z][1],TASKAI[z][2]);

            }
        }
        output->InsertNextCell(VTK_POLYHEDRON,
                               pointsIDS.size(), pointsIDS.data(),
                               face_count, dodechedronFaces->GetPointer());
    }
    output->SetPoints(points);
    points->Delete();
    vtkCellData* cdata =output->GetCellData();
    vtkFieldData* fdata = output->GetFieldData();


    AddToVTKScalar(fdata,STEP,"STEP");
    AddToVTKScalar(fdata,TIME,"TIME");
    AddToVTKScalar(fdata,MAX_FACE_COUNT,"MAX_FACE_COUNT");


    AddToVTKVector(cdata,VELOCITY,"VELOCITY");
    AddToVTKVector(cdata,ANGULAR_VELOCITY,"ANGULAR_VELOCITY");
    AddToVTKScalar(cdata,FIX,"FIX");
    AddToVTKScalar(cdata,MATERIAL,"MATERIAL");
    AddToVTKScalar(cdata,WALL,"WALL");
    AddToVTKScalar(cdata,RADIUS,"RADIUS");




    if(gr.exists("FORCE"))
    {
        std::vector<REAL4> FORCE;
        FORCE.resize(PARTICLE_COUNT);
        readDataSetREAL4(gr,FORCE,"FORCE",PARTICLE_COUNT);
        AddToVTKVector(cdata,FORCE,"FORCE");
    }
    if(gr.exists("TORQUE"))
    {
        std::vector<REAL4> TORQUE;
        TORQUE.resize(PARTICLE_COUNT);
        readDataSetREAL4(gr,TORQUE,"TORQUE",PARTICLE_COUNT);
        AddToVTKVector(cdata,TORQUE,"TORQUE");
    }

    if(gr.exists("NN_COUNT"))
    {
        std::vector<int> NN_COUNT;
        NN_COUNT.resize(PARTICLE_COUNT);
        readDataSetINT(gr,NN_COUNT,"NN_COUNT",PARTICLE_COUNT);
        AddToVTKScalar(cdata,NN_COUNT,"NN_COUNT");
    }
    if(gr.exists("PP_COUNT"))
    {
        std::vector<int> PP_COUNT;
        PP_COUNT.resize(PARTICLE_COUNT);
        readDataSetINT(gr,PP_COUNT,"PP_COUNT",PARTICLE_COUNT);
        AddToVTKScalar(cdata,PP_COUNT,"PP_COUNT");
    }
    if(gr.exists("MASS"))
    {
        std::vector<REAL> MASS;
        MASS.resize(PARTICLE_COUNT);
        readDataSetREAL(gr,MASS,"MASS",PARTICLE_COUNT);
        AddToVTKScalar(cdata,MASS,"MASS");
    }

    if(gr.exists("VOLUME"))
    {
        std::vector<REAL> VOLUME;
        VOLUME.resize(PARTICLE_COUNT);
        readDataSetREAL(gr,VOLUME,"VOLUME",PARTICLE_COUNT);
        AddToVTKScalar(cdata,VOLUME,"VOLUME");
    }


    if(gr.exists("PARTICLE_ID"))
    {
        std::vector<int> PARTICLE_ID;
        PARTICLE_ID.resize(PARTICLE_COUNT);
        readDataSetINT(gr,PARTICLE_ID,"PARTICLE_ID",PARTICLE_COUNT);
        AddToVTKScalar(cdata,PARTICLE_ID,"PARTICLE_ID");
    }
    if(gr.exists("PARTICLE_CPU"))
    {
        std::vector<int> PARTICLE_CPU;
        PARTICLE_CPU.resize(PARTICLE_COUNT);
        readDataSetINT(gr,PARTICLE_CPU,"PARTICLE_CPU",PARTICLE_COUNT);
        AddToVTKScalar(cdata,PARTICLE_CPU,"PARTICLE_CPU");
    }
    cdata->SetActiveVectors("VELOCITY");
    cdata->SetActiveScalars("FIX");




    output->Modified();
    return 1;
}

void McdemReader::PrintSelf(ostream& os, vtkIndent indent)
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

int McdemReader::RequestInformation(vtkInformation *vtkNotUsed(request), vtkInformationVector **vtkNotUsed(inputVector), vtkInformationVector *outputVector)
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
        if (p.path().string().find(".mcdem") != std::string::npos)
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


    delete[] trange;

    return 1;
}
