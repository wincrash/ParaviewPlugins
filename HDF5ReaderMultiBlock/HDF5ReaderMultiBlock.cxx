#include "HDF5ReaderMultiBlock.h"

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
#include <vtkMultiBlockDataSet.h>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/lexical_cast.hpp>
#include <vtkPolyData.h>
#include <vtkPolyDataAlgorithm.h>
namespace fs = std::experimental::filesystem;
vtkIdType HDF5ReaderMultiBlock::findClosestSolutionIndex(double Time)
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

void AddToVTKVector(vtkFieldData *location,std::vector<double>&array,std::string name)
{
    if(array.size()>0)
    {
        //vtkSmartPointer<vtkDoubleArray>data = vtkSmartPointer<vtkDoubleArray>::New();
        vtkDoubleArray*data=vtkDoubleArray::New();
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








vtkStandardNewMacro(HDF5ReaderMultiBlock);



HDF5ReaderMultiBlock::HDF5ReaderMultiBlock()
{
    this->FileName = NULL;
    this->DirectoryName = NULL;
    this->SetNumberOfInputPorts(0);
    this->SetNumberOfOutputPorts(2);
    this->times = vtkSmartPointer<vtkDoubleArray>::New();
    this->timesNames = vtkSmartPointer<vtkVariantArray>::New();
}

int HDF5ReaderMultiBlock::RequestData(
        vtkInformation *vtkNotUsed(request),
        vtkInformationVector **vtkNotUsed(inputVector),
        vtkInformationVector *outputVector)
{

    vtkInformation *outInfo = outputVector->GetInformationObject(0);
    vtkInformation *outInfo1 = outputVector->GetInformationObject(1);


    double requestedTime = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());
    std::cout<<"  "<<findClosestSolutionIndex(requestedTime)<<"\n";

    double requestedTime1 = outInfo1->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());
    std::cout<<"  "<<findClosestSolutionIndex(requestedTime1)<<"\n";


    vtkPolyData *output = vtkPolyData::SafeDownCast(
                outInfo->Get(vtkDataObject::DATA_OBJECT()));
    vtkPolyData *outputMesh= vtkPolyData::SafeDownCast(
                outInfo1->Get(vtkDataObject::DATA_OBJECT()));
    outputMesh->Modified();

    // Here is where you would read the data from the file. In this example,
    // we simply create a point.
    std::cout<<"iddd "<<findClosestSolutionIndex(requestedTime)<<"\n";
    std::string tempFilename=filenames[findClosestSolutionIndex(requestedTime)];
    std::cout<<"filename "<<tempFilename<<"\n";

    //vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();


    //For displacement
    h5cpp::File file1(filenames[findClosestSolutionIndex(0)], "r");
    auto gr0 = file1.root();
    std::vector<double>oldPosition=GetData(gr0,"POSITIONS",4);
    file1.close();





    h5cpp::File file(tempFilename, "r");
    auto gr = file.root();

    h5cpp::iterator it=gr.begin();
    std::vector<std::string> pavadinimai;
    std::vector<int> rankas;
    while(it!=gr.end())
    {
        if(it.dereference().compare("POSITIONS")==0)
        {
            it++;
            continue;
        }
        if(it.dereference().compare("UNIQUE_RADIUS")==0)
        {
            it++;
            continue;
        }
        pavadinimai.push_back(it.dereference());

        auto dt = gr.open_dataset(it.dereference());
        rankas.push_back(dt.get_dataspace().get_rank());

        it++;
    }

    std::vector<double> STEP;
    STEP.push_back(gr.attrs().get<int>("STEP"));
    std::vector<double> TIME;    
    TIME.push_back(gr.attrs().get<double>("TIME"));
    std::vector<double> TIME_TEMPERATURE;
    if(gr.attrs().exists("TIME_TEMPERATURE"))
    {
        TIME_TEMPERATURE.push_back(gr.attrs().get<double>("TIME_TEMPERATURE"));
    }
    //vtkPolyData*output = vtkPolyData::New();
    vtkPoints*points = vtkPoints::New();
    std::vector<double> POSITIONS=GetData(gr,"POSITIONS",4);
    std::vector<double> RADIUSAS=GetData(gr,"RADIUS",1);
    points->SetNumberOfPoints(POSITIONS.size()/4);
    std::vector<double> DEFORM;
    DEFORM.resize(POSITIONS.size());
    for(int i=0;i<POSITIONS.size()/4;i++)
    {
        int idd=4*i;
        points->SetPoint(i,POSITIONS[idd],POSITIONS[idd+1],POSITIONS[idd+2]);
        DEFORM[idd]=-oldPosition[idd]+POSITIONS[idd];
        DEFORM[idd+1]=-oldPosition[idd+1]+POSITIONS[idd+1];
        DEFORM[idd+2]=-oldPosition[idd+2]+POSITIONS[idd+2];
        DEFORM[idd+3]=0;
    }

    output->SetPoints(points);
    vtkPointData* pdata = output->GetPointData();
    vtkCellData* cdata = output->GetCellData();
    vtkFieldData* fdata = output->GetFieldData();

    AddToVTKVector(pdata,DEFORM,"DISPLACEMENT");


    AddToVTKScalar(fdata,STEP,"STEP");
    AddToVTKScalar(fdata,TIME,"TIME");
    AddToVTKScalar(fdata,TIME_TEMPERATURE,"TIME_TEMPERATURE");




    std::vector<double> data;
    data=GetData(gr,"UNIQUE_RADIUS",1);
    AddToVTKScalar(fdata,data,"UNIQUE_RADIUS");


    for(int i=0;i<pavadinimai.size();i++)
    {
        try {
            if (boost::starts_with(pavadinimai[i], "BOUNDARY"))
                continue;
            //foo_value = boost::lexical_cast<int>(argv[1]+6);
        } catch (boost::bad_lexical_cast) {
            // bad parameter
        }
        try {
            if (boost::starts_with(pavadinimai[i], "BOND_"))
                continue;
            //foo_value = boost::lexical_cast<int>(argv[1]+6);
        } catch (boost::bad_lexical_cast) {
            // bad parameter
        }




        if(rankas[i]==2)
        {
            data=GetData(gr,pavadinimai[i],4);
            if(data.size()==POSITIONS.size())
                AddToVTKVector(pdata,data,pavadinimai[i]);
            else
                AddToVTKVector(cdata,data,pavadinimai[i]);
            std::cout<<"Added vector array "<<pavadinimai[i]<<" \n";
        }
        if(rankas[i]==1)
        {
            data=GetData(gr,pavadinimai[i],1);
            if(data.size()==(POSITIONS.size()/4))
                AddToVTKScalar(pdata,data,pavadinimai[i]);
            else
                AddToVTKScalar(cdata,data,pavadinimai[i]);

            std::cout<<"Added scalar array "<<pavadinimai[i]<<" \n";
        }
    }

    if(gr.exists("BOND_PARTICLE_1") &&gr.exists("BOND_PARTICLE_2")){

        std::vector<double> particle1=GetData(gr,"BOND_PARTICLE_1",1);
        std::vector<double> particle2=GetData(gr,"BOND_PARTICLE_2",1);

        for(int i=0;i<pavadinimai.size();i++)
        {

            try {
                if (!boost::starts_with(pavadinimai[i], "BOND_"))
                    continue;
                //foo_value = boost::lexical_cast<int>(argv[1]+6);
            } catch (boost::bad_lexical_cast) {
                // bad parameter
            }

            try {
                if (boost::starts_with(pavadinimai[i], "BOND_PARTICLE_1"))
                    continue;
                //foo_value = boost::lexical_cast<int>(argv[1]+6);
            } catch (boost::bad_lexical_cast) {
                // bad parameter
            }
            try {
                if (boost::starts_with(pavadinimai[i], "BOND_PARTICLE_2"))
                    continue;
                //foo_value = boost::lexical_cast<int>(argv[1]+6);
            } catch (boost::bad_lexical_cast) {
                // bad parameter
            }


            if(rankas[i]==2)
            {
                data=GetData(gr,pavadinimai[i],4);
                if(data.size()==particle1.size()*4)
                    AddToVTKVector(cdata,data,pavadinimai[i]);
                std::cout<<"Added vector array "<<pavadinimai[i]<<" \n";
            }
            if(rankas[i]==1)
            {
                data=GetData(gr,pavadinimai[i],1);
                if(data.size()==particle1.size())
                    AddToVTKScalar(cdata,data,pavadinimai[i]);

                std::cout<<"Added scalar array "<<pavadinimai[i]<<" \n";
            }
        }




       /*

        data=GetData(gr,"BOND_STATE",1);
        AddToVTKScalar(cdata,data,"BOND_STATE");

        data=GetData(gr,"BOND_F_LIMIT_N",1);
        AddToVTKScalar(cdata,data,"BOND_F_LIMIT_N");

        data=GetData(gr,"BOND_F_LIMIT_T",1);
        AddToVTKScalar(cdata,data,"BOND_F_LIMIT_T");

        data=GetData(gr,"BOND_COUNT",1);
        AddToVTKScalar(pdata,data,"BOND_COUNT");

        data=GetData(gr,"BOND_TYPE",1);
        AddToVTKScalar(cdata,data,"BOND_TYPE");
        data=GetData(gr,"BOND_AREA",1);
        AddToVTKScalar(cdata,data,"BOND_AREA");
*/


        vtkCellArray*lines=vtkCellArray::New();
        std::vector<int> CellsInParticles;
        CellsInParticles.resize(points->GetNumberOfPoints(),0);
        vtkDoubleArray*lineLengthProc=vtkDoubleArray::New();
        lineLengthProc->SetName("LENGTH_PROC");
        lineLengthProc->SetNumberOfComponents(1);
        lineLengthProc->SetNumberOfTuples(particle1.size());


        for(int i=0;i<particle1.size();i++)
        {
            lines->InsertNextCell(2);
            int id1=particle1[i];
            int id2=particle2[i];


            //std::vector<double> POSITIONS=GetData(gr,"POSITIONS",4);
            //std::vector<double> RADIUSAS=GetData(gr,"RADIUS",1);
            double radius1=RADIUSAS[id1];
            double radius2=RADIUSAS[id2];
            double ilgis=std::sqrt((POSITIONS[id1*4+0]-POSITIONS[id2*4+0])*(POSITIONS[id1*4+0]-POSITIONS[id2*4+0])+
                    (POSITIONS[id1*4+1]-POSITIONS[id2*4+1])*(POSITIONS[id1*4+1]-POSITIONS[id2*4+1])+
                    (POSITIONS[id1*4+2]-POSITIONS[id2*4+2])*(POSITIONS[id1*4+2]-POSITIONS[id2*4+2]));
            lineLengthProc->SetTuple1(i,ilgis/(radius1+radius2));
            lines->InsertCellPoint(id1);
            lines->InsertCellPoint(id2);
            CellsInParticles[id1]++;
            CellsInParticles[id2]++;
        }
        cdata->AddArray(lineLengthProc);

        for(int i=0;i<CellsInParticles.size();i++)
        {
            if(CellsInParticles[i]==0)
            {
                lines->InsertNextCell(2);
                lines->InsertCellPoint(i);
                lines->InsertCellPoint(i);
                for(int k=0;k<output->GetCellData()->GetNumberOfArrays();k++)
                {
                    if(output->GetCellData()->GetArray(k)->GetNumberOfComponents()==1)
                    output->GetCellData()->GetArray(k)->InsertNextTuple1(0);
                    if(output->GetCellData()->GetArray(k)->GetNumberOfComponents()==3)
                    output->GetCellData()->GetArray(k)->InsertNextTuple3(0,0,0);
                }

            }
        }
        output->SetLines(lines);
        /*output->SetLines(lines);
        vtkCellArray*verts=vtkCellArray::New();

        output->SetVerts(verts);*/
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
    if(output->GetPointData()->HasArray("RADIUS"))
    {
        output->GetPointData()->SetActiveScalars("RADIUS");
    }
    if(output->GetPointData()->HasArray("VELOCITY"))
    {
        output->GetPointData()->SetActiveVectors("VELOCITY");
    }
    points->Delete();

    if(gr.exists("BOUNDARY_IDS") &&gr.exists("BOUNDARY_POINTS")){

        vtkPoints*bpoints = vtkPoints::New();
        std::vector<double> BPOSITIONS=GetData(gr,"BOUNDARY_POINTS",4);
        for(int i=0;i<BPOSITIONS.size()/4;i++)
        {
            bpoints->InsertNextPoint(BPOSITIONS[i*4],BPOSITIONS[i*4+1],BPOSITIONS[i*4+2]);
            std::cout<<BPOSITIONS[i*4]<<" "<<BPOSITIONS[i*4+1]<<" "<<BPOSITIONS[i*4+2]<<" \n";
        }
        std::vector<double> IDS=GetData(gr,"BOUNDARY_IDS",4);
        vtkCellArray* bcells=vtkCellArray::New();


        for(int i=0;i<IDS.size()/4;i++)
        {
           /* if(IDS[i*4+3]<1000000)
            {
                bcells->InsertNextCell(4);
                bcells->InsertCellPoint(IDS[i*4+0]);
                bcells->InsertCellPoint(IDS[i*4+1]);
                bcells->InsertCellPoint(IDS[i*4+2]);
                bcells->InsertCellPoint(IDS[i*4+3]);
            }
            else*/
            {
                bcells->InsertNextCell(3);
                bcells->InsertCellPoint(IDS[i*4+0]);
                bcells->InsertCellPoint(IDS[i*4+1]);
                bcells->InsertCellPoint(IDS[i*4+2]);
            }
        }

        outputMesh->SetPoints(bpoints);
        outputMesh->SetPolys(bcells);
        bpoints->Delete();
        bcells->Delete();

        vtkCellData* bcdata = outputMesh->GetCellData();
        data=GetData(gr,"BOUNDARY_TEMPERATURE",1);
        AddToVTKScalar(bcdata,data,"BOUNDARY_TEMPERATURE");
        data=GetData(gr,"BOUNDARY_FORCE",1);
        AddToVTKScalar(bcdata,data,"BOUNDARY_FORCE");


        data=GetData(gr,"BOUNDARY_VELOCITY",4);
        AddToVTKVector(bcdata,data,"BOUNDARY_VELOCITY");

    }


    file.close();
    return 1;
}

void HDF5ReaderMultiBlock::PrintSelf(ostream& os, vtkIndent indent)
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

int HDF5ReaderMultiBlock::RequestInformation(vtkInformation *vtkNotUsed(request), vtkInformationVector **vtkNotUsed(inputVector), vtkInformationVector *outputVector)
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
        if (p.path().string().find(".h5") != std::string::npos)
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

    vtkInformation * outInfo1 = outputVector->GetInformationObject(1);
    outInfo1->Set(
                vtkStreamingDemandDrivenPipeline::TIME_STEPS(),
                times->GetPointer(0), nsteps
                );
    outInfo1->Set(
                vtkStreamingDemandDrivenPipeline::TIME_RANGE(),
                trange,2);




    return 1;
}
