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
void AddToVTKScalar(auto *location,std::vector<double>&array,std::string name)
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

void AddToVTKVector(auto *location,std::vector<double>&array,std::string name)
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
    this->SetNumberOfOutputPorts(1);
    this->times = vtkSmartPointer<vtkDoubleArray>::New();
    this->timesNames = vtkSmartPointer<vtkVariantArray>::New();
}

int HDF5ReaderMultiBlock::RequestData(
        vtkInformation *vtkNotUsed(request),
        vtkInformationVector **vtkNotUsed(inputVector),
        vtkInformationVector *outputVector)
{

    vtkInformation *outInfo = outputVector->GetInformationObject(0);
    double requestedTime = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());
    std::cout<<"  "<<findClosestSolutionIndex(requestedTime)<<"\n";



    // get the ouptut
    vtkMultiBlockDataSet *outputBlock = vtkMultiBlockDataSet::SafeDownCast(
                outInfo->Get(vtkDataObject::DATA_OBJECT()));
    vtkPolyData* output=vtkPolyData::New();

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
    //vtkPolyData*output = vtkPolyData::New();
    vtkPoints*points = vtkPoints::New();
    std::vector<double> POSITIONS=GetData(gr,"POSITIONS",4);
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
        data=GetData(gr,"BOND_STATE",1);
        AddToVTKScalar(cdata,data,"BOND_STATE");

        data=GetData(gr,"BOND_F_LIMIT_N",1);
        AddToVTKScalar(cdata,data,"BOND_F_LIMIT_N");

        data=GetData(gr,"BOND_F_LIMIT_T",1);
        AddToVTKScalar(cdata,data,"BOND_F_LIMIT_T");

        data=GetData(gr,"BOND_COUNT",1);
        AddToVTKScalar(pdata,data,"BOND_COUNT");

        std::vector<double> particle1=GetData(gr,"BOND_PARTICLE_1",1);
        std::vector<double> particle2=GetData(gr,"BOND_PARTICLE_2",1);
        vtkCellArray*lines=vtkCellArray::New();
        for(int i=0;i<particle1.size();i++)
        {
            lines->InsertNextCell(2);
            int id1=particle1[i];
            int id2=particle2[i];

            lines->InsertCellPoint(id1);
            lines->InsertCellPoint(id2);
        }
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
    if(output->GetPointData()->HasArray("RADIUS"))
    {
        output->GetPointData()->SetActiveScalars("RADIUS");
    }
    if(output->GetPointData()->HasArray("VELOCITY"))
    {
        output->GetPointData()->SetActiveVectors("VELOCITY");
    }
    points->Delete();
    outputBlock->SetBlock(0,output);

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
            /*if(IDS[i*4+3]>=0)
            {
                bcells->InsertNextCell(4);
                bcells->InsertCellPoint(IDS[i*4+0]);
                bcells->InsertCellPoint(IDS[i*4+1]);
                bcells->InsertCellPoint(IDS[i*4+2]);
                bcells->InsertCellPoint(IDS[i*4+3]);
            }*/
          //  else
            {
                bcells->InsertNextCell(3);
                bcells->InsertCellPoint(IDS[i*4+0]);
                bcells->InsertCellPoint(IDS[i*4+1]);
                bcells->InsertCellPoint(IDS[i*4+2]);
            }
        }
        vtkPolyData* bpoly=vtkPolyData::New();
        bpoly->SetPoints(bpoints);
        bpoly->SetPolys(bcells);

        vtkCellData* bcdata = bpoly->GetCellData();
        vtkFieldData* bfdata = bpoly->GetFieldData();
        data=GetData(gr,"BOUNDARY_TEMPERATURE",1);
        AddToVTKScalar(bcdata,data,"BOUNDARY_TEMPERATURE");
        outputBlock->SetBlock(1,bpoly);

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




    return 1;
}