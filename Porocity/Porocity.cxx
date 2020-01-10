#include "vtkObjectFactory.h" //for new() macro
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkCellArray.h"
#include "vtkUnstructuredGrid.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkDataArray.h"
#include "vtkDoubleArray.h"
#include "Porocity.h"
#include "vtkPolygon.h"
#include "vtkImageAlgorithm.h"
#include "vtkImageData.h"
#include "vtkUnstructuredGrid.h"
#include "vtkThreshold.h"
#include <vector>
struct Taskas {
    double x;
    double y;
    double z;
    double r;
};


vtkStandardNewMacro(Porocity);


//-----------------------------------------------------------------------------
Porocity::Porocity()
{
    CELL_SIZE=0.1;
}


//-----------------------------------------------------------------------------
Porocity::~Porocity()
{
}


//----------------------------------------------------------------------------
void Porocity::AddSourceConnection(vtkAlgorithmOutput* input)
{
  this->AddInputConnection(1, input);
}


//----------------------------------------------------------------------------
void Porocity::RemoveAllSources()
{
  this->SetInputConnection(1, 0);
}


//----------------------------------------------------------------------------
int Porocity::FillInputPortInformation( int port, vtkInformation* info )
{
  if (!this->Superclass::FillInputPortInformation(port, info))
    {
    return 0;
    }
  if ( port == 0 )
    {
    //info->Set( vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData" );
    info->Set( vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet" );
    return 1;
    }
  return 0;
}


//----------------------------------------------------------------------------
int Porocity::RequestData(vtkInformation *vtkNotUsed(request),
                              vtkInformationVector **inputVector,
                              vtkInformationVector *outputVector)
{
  // get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // get the input and output
//  vtkPolyData *input = vtkPolyData::SafeDownCast(
  vtkDataSet *input = vtkDataSet::SafeDownCast(
                                               inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkUnstructuredGrid *output = vtkUnstructuredGrid::SafeDownCast(
                                                  outInfo->Get(vtkDataObject::DATA_OBJECT()));

 vtkSmartPointer<vtkImageData> imageData =vtkSmartPointer< vtkImageData>::New();
  //vtkImageData* imageData =vtkImageData::New();


  std::vector<Taskas> taskai;


  double bounds[6]={100000000,-100000000,100000000,-100000000,100000000,-100000000};

  for(int i=0;i<input->GetNumberOfPoints();i++)
  {
      double T[3];
      input->GetPoint(i,T);
      Taskas t;
      t.x=T[0];
      t.y=T[1];
      t.z=T[2];
      t.r=input->GetPointData()->GetArray("RADIUS")->GetTuple1(i);
      taskai.push_back(t);
      if((t.x-t.r)<bounds[0]) bounds[0]=t.x-t.r;
      if((t.y-t.r)<bounds[2]) bounds[2]=t.y-t.r;
      if((t.z-t.r)<bounds[4]) bounds[4]=t.z-t.r;

      if((t.x+t.r)>bounds[1]) bounds[1]=t.x+t.r;
      if((t.y+t.r)>bounds[3]) bounds[3]=t.y+t.r;
      if((t.z+t.r)>bounds[5]) bounds[5]=t.z+t.r;
  }
  //0.007
 // double CELL_SIZE=0.0025*10;
  int NX=std::floor((bounds[1]-bounds[0])/CELL_SIZE)+1;
  int NY=std::floor((bounds[3]-bounds[2])/CELL_SIZE)+1;
  int NZ=std::floor((bounds[5]-bounds[4])/CELL_SIZE)+1;
  // Specify the size of the image data
  imageData->SetDimensions(NX+2,NY+2,NZ+2);
  imageData->SetSpacing(CELL_SIZE,CELL_SIZE,CELL_SIZE);
  imageData->SetOrigin(bounds[0],bounds[2],bounds[4]);
  vtkSmartPointer<vtkDoubleArray>count=vtkSmartPointer<vtkDoubleArray>::New();
  count->SetName("PARTICLES_COUNT");
  count->SetNumberOfComponents(1);
  count->SetNumberOfTuples(imageData->GetNumberOfCells());

  vtkSmartPointer<vtkDoubleArray> vol=vtkSmartPointer<vtkDoubleArray>::New();
  vol->SetName("PARTICLES_VOL");
  vol->SetNumberOfComponents(1);
  vol->SetNumberOfTuples(imageData->GetNumberOfCells());

  vtkSmartPointer<vtkDoubleArray> porocity=vtkSmartPointer<vtkDoubleArray>::New();
  porocity->SetName("POROCITY");
  porocity->SetNumberOfComponents(1);
  porocity->SetNumberOfTuples(imageData->GetNumberOfCells());

  for(int i=0;i<imageData->GetNumberOfCells();i++)
  {
      porocity->SetTuple1(i,0);
      vol->SetTuple1(i,0);
      count->SetTuple1(i,0);
  }
  int ijk[3]={0,0,0};
  double CELES_TURIS=CELL_SIZE*CELL_SIZE*CELL_SIZE;
  for(size_t i=0;i<taskai.size();i++)
  {

      Taskas t;
      t=taskai[i];
      t.x=t.x-bounds[0];
      t.y=t.y-bounds[2];
      t.z=t.z-bounds[4];
      ijk[0]=std::floor(t.x/CELL_SIZE)+1;
      ijk[1]=std::floor(t.y/CELL_SIZE)+1;
      ijk[2]=std::floor(t.z/CELL_SIZE)+1;

      int id=imageData->ComputeCellId (ijk);
     // std::cout<<id<<" "<<ijk[0]<<" "<<ijk[1]<<" "<<ijk[2]<<"  "<<t.r<<"\n";
      double kiekis=count->GetTuple1(id);
      double turis=vol->GetTuple1(id);
      double por=porocity->GetTuple1(id);
      kiekis=kiekis+1.0;
      turis=turis+4.0/3.0*t.r*t.r*t.r*3.14159265359;
      por=turis/CELES_TURIS;
      count->SetTuple1(id,kiekis);
      vol->SetTuple1(id,turis);
      porocity->SetTuple1(id,por);

  }
  imageData->GetCellData()->SetScalars(porocity);
  imageData->GetCellData()->AddArray(vol);
  imageData->GetCellData()->AddArray(count);

 // imageData->Print(std::cout);
  vtkSmartPointer<vtkThreshold>tr=vtkSmartPointer<vtkThreshold>::New();
  tr->AllScalarsOn();
  tr->SetInputData(imageData);
  tr->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, "POROCITY");
  tr->ThresholdBetween(-1000000000,1000000000);
  tr->Update();
 // tr->GetOutput()->Print(std::cout);
  output->DeepCopy(tr->GetOutput());

  return 1;
}



////////// External Operators /////////////

void Porocity::PrintSelf(ostream &os, vtkIndent indent)
{
}
