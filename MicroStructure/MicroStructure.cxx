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
#include "MicroStructure.h"
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


vtkStandardNewMacro(MicroStructure);


//-----------------------------------------------------------------------------
MicroStructure::MicroStructure()
{
    CELL_SIZE=1;
}


//-----------------------------------------------------------------------------
MicroStructure::~MicroStructure()
{
}


//----------------------------------------------------------------------------
void MicroStructure::AddSourceConnection(vtkAlgorithmOutput* input)
{
  this->AddInputConnection(1, input);
}


//----------------------------------------------------------------------------
void MicroStructure::RemoveAllSources()
{
  this->SetInputConnection(1, 0);
}


//----------------------------------------------------------------------------
int MicroStructure::FillInputPortInformation( int port, vtkInformation* info )
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
int MicroStructure::RequestData(vtkInformation *vtkNotUsed(request),
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


  double bounds[6];
  input->GetBounds(bounds);
    vtkDoubleArray*f=vtkDoubleArray::New();
    f->SetName("FILTRAS");
    f->SetNumberOfComponents(1);
    f->SetNumberOfTuples(input->GetNumberOfPoints());


  for(int i=0;i<input->GetNumberOfPoints();i++)
  {
      double T[3];
      input->GetPoint(i,T);

      Taskas t;
      t.x=T[0];
      t.y=T[1];
      t.z=T[2];
      t.r=input->GetPointData()->GetArray("RADIUS")->GetTuple1(i);
      int id=0;
      int cx=std::ceil(T[0]/CELL_SIZE);
      int cy=std::ceil(T[1]/CELL_SIZE);
      int cz=std::ceil(T[2]/CELL_SIZE);
/*
      double rx=cx*CELL_SIZE+CELL_SIZE/2-T[0];
      double ry=cy*CELL_SIZE+CELL_SIZE/2-T[1];
      double rz=cz*CELL_SIZE+CELL_SIZE/2-T[2];
      double ilgis=std::sqrt(rx*rx+ry*ry+rz*rz);
      if(ilgis<RADIUS)
          id=1;*/


      for( int x=cx-1;x<cx+2;x++)
          for( int y=cy-1;y<cy+2;y++)
              for( int z=cz-1;z<cz+2;z++)
              {
                  double rx=x*CELL_SIZE-T[0];
                  double ry=y*CELL_SIZE-T[1];
                  double rz=z*CELL_SIZE-T[2];
                  double ilgis=std::sqrt(rx*rx+ry*ry+rz*rz);
                  if(ilgis<RADIUS)
                      id=1;


              }



      f->SetTuple1(i,id);
  }
  input->GetPointData()->AddArray(f);
//input->Modified();

 // imageData->Print(std::cout);
  vtkSmartPointer<vtkThreshold>tr=vtkSmartPointer<vtkThreshold>::New();
  tr->AllScalarsOn();
  tr->SetInputData(input);
  tr->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, "FILTRAS");
  tr->ThresholdBetween(-1,0);
  tr->Update();
  tr->GetOutput()->GetPointData()->SetActiveScalars("RADIUS");
 // tr->GetOutput()->Print(std::cout);
  output->DeepCopy(tr->GetOutput());
output->GetPointData()->SetActiveScalars("RADIUS");
  return 1;
}



////////// External Operators /////////////

void MicroStructure::PrintSelf(ostream &os, vtkIndent indent)
{
}
