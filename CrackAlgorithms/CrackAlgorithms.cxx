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
#include "CrackAlgorithms.h"
#include "vtkCellCenter3D.h"
#include "vtkVoronoi3D.h"

vtkStandardNewMacro(CrackAlgorithms);

//-----------------------------------------------------------------------------
CrackAlgorithms::CrackAlgorithms()
{
}


//-----------------------------------------------------------------------------
CrackAlgorithms::~CrackAlgorithms()
{
}


//----------------------------------------------------------------------------
void CrackAlgorithms::AddSourceConnection(vtkAlgorithmOutput* input)
{
  this->AddInputConnection(1, input);
}


//----------------------------------------------------------------------------
void CrackAlgorithms::RemoveAllSources()
{
  this->SetInputConnection(1, 0);
}


//----------------------------------------------------------------------------
int CrackAlgorithms::FillInputPortInformation( int port, vtkInformation* info )
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
int CrackAlgorithms::RequestData(vtkInformation *vtkNotUsed(request),
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
  vtkPolyData *output = vtkPolyData::SafeDownCast(
                                                  outInfo->Get(vtkDataObject::DATA_OBJECT()));

/*
vispartdem::vtkCellCenter3D::vtkCellCenter3D *cellCenters=vispartdem::vtkCellCenter3D::vtkCellCenter3D::New();
cellCenters->SetInputData(input);
cellCenters->SetStateArray("BOND_STATE");
cellCenters->SetResultArrayName("RESULT");
cellCenters->Update();
output->ShallowCopy(cellCenters->GetOutput());
cellCenters->Delete();
*/
  vispartdem::voronoi3D_LOCAL::vtkVoronoi3D* vor = vispartdem::voronoi3D_LOCAL::vtkVoronoi3D::New();
  vor->SetInputData(input);
  vor->SetStateArray_("BOND_STATE");
  vor->SetResultArrayName("BUSENA");
  vor->SetDeformations_(true);
  vor->Setcheckpailgejima_(false);
  vor->Update();

output->ShallowCopy(vor->GetOutput());
vor->Delete();





  return 1;
}



////////// External Operators /////////////

void CrackAlgorithms::PrintSelf(ostream &os, vtkIndent indent)
{
}
