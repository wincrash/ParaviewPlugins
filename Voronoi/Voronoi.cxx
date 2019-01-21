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
#include "Voronoi.h"


vtkStandardNewMacro(Voronoi);


struct cele3d {
    int id;
    int nrvertext;
    int nrfaces;
    int* facevertextnr;
    int * kaimynai;
    int **faceorder;
    float** vertexes;
    ~cele3d() {

        int k = 0;

        for (k = 0; k < nrfaces; k++) {
            delete[] faceorder[k];

        }

        for (k = 0; k < nrvertext; k++) {

            delete[] vertexes[k];

        }
        if (nrfaces != 0 && nrvertext != 0) {
            delete[] vertexes;
            delete[] kaimynai;
            delete[] facevertextnr;
            delete[] faceorder;
        }

    }
    ;
    cele3d() {
    }
    ;

};
#include "../voro/voro++.cpp"
//-----------------------------------------------------------------------------
Voronoi::Voronoi()
{
}


//-----------------------------------------------------------------------------
Voronoi::~Voronoi()
{
}


//----------------------------------------------------------------------------
void Voronoi::AddSourceConnection(vtkAlgorithmOutput* input)
{
  this->AddInputConnection(1, input);
}


//----------------------------------------------------------------------------
void Voronoi::RemoveAllSources()
{
  this->SetInputConnection(1, 0);
}


//----------------------------------------------------------------------------
int Voronoi::FillInputPortInformation( int port, vtkInformation* info )
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
int Voronoi::RequestData(vtkInformation *vtkNotUsed(request),
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

vtkSmartPointer<vtkUnstructuredGrid> OutputPolydata = vtkSmartPointer<vtkUnstructuredGrid>::New();


   double * bounds=input->GetBounds();

   double *maxr=input->GetPointData()->GetArray("RADIUS")->GetRange();
   double R=maxr[1];
    bounds[0]=bounds[0]-R;
    bounds[1]=bounds[1]+R;

    bounds[2]=bounds[2]-R;
    bounds[3]=bounds[3]+R;

    bounds[4]=bounds[4]-R;
    bounds[5]=bounds[5]+R;

 container_poly c(bounds[0], bounds[1], bounds[2], bounds[3], bounds[4], bounds[5], 8, 8, 8, false, false, false, 8);
  cele3d *faces = new cele3d[input->GetNumberOfPoints()];
  for (int i = 0; i < input->GetNumberOfPoints(); ++i) {
      c.put(i,input->GetPoint(i)[0],input->GetPoint(i)[1],input->GetPoint(i)[2],input->GetPointData()->GetArray("RADIUS")->GetTuple1(i));
      faces[i].nrfaces = 0;
      faces[i].id = i;
      faces[i].nrvertext = 0;
  }
  c.print_all_custom("%i %w %s %a %n %t %P", faces);




std::vector<vtkDoubleArray*> arrays;
for(int i=0;i<input->GetPointData()->GetNumberOfArrays();i++)
{
    arrays.push_back(vtkDoubleArray::New());
    arrays[i]->SetName(input->GetPointData()->GetArray(i)->GetName());
    arrays[i]->SetNumberOfComponents(input->GetPointData()->GetArray(i)->GetNumberOfComponents());
    arrays[i]->SetNumberOfTuples(input->GetPointData()->GetArray(i)->GetNumberOfTuples());
}



OutputPolydata->Allocate(1,1);


  vtkSmartPointer<vtkPoints> points=vtkSmartPointer<vtkPoints>::New();

  for(int i=0;i<input->GetNumberOfPoints();i++)
  {
      for(int h=0;h<arrays.size();h++)
      {
          if(arrays[h]->GetNumberOfComponents()==1)
          {
              arrays[h]->SetTuple1(i,input->GetPointData()->GetArray(h)->GetTuple1(i));
            //  double val=input->GetPointData()->GetArray(h)->GetTuple1(i);
             // std::cout<<val<<"\n";
            //arrays[h]->SetTuple1(i,val);
          }else
          {
              double*vals=input->GetPointData()->GetArray(h)->GetTuple3(i);
              arrays[h]->SetTuple3(i,vals[0],vals[1],vals[2]);
          }

      }

       int start=points->GetNumberOfPoints();
      vtkIdType *pointIds=new vtkIdType[faces[i].nrvertext];
      for(int z=0;z<faces[i].nrvertext;z++)
      {
          points->InsertNextPoint(faces[i].vertexes[z][0],faces[i].vertexes[z][1],faces[i].vertexes[z][2]);
          pointIds[z]=points->GetNumberOfPoints()-1;
      }

      std::cout<<faces[i].id<<" "<<faces[i].nrvertext<<" "<<faces[i].nrfaces<<"   "<<faces[i].facevertextnr<<"\n";
      vtkSmartPointer<vtkCellArray> FF = vtkSmartPointer<vtkCellArray>::New();
        FF->Allocate(1,1);



      for(int k=0;k<faces[i].nrfaces;k++)
      {
          FF->InsertNextCell(faces[i].facevertextnr[k]);
          for(int j=0;j<faces[i].facevertextnr[k];j++)
          {
              FF->InsertCellPoint(start+faces[i].faceorder[k][j]);
          }
      }
      OutputPolydata->InsertNextCell(VTK_POLYHEDRON, faces[i].nrvertext, pointIds,
          faces[i].nrfaces, FF->GetPointer());
        delete[]pointIds;

  }
delete[ ]faces;


OutputPolydata->SetPoints(points);
for(int h=0;h<input->GetPointData()->GetNumberOfArrays();h++)
{
    OutputPolydata->GetCellData()->AddArray(arrays[h]);
    //arrays[h]->Delete();
}
output->ShallowCopy(OutputPolydata);

  return 1;
}



////////// External Operators /////////////

void Voronoi::PrintSelf(ostream &os, vtkIndent indent)
{
}
