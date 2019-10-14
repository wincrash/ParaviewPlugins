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
#include "vtkPolygon.h"


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
  vtkPolyData *output = vtkPolyData::SafeDownCast(
                                                  outInfo->Get(vtkDataObject::DATA_OBJECT()));



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






  vtkSmartPointer<vtkPoints> points=vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> cells=vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkDoubleArray> area=vtkSmartPointer<vtkDoubleArray>::New();
    area->SetName("AREA");
    area->SetNumberOfComponents(1);
    vtkSmartPointer<vtkDoubleArray> L1=vtkSmartPointer<vtkDoubleArray>::New();
    L1->SetName("L1");
    L1->SetNumberOfComponents(1);
    vtkSmartPointer<vtkDoubleArray> L2=vtkSmartPointer<vtkDoubleArray>::New();
    L2->SetName("L2");
    L2->SetNumberOfComponents(1);



    for (int i = 0; i < input->GetNumberOfPoints(); ++i) {
        int offset=points->GetNumberOfPoints();

        for(int k=0;k<faces[i].nrvertext;k++)
        {
         points->InsertNextPoint(faces[i].vertexes[k][0],faces[i].vertexes[k][1],faces[i].vertexes[k][2]);

        }
           for(int k=0;k<faces[i].nrfaces;k++)
           {
               if(faces[i].kaimynai[k]<0) continue;
           //    if(faces[i].facevertextnr[k]!=4)continue;

               cells->InsertNextCell(faces[i].facevertextnr[k]);
                              vtkSmartPointer<vtkPoints> points1=vtkSmartPointer<vtkPoints>::New();
                  for(int z=0;z<faces[i].facevertextnr[k];z++)
                  {
                    cells->InsertCellPoint(offset+faces[i].faceorder[k][z]);
                    int hh=faces[i].faceorder[k][z];
                     points1->InsertNextPoint(faces[i].vertexes[hh][0],faces[i].vertexes[hh][1],faces[i].vertexes[hh][2]);
                  }
                  vtkPolygon*poly=vtkPolygon::New();
                  poly->Initialize(faces[i].facevertextnr[k],points1);
                  //std::cout<<poly->ComputeArea()<<"\n";
                  area->InsertNextTuple1(poly->ComputeArea());
                  poly->Delete();
/*
                  double p1[3];
                  double p2[3];
                  double p3[3];
                  double p4[3];
                  points1->GetPoint(0,p1);
                  points1->GetPoint(1,p2);
                  points1->GetPoint(2,p3);
                  points1->GetPoint(3,p4);
                  L1->InsertNextTuple1(sqrt((p1[0]-p3[0])*(p1[0]-p3[0])+(p1[1]-p3[1])*(p1[1]-p3[1])+(p1[2]-p3[2])*(p1[2]-p3[2])));
                  L2->InsertNextTuple1(sqrt((p2[0]-p4[0])*(p2[0]-p4[0])+(p2[1]-p4[1])*(p2[1]-p4[1])+(p2[2]-p4[2])*(p2[2]-p4[2])));*/

           }
    }
  output->SetPoints(points);
  output->SetPolys(cells);
  output->GetCellData()->SetScalars(area);
  //output->GetCellData()->AddArray(L1);
  //output->GetCellData()->AddArray(L2);
delete[ ]faces;

  return 1;
}



////////// External Operators /////////////

void Voronoi::PrintSelf(ostream &os, vtkIndent indent)
{
}
