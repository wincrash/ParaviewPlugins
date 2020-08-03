#include "vtkObjectFactory.h" //for new() macro
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkUnstructuredGrid.h"
#include "vtkSmartPointer.h"
#include "vtkCellArray.h"
#include "vtkUnstructuredGrid.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkDataArray.h"
#include "vtkDoubleArray.h"
#include "Voronoi.h"
#include "vtkPolygon.h"
#include "vtkIntArray.h"

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
        //info->Set( vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid" );
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
    //  vtkUnstructuredGrid *input = vtkUnstructuredGrid::SafeDownCast(
    vtkDataSet *input = vtkDataSet::SafeDownCast(
                inInfo->Get(vtkDataObject::DATA_OBJECT()));
    vtkUnstructuredGrid *output = vtkUnstructuredGrid::SafeDownCast(
                outInfo->Get(vtkDataObject::DATA_OBJECT()));
    output->Allocate(1000,1000);


    double * bounds=input->GetBounds();

    double *maxr=input->GetPointData()->GetArray("RADIUS")->GetRange();
    double R=maxr[1]*2;
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
   /* vtkSmartPointer<vtkDoubleArray> area=vtkSmartPointer<vtkDoubleArray>::New();
    area->SetName("AREA");
    area->SetNumberOfComponents(1);
    vtkSmartPointer<vtkDoubleArray> L1=vtkSmartPointer<vtkDoubleArray>::New();
    L1->SetName("L1");
    L1->SetNumberOfComponents(1);
    vtkSmartPointer<vtkDoubleArray> L2=vtkSmartPointer<vtkDoubleArray>::New();
    L2->SetName("L2");
    L2->SetNumberOfComponents(1);
    */

    vtkSmartPointer<vtkIntArray> FACE_COUNT=vtkSmartPointer<vtkIntArray>::New();
    FACE_COUNT->SetName("FACE_COUNT");
    FACE_COUNT->SetNumberOfComponents(1);

    vtkSmartPointer<vtkIntArray> VERTICES_COUNT=vtkSmartPointer<vtkIntArray>::New();
    VERTICES_COUNT->SetName("VERTICES_COUNT");
    VERTICES_COUNT->SetNumberOfComponents(1);

    vtkSmartPointer<vtkIntArray> MAX_VERTICES_PER_FACE=vtkSmartPointer<vtkIntArray>::New();
    MAX_VERTICES_PER_FACE->SetName("MAX_VERTICES_PER_FACE");
    MAX_VERTICES_PER_FACE->SetNumberOfComponents(1);







    for(size_t i=0;i<input->GetNumberOfPoints();i++)
    {

//std::cout<<"particle i "<<i<<"\n";
        int pointStart=points->GetNumberOfPoints();
        int maxas=0;
        vtkIdType dodechedronPointsIds[faces[i].nrvertext];
        for(size_t k=0;k<faces[i].nrvertext;k++)
        {
            dodechedronPointsIds[k]=points->GetNumberOfPoints();
            points->InsertNextPoint(faces[i].vertexes[k][0],faces[i].vertexes[k][1],faces[i].vertexes[k][2]);
        }
        vtkSmartPointer<vtkCellArray> dodechedronFaces =vtkSmartPointer<vtkCellArray>::New();

        for (size_t z = 0; z < faces[i].nrfaces; z++)
        {
           // if(faces[i].kaimynai[z]<0)continue;
            if(maxas<faces[i].facevertextnr[z]) maxas=faces[i].facevertextnr[z];
            dodechedronFaces->InsertNextCell(faces[i].facevertextnr[z]);
            for(size_t k=0;k<faces[i].facevertextnr[z];k++)
                dodechedronFaces->InsertCellPoint(faces[i].faceorder[z][k]+pointStart);
        }

        output->InsertNextCell(VTK_POLYHEDRON,
                               faces[i].nrvertext, dodechedronPointsIds,
                               faces[i].nrfaces, dodechedronFaces->GetPointer());
         MAX_VERTICES_PER_FACE->InsertNextTuple1(maxas);
         VERTICES_COUNT->InsertNextTuple1(faces[i].nrvertext);
         FACE_COUNT->InsertNextTuple1(faces[i].nrfaces);
    }














    output->SetPoints(points);
    // output->SetPolys(cells);
    //output->GetCellData()->SetScalars(area);
    output->GetCellData()->AddArray(MAX_VERTICES_PER_FACE);
    output->GetCellData()->AddArray(VERTICES_COUNT);
    output->GetCellData()->AddArray(FACE_COUNT);
    //output->GetCellData()->AddArray(L1);
    //output->GetCellData()->AddArray(L2);
    delete[ ]faces;

    return 1;
}



////////// External Operators /////////////

void Voronoi::PrintSelf(ostream &os, vtkIndent indent)
{
}
