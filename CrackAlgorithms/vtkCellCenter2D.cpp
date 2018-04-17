/*

 #     # ###  #####  ######     #    ######  ####### ######  ####### #     #
 #     #  #  #     # #     #   # #   #     #    #    #     # #       ##   ##
 #     #  #  #       #     #  #   #  #     #    #    #     # #       # # # #
 #     #  #   #####  ######  #     # ######     #    #     # #####   #  #  #
 #   #   #        # #       ####### #   #      #    #     # #       #     #
 # #    #  #     # #       #     # #    #     #    #     # #       #     #
 #    ###  #####  #       #     # #     #    #    ######  ####### #     #


 * Author       : Ruslan Pacevic
 * E-mail       : rpa@sc.vgtu.lt
 * Vendor       : VGTU
 * Home page    : http://lsl.vgtu.lt/vispartdem
 */

#include "vtkCellCenter2D.h"
namespace vispartdem {
namespace vtkCellCenter2D {
vtkCxxRevisionMacro(vtkCellCenter2D, "$Revision: 1.42 $")
vtkStandardNewMacro(vtkCellCenter2D)

vtkCellCenter2D::vtkCellCenter2D() {
    this->StateArray = "";
    onetime = true;
    max_id = 0;
}
vtkCellCenter2D::~vtkCellCenter2D() {

}

int vtkCellCenter2D::RequestData(vtkInformation *vtkNotUsed(request), vtkInformationVector **inputVector, vtkInformationVector *outputVector) {
    // get the info objects
    Timeris t0, t1, t2, t3, t4;
    t4.start();
    t0.start();
    vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
    vtkInformation *outInfo = outputVector->GetInformationObject(0);
    // get the input and ouptut
    vtkPolyData *input = vtkPolyData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
    vtkPolyData *output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
    vtkDataArray*state = input->GetCellData()->GetArray(this->StateArray.c_str()); //gauna nurodyta state masyva
    if (state == NULL) {
        vtkErrorMacro( << "Nenurodytas state masyvas (turi buti idetas i celldata)!\n");
        return 0;
    }
    vtkPoints*points = vtkPoints::New(); //tie patys taskai
    vtkCellArray*faces = vtkCellArray::New(); // linijos bus dedami cia
    vtkDoubleArray**arrays = new vtkDoubleArray*[input->GetCellData()->GetNumberOfArrays()];
    for (int i = 0; i < input->GetCellData()->GetNumberOfArrays(); ++i) //sukuriam arrays i kuriuos kopijuosim celldata reiksmes
    {
        arrays[i] = vtkDoubleArray::New();
        arrays[i]->SetName(input->GetCellData()->GetArray(i)->GetName());
        output->GetCellData()->AddArray(arrays[i]);
        arrays[i]->Delete();
    }
    delete[] arrays;
    arrays = NULL;

    vispartdem::INT numberofcells = input->GetNumberOfCells();
    vispartdem::INT numberofpoints = input->GetNumberOfPoints();
    vispartdem::INT* connections = (vispartdem::INT*) (((input->GetLines()->GetPointer())));
    if (onetime == true) {
        indexing.resize(numberofpoints, 0);
        lIndex.first=-1;
        lIndex.second=-1;
        Ltotriag.resize(numberofcells,lIndex);
        onetime = false;
        cellSearch = new vispartdem::CellSearch(numberofpoints, numberofcells, connections);        
    }
    vispartdem::INT id1, id2, id,cellid1,cellid2,cellid3;
    //////////////////////////////////////////////////////////////////
    t0.stop();
    t1.start();
    vector<vispartdem::INT> tempIndex1;
    vector<vispartdem::INT> tempCellIndex;
    tempCellIndex.resize(numberofcells, 0);
    tempIndex1.resize(numberofpoints, 0);
    for (id = 0; id < numberofcells; ++id) {
        if (state->GetTuple1(id) == 0) {
            id1 = connections[3 * id + 1];
            id2 = connections[3 * id + 2];
            if (tempIndex1[id1] == 0) {
                tempIndex1[id1] = 1;
                vector<vispartdem::INT> *pointCells1 = cellSearch->getPointNeighboursCells(id1);
                for (vector<vispartdem::INT>::iterator iterator1 = pointCells1->begin(); iterator1 != pointCells1->end(); iterator1++) {
                    tempCellIndex[*iterator1] = 1;
                }
            }
            if (tempIndex1[id2] == 0) {
                tempIndex1[id2] = 1;
                vector<vispartdem::INT> *pointCells2 = cellSearch->getPointNeighboursCells(id2);
                for (vector<vispartdem::INT>::iterator iterator1 = pointCells2->begin(); iterator1 != pointCells2->end(); iterator1++) {
                    tempCellIndex[*iterator1] = 1;
                }
            }
        }
    }
    t1.stop();
    t2.start();
    //////////////////////////////////////////////////////////////////
    vispartdem::CombinationsUtils* comb2 = new vispartdem::CombinationsUtils();
    for (id = 0; id < numberofpoints; ++id) {
        if(indexing[id]==0 &&tempIndex1[id]!=5)
        {
            indexing[id]=1;
             vector<vispartdem::INT> *pointCells = cellSearch->getPointNeighbours(id);
            vector<vispartdem::INT> *cellspoints = cellSearch->getPointNeighboursCells(id);
            comb2->calcCombinations(pointCells->size(), 2);
            for (unsigned int k = 0; k < comb2->combination.size(); ++k) {
                id1 = (*pointCells)[comb2->combination[k]];
                id2 = (*pointCells)[comb2->combination[k + 1]];
                cellid1 = (*cellspoints)[comb2->combination[k]];
                cellid2 = (*cellspoints)[comb2->combination[k + 1]];
                k = k + 1;
                cellid3=cellSearch->cellID(id1,id2);
                if(cellid3==-1)
                {
                    continue;
                }
                triag.first.first = id1;
                triag.first.second = id2;
                triag.second = id;
                sort3(triag);
                if(triangMap.find(triag)!=triangMap.end())
                {
                    continue;
                }
                max_id++;
                triangMap[triag]=max_id;
                revtriangMap[max_id]=triag;
                if(Ltotriag[cellid1].first==-1)
                {
                    Ltotriag[cellid1].first=max_id;
                }else
                {
                    Ltotriag[cellid1].second=max_id;
                }
                if(Ltotriag[cellid2].first==-1)
                {
                    Ltotriag[cellid2].first=max_id;
                }else
                {
                    Ltotriag[cellid2].second=max_id;
                }
                if(Ltotriag[cellid3].first==-1)
                {
                    Ltotriag[cellid3].first=max_id;
                }else
                {
                    Ltotriag[cellid3].second=max_id;
                }
            }

        }
    }
    t2.stop();
    t3.start();
    //////////////////////////////////////////////////////////////////
    vtkIntArray*errorarray = vtkIntArray::New();
    errorarray->SetName(ResultArrayName.c_str());
    errorarray->SetNumberOfComponents(1);
    simpleMap pointsIndex;
    double p1[3];
    double p2[3];
    double p3[3];
    double center[3];
    simpleMap::iterator lb1, lb2;
    LineIndexVal line;
    double point1[3];
    double point2[3];
    for (id = 0; id < numberofcells; ++id) {
        if (tempCellIndex[id] != 0) {
            line = Ltotriag[id];
            if (line.first != -1 && line.second != -1) {
                faces->InsertNextCell(2);
                copyCellsValues(id, input->GetCellData(), output->GetCellData());
                errorarray->InsertNextTuple1( 0);
                lb1 = pointsIndex.lower_bound(line.first);
                if (lb1 != pointsIndex.end() && !(pointsIndex.key_comp()(line.first, lb1->first))) {
                    faces->InsertCellPoint(lb1->second);
                } else {
                    pointsIndex.insert(lb1, simpleMap::value_type(line.first, points->GetNumberOfPoints())); // Use lb as a hint to insert,
                    faces->InsertCellPoint(points->GetNumberOfPoints());
                    triangle_index ids = revtriangMap[line.first];
                    input->GetPoint(ids.first.first, p1);
                    input->GetPoint(ids.first.second, p2);
                    input->GetPoint(ids.second, p3);
                    TriangleCenter(p1, p2, p3, center);
                    points->InsertNextPoint(center);
                }
                lb2 = pointsIndex.lower_bound(line.second);
                if (lb2 != pointsIndex.end() && !(pointsIndex.key_comp()(line.second, lb2->first))) {
                    faces->InsertCellPoint(lb2->second);
                } else {
                    pointsIndex.insert(lb2, simpleMap::value_type(line.second, points->GetNumberOfPoints())); // Use lb as a hint to insert,
                    faces->InsertCellPoint(points->GetNumberOfPoints());
                    triangle_index ids = revtriangMap[line.second];
                    input->GetPoint(ids.first.first, p1);
                    input->GetPoint(ids.first.second, p2);
                    input->GetPoint(ids.second, p3);
                    TriangleCenter(p1, p2, p3, center);
                    points->InsertNextPoint(center);
                }

            }
        }
    }


    t3.stop();
    output->GetCellData()->AddArray(errorarray);
    errorarray->Delete();
    output->SetPoints(points);
    output->SetLines(faces);
    faces->Delete();
    points->Delete();
    output->Squeeze();
    t4.stop();
    //////////////////////////////////////////////////////////////////
    t0.print("incircle_2D_paruosimas");
    t1.print("incircle_2D_surinkimas_nagrinejamu_jungciu");
    t2.print("incircle_2D_suradimas_trukstamu_trikampiu");
    t3.print("incircle_2D_generavimas_voronoi_edge");
    t4.print("incircle_2D_bendras_laikas");
    return 1;
}

void vtkCellCenter2D::TriangleCenter(double*p1, double *p2, double*p3, double *center) {
    /*double a=vtkMath::Distance2BetweenPoints(p2,p3);
    double b=vtkMath::Distance2BetweenPoints(p1,p3);
    double c=vtkMath::Distance2BetweenPoints(p1,p2);
    double P=a+b+c;
    center[0]=a*p1[0]+b*p2[0]+c*p3[0];
    center[1]=a*p1[1]+b*p2[1]+c*p3[1];
    center[2]=a*p1[2]+b*p2[2]+c*p3[2];
    center[0]=center[0]/P;
    center[1]=center[1]/P;
    center[2]=center[2]/P;*/
    /*center[0]=(p1[0]+p2[0]+p3[0])/3;
    center[1]=(p1[1]+p2[1]+p3[1])/3;
   center[2]=(p1[2]+p2[2]+p3[2])/3;*/
    vtkTriangle::TriangleCenter(p1, p2, p3, center);
}
int vtkCellCenter2D::FillInputPortInformation(int, vtkInformation *info) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
    return 1;
}

void vtkCellCenter2D::PrintSelf(ostream& os, vtkIndent indent) {
    this->Superclass::PrintSelf(os, indent);
    os << "State array name " << StateArray << endl;
}
inline void vtkCellCenter2D::copyCellsValues(vispartdem::INT &i, vtkCellData*icelldata, vtkCellData * ocelldata) //kopijuouja celldatoje esancias reiksmes
{
    for (int m = 0; m < icelldata->GetNumberOfArrays(); ++m) {
        ocelldata->GetArray(m)->InsertNextTuple(icelldata->GetArray(m)->GetTuple(i));
    }
}

}
}
