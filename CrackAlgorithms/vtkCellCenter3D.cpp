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

#include "vtkCellCenter3D.h"
namespace vispartdem {
namespace vtkCellCenter3D {

vtkStandardNewMacro(vtkCellCenter3D)

vtkCellCenter3D::vtkCellCenter3D() {
    this->StateArray = "";
    onetime = true;
}
vtkCellCenter3D::~vtkCellCenter3D() {

}

void vtkCellCenter3D::GetPyramidID(_tetrahedron & piramid, vispartdem::INT & id1, vispartdem::INT & id2, vispartdem::INT & id3, vispartdem::INT & id4) {
    int a[4];
    a[0] = id1;
    a[1] = id2;
    a[2] = id3;
    a[3] = id4;
    sort(a, a + 4);
    piramid.first.first.first = a[0];
    piramid.first.first.second = a[1];
    piramid.first.second = a[2];
    piramid.second = a[3];
}

void vtkCellCenter3D::GetOctahedronID(_octahedron & ket, vispartdem::INT & id1, vispartdem::INT & id2, vispartdem::INT & id3, vispartdem::INT & id4, vispartdem::INT & id5, vispartdem::INT & id6) {
    int a[6];
    a[0] = id1;
    a[1] = id2;
    a[2] = id3;
    a[3] = id4;
    a[4] = id5;
    a[5] = id6;
    sort(a, a + 6);
    ket.first.first.first.first.first = a[0];
    ket.first.first.first.first.second = a[1];
    ket.first.first.first.second = a[2];
    ket.first.first.second = a[3];
    ket.first.second = a[4];
    ket.second = a[5];
}
bool sort_pred(const pair<int, int>& left, const pair<int, int>& right) {
    return left.second < right.second;
}
int vtkCellCenter3D::RequestData(vtkInformation *vtkNotUsed(request), vtkInformationVector **inputVector, vtkInformationVector *outputVector) {
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
    vtkCellArray*faces = vtkCellArray::New(); /// trikiampiai bus dedami cia

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
        indexingP.resize(numberofpoints, 0);
        onetime = false;
        cellSearch = new vispartdem::CellSearch(numberofpoints, numberofcells, connections);
    }

    vispartdem::INT id1, id2, id3, id4, id, d1;

    t0.stop();
    t1.start();

    vector<vispartdem::INT> tempIndex1;
    vector<vispartdem::INT> tempCellIndex;
    tempCellIndex.resize(numberofcells, 0);
    tempIndex1.resize(numberofpoints, 0);
    for (id = 0; id < numberofcells; ++id) {
        if (state->GetTuple1(id) != 0) {
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

    //-----get numbers of cells and points

    _octahedron ket;
    _tetrahedron pyra;

    MyStruct pyra_struct; //struktura kurioje saugomi tetraedrai ir octaedrai
    pyra_struct.ids.resize(6, -1);

    MyStruct tetra_struct; //struktura kurioje saugomi tetraedrai ir octaedrai
    tetra_struct.ids.resize(4, -1);
    int istrizainiu_kiekis;
    vector<vispartdem::INT> keturkapioIstrizaines; //ieskosim istrizainiu kuriu neturi buti
    keturkapioIstrizaines.resize(4, -1);
    vispartdem::CombinationsUtils* comb3 = new vispartdem::CombinationsUtils();
    vispartdem::CombinationsUtils* comb4 = new vispartdem::CombinationsUtils();
    vispartdem::INT cellID1, cellID2, cellID3;
    for (id = 0; id < numberofpoints; ++id) {
        if (tempIndex1[id] != 0) {
            if (indexingP[id] == 0) {
                indexingP[id] = 1;
                //sakninis pointas
                vector<vispartdem::INT> *pointCells = cellSearch->getPointNeighbours(id);
                comb3->calcCombinations(pointCells->size(), 3);
                comb4->calcCombinations(pointCells->size(), 4);

                for (unsigned int k = 0; k < comb3->combination.size(); ++k) {
                    id1 = (*pointCells)[comb3->combination[k]];
                    id2 = (*pointCells)[comb3->combination[k + 1]];
                    id3 = (*pointCells)[comb3->combination[k + 2]];
                    k = k + 2;
                    //darom tetraedra
                    GetPyramidID(pyra, id, id1, id2, id3);
                    if (pyramid.find(pyra) == pyramid.end()) {
                        //tikrinam ar buvo padarytas toks tetraedras
                        if ((cellID1 = cellSearch->cellID(id1, id2)) != -1) { //sutikrinam visas junktis
                            if ((cellID2 = cellSearch->cellID(id2, id3)) != -1) {
                                if ((cellID3 = cellSearch->cellID(id3, id1)) != -1) {
                                    pyramid[pyra] = 1; //pazymime kad buvo bus toks tetraedras
                                    id4 = -2; // pazymime kad buvo sudarytas tetraedras ir neieskosime tetraedro
                                    tetra_struct.ids[0] = id1; //pagrindas A
                                    tetra_struct.ids[1] = id2; //pagrindas	B
                                    tetra_struct.ids[2] = id3; //pagrindas	C
                                    tetra_struct.ids[3] = id; //virsune	D
                                    realationByCell[cellID1].push_back(struktura.size());
                                    realationByCell[cellID2].push_back(struktura.size());
                                    realationByCell[cellID3].push_back(struktura.size());
                                    realationByCell[cellSearch->cellID(id, id1)].push_back(struktura.size());
                                    realationByCell[cellSearch->cellID(id, id2)].push_back(struktura.size());
                                    realationByCell[cellSearch->cellID(id, id3)].push_back(struktura.size());
                                    struktura.push_back(tetra_struct); //sudedam i struktura, kurioje saugomi tetraedrai ir octaedrai
                                }
                            }
                        }
                    }
                }
                for (unsigned int k = 0; k < comb4->combination.size(); ++k) {
                    id1 = (*pointCells)[comb4->combination[k]];
                    id2 = (*pointCells)[comb4->combination[k + 1]];
                    id3 = (*pointCells)[comb4->combination[k + 2]];
                    id4 = (*pointCells)[comb4->combination[k + 3]];
                    k = k + 3;
                    istrizainiu_kiekis = 0;

                    if (!cellSearch->cellExist(id1, id2)) {
                        keturkapioIstrizaines[istrizainiu_kiekis++] = id1;
                        keturkapioIstrizaines[istrizainiu_kiekis++] = id2;
                    }
                    if (!cellSearch->cellExist(id2, id3)) {
                        keturkapioIstrizaines[istrizainiu_kiekis++] = id2;
                        keturkapioIstrizaines[istrizainiu_kiekis++] = id3;
                    }

                    if (!cellSearch->cellExist(id3, id4)) {
                        if (istrizainiu_kiekis == 4) {
                            continue;
                        }
                        keturkapioIstrizaines[istrizainiu_kiekis++] = id3;
                        keturkapioIstrizaines[istrizainiu_kiekis++] = id4;
                    }

                    if (!cellSearch->cellExist(id4, id1)) {
                        if (istrizainiu_kiekis == 4) {
                            continue;
                        }
                        keturkapioIstrizaines[istrizainiu_kiekis++] = id4;
                        keturkapioIstrizaines[istrizainiu_kiekis++] = id1;
                    }

                    if (!cellSearch->cellExist(id1, id3)) {
                        if (istrizainiu_kiekis == 4) {
                            continue;
                        }
                        keturkapioIstrizaines[istrizainiu_kiekis++] = id1;
                        keturkapioIstrizaines[istrizainiu_kiekis++] = id3;
                    }

                    if (!cellSearch->cellExist(id2, id4)) {
                        if (istrizainiu_kiekis == 4) {
                            continue;
                        }
                        keturkapioIstrizaines[istrizainiu_kiekis++] = id2;
                        keturkapioIstrizaines[istrizainiu_kiekis++] = id4;
                    }
                    if (istrizainiu_kiekis != 4) {
                        continue;
                    }

                    vector<vispartdem::INT>* pointCells1 = cellSearch->getPointNeighbours(id1);
                    vector<vispartdem::INT>::iterator end = pointCells1->end();
                    for (vector<vispartdem::INT>::iterator iterator = pointCells1->begin(); iterator != end; iterator++) {
                        d1 = *iterator;
                        if (d1 == id) {
                            //svarbu kad nebutu jau pradinis mazgo id
                            continue;
                        }
                        if (!cellSearch->cellExist(id2, d1)) {
                            continue;
                        }
                        if (!cellSearch->cellExist(id3, d1)) {
                            continue;
                        }
                        if (!cellSearch->cellExist(id4, d1)) {
                            continue;
                        }
                        if (d1 != id && d1 != id1 && d1 != id2 && d1 != id3 && d1 != id4) {
                            //tikrinam ar tikrai galima kurti octahedrona
                            GetOctahedronID(ket, id, id1, id2, id3, id4, d1);
                            if (keturkampiai.find(ket) == keturkampiai.end()) {
                                //tikrinam ar nebuvo sukurtas oc
                                keturkampiai[ket] = 1; //pazymime kad toki jau radome (tai yra 2 piramides su bendru pagrindu keturkampiu)
                                //pagal istrizainiu nebuvimo indeksus sukonstrojame keturkampi kuri apiesimime tokia tvarka (0,2,1,3), bet mes dar nezinome ar apeiname desines rankos ar kaires rankos taisykle

                                pyra_struct.ids[0] = keturkapioIstrizaines[0]; //b
                                pyra_struct.ids[1] = keturkapioIstrizaines[2]; //c
                                pyra_struct.ids[2] = keturkapioIstrizaines[1]; //d
                                pyra_struct.ids[3] = keturkapioIstrizaines[3]; //e
                                pyra_struct.ids[4] = id; //a
                                pyra_struct.ids[5] = d1; //f

                                realationByCell[cellSearch->cellID(pyra_struct.ids[0], pyra_struct.ids[1])].push_back(struktura.size());
                                realationByCell[cellSearch->cellID(pyra_struct.ids[1], pyra_struct.ids[2])].push_back(struktura.size());
                                realationByCell[cellSearch->cellID(pyra_struct.ids[2], pyra_struct.ids[3])].push_back(struktura.size());
                                realationByCell[cellSearch->cellID(pyra_struct.ids[3], pyra_struct.ids[0])].push_back(struktura.size());
                                realationByCell[cellSearch->cellID(id, pyra_struct.ids[0])].push_back(struktura.size());
                                realationByCell[cellSearch->cellID(id, pyra_struct.ids[1])].push_back(struktura.size());
                                realationByCell[cellSearch->cellID(id, pyra_struct.ids[2])].push_back(struktura.size());
                                realationByCell[cellSearch->cellID(id, pyra_struct.ids[3])].push_back(struktura.size());
                                realationByCell[cellSearch->cellID(d1, pyra_struct.ids[0])].push_back(struktura.size());
                                realationByCell[cellSearch->cellID(d1, pyra_struct.ids[1])].push_back(struktura.size());
                                realationByCell[cellSearch->cellID(d1, pyra_struct.ids[2])].push_back(struktura.size());
                                realationByCell[cellSearch->cellID(d1, pyra_struct.ids[3])].push_back(struktura.size());
                                struktura.push_back(pyra_struct); //idedam i bendra array kur yra ir tetraedrai ir octaedrai
                            }
                        }
                    }
                }
            }
        }
    }

    ////
    t2.stop();
    t3.start();


    vtkIntArray*errorarray = vtkIntArray::New();
    errorarray->SetName(ResultArrayName.c_str());
    errorarray->SetNumberOfComponents(1);
    double center[3];
    vector<vispartdem::INT> IndexForNewPoints;
    IndexForNewPoints.resize(struktura.size(), -1);
    vector<pair<int, int> > tt;
    for (id = 0; id < numberofcells; ++id) {
        if (tempCellIndex[id] != 0) {
            tt.resize(realationByCell[id].size(), make_pair(-1, -1));
            for (unsigned int k = 0; k < realationByCell[id].size(); k++) {
                if (IndexForNewPoints[realationByCell[id][k]] == -1) {//id centro tos structuros
                    IndexForNewPoints[realationByCell[id][k]] = points->GetNumberOfPoints();//priskiriam id
                    CalcCenter(struktura[realationByCell[id][k]].ids, input, center);//suskaiciuojam centra
                    points->InsertNextPoint(center);//idedam centra
                }
                tt[k] = make_pair(k, struktura[realationByCell[id][k]].ids.size());
            }
            if (tt.size() == 4) {
                sort(tt.begin(), tt.end(), sort_pred);
                swap(tt[1], tt[2]);
            }

         /*   if(tt.size()==4)
            {
                points->GetPoint(IndexForNewPoints[realationByCell[id][tt[0].first]],point1);
                points->GetPoint(IndexForNewPoints[realationByCell[id][tt[1].first]],point2);
                points->GetPoint(IndexForNewPoints[realationByCell[id][tt[2].first]],point3);
                points->GetPoint(IndexForNewPoints[realationByCell[id][tt[3].first]],point4);

                vtkTriangle::ComputeNormal(point1, point2, point3, normal);

                double d = point1[0] * normal[0] + point1[1] * normal[1] + point1[2] * normal[2];
                cout<<"dddd "<<d<<endl;
             //   double DD = (normal[0] * virsune1[0] + normal[1] * virsune1[1] + normal[2] * virsune1[2] + d) / vtkMath::Norm(normal);
                if (d < 0) {
                  //  swap(tt[1], tt[3]);
                }
                if (d> 0) {
                //    swap(tt[1], tt[3]);
                }

            }
            if(tt.size()==3)
            {
                points->GetPoint(IndexForNewPoints[realationByCell[id][tt[0].first]],point1);
                points->GetPoint(IndexForNewPoints[realationByCell[id][tt[1].first]],point2);
                points->GetPoint(IndexForNewPoints[realationByCell[id][tt[2].first]],point3);

                vtkTriangle::ComputeNormal(point1, point2, point3, normal);

                double d = point1[0] * normal[0] + point1[1] * normal[1] + point1[2] * normal[2];
                cout<<"dddd "<<d<<endl;
             //   double DD = (normal[0] * virsune1[0] + normal[1] * virsune1[1] + normal[2] * virsune1[2] + d) / vtkMath::Norm(normal);
                if (d < 0) {
                  //  swap(tt[1], tt[2]);
                }

            }*/
            if(tt.size()>2)
            {
            faces->InsertNextCell(tt.size());
            for (unsigned int z = 0; z < tt.size(); z++) {
                 faces->InsertCellPoint(IndexForNewPoints[realationByCell[id][tt[z].first]]);
            }
            copyCellsValues(id, input->GetCellData(), output->GetCellData());
            errorarray->InsertNextTuple1(I_ATITIKT);
            }
        }

    }
    t3.stop();
    output->GetCellData()->AddArray(errorarray);
    output->SetPoints(points);
    output->SetPolys(faces);
    output->Squeeze();
    errorarray->Delete();
    points->Delete();
    faces->Delete();

    t4.stop();
    t0.print("insphere_3D_paruosimas");
    t1.print("insphere_3D_surinkimas_nagrinejamu_jungciu");
    t2.print("insphere_3D_suradimas_trukstamu_tetrahedru_octahedru");
    t3.print("insphere_3D_generavimas_voronoi_faces");
    t4.print("insphere_3D_bendras_laikas");
    return 1;
}
void vtkCellCenter3D::CalcCenter(vector<vispartdem::INT> &ids, vtkPolyData*input, double *center) {

    double p1[3];
    double p2[3];
    double p3[3];
    double p4[3];

    center[0] = 0;
    center[1] = 0;
    center[2] = 0;

    if (ids.size() == 4) {
        input->GetPoint(ids[0], p1);
        input->GetPoint(ids[1], p2);
        input->GetPoint(ids[2], p3);
        input->GetPoint(ids[3], p4);


        center[0]=(p1[0]+p2[0]+p3[0]+p4[0])/4;
        center[1]=(p1[1]+p2[1]+p3[1]+p4[1])/4;
        center[2]=(p1[2]+p2[2]+p3[2]+p4[2])/4;

        // vtkTetra::Insphere(p1, p2, p3, p4, center);
    }
    if (ids.size() == 6) {



        input->GetPoint(ids[0], p1);
        input->GetPoint(ids[1], p2);
        input->GetPoint(ids[2], p3);
        input->GetPoint(ids[3], p4);

        double p5[3],p6[3];
        input->GetPoint(ids[4], p5);
        input->GetPoint(ids[5], p6);
        center[0]=(p1[0]+p2[0]+p3[0]+p4[0]+p5[0]+p6[0])/6;
        center[1]=(p1[1]+p2[1]+p3[1]+p4[1]+p5[1]+p6[1])/6;
        center[2]=(p1[2]+p2[2]+p3[2]+p4[2]+p5[2]+p6[2])/6;

        /*
                 double midpoint1[3];
        double midpoint2[3];
         vtkMath::Add(p1, p3, midpoint1);
        vtkMath::Add(p2, p4, midpoint2);
        vtkMath::MultiplyScalar(midpoint1, 0.5);
        vtkMath::MultiplyScalar(midpoint2, 0.5);
        vtkMath::Add(midpoint1, midpoint2, center);
        vtkMath::MultiplyScalar(center, 0.5);*/

        /*
         double p5[3];
         double p6[3];
         double center1[3];
         double center2[3];
         double center3[3];
         double center4[3];
         input->GetPoint(ids[4], p5);
         input->GetPoint(ids[5], p6);
         vtkTetra::Insphere(p1, p2, p3, p5, center1);
         vtkTetra::Insphere(p1, p2, p4, p5, center2);
         vtkTetra::Insphere(p1, p2, p3, p6, center3);
         vtkTetra::Insphere(p1, p2, p4, p6, center4);
         vtkMath::Add(center1, center2, midpoint1);
         vtkMath::Add(center3, center4, midpoint2);
         vtkMath::MultiplyScalar(midpoint1, 0.5);
         vtkMath::MultiplyScalar(midpoint2, 0.5);
         vtkMath::Add(midpoint1, midpoint2, center);
         vtkMath::MultiplyScalar(center, 0.5);*/
    }

}
int vtkCellCenter3D::FillInputPortInformation(int, vtkInformation *info) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
    return 1;
}

void vtkCellCenter3D::PrintSelf(ostream& os, vtkIndent indent) {
    this->Superclass::PrintSelf(os, indent);
    os << "State array name " << StateArray << endl;
}
inline void vtkCellCenter3D::copyCellsValues(vispartdem::INT &i, vtkCellData*icelldata, vtkCellData * ocelldata) //kopijuouja celldatoje esancias reiksmes
{
    for (int m = 0; m < icelldata->GetNumberOfArrays(); ++m) {
        ocelldata->GetArray(m)->InsertNextTuple(icelldata->GetArray(m)->GetTuple(i));
    }
}
}
}
