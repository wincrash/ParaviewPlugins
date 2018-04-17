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

#include "vtkVoronoi3D.h"
namespace vispartdem {
namespace voronoi3D_LOCAL {

vtkStandardNewMacro(vtkVoronoi3D);

vtkVoronoi3D::vtkVoronoi3D() {
    onetime_ = true;
    v = new VoroCPPManager();
}

vtkVoronoi3D::~vtkVoronoi3D() {
}

int vtkVoronoi3D::LocalRun3D() {
    Timeris* timer = new Timeris[5];
    if (onetime_) {
        timer[2].start();
        UzpildytiKaimynus();
        timer[2].stop();
        tree_=vtkKdTree::New();
    }

    timer[1].start();
    if(Deformations_)
    {
        tree_=vtkKdTree::New();
        tree_->BuildLocatorFromPoints(input->GetPoints());
    }
    timer[1].stop();
    timer[0].start();
    ParasuomiejiDarbai();

    int * check_points = new int[input->GetNumberOfPoints()];
    for (int i = 0; i < input->GetNumberOfPoints(); ++i) {
        check_points[i] = 0;
    }
    vtkIdType* connections = input->GetLines()->GetPointer();
    int id;
    timer[0].stop();
    timer[3].start();
    for (int i = 0; i < input->GetNumberOfLines(); ++i) {
        if (state->GetTuple1(i) != 0) {
            id = connections[i * 3 + 1];
            if (check_points[id] == 0) {
                CreateVoronoiCell3D(id);
                check_points[id] = 1;
            }
            id = connections[i * 3 + 2];
            if (check_points[id] == 0) {
                CreateVoronoiCell3D(id);
                check_points[id] = 1;
            }
        }
    }
    timer[3].stop();
    //int kk = 0;
    //int kkk = 0;
    timer[4].start();
    if (Deformations_) {
        double point1[3];
        double point2[3];
        for (int i = 0; i < input->GetNumberOfCells(); ++i) {
           // if (state->GetTuple1(i) == 0) {
           //     kkk++;
           // }
            if (indexavimas_[i] <= -1) {
                //kk++;
                input->GetPoint(connections[i * 3 + 1], point1);
                input->GetPoint(connections[i * 3 + 2], point2);
                if (checkpailgejima_) {
                    if (isElongated(i, DISTANCE(point1, point2))) {
                        IdetiLinija(point1, point2, i, I_PAILG);
                    }
                }
                IdetiLinija(point1, point2, i, I_NEATITIKT);
            }
        }

    }
    timer[4].stop();


    timer[0].start();
    CleanUp();
    delete[] check_points;
    timer[0].stop();
    //cout<< "\n--------========================------------------=================----------\n";
   // cout<<"Deformaciju_kiekis "<<kk<<endl;
    timer[2].print("Local_3D_Kaimynu_surinkimas");
    timer[0].print("Local_3D_Parasoumieji_darbai");
    timer[1].print("Local_3D_KDTree_medzio_sukurimas");
    timer[3].print("Local_3D_VoronoiGeneration");
    timer[4].print("Local_3D_Deformacijos");
    //cout<< "\n--------========================------------------=================----------\n";
    return 1;
}

inline bool vtkVoronoi3D::TestBounds(double*bounds, double*p) {
    return (bounds[0] <= p[0] && bounds[1] >= p[0] && bounds[2] <= p[1] && bounds[3] >= p[1] && bounds[4] <= p[2] && bounds[5] >= p[2]);
}

int vtkVoronoi3D::RequestData(vtkInformation * vtkNotUsed(request), vtkInformationVector **inputVector, vtkInformationVector * outputVector) {
    vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
    vtkInformation *outInfo = outputVector->GetInformationObject(0);
    input = vtkPolyData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
    output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

    return LocalRun3D();
}

inline void vtkVoronoi3D::copyCellsValues(int &i, vtkCellData*icelldata, vtkCellData * ocelldata) //kopijuouja celldatoje esancias reiksmes
{
    for (int m = 0; m < icelldata->GetNumberOfArrays(); ++m) {
        ocelldata->GetArray(m)->InsertNextTuple(icelldata->GetArray(m)->GetTuple(i));
    }
}

void vtkVoronoi3D::SetResultArrayName(std::string ArrayName) {
    this->name_ = ArrayName;
    this->Modified();
}

inline bool vtkVoronoi3D::isElongated(int cellid, double dist) {
    double temp = dist - connectionlengh_->GetTuple1(cellid);
    return (temp > skirtumas_->GetTuple1(cellid));

}
void vtkVoronoi3D::CreateVoronoiCell3D(int &mano_id) {
    set<int>::iterator it;
    vector<int> kaimynai;
    vector<int> face_kiekis;
    double point1[3];
    double point2[3];
    input->GetPoint(mano_id, point1);
    v->Voro_Cell_init(bounds_[0] - point1[0], bounds_[1] - point1[0], bounds_[2] - point1[1], bounds_[3] - point1[1], bounds_[4] - point1[2], bounds_[5] - point1[2]);
    for (unsigned int i = 0; i < neighbours[mano_id]->kaimynai.size(); ++i) {
        if (indexavimas_[neighbours[mano_id]->kaimynai[i]->cellid_] == 0) {
            indexavimas_[neighbours[mano_id]->kaimynai[i]->cellid_] = -1;
        }
        input->GetPoint(neighbours[mano_id]->kaimynai[i]->kaimyn_id_, point2);
        if (radius_array == NULL) {
            v->nplane(point2[0] - point1[0], point2[1] - point1[1], point2[2] - point1[2], i);
        } else {
            /*     v->nplane(
             point2[0] - point1[0],
             point2[1] - point1[1],
             point2[2] - point1[2],
             (point2[0] - point1[0]) * (point2[0] - point1[0])
             + (point2[1] - point1[1]) * (point2[1] - point1[1])
             + (point2[2] - point1[2]) * (point2[2] - point1[2])
             + radius_array->GetTuple1(
             neighbours[mano_id]->kaimynai[i]->kaimyn_id_)
             * radius_array->GetTuple1(
             neighbours[mano_id]->kaimynai[i]->kaimyn_id_)
             - radius_array->GetTuple1(
             neighbours[mano_id]->kaimynai[i]->kaimyn_id_)
             * radius_array->GetTuple1(
             neighbours[mano_id]->kaimynai[i]->kaimyn_id_),
             i);*/
        }
    }
    if(Deformations_){

        for (it = neighbours[mano_id]->pirmo_lygio_kaimynai.begin(); it != neighbours[mano_id]->pirmo_lygio_kaimynai.end(); it++) {
            input->GetPoint(*it, point2);
            if (radius_array == NULL) {
                v->nplane(point2[0] - point1[0], point2[1] - point1[1], point2[2] - point1[2], -500);
            } else {

            }
        }

    }
    v->output_neighbors(kaimynai);
    v->output_face_orders(face_kiekis);
    int **indexes = new int*[face_kiekis.size()];
    for (unsigned int i = 0; i < face_kiekis.size(); ++i) {
        indexes[i] = new int[face_kiekis[i]];
    }
    v->output_face_vertices(indexes);
    int *geros_virsunes = new int[v->Voro_cell_p()];
    if(!Deformations_){
        for (int t = 0; t < v->Voro_cell_p(); t++) {
            point2[0] = point1[0] + v->Voro_pts()[t * 3] * 0.5;
            point2[1] = point1[1] + v->Voro_pts()[t * 3 + 1] * 0.5;
            point2[2] = point1[2] + v->Voro_pts()[t * 3 + 2] * 0.5;
            geros_virsunes[t] = points->GetNumberOfPoints();
            points->InsertNextPoint(point2);
        }
    }
    else
    {
        int tempID;
        double x,y,z;
        double point3[3];
        bool yra;
        vtkIdList*list1 = vtkIdList::New();
        for (int t = 0; t < v->Voro_cell_p(); t++) {
            tempID = t * 3;
            x = v->Voro_pts()[tempID] * 0.5;
            y = v->Voro_pts()[tempID + 1] * 0.5;
            z = v->Voro_pts()[tempID + 2] * 0.5;            
            point3[0] = point1[0] + x;
            point3[1] = point1[1] + y;
            point3[2] = point1[2] + z;
            yra = false;
            double dist = sqrt(x * x + y * y + z * z);
            tree_->FindPointsWithinRadius(dist, point3, list1);
            if (list1->GetNumberOfIds() > 5) {
                yra = true;
            }
            if (yra == false) {
                geros_virsunes[t]=points->GetNumberOfPoints();
                points->InsertNextPoint(point3);
            } else {
                geros_virsunes[t]=-1;
            }

        }
    }


    for (unsigned int i = 0; i < face_kiekis.size(); ++i) {
        if (kaimynai[i] < 0) {
            continue;
        }
        if (indexavimas_[neighbours[mano_id]->kaimynai[kaimynai[i]]->cellid_] == -1) {
            int blogai=0;
            if(Deformations_){
                for (int k = 0; (k < face_kiekis[i])&&blogai==0; ++k) {
                    if(geros_virsunes[indexes[i][k]]==-1){
                        blogai=1;

                    }
                }
            }
            if(blogai==0)
            {
                voronoij->InsertNextCell(face_kiekis[i]);
                for (int k = 0; k < face_kiekis[i]; ++k) {
                    voronoij->InsertCellPoint(geros_virsunes[indexes[i][k]]);
                }
                this->copyCellsValues(neighbours[mano_id]->kaimynai[kaimynai[i]]->cellid_, icelldata, ocelldata);
                errorarray->InsertNextTuple1(I_ATITIKT);
                indexavimas_[neighbours[mano_id]->kaimynai[kaimynai[i]]->cellid_] = 1;
                if (Deformations_ && checkpailgejima_) {
                    input->GetPoint(mano_id, point1);
                    input->GetPoint(neighbours[mano_id]->kaimynai[kaimynai[i]]->kaimyn_id_, point2);
                    if (isElongated(neighbours[mano_id]->kaimynai[kaimynai[i]]->cellid_, DISTANCE(point1, point2))) {
                        IdetiLinija(point1, point2, neighbours[mano_id]->kaimynai[kaimynai[i]]->cellid_, I_PAILG);
                    }
                }
            }else
            {
                indexavimas_[neighbours[mano_id]->kaimynai[kaimynai[i]]->cellid_] = -2;
            }
        }
        delete[] indexes[i];
    }

    delete[] geros_virsunes;

    delete[] indexes;
    while (!face_kiekis.empty()) {
        face_kiekis.pop_back();
    }
    while (!kaimynai.empty()) {
        kaimynai.pop_back();
    }
    face_kiekis.clear();
    kaimynai.clear();
}
void vtkVoronoi3D::UzpildytiKaimynus() {
    neighbours = new Neighbour*[input->GetNumberOfPoints()];
    for (int i = 0; i < input->GetNumberOfPoints(); ++i) {
        neighbours[i] = new Neighbour;
    }

    vtkIdType* connections = input->GetLines()->GetPointer();
    int id1, id2;
    for (int i = 0; i < input->GetNumberOfLines(); ++i) {
        id1 = connections[i * 3 + 1];
        id2 = connections[i * 3 + 2];
        neighbours[id1]->kaimynai.push_back(new kaimynukas);
        neighbours[id1]->kaimynai.back()->cellid_ = i;
        neighbours[id1]->kaimynai.back()->kaimyn_id_ = id2;
        neighbours[id2]->kaimynai.push_back(new kaimynukas);
        neighbours[id2]->kaimynai.back()->cellid_ = i;
        neighbours[id2]->kaimynai.back()->kaimyn_id_ = id1;
    }

    for (int i = 0; i < input->GetNumberOfPoints(); ++i) {
        for (unsigned int k = 0; k < neighbours[i]->kaimynai.size(); ++k) {
            for (unsigned int z = 0; z < neighbours[neighbours[i]->kaimynai[k]->kaimyn_id_]->kaimynai.size(); ++z) {
                neighbours[i]->pirmo_lygio_kaimynai.insert(neighbours[neighbours[i]->kaimynai[k]->kaimyn_id_]->kaimynai[z]->kaimyn_id_);
            }
        }
        for (unsigned int k = 0; k < neighbours[i]->kaimynai.size(); ++k) {
            neighbours[i]->pirmo_lygio_kaimynai.erase(neighbours[i]->kaimynai[k]->kaimyn_id_);
        }
        neighbours[i]->pirmo_lygio_kaimynai.erase(i);
    }

}

void vtkVoronoi3D::IdetiLinija(double*point1, double*point2, int &cell_id, int error_type) {
    points->InsertNextPoint(point1);
    points->InsertNextPoint(point2);
    voronoij->InsertNextCell(2);
    voronoij->InsertCellPoint(points->GetNumberOfPoints() - 2);
    voronoij->InsertCellPoint(points->GetNumberOfPoints() - 1);
    copyCellsValues(cell_id, icelldata, ocelldata);
    errorarray->InsertNextTuple1(error_type);
}

void vtkVoronoi3D::ParasuomiejiDarbai() {
    points = vtkPoints::New();
    bounds_ = input->GetBounds();
    bounds_[0]=bounds_[0]+0.0001;
    bounds_[1]=bounds_[1]+0.0001;
    bounds_[2]=bounds_[2]+0.0001;

    bounds_[3]=bounds_[3]-0.0001;
    bounds_[4]=bounds_[4]-0.0001;
    bounds_[5]=bounds_[5]-0.0001;
    voronoij = vtkCellArray::New();
    icelldata = input->GetCellData();
    ocelldata = output->GetCellData();
    state = input->GetCellData()->GetArray(StateArray_.c_str()); //sclaras pagal kuri ziuresim ar nutruke ar ne
    errorarray = vtkIntArray::New();
    errorarray->SetNumberOfComponents(1);
    errorarray->SetName(name_.c_str());
    if (onetime_ && checkpailgejima_) {
        connectionlengh_ = vtkDoubleArray::New();
        connectionlengh_->SetNumberOfComponents(1);
        skirtumas_ = vtkDoubleArray::New();
        skirtumas_->SetNumberOfComponents(1);
        skirtumas_->SetNumberOfTuples(praddata_->GetNumberOfCells());
        connectionlengh_->SetNumberOfTuples(praddata_->GetNumberOfCells());
        double p1[3];
        double p2[3];
        for (int i = 0; i < praddata_->GetNumberOfCells(); ++i) {
            praddata_->GetPoint(praddata_->GetCell(i)->GetPointId(0), p1);
            praddata_->GetPoint(praddata_->GetCell(i)->GetPointId(1), p2);
            connectionlengh_->SetTuple1(i, DISTANCE(p1, p2));
            skirtumas_->SetTuple1(i, (nuokrypis_) * connectionlengh_->GetTuple1(i));
        }
    }
    onetime_ = false;
    vtkDoubleArray**arrays = new vtkDoubleArray*[icelldata->GetNumberOfArrays()];
    for (int i = 0; i < icelldata->GetNumberOfArrays(); ++i) //sukuriam arrays i kuriuos kopijuosim celldata reiksmes
    {
        arrays[i] = vtkDoubleArray::New();
        arrays[i]->SetName(icelldata->GetArray(i)->GetName());
        ocelldata->AddArray(arrays[i]);
        arrays[i]->Delete();
    }
    delete[] arrays;
    arrays = NULL;
    indexavimas_ = new int[input->GetNumberOfCells()];
    for (int i = 0; i < input->GetNumberOfCells(); ++i) {
        indexavimas_[i] = 0;
    }
    radius_array = NULL;
    if (strcmp(this->radius_.c_str(), "") != 0) {
        radius_array = input->GetPointData()->GetArray(this->radius_.c_str());
    }
    ocelldata->AddArray(errorarray);
    output->SetPoints(points);
    output->SetPolys(voronoij);

}

void vtkVoronoi3D::CleanUp() {
    errorarray->Delete();
    points->Delete();
    voronoij->Delete();
    output->Squeeze();
    delete[] indexavimas_;
    indexavimas_ = NULL;

}

inline double vtkVoronoi3D::DISTANCE(double*p1, double*p2) //grazina atstuma tarp dvieju tasku
{
    double deltax = p1[0] - p2[0];
    double deltay = p1[1] - p2[1];
    double deltaz = p1[2] - p2[2];
    return sqrt(deltax * deltax + deltay * deltay + deltaz * deltaz);
}
}
}
