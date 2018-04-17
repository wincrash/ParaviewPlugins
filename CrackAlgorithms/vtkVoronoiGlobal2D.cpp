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

#include "vtkVoronoiGlobal2D.h"
namespace vispartdem {
namespace voronoi2D_GLOBAL {

vtkCxxRevisionMacro(vtkVoronoiGlobal2D, "$Revision: 1.42 $");
vtkStandardNewMacro(vtkVoronoiGlobal2D);

vtkVoronoiGlobal2D::vtkVoronoiGlobal2D() {
    onetime_ = true;
    v = new VoroCPPManager();
}

vtkVoronoiGlobal2D::~vtkVoronoiGlobal2D() {
}

int vtkVoronoiGlobal2D::GlobalRun2D() {
    Timeris* timer = new Timeris[5];
    if (onetime_) {
        timer[2].start();
        UzpildytiKaimynus ();
        timer[2].stop();
    }
    timer[0].start();
    ParasuomiejiDarbai();
    double kart = 0.000001;
    double x_min = bounds_[0] - kart, x_max = bounds_[1] + kart;
    double y_min = bounds_[2] - kart, y_max = bounds_[3] + kart;
    double z_min = bounds_[4] - kart, z_max = bounds_[5] + kart;
    /* cia konkreciai dvimaciui atvejui,
     *  kadangi voro++ generuoja tik 3d atveji,
     *  tai z asi turim praplesti truput,
     *  o veliau kai bus dedami taskai
     *  bus z kordinate is input->GetPoint(0)[2] pirmo tasko z koordinate
     */
    double z_koord = input->GetPoint(0)[2];
    if (bounds_[4] == bounds_[5]) {
        z_min = -0.001;
        z_max = 0.001;
    }
    // Set up the number of blocks that the container is divided into

    int numberofpoints = input->GetNumberOfPoints();
    timer[0].stop();
    // Create a container with the geometry given above, and make it
    // non-periodic in each of the three coordinates. Allocate space for
    // eight particles within each computational block
    timer[3].start();
    cele3d* faces = NULL;
    if (strcmp(this->radius_.c_str(), "") == 0) {
        faces= v->Voro_Cont_getInfo2D(x_min, x_max, y_min, y_max, z_min, z_max,input);
    } else {
        /*container_poly con(x_min, x_max, y_min, y_max, z_min, z_max, n_x, n_y, n_z, false, false, false, 8);
        vtkDataArray*array = input->GetPointData()->GetArray(this->radius_.c_str());
        faces = new cele3d[numberofpoints];
        for (int i = 0; i < numberofpoints; ++i) {
            con.put(i, input->GetPoint(i)[0], input->GetPoint(i)[1], input->GetPoint(i)[2], array->GetTuple1(i)); //uzpildom voro++ pozicijom
            faces[i].nrfaces = 0;
            faces[i].id = -1;
            faces[i].nrvertext = 0;
        }
        //----------------- padarom outputa----------------//
        con.print_all_custom("%i %w %s %a %n %t %P", faces); //paduodam structura, kuria uzpildo visa info reikalinga voronoi
        */
    }
    //-------- pradedama deliojimus
    int id = 0; //keik yra idetu poziciju
    int manoid = 0; //einamo mazgo id
    double point1[3], point2[3];
    int cell_id=0;
    for (manoid = 0; manoid < numberofpoints; ++manoid) //sukam per visus mazgus
    {
        for (int n = 0; n < faces[manoid].nrfaces; ++n) //sukam per voronoi surastu kaimynus
        {
            for (unsigned int m = 0; m < neighbours[manoid]->kaimynai.size(); ++m) //sukam per to mazgo kaimynus
            {
                cell_id=neighbours[manoid]->kaimynai[m]->cellid_;
                if (neighbours[manoid]->kaimynai[m]->kaimyn_id_ == faces[manoid].kaimynai[n] && indexavimas_[cell_id] == 0) //tikrinam ar kaimynai sutampa, ir ar jau nebevo nubrezta ta tiese
                {
                    indexavimas_[cell_id] = 1;
                    for (int j = 0; j < faces[manoid].facevertextnr[n]; ++j) //sukam per daugiakampio virsunes
                    {
                        if (faces[manoid].vertexes[faces[manoid].faceorder[n][j]][2] > 0) //2d atveju daugiakampis yra keturkampis, tai is jo pasiimame tik virsutine linija( tai yra 2 taskai), o visa kita atmetame
                        {
                            points->InsertNextPoint(faces[manoid].vertexes[faces[manoid].faceorder[n][j]][0], faces[manoid].vertexes[faces[manoid].faceorder[n][j]][1], z_koord); // z koordinate pakeiciame tikraja is 0 realaus tasko
                            id++; //skaiciuojame kiek ideta koordinaciu
                        }
                    }
                    voronoij->InsertNextCell(2); //kuriame linija imdami du paskutinius taskus
                    voronoij->InsertCellPoint(id - 1);
                    voronoij->InsertCellPoint(id - 2);
                    copyCellsValues(cell_id, icelldata, ocelldata);
                    errorarray->InsertNextTuple1(I_ATITIKT);
                }
            }

        }

    }

    timer[3].stop();
    if (Deformations_) {
        timer[4].start();
        vtkIdType* connections = input->GetLines()->GetPointer();
        for(int i=0;i<input->GetNumberOfCells();i++){
            if (indexavimas_[i] == 0) {
                indexavimas_[i] = 2;
                input->GetPoint(connections[i*3+1], point1);
                input->GetPoint(connections[i*3+2], point2);
                IdetiLinija(point1, point2, i, I_NEATITIKT);
            }
            if (checkpailgejima_) {
                input->GetPoint(connections[i*3+1], point1);
                input->GetPoint(connections[i*3+2], point2);
                IdetiLinija(point1, point2, i, I_NEATITIKT);
                if (isElongated(i, DISTANCE(point1, point2))) {
                    IdetiLinija(point1, point2, i, I_PAILG);
                }
            }
        }
        timer[4].stop();
    }
    timer[0].start();
    CleanUp();
    delete[] faces;
    timer[0].stop();
    //cout<< "\n--------========================------------------=================----------\n";
    timer[2].print("Global_2D_Kaimynu_surinkimas");
    timer[0].print("Global_2D_Parasoumieji_darbai");
    timer[1].print("Global_2D_KDTree_medzio_sukurimas");
    timer[3].print("Global_2D_VoronoiGeneration");
    timer[4].print("Global_2D_Deformacijos");
    //cout<< "\n--------========================------------------=================----------\n";
    return 1;
}

inline bool vtkVoronoiGlobal2D::TestBounds(double*bounds, double*p) {
    return (bounds[0] <= p[0] && bounds[1] >= p[0] && bounds[2] <= p[1] && bounds[3] >= p[1] && bounds[4] <= p[2] && bounds[5] >= p[2]);
}

int vtkVoronoiGlobal2D::RequestData(vtkInformation * vtkNotUsed(request), vtkInformationVector **inputVector, vtkInformationVector * outputVector) {
    vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
    vtkInformation *outInfo = outputVector->GetInformationObject(0);
    input = vtkPolyData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
    output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

    return GlobalRun2D();
}

inline void vtkVoronoiGlobal2D::copyCellsValues(int &i, vtkCellData*icelldata, vtkCellData * ocelldata) //kopijuouja celldatoje esancias reiksmes
{
    for (int m = 0; m < icelldata->GetNumberOfArrays(); ++m) {
        ocelldata->GetArray(m)->InsertNextTuple(icelldata->GetArray(m)->GetTuple(i));
    }
}

void vtkVoronoiGlobal2D::SetResultArrayName(std::string ArrayName) {
    this->name_ = ArrayName;
    this->Modified();
}

inline bool vtkVoronoiGlobal2D::isElongated(int cellid, double dist) {
    double temp = dist - connectionlengh_->GetTuple1(cellid);
    return (temp > skirtumas_->GetTuple1(cellid));

}
void vtkVoronoiGlobal2D::UzpildytiKaimynus() {
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

}

void vtkVoronoiGlobal2D::IdetiLinija(double*point1, double*point2, int &cell_id, int error_type) {
    points->InsertNextPoint(point1);
    points->InsertNextPoint(point2);
    voronoij->InsertNextCell(2);
    voronoij->InsertCellPoint(points->GetNumberOfPoints() - 2);
    voronoij->InsertCellPoint(points->GetNumberOfPoints() - 1);
    copyCellsValues(cell_id, icelldata, ocelldata);
    errorarray->InsertNextTuple1(error_type);
}

void vtkVoronoiGlobal2D::ParasuomiejiDarbai() {
    points = vtkPoints::New();
    bounds_ = input->GetBounds();
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

void vtkVoronoiGlobal2D::CleanUp() {
    errorarray->Delete();
    points->Delete();
    voronoij->Delete();
    output->Squeeze();
    delete[] indexavimas_;
    indexavimas_ = NULL;

}

inline double vtkVoronoiGlobal2D::DISTANCE(double*p1, double*p2) //grazina atstuma tarp dvieju tasku
{
    double deltax = p1[0] - p2[0];
    double deltay = p1[1] - p2[1];
    double deltaz = p1[2] - p2[2];
    return sqrt(deltax * deltax + deltay * deltay + deltaz * deltaz);
}
}
}
