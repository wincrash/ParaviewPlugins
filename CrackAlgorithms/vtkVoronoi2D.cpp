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

#include "vtkVoronoi2D.h"

namespace vispartdem {
namespace voronoi2D_LOCAL {
vtkCxxRevisionMacro(vtkVoronoi2D, "$Revision: 1.42 $");
vtkStandardNewMacro(vtkVoronoi2D);

vtkVoronoi2D::vtkVoronoi2D() {
    onetime_ = true;
}

vtkVoronoi2D::~vtkVoronoi2D() {
}

int vtkVoronoi2D::LocalRun2D() {
    Timeris* timer = new Timeris[5];
    if (onetime_) {
        timer[2].start();
        UzpildytiKaimynus();
        timer[2].stop();
    }
    timer[0].start();
    ParasuomiejiDarbai();
    vector<int> check_points(input->GetNumberOfPoints(), 0);
    timer[0].stop();
    timer[1].start();
    if(Deformations_){
        tree_ = vtkKdTree::New();
        tree_->BuildLocatorFromPoints(input->GetPoints());
    }
    timer[1].stop();
    timer[3].start();
    vtkIdType* connections = input->GetLines()->GetPointer();
    int id, id1;
    for (int i = 0; i < input->GetNumberOfLines(); ++i) {
        if (state->GetTuple1(i) == 0) {
            id = connections[i * 3 + 1];
            if (check_points[id] == 0) {
#ifdef VORONOI2D_NEW
                GenerateVoronoiCell2D(id);
#else
                CreateVoronoiCell(id);
#endif
                check_points[id] = 1;
            }
            id = connections[i * 3 + 2];
            if (check_points[id] == 0) {
#ifdef VORONOI2D_NEW
                GenerateVoronoiCell2D(id);
#else
                CreateVoronoiCell(id);
#endif
                check_points[id] = 1;
            }
        }
    }
    int kiekdef=0;
    timer[3].stop();
    if (Deformations_) {
        timer[4].start();
        double point1[3];
        double point2[3];
        double point3[3];
        double point4[3];
        int kkk = -1;
        for (int i = 0; i < input->GetNumberOfCells(); ++i) {
            if (indexavimas_[i] == 2) {
                id = connections[i * 3 + 1];
                id1 = connections[i * 3 + 2];
                kkk = -1;
                int kiek = 0;
                for (unsigned int k = 0; k < neighbours[id]->kaimynai.size(); k++) {
                    for (unsigned int z = 0; z < neighbours[id1]->kaimynai.size(); z++) {
                        if (neighbours[id1]->kaimynai[z]->kaimyn_id_ == neighbours[id]->kaimynai[k]->kaimyn_id_) {
                            kkk = neighbours[id1]->kaimynai[z]->kaimyn_id_;
                            kiek++;
                        }
                    }
                }

                if (kkk != -1 && kiek == 1) {

                    input->GetPoint(connections[i * 3 + 1], point1);
                    input->GetPoint(connections[i * 3 + 2], point2);
                    input->GetPoint(kkk, point3);
                    bool r = CircleCenter(point1, point2, point3, point4);
                    if (r) {
                        if (TestBounds(bounds_, point4)) {
                            point3[0] = (point1[0] + point2[0]) / 2;
                            point3[1] = (point1[1] + point2[1]) / 2;
                            point3[2] = (point1[2] + point2[2]) / 2;
                            IdetiLinija(point3, point4, i, I_ATITIKT);
                        }
                    }
                }
            }
            if (checkpailgejima_) {
                input->GetPoint(connections[i * 3 + 1], point1);
                input->GetPoint(connections[i * 3 + 2], point2);
                if (indexavimas_[i] != 0) {
                    if (isElongated(i, DISTANCE(point1, point2))) {
                       IdetiLinija(point1, point2, i, I_PAILG);
                    }
                }
            }
            if (indexavimas_[i] == -1) {
                input->GetPoint(connections[i * 3 + 1], point1);
                input->GetPoint(connections[i * 3 + 2], point2);
                IdetiLinija(point1, point2, i, I_NEATITIKT);
                kiekdef++;
            }
        }
        timer[4].stop();
    }
    timer[0].start();
    CleanUp();
    if(Deformations_){
        tree_->Delete();
    }
    timer[0].stop();
#ifdef VORONOI2D_NEW
    timer[2].print("LocalNEW_2D_Kaimynu_surinkimas");
    timer[0].print("LocalNEW_2D_Parasoumieji_darbai");
    timer[1].print("LocalNEW_2D_KDTree_medzio_sukurimas");
    timer[3].print("LocalNEW_2D_VoronoiGeneration");
    timer[4].print("LocalNEW_2D_Deformacijos");
#else
    timer[2].print("Local_2D_Kaimynu_surinkimas");
    timer[0].print("Local_2D_Parasoumieji_darbai");
    timer[1].print("Local_2D_KDTree_medzio_sukurimas");
    timer[3].print("Local_2D_VoronoiGeneration");
    timer[4].print("Local_2D_Deformacijos");
#endif
    delete[] timer;
    return 1;
}

inline bool vtkVoronoi2D::TestBounds(double*bounds, double*p) {
    return (bounds[0] <= p[0] && bounds[1] >= p[0] && bounds[2] <= p[1] && bounds[3] >= p[1] && bounds[4] <= p[2] && bounds[5] >= p[2]);
}

bool vtkVoronoi2D::CircleCenter(double*p1, double*p2, double*p3, double*center) {
    double s = 0.5 * ((p2[0] - p3[0]) * (p1[0] - p3[0]) - (p2[1] - p3[1]) * (p3[1] - p1[1]));
    double sUnder = (p1[0] - p2[0]) * (p3[1] - p1[1]) - (p2[1] - p1[1]) * (p1[0] - p3[0]);

    if (sUnder == 0)
        return false; //insufficient data to calculate center

    s /= sUnder;
    center[0] = 0.5 * (p1[0] + p2[0]) + s * (p2[1] - p1[1]); // center x coordinate
    center[1] = 0.5 * (p1[1] + p2[1]) + s * (p1[0] - p2[0]); // center y coordinate
    center[2] = (p1[2] + p2[2] + p3[2]) / 3;
    return true;

}





void vtkVoronoi2D::initBounds(double x1,double x2,double y1,double y2,vector<node> &vorEdges,vector<double*>&vertexes){
    double*p1=new double[3];
    double*p2=new double[3];
    double*p3=new double[3];
    double*p4=new double[3];
    p1[0]=x1;
    p1[1]=y1;
    p1[2]=0;

    p2[0]=x2;
    p2[1]=y1;
    p2[2]=0;

    p3[0]=x2;
    p3[1]=y2;
    p3[2]=0;

    p4[0]=x1;
    p4[1]=y2;
    p4[2]=0;

    vertexes.push_back(p1);
    vertexes.push_back(p2);
    vertexes.push_back(p3);
    vertexes.push_back(p4);
    node n;
    n.line_id=-1;
    n.p_prev=0;
    n.p_next=1;
    vorEdges.push_back(n);
    n.line_id=-2;
    n.p_prev=1;
    n.p_next=2;
    vorEdges.push_back(n);
    n.line_id=-3;
    n.p_prev=2;
    n.p_next=3;
    vorEdges.push_back(n);
    n.line_id=-4;
    n.p_prev=3;
    n.p_next=0;
    vorEdges.push_back(n);
}



void vtkVoronoi2D::GenerateVoronoiCell2D(int &manoid)
{
    Neighbour* mano_kaimynas = neighbours[manoid];
    int kiekis = mano_kaimynas->kaimynai.size(); //celiu skaicius
    double origin[3];
    input->GetPoint(manoid,origin);
    vector<double*> vertexes;
    vector<node> vorEdges;

    vector<node> vorEdgesCopy;
    myline lines_main;

    double p[3];

    initBounds(bounds_[0],bounds_[1],bounds_[2],bounds_[3],vorEdges,vertexes);
    int prev_state,next_state;
    node n;
    int buvo1=0;
    int buvo2=0;
    for(int i=0;i<kiekis;i++)
    {
        input->GetPoint(mano_kaimynas->kaimynai[i]->kaimyn_id_,p);
        lines_main.line_id=mano_kaimynas->kaimynai[i]->cellid_;
        lines_main.init(origin,p);
        n.line_id=-1;
        n.p_next=-1;
        n.p_prev=-1;
        buvo1=0;
        buvo2=0;
        for(int current=0;current<vorEdges.size();current++)
        {
            if(vorEdges[current].p_prev==-1||vorEdges[current].p_next==-1)
            {
                //                   continue;
                prev_state=OUTSIDE;
                next_state=OUTSIDE;
            }else
            {
                prev_state=lines_main.OutsidePlane(vertexes[vorEdges[current].p_prev]);
                next_state=lines_main.OutsidePlane(vertexes[vorEdges[current].p_next]);
            }
            if(prev_state==OUTSIDE&&next_state==INSIDE)
            {
                buvo1++;
            }
            if(prev_state==INSIDE&&next_state==OUTSIDE)
            {
                buvo2++;
            }
            if(buvo1>1 || buvo2>1)
            {
                swap(vorEdges[current].p_next,vorEdges[current].p_prev);
                swap(prev_state,next_state);
            }
            if(prev_state==OUTSIDE&&next_state==OUTSIDE)
            {
                continue;
            }else
            {
                n.line_id=lines_main.line_id;
                vorEdgesCopy.push_back(vorEdges[current]);
                if(prev_state==OUTSIDE&&next_state==INSIDE){
                    double *X=new double[3];
                    X[2]=z_val;
                    int ret=lines_main.intesection(vertexes[vorEdges[current].p_prev],vertexes[vorEdges[current].p_next],X);
                    vertexes.push_back(X);
                    n.p_prev=vertexes.size()-1;
                    vorEdgesCopy.back().p_prev=n.p_prev;
                }
                if(prev_state==INSIDE&&next_state==OUTSIDE){
                    double *X=new double[3];
                    X[2]=z_val;
                    int ret=lines_main.intesection(vertexes[vorEdges[current].p_prev],vertexes[vorEdges[current].p_next],X);
                    vertexes.push_back(X);
                    n.p_next=vertexes.size()-1;
                    vorEdgesCopy.back().p_next=n.p_next;
                }
            }
        }
        if(n.line_id!=-1)
        {
            vorEdgesCopy.push_back(n);
        }
        vorEdges.clear();
        for(int k=0;k<vorEdgesCopy.size();k++)
        {
            node n1;
            n1.line_id=vorEdgesCopy[k].line_id;
            n1.p_prev=vorEdgesCopy[k].p_prev;
            n1.p_next=vorEdgesCopy[k].p_next;
            vorEdges.push_back(n1);
        }
        if(!vorEdgesCopy.empty())
        {
            vorEdgesCopy.erase(vorEdgesCopy.begin(),vorEdgesCopy.end());
        }

    }
    vector<int> checkedPoints(vertexes.size(),-1);
    for(int k=0;k<vorEdges.size();k++)
    {
        if(vorEdges[k].line_id>=0)
        {
            int g=vorEdges[k].p_prev;
            int gn=vorEdges[k].p_next;
            if(indexavimas_[vorEdges[k].line_id]==0)
            {
                if(Deformations_)
                {
                    if(checkedPoints[g]==-1)
                    {
                        checkedPoints[g]=CheckPoint(origin,vertexes[g]);
                    }
                    if(checkedPoints[gn]==-1)
                    {
                        checkedPoints[gn]=CheckPoint(origin,vertexes[gn]);
                    }
                    if(checkedPoints[g]==2 &&checkedPoints[gn]==2)
                    {
                        IdetiLinija(vertexes[g], vertexes[gn],vorEdges[k].line_id, I_ATITIKT);
                        indexavimas_[vorEdges[k].line_id] = 3;
                    }else
                    {
                        indexavimas_[vorEdges[k].line_id] = -1;
                    }
                }else
                {
                    if(vorEdges[k].line_id>=0&&g>=0&&gn>=0)
                    {
                        IdetiLinija(vertexes[g], vertexes[gn],vorEdges[k].line_id, I_ATITIKT);
                        indexavimas_[vorEdges[k].line_id] = 3;
                    }
                }
            }
        }
    }
    for(int i=0;i<vertexes.size();i++)
    {
        delete[] vertexes[i];
    }
    vertexes.clear();
    vorEdges.clear();
}
int  vtkVoronoi2D::CheckPoint(double *origin,double*point){
    double radius = DISTANCE(origin, point);
    radius=radius+radius*0.0001;
    vtkIdList*nrpoints=vtkIdList::New();
    tree_->FindPointsWithinRadius(radius, point, nrpoints);


    if (nrpoints->GetNumberOfIds() <= 3) {
        nrpoints->Delete();
        return 2;
    }else{
        nrpoints->Delete();
        return 0;
    }

}

void vtkVoronoi2D::CreateVoronoiCell(int &manoid) {
    Neighbour* mano_kaimynas = neighbours[manoid];
    int kiekis = mano_kaimynas->kaimynai.size(); //celiu skaicius
    Calc2DAngles(manoid, input);
    sort(mano_kaimynas->kaimynai.begin(), mano_kaimynas->kaimynai.end(), kaimynukas::cmp);

    vector<VoroLine*> lines;
    if (kiekis >= 3) {
        for (int i = 1; i < kiekis - 1; ++i) {
            lines.push_back(new VoroLine);
            lines.back()->left_ = i - 1;
            lines.back()->main_ = i;
            lines.back()->right_ = i + 1;
        }

        lines.push_back(new VoroLine);
        lines.back()->left_ = kiekis - 1;
        lines.back()->main_ = 0;
        lines.back()->right_ = 1;

        lines.push_back(new VoroLine);
        lines.back()->left_ = 0;
        lines.back()->main_ = kiekis - 1;
        lines.back()->right_ = kiekis - 2;
    }


    int temp_id;
    bool flag = true;
    bool flag1 = true;
    vtkIdList *nrpoints = vtkIdList::New();
    double point1_[3];
    double point2_[3];
    double point3_[3];
    double point4_[3];
    double center_left[3];
    double center_right[3];
    input->GetPoint(manoid, point1_);
    for (unsigned i = 0; i < lines.size(); ++i) {
        if (indexavimas_[mano_kaimynas->kaimynai[lines.at(i)->main_]->cellid_] == 1) {
            continue;
        }
        if (indexavimas_[mano_kaimynas->kaimynai[lines.at(i)->main_]->cellid_] == -1) {
            continue;
        }
        if (indexavimas_[mano_kaimynas->kaimynai[lines.at(i)->main_]->cellid_] == 0) {
            indexavimas_[mano_kaimynas->kaimynai[lines.at(i)->main_]->cellid_] = 2;
        }
        input->GetPoint(mano_kaimynas->kaimynai[lines.at(i)->main_]->kaimyn_id_, point2_);
        input->GetPoint(mano_kaimynas->kaimynai[lines.at(i)->left_]->kaimyn_id_, point3_);
        input->GetPoint(mano_kaimynas->kaimynai[lines.at(i)->right_]->kaimyn_id_, point4_);

        bool r1 = CircleCenter(point1_, point2_, point3_, center_left);
        bool r2 = CircleCenter(point1_, point2_, point4_, center_right);

        if (r1 == false || r2 == false) {
            continue;
        }

        flag = true;
        flag1 = true;
        if(Deformations_){
            double radius = DISTANCE(center_left, point1_);
            tree_->FindPointsWithinRadius(radius, center_left, nrpoints);

            if (nrpoints->GetNumberOfIds() <= 3) {
                for (int j = 0; j < nrpoints->GetNumberOfIds(); ++j) {
                    temp_id = nrpoints->GetId(j);
                    if ((manoid != temp_id) && (mano_kaimynas->kaimynai[lines.at(i)->main_]->kaimyn_id_ != temp_id) && (mano_kaimynas->kaimynai[lines.at(i)->left_]->kaimyn_id_ != temp_id)) {
                        flag = false;

                        break;
                    }
                }
            } else {
                flag = false;
            }
            radius = DISTANCE(center_right, point1_);
            tree_->FindPointsWithinRadius(radius, center_right, nrpoints);
            flag1 = true;
            if (nrpoints->GetNumberOfIds() <= 3) {
                for (int j = 0; j < nrpoints->GetNumberOfIds(); ++j) {
                    temp_id = nrpoints->GetId(j);
                    if ((manoid != temp_id) && (mano_kaimynas->kaimynai[lines.at(i)->main_]->kaimyn_id_ != temp_id) && (mano_kaimynas->kaimynai[lines.at(i)->right_]->kaimyn_id_ != temp_id)) {
                        flag1 = false;
                        break;
                    }
                }
            } else {
                flag1 = false;
            }

            if (flag == false) {

                int id2 = lines.at(i)->left_;
                if (id2 == 0) {
                    id2 = kiekis - 1;
                } else {
                    --id2;
                }

                input->GetPoint(mano_kaimynas->kaimynai[id2]->kaimyn_id_, point3_);
                bool r1 = CircleCenter(point1_, point2_, point3_, center_left);
                if (r1) {
                    double radius = DISTANCE(center_left, point1_);
                    tree_->FindPointsWithinRadius(radius, center_left, nrpoints);
                    flag = true;
                    if (nrpoints->GetNumberOfIds() <= 3) {
                        for (int j = 0; j < nrpoints->GetNumberOfIds(); ++j) {
                            temp_id = nrpoints->GetId(j);
                            if ((manoid != temp_id) && (mano_kaimynas->kaimynai[lines.at(i)->main_]->kaimyn_id_ != temp_id) && (mano_kaimynas->kaimynai[id2]->kaimyn_id_ != temp_id)) {
                                flag = false;

                                break;
                            }
                        }

                    } else {
                        flag = false;
                    }

                }
            }

            if (flag1 == false) {
                int id2 = lines.at(i)->right_;
                if (id2 == kiekis - 1) {
                    id2 = 0;
                } else {
                    ++id2;
                }
                input->GetPoint(mano_kaimynas->kaimynai[id2]->kaimyn_id_, point4_);
                bool r1 = CircleCenter(point1_, point2_, point4_, center_right);
                if (r1) {
                    double radius = DISTANCE(center_right, point1_);
                    tree_->FindPointsWithinRadius(radius, center_right, nrpoints);
                    flag1 = true;
                    if (nrpoints->GetNumberOfIds() <= 3) {
                        for (int j = 0; j < nrpoints->GetNumberOfIds(); ++j) {
                            temp_id = nrpoints->GetId(j);
                            if ((manoid != temp_id) && (mano_kaimynas->kaimynai[lines.at(i)->main_]->kaimyn_id_ != temp_id) && (mano_kaimynas->kaimynai[id2]->kaimyn_id_ != temp_id)) {
                                flag1 = false;
                                break;
                            }
                        }
                    } else {
                        flag1 = false;
                    }

                }
            }
        }
        if (flag && flag1) {
            if (TestBounds(bounds_, center_left) && TestBounds(bounds_, center_right)) {
                IdetiLinija(center_left, center_right, mano_kaimynas->kaimynai[lines.at(i)->main_]->cellid_, I_ATITIKT);
                indexavimas_[mano_kaimynas->kaimynai[lines.at(i)->main_]->cellid_] = 1;

            } else {

                if (TestBounds(bounds_, center_left)) {
                    point3_[0] = (point1_[0] + point2_[0]) / 2;
                    point3_[1] = (point1_[1] + point2_[1]) / 2;
                    point3_[2] = (point1_[2] + point2_[2]) / 2;
                    IdetiLinija(point3_, center_left, mano_kaimynas->kaimynai[lines.at(i)->main_]->cellid_, I_ATITIKT);
                    indexavimas_[mano_kaimynas->kaimynai[lines.at(i)->main_]->cellid_] = 1;
                }
                if (TestBounds(bounds_, center_right)) {
                    point3_[0] = (point1_[0] + point2_[0]) / 2;
                    point3_[1] = (point1_[1] + point2_[1]) / 2;
                    point3_[2] = (point1_[2] + point2_[2]) / 2;
                    IdetiLinija(point3_, center_right, mano_kaimynas->kaimynai[lines.at(i)->main_]->cellid_, I_ATITIKT);
                    indexavimas_[mano_kaimynas->kaimynai[lines.at(i)->main_]->cellid_] = 1;
                }
            }
        } else {
            if ((!flag && !flag1)) {
                indexavimas_[mano_kaimynas->kaimynai[lines.at(i)->main_]->cellid_] = -1;
            }
        }
    }

    while (!lines.empty()) {
        VoroLine* l = lines.back();
        delete l;
        lines.pop_back();
    }
    nrpoints->Delete();
}

int vtkVoronoi2D::RequestData(vtkInformation * vtkNotUsed(request), vtkInformationVector **inputVector, vtkInformationVector * outputVector) {
    vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
    vtkInformation *outInfo = outputVector->GetInformationObject(0);
    input = vtkPolyData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
    output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

    return LocalRun2D();
}

inline void vtkVoronoi2D::copyCellsValues(int &i, vtkCellData*icelldata, vtkCellData * ocelldata) //kopijuouja celldatoje esancias reiksmes
{
    for (int m = 0; m < icelldata->GetNumberOfArrays(); ++m) {
        ocelldata->GetArray(m)->InsertNextTuple(icelldata->GetArray(m)->GetTuple(i));
    }
}

void vtkVoronoi2D::SetResultArrayName(std::string ArrayName) {
    this->name_ = ArrayName;
    this->Modified();
}

inline bool vtkVoronoi2D::isElongated(int cellid, double dist) {
    double temp = dist - connectionlengh_->GetTuple1(cellid);
    return (temp > skirtumas_->GetTuple1(cellid));

}

void vtkVoronoi2D::UzpildytiKaimynus() {
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

void vtkVoronoi2D::Calc2DAngles(int &id, vtkPolyData * input) {
    double localpoint_[2];
    double point1_[3];
    double point2_[3];
    input->GetPoint(id, point1_);

    int kiekis = neighbours[id]->kaimynai.size();
    for (int i = 0; i < kiekis; ++i) { //suskaiciuoja kampus, kampai nuo 0 iki 360 naudojant polines koordinates
        input->GetPoint(neighbours[id]->kaimynai[i]->kaimyn_id_, point2_);
        localpoint_[0] = point2_[0] - point1_[0];
        localpoint_[1] = point2_[1] - point1_[1];
        if (localpoint_[0] == 0 && localpoint_[1] == 0) {
            neighbours[id]->kaimynai[i]->angle_ = 0;

        }
        if (localpoint_[0] >= 0 && localpoint_[1] >= 0) {
            neighbours[id]->kaimynai[i]->angle_ = atan(localpoint_[1] / localpoint_[0]);
        }

        if (localpoint_[0] >= 0 && localpoint_[1] < 0) {
            neighbours[id]->kaimynai[i]->angle_ = atan(localpoint_[1] / localpoint_[0]) + TWOPI;
        }

        if (localpoint_[0] < 0) {
            neighbours[id]->kaimynai[i]->angle_ = atan(localpoint_[1] / localpoint_[0]) + PI;
        }

        if (localpoint_[0] == 0 && localpoint_[1] > 0) {
            neighbours[id]->kaimynai[i]->angle_ = HALFPI;
        }
        if (localpoint_[0] == 0 && localpoint_[1] < 0) {
            neighbours[id]->kaimynai[i]->angle_ = PI32;
        }
        neighbours[id]->kaimynai[i]->angle_ = neighbours[id]->kaimynai[i]->angle_ * PI180;
    }
}

void vtkVoronoi2D::IdetiLinija(double*point1, double*point2, int &cell_id, int error_type) {
    points->InsertNextPoint(point1);
    points->InsertNextPoint(point2);
    voronoij->InsertNextCell(2);
    voronoij->InsertCellPoint(points->GetNumberOfPoints() - 2);
    voronoij->InsertCellPoint(points->GetNumberOfPoints() - 1);
    copyCellsValues(cell_id, icelldata, ocelldata);
    errorarray->InsertNextTuple1(error_type);
}

void vtkVoronoi2D::ParasuomiejiDarbai() {
    points = vtkPoints::New();
    bounds_ = input->GetBounds();
    z_val=bounds_[4];
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
    output->SetLines(voronoij);

}

void vtkVoronoi2D::CleanUp() {
    errorarray->Delete();
    points->Delete();
    voronoij->Delete();
    output->Squeeze();
    delete[] indexavimas_;
    indexavimas_ = NULL;

}

inline double vtkVoronoi2D::DISTANCE(double*p1, double*p2) //grazina atstuma tarp dvieju tasku
{
    double deltax = p1[0] - p2[0];
    double deltay = p1[1] - p2[1];
    double deltaz = p1[2] - p2[2];
    return sqrt(deltax * deltax + deltay * deltay + deltaz * deltaz);
}
}
}
