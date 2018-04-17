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

#include "vtkVoronoi2DWithIncircle.h"

namespace vispartdem {
namespace voronoi2D_LOCAL_INCIRCLE {
vtkCxxRevisionMacro(vtkVoronoi2DWithIncircle, "$Revision: 1.42 $");
vtkStandardNewMacro(vtkVoronoi2DWithIncircle);

vtkVoronoi2DWithIncircle::vtkVoronoi2DWithIncircle() {
    onetime_ = true;
}

vtkVoronoi2DWithIncircle::~vtkVoronoi2DWithIncircle() {
}

int vtkVoronoi2DWithIncircle::LocalRun2D() {
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

                GenerateVoronoiCell2D(id);

                check_points[id] = 1;
            }
            id = connections[i * 3 + 2];
            if (check_points[id] == 0) {

                GenerateVoronoiCell2D(id);

                check_points[id] = 1;
            }
        }
    }
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
                } else {
                    //    indexavimas_[i] = -1;
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
    //cout<< "\n--------========================------------------=================----------\n";

    timer[2].print("LocalNEW_2D_Kaimynu_surinkimas");
    timer[0].print("LocalNEW_2D_Parasoumieji_darbai");
    timer[1].print("LocalNEW_2D_KDTree_medzio_sukurimas");
    timer[3].print("LocalNEW_2D_VoronoiGeneration");
    timer[4].print("LocalNEW_2D_Deformacijos");

    //cout<< "\n--------========================------------------=================----------\n";
    delete[] timer;
    return 1;
}

void vtkVoronoi2DWithIncircle::GetIncircleStruct(int main_point,vector<InCircleStruct> &struktura)
{
    for(int ii=0;ii<neighbours[main_point]->kaimynai.size();ii++){
        InCircleStruct st;
        st.cell_id=neighbours[main_point]->kaimynai[ii]->cellid_;
        st.id_main=main_point;
        st.id1=neighbours[main_point]->kaimynai[ii]->kaimyn_id_;
        st.id2=-1;
        for(int kk=0;kk<neighbours[main_point]->kaimynai.size();kk++){
            for(int zz=0;zz<neighbours[st.id1]->kaimynai.size();zz++){
                if(neighbours[st.id1]->kaimynai[zz]->kaimyn_id_!=main_point)
                {
                    if(neighbours[st.id1]->kaimynai[zz]->kaimyn_id_==neighbours[main_point]->kaimynai[kk]->kaimyn_id_)
                    {
                        if(st.id2==-1){
                            st.id2=neighbours[st.id1]->kaimynai[zz]->kaimyn_id_;
                        }else
                        {
                            st.id3=neighbours[st.id1]->kaimynai[zz]->kaimyn_id_;
                        }
                    }
                }
            }
        }
        struktura.push_back(st);
    }
}

inline bool vtkVoronoi2DWithIncircle::TestBounds(double*bounds, double*p) {
    return (bounds[0] <= p[0] && bounds[1] >= p[0] && bounds[2] <= p[1] && bounds[3] >= p[1] && bounds[4] <= p[2] && bounds[5] >= p[2]);
}


bool vtkVoronoi2DWithIncircle::CircleCenter(double*p1, double*p2, double*p3, double*center) {
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

void vtkVoronoi2DWithIncircle::initBounds(double x1,double x2,double y1,double y2,vector<node> &vorEdges,vector<double*>&vertexes){
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


void vtkVoronoi2DWithIncircle::GenerateVoronoiCell2D(int &manoid)
{
    vtkIdType* connections = input->GetLines()->GetPointer();
    vector<InCircleStruct> incirclest;
    GetIncircleStruct(manoid,incirclest);
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

    cout<<"================================================\n";
    //int p_prev;
    //int p_next;
    //int line_id;
    //  int line_id_prev;
    //  int line_id_prev_vect_edges;
    //  int line_id_next;
    // int line_id_next_vect_edges;
    for(int t=0;t<vorEdges.size();t++)
    {
        vorEdges[t].line_id_next_vect_edges=-1;
        vorEdges[t].line_id_prev_vect_edges=-1;
        for(int y=0;y<vorEdges.size();y++)
        {
            if(y==t)
            {
                continue;
            }
            if(vorEdges[t].p_prev==vorEdges[y].p_prev ||vorEdges[t].p_prev==vorEdges[y].p_next)
            {
                vorEdges[t].line_id_prev=vorEdges[y].line_id;
                vorEdges[t].line_id_prev_vect_edges=y;
            }
            if(vorEdges[t].p_next==vorEdges[y].p_prev ||vorEdges[t].p_next==vorEdges[y].p_next)
            {
                vorEdges[t].line_id_next=vorEdges[y].line_id;
                vorEdges[t].line_id_next_vect_edges=y;
            }

        }
    }
    for(int t=0;t<vorEdges.size();t++)
    {
        cout<<vorEdges[t].p_prev<<" aa "<<vorEdges[t].p_next<<" ============= "<<vorEdges[t].line_id_prev_vect_edges<<"  n "<< vorEdges[t].line_id_next_vect_edges<<endl;
    }


    cout<<"================================================222222222\n";


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
                        indexavimas_[vorEdges[k].line_id] = 3;//////////////neat
                        for(int h=0;h<incirclest.size();h++)
                        {
                            if(incirclest[h].cell_id==vorEdges[k].line_id)
                            {
                                double center1[3],center2[3],pp1[3],pp2[3],pp3[3],pp4[3];
                                input->GetPoint(incirclest[h].id_main,pp1);
                                input->GetPoint(incirclest[h].id1,pp2);
                                input->GetPoint(incirclest[h].id2,pp3);
                                input->GetPoint(incirclest[h].id3,pp4);
                                TriangleCenter(pp1,pp2,pp3,center1);
                                TriangleCenter(pp1,pp2,pp4,center2);
                                IdetiLinija(center1, center2,vorEdges[k].line_id, I_NEATITIKT);
                            }
                        }
                    }
                }else
                {
                    if(vorEdges[k].line_id>=0&&g>=0&&gn>=0)//cia kai isjungtos deformacijos
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
int  vtkVoronoi2DWithIncircle::CheckPoint(double *origin,double*point){
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
void vtkVoronoi2DWithIncircle::TriangleCenter(double*p1, double *p2, double*p3, double *center) {
    double a=vtkMath::Distance2BetweenPoints(p2,p3);
    double b=vtkMath::Distance2BetweenPoints(p1,p3);
    double c=vtkMath::Distance2BetweenPoints(p1,p2);
    double P=a+b+c;
    center[0]=a*p1[0]+b*p2[0]+c*p3[0];
    center[1]=a*p1[1]+b*p2[1]+c*p3[1];
    center[2]=a*p1[2]+b*p2[2]+c*p3[2];
    center[0]=center[0]/P;
    center[1]=center[1]/P;
    center[2]=center[2]/P;
    /*center[0]=(p1[0]+p2[0]+p3[0])/3;
    center[1]=(p1[1]+p2[1]+p3[1])/3;
   center[2]=(p1[2]+p2[2]+p3[2])/3;*/
    //    vtkTriangle::TriangleCenter(p1, p2, p3, center);
}
int vtkVoronoi2DWithIncircle::RequestData(vtkInformation * vtkNotUsed(request), vtkInformationVector **inputVector, vtkInformationVector * outputVector) {
    vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
    vtkInformation *outInfo = outputVector->GetInformationObject(0);
    input = vtkPolyData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
    output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

    return LocalRun2D();
}

inline void vtkVoronoi2DWithIncircle::copyCellsValues(int &i, vtkCellData*icelldata, vtkCellData * ocelldata) //kopijuouja celldatoje esancias reiksmes
{
    for (int m = 0; m < icelldata->GetNumberOfArrays(); ++m) {
        ocelldata->GetArray(m)->InsertNextTuple(icelldata->GetArray(m)->GetTuple(i));
    }
}

void vtkVoronoi2DWithIncircle::SetResultArrayName(std::string ArrayName) {
    this->name_ = ArrayName;
    this->Modified();
}

inline bool vtkVoronoi2DWithIncircle::isElongated(int cellid, double dist) {
    double temp = dist - connectionlengh_->GetTuple1(cellid);
    return (temp > skirtumas_->GetTuple1(cellid));

}

void vtkVoronoi2DWithIncircle::UzpildytiKaimynus() {
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

void vtkVoronoi2DWithIncircle::IdetiLinija(double*point1, double*point2, int &cell_id, int error_type) {
    points->InsertNextPoint(point1);
    points->InsertNextPoint(point2);
    voronoij->InsertNextCell(2);
    voronoij->InsertCellPoint(points->GetNumberOfPoints() - 2);
    voronoij->InsertCellPoint(points->GetNumberOfPoints() - 1);
    copyCellsValues(cell_id, icelldata, ocelldata);
    errorarray->InsertNextTuple1(error_type);
}

void vtkVoronoi2DWithIncircle::ParasuomiejiDarbai() {
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

void vtkVoronoi2DWithIncircle::CleanUp() {
    errorarray->Delete();
    points->Delete();
    voronoij->Delete();
    output->Squeeze();
    delete[] indexavimas_;
    indexavimas_ = NULL;

}

inline double vtkVoronoi2DWithIncircle::DISTANCE(double*p1, double*p2) //grazina atstuma tarp dvieju tasku
{
    double deltax = p1[0] - p2[0];
    double deltay = p1[1] - p2[1];
    double deltaz = p1[2] - p2[2];
    return sqrt(deltax * deltax + deltay * deltay + deltaz * deltaz);
}
}
}
