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

#ifndef _vtkVoronoi2D_H
#define _vtkVoronoi2D_H

#include "vtkObjectFactory.h"
#include "vtkPolyDataAlgorithm.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkPolyData.h"
#include "vtkCellData.h"
#include "vtkPointData.h"
#include "vtkDataArray.h"
#include "vtkDoubleArray.h"
#include "vtkCellArray.h"
#include "vtkIdList.h"
#include "vtkKdTree.h"
#include "vtkMath.h"
#include "Timeris.h"
#include <algorithm>
#include <cmath>
#include <set>
#include <map>
#include <string>
#include "vtkTriangle.h"

namespace vispartdem {
namespace voronoi2D_LOCAL {

#define PI 3.14159265358979323846;
#define PI180 57.295779513;
#define HALFPI 1.570796327
#define PI32 4.71238898
#define TWOPI 6.283185307
#define R_DELTA 0.0000001
#define I_PAILG 2
#define I_NEATITIKT 1
#define I_ATITIKT 0

#define OUTSIDE 1
#define INSIDE 0
#define ONLINE 2
#define NO_INTERSECTION 0
#define EXIST_INTERSECTION 1
#define ELPS 0.00000001


using namespace std;

struct kaimynukas {
    int kaimyn_id_;
    int cellid_;
    double angle_;
    static bool cmp(const kaimynukas *a, const kaimynukas * b) {
        return a->angle_ < b->angle_;
    }
};

typedef struct MazgoKaimynai {
    vector<kaimynukas*> kaimynai;
} Neighbour;

struct VoroLine {
    int left_;
    int main_;
    int right_;
};

struct ParamLine{
    double p1[3];
    double p2[3];
    double pred1[3];
    double pred2[3];
    int line_id;
    double deltay;
    double deltax;
    double deltapred;
    ParamLine(){
        line_id=-1;
    }

    void init( double *sp,double *dp)
    {
        for(int i=0;i<3;i++)
        {
            p1[i]=sp[i];
            p2[i]=dp[i];
            pred1[i]=(p1[i]+p2[i])*0.5;
            pred2[i]=pred1[i]-p1[i];
        }

        swap(pred2[0],pred2[1]);
        pred2[0]=-pred2[0];
        for(int i=0;i<3;i++){
            pred2[i]=pred1[i]+pred2[i];
        }
        deltay=pred1[1]-pred2[1];
        deltax=pred2[0]-pred1[0];
        deltapred=pred1[0]*pred2[1]-pred2[0]*pred1[1];

    }

    ParamLine(const ParamLine& c)
    {
        for(int i=0;i<3;i++)
        {
            p1[i]=c.p1[i];
            p2[i]=c.p2[i];
            pred1[i]=c.pred1[i];
            pred2[i]=c.pred2[i];
        }
        line_id=c.line_id;
        deltax=c.deltax;
        deltay=c.deltay;
        deltapred=c.deltapred;

    }
    int OutsidePlane(double *p)
    {
        double a=deltay*p[0]+deltax*p[1]+deltapred;
        // double a=(pred1[1]-pred2[1])*p[0]+(pred2[0]-pred1[0])*p[1]+(pred1[0]*pred2[1]-pred2[0]*pred1[1]);
        if(a==0)
        {
            return ONLINE;
        }
        if(a>0)
        {
            return INSIDE;
        }else
        {
            return OUTSIDE;
        }
    }

    int intesection(double *p1,double *p2,double *X)
    {
        double x1,x2,x3,x4;
        double y1,y2,y3,y4;
        x1=pred1[0];
        y1=pred1[1];
        x2=pred2[0];
        y2=pred2[1];
        x3=p1[0];
        y3=p1[1];
        x4=p2[0];
        y4=p2[1];
        double dalyba=(x1-x2)*(y3-y4)-(y1-y2)*(x3-x4);
        if(dalyba<ELPS&&dalyba>-ELPS)
        {
            return NO_INTERSECTION;
        }
        X[0]=((x1*y2-y1*x2)*(x3-x4)-(x1-x2)*(x3*y4-y3*x4))/dalyba;
        X[1]=((x1*y2-y1*x2)*(y3-y4)-(y1-y2)*(x3*y4-y3*x4))/dalyba;
        return EXIST_INTERSECTION;
    }
};
typedef ParamLine myline;


struct node
{
    int p_prev;
    int p_next;
    int line_id;
};
class VTK_GRAPHICS_EXPORT vtkVoronoi2D: public vtkPolyDataAlgorithm {
public:
    vtkTypeRevisionMacro(vtkVoronoi2D, vtkPolyDataAlgorithm)
    static vtkVoronoi2D *New();

    vtkVoronoi2D();
    virtual ~vtkVoronoi2D();vtkSetMacro(radius_,string)
    ;
    ;vtkGetMacro(radius_,string)
    ;vtkSetMacro(StateArray_,string)
    ;vtkGetMacro(StateArray_,string)
    ;vtkSetMacro(checkpailgejima_, bool)
    ;vtkGetMacro(checkpailgejima_, bool)
    ;vtkSetMacro(nuokrypis_, double)
    ;vtkGetMacro(nuokrypis_, double)
    ;vtkSetMacro(Deformations_, bool)
    ;vtkGetMacro(Deformations_, bool)
    ;vtkSetMacro(praddata_,vtkDataSet*)
    ;vtkGetMacro(praddata_,vtkDataSet*)

    void SetResultArrayName(std::string ArrayName);

protected:
    int LocalRun2D();
    inline bool TestBounds(double*bounds, double*p);
    void ParasuomiejiDarbai();
    void CleanUp();
    inline double DISTANCE(double*p1, double*p2);
    bool CircleCenter(double*p1, double*p2, double*p3, double*center);
    void CreateVoronoiCell(int &);
    void UzpildytiKaimynus();
    void IdetiLinija(double*point1, double*point2, int &cell_id, int error_type);
    void Calc2DAngles(int &id, vtkPolyData*input);
    inline bool isElongated(int cellid, double dist);

    inline void copyCellsValues(int &i, vtkCellData*, vtkCellData*);
    virtual int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

    virtual int FillInputPortInformation(int port, vtkInformation *info) {
        info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
        return 1;
    }
    void initBounds(double x1,double x2,double y1,double y2,vector<node> &vorEdges,vector<double*>&vertexes);
    //void GenerateVoronoiCell2D(double*origin, vector<double*> &points,vector<double*> &vertexes,vector<node> &vorEdges);
    void GenerateVoronoiCell2D(int &manoid);
    int CheckPoint(double *origin, double*point);
    std::string radius_;
    bool onetime_;
    bool checkpailgejima_;
    double nuokrypis_;
    bool Deformations_;
    string StateArray_;
    vtkDataSet*praddata_;
    vtkDoubleArray* connectionlengh_;
    vtkDoubleArray* skirtumas_;
    std::string name_;
    vtkKdTree *tree_;
    double *bounds_;
    int*indexavimas_;
    double z_val;
    Neighbour**neighbours;
    vtkIntArray*errorarray;
    vtkCellData*icelldata;
    vtkCellData*ocelldata;
    vtkPolyData*input;
    vtkPoints*points;
    vtkCellArray*voronoij;
    vtkPolyData*output;
    vtkDataArray*state;
    vtkDataArray*radius_array;
    //vector<int> nagrinejamu;
private:

};
}
}
#endif  /* _vtkVoronoi2D_H */
