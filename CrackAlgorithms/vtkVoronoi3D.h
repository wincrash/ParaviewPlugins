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

#ifndef _vtkVoronoi3D_H
#define _vtkVoronoi3D_H

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
#include "VoroCPPManager.h"


namespace vispartdem {
namespace voronoi3D_LOCAL {
#define I_PAILG 2
#define I_NEATITIKT 1
#define I_ATITIKT 0
#define PI 3.14159265358979323846;
#define PI180 57.295779513;
#define HALFPI 1.570796327
#define PI32 4.71238898
#define TWOPI 6.283185307

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
	set<int> pirmo_lygio_kaimynai;
} Neighbour;

struct VoroLine {
	int left_;
	int main_;
	int right_;
};

class vtkVoronoi3D: public vtkPolyDataAlgorithm {
public:
	static vtkVoronoi3D *New();

	vtkVoronoi3D();
	virtual ~vtkVoronoi3D();vtkSetMacro(radius_,string)
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
	int LocalRun3D();
	inline bool TestBounds(double*bounds, double*p);
	void ParasuomiejiDarbai();
	void CleanUp();
	inline double DISTANCE(double*p1, double*p2);

	void UzpildytiKaimynus();
	void IdetiLinija(double*point1, double*point2, int &cell_id, int error_type);
	inline bool isElongated(int cellid, double dist);

	inline void copyCellsValues(int &i, vtkCellData*, vtkCellData*);
	virtual int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
	void CreateVoronoiCell3D(int &mano_id);
	virtual int FillInputPortInformation(int port, vtkInformation *info) {
		info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
		return 1;
	}

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
	VoroCPPManager *v;

private:

};
}
}
#endif  /* _vtkVoronoi3D_H */
