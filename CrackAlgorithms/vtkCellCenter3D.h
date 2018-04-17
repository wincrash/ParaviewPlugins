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

#ifndef vtkCellCenter3D_H_
#define vtkCellCenter3D_H_

#include "vtkPolyDataAlgorithm.h"
#include "vtkCellArray.h"
#include "vtkDataArray.h"
#include "vtkDataSet.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkDoubleArray.h"
#include "vtkCellData.h"
#include "vtkMath.h"
#include <vector>
#include <set>

#include "CombinationsUtils.h"
#include <map>
#include "CellSearch.h"
#include "Timeris.h"
#include "vtkTetra.h"
#include "vtkTriangle.h"
#define I_PAILG 2
#define I_NEATITIKT 1
#define I_ATITIKT 0
using namespace std;
namespace vispartdem {
namespace vtkCellCenter3D {
typedef vector<vispartdem::INT> myvector;
typedef pair<pair<pair<vispartdem::INT, vispartdem::INT>, vispartdem::INT>, vispartdem::INT> _tetrahedron;
typedef pair<pair<pair<pair<pair<vispartdem::INT, vispartdem::INT>, vispartdem::INT>, vispartdem::INT>, vispartdem::INT>, vispartdem::INT> _octahedron;
typedef  vector<vispartdem::INT>  VectorArray;
typedef map<vispartdem::INT, vector<vispartdem::INT> > CellToStructure;
typedef struct Structure3D {
	vector<vispartdem::INT> ids;
	Structure3D() {
		ids.reserve(6);
	}
	Structure3D(const Structure3D &c) {
		ids.resize(c.ids.size());
		copy(c.ids.begin(), c.ids.end(), ids.begin());
	}

} MyStruct;

class  vtkCellCenter3D: public vtkPolyDataAlgorithm {
public:
    static vtkCellCenter3D *New();
	void PrintSelf(ostream& os, vtkIndent indent);
	vtkSetMacro(ResultArrayName,string)
	;vtkGetMacro(ResultArrayName,string);
	vtkSetMacro(StateArray,string)
	;vtkGetMacro(StateArray,string)
	;


protected:
    vtkCellCenter3D();
    ~vtkCellCenter3D();
	// Usual data generation method
	virtual int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

	virtual int FillInputPortInformation(int port, vtkInformation *info);

	void GetPyramidID(_tetrahedron &piramid, vispartdem::INT &id1, vispartdem::INT &id2, vispartdem::INT &id3, vispartdem::INT &id4);
	void GetOctahedronID(_octahedron &ket, vispartdem::INT &id1, vispartdem::INT &id2, vispartdem::INT &id3, vispartdem::INT &id4, vispartdem::INT &id5, vispartdem::INT &id6);
private:
	inline void copyCellsValues(vispartdem::INT &i, vtkCellData*, vtkCellData*);
    vtkCellCenter3D(const vtkCellCenter3D&); // Not implemented.
    void operator=(const vtkCellCenter3D&); // Not implemented.
	void CalcCenter(vector<vispartdem::INT> &ids, vtkPolyData*input, double *center);
	string StateArray;
		string ResultArrayName;
	bool onetime;
	vispartdem::CellSearch* cellSearch;
	vector<MyStruct> struktura;

//reikalingi
	vector<vispartdem::INT> indexingP;
	map<_tetrahedron, vispartdem::INT> pyramid;
	map<_octahedron, vispartdem::INT> keturkampiai;
	CellToStructure realationByCell;

};
}
}
#endif /* vtkPolydataDelaunay2D_H_ */
