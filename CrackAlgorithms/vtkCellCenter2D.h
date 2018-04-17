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

#ifndef vtkCellCenter2D_H_
#define vtkCellCenter2D_H_

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
#include <sstream>
#include "vtkCellData.h"
#include <vector>
#include <set>
#include <map>
#include "CellSearch.h"
#include "vtkTriangle.h"
#include "Timeris.h"
#include "CombinationsUtils.h"


#define I_PAILG 2
#define I_NEATITIKT 1
#define I_ATITIKT 0
using namespace std;
namespace vispartdem {
namespace vtkCellCenter2D {
typedef pair<pair<vispartdem::INT, vispartdem::INT>, vispartdem::INT> triangle_index;
typedef map<triangle_index, vispartdem::INT> inCircleMapType;
typedef map<vispartdem::INT, triangle_index> rev_inCircleMapType;
typedef pair<vispartdem::INT, vispartdem::INT> LineIndexVal;
typedef map<vispartdem::INT, LineIndexVal> inCircleLineIndex;
typedef map<vispartdem::INT, vispartdem::INT> simpleMap;


typedef vector<LineIndexVal> linesToTriangles;

class VTK_GRAPHICS_EXPORT vtkCellCenter2D: public vtkPolyDataAlgorithm {
public:
vtkTypeRevisionMacro(vtkCellCenter2D, vtkPolyDataAlgorithm)

	;vtkGetMacro(StateArray,string)
	;vtkSetMacro(StateArray,string)
	;vtkSetMacro(ResultArrayName,string)
	;vtkGetMacro(ResultArrayName,string)
    static vtkCellCenter2D *New();
	void PrintSelf(ostream& os, vtkIndent indent);

protected:
    vtkCellCenter2D();
    ~vtkCellCenter2D();
	// Usual data generation method
	virtual int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
	virtual int FillInputPortInformation(int port, vtkInformation *info);
	inline void sort3(triangle_index& t) {
		if (t.first.first > t.first.second)
			swap(t.first.first, t.first.second);
		if (t.first.first > t.second)
			swap(t.first.first, t.second);
		if (t.first.second > t.second)
			swap(t.first.second, t.second);
	}

private:
	void TriangleCenter(double*p1, double *p2, double*p3, double *center);
	inline void copyCellsValues(vispartdem::INT &i, vtkCellData*, vtkCellData*);
    vtkCellCenter2D(const vtkCellCenter2D&); // Not implemented.
    void operator=(const vtkCellCenter2D&); // Not implemented.
	string StateArray;
	string ResultArrayName;
	bool onetime;
	vispartdem::INT max_id;
	///
	triangle_index triag;
	inCircleMapType triangMap;
	rev_inCircleMapType revtriangMap;
	inCircleLineIndex linesIndex;
	LineIndexVal lIndex;
    linesToTriangles Ltotriag;
	//
	vector<vispartdem::INT> indexing;
	vispartdem::CellSearch* cellSearch;





};
}
}
#endif /* vtkPolydataDelaunay2D_H_ */
