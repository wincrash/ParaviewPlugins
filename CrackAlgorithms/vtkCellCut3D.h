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
#ifndef vtkCellCut3D_H_
#define vtkCellCut3D_H_

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
#include <map>
#include <algorithm>
#include "vtkPolyDataAlgorithm.h"
#include "vtkPolyData.h"
#include "vtkPyramid.h"
#include "vtkTriangle.h"
#include "vtkTetra.h"
#include "vtkAppendFilter.h"
#include <bitset>
#include "ReformCellStructure.h"
#include "Timeris.h"
#define MIDPOINT(a,b,ab)\
{\
    ab[0]=(a[0]+b[0])/2;\
    ab[1]=(a[1]+b[1])/2;\
    ab[2]=(a[2]+b[2])/2;\
    }
using namespace std;
namespace vispartdem {
namespace CellCut3D {

#define TETRA_INDEX_COUNT 63
#define TETRA_INDEX_COUNT_IN_CASE 32

#define PYRAMID_INDEX_COUNT 255
#define PYRAMID_INDEX_COUNT_IN_CASE 48

class VTK_GRAPHICS_EXPORT vtkCellCut3D: public vtkPolyDataAlgorithm {
public:
    static vtkCellCut3D *New();vtkTypeRevisionMacro(vtkCellCut3D, vtkPolyDataAlgorithm)
	void PrintSelf(ostream& os, vtkIndent indent);
	void SetStateArray(std::string ResultArrayName);
	const char *getStateArray();

protected:
	void CreateTetraMarching(vtkPolyData* input, vtkPoints*points, vtkCellArray*cellArray, vtkCellArray*cellArrayLines, vtkCellArray*cellArrayVertex, vtkIntArray*newstate, vtkDataArray*state);
	void CreatePyramidMarching(vtkPolyData* input, vtkPoints*points, vtkCellArray*cellArray, vtkCellArray*cellArrayLines, vtkCellArray*cellArrayVertex, vtkIntArray*newstate, vtkDataArray*state);
    vtkCellCut3D();
    ~vtkCellCut3D();
	// Usual data generation method
	virtual int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
	virtual int FillInputPortInformation(int port, vtkInformation *info);
    void sort1(INT& id1,INT &id2){
        if(id1>id2)
        {
            swap(id1,id2);
        }
    }

private:
    int getPointID(vtkPolyData *input,INT* connections,vtkPoints*points,int cellID);
    vtkCellCut3D(const vtkCellCut3D&); // Not implemented.
    void operator=(const vtkCellCut3D&); // Not implemented.
	string StateArray;
	bool onetime;
	Timeris tt;
	vispartdem::ReformCellStructure structura3D;
	static const int iTetraIndexes[TETRA_INDEX_COUNT][TETRA_INDEX_COUNT_IN_CASE];
	static const int iTetraCases[TETRA_INDEX_COUNT];

	static const int iPyramidIndexes[PYRAMID_INDEX_COUNT][PYRAMID_INDEX_COUNT_IN_CASE];
	static const int iPyramidCases[PYRAMID_INDEX_COUNT];

    vector<int> IndexuotiPointai;



    map<pair<int,int>, int> buvo2;
    pair<int,int> du;
    map<int,int> buvo1;
    map<pair<pair<int,int>,int>, int> buvo3;
    pair<pair<int,int>,int> trys;



};
}
}
#endif
