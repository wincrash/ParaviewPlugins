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
#ifndef vtkCellBased3D_H_
#define vtkCellBased3D_H_

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
#include <map>
#include <algorithm>
#include "vtkUnstructuredGridAlgorithm.h"
#include "vtkUnstructuredGrid.h"
#include "vtkPyramid.h"
#include "vtkMath.h"
#include "vtkTriangle.h"
#include "vtkTetra.h"
#include "vtkAppendFilter.h"
#include "ReformCellStructure.h"
#include "Timeris.h"
#include "vtkExtractCells.h"
using namespace std;
namespace vispartdem {
namespace cellbased3D {

class VTK_GRAPHICS_EXPORT vtkCellBased3D: public vtkUnstructuredGridAlgorithm {
public:
    static vtkCellBased3D *New();vtkTypeRevisionMacro(vtkCellBased3D, vtkUnstructuredGridAlgorithm)
	void PrintSelf(ostream& os, vtkIndent indent);
	void SetStateArray(std::string ResultArrayName);
	const char *getStateArray();
protected:
    vtkCellBased3D();
    ~vtkCellBased3D();
	void CreateTetraUnstructuredGrid(vtkUnstructuredGrid *output, vtkPolyData* input, vtkDataArray*state);
	void CreatePyramidUnstructuredGrid(vtkUnstructuredGrid *output, vtkPolyData* input, vtkDataArray*state);
	// Usual data generation method
	virtual int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
	virtual int FillInputPortInformation(int port, vtkInformation *info);

private:
    vtkCellBased3D(const vtkCellBased3D&); // Not implemented.
    void operator=(const vtkCellBased3D&); // Not implemented.
	string StateArray;
	bool onetime;
	vispartdem::ReformCellStructure structura3D;
	vector<vector<vispartdem::INT> > cellsIDTetra;
	vector<vector<vispartdem::INT> > cellsIDPyra;
	vector<vector<vispartdem::INT> > IDTetra;
	vector<vector<vispartdem::INT> > IDPyra;
    vtkUnstructuredGrid*cells;
};
}
}
#endif
