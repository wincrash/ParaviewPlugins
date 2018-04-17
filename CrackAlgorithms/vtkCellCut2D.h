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

#ifndef vtkCellCut2D_H_
#define vtkCellCut2D_H_

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
#include "ReformCellStructure.h"
#include "vtkMath.h"
#include "Timeris.h"
using namespace std;
namespace vispartdem {
namespace CellCut3D {
class VTK_GRAPHICS_EXPORT vtkCellCut2D: public vtkPolyDataAlgorithm {
public:
vtkTypeRevisionMacro(vtkCellCut2D, vtkPolyDataAlgorithm)
    static vtkCellCut2D *New();
	void PrintSelf(ostream& os, vtkIndent indent);
	void SetStateArray(std::string ResultArrayName);
	const char *getStateArray();
protected:
    vtkCellCut2D();
    ~vtkCellCut2D();
	virtual int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
	virtual int FillInputPortInformation(int port, vtkInformation *info);
private:
    vtkCellCut2D(const vtkCellCut2D&); // Not implemented.
    void operator=(const vtkCellCut2D&); // Not implemented.
	string StateArray;
	bool onetime;
	vispartdem::ReformCellStructure structura;
};
}
}
#endif /* vtkPolydataDelaunay2D_H_ */
