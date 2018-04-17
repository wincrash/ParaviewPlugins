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

#ifndef vtkCellBased2D_H_
#define vtkCellBased2D_H_

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
#include "Timeris.h"
using namespace std;
namespace vispartdem {
namespace cellbased2D {

class VTK_GRAPHICS_EXPORT vtkCellBased2D: public vtkPolyDataAlgorithm {
public:
vtkTypeRevisionMacro(vtkCellBased2D, vtkPolyDataAlgorithm)
    static vtkCellBased2D *New();
	void PrintSelf(ostream& os, vtkIndent indent);
	void SetStateArray(std::string ResultArrayName);
	const char *getStateArray();
protected:
    vtkCellBased2D();
    ~vtkCellBased2D();
	// Usual data generation method
	virtual int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

	virtual int FillInputPortInformation(int port, vtkInformation *info);

private:
    vtkCellBased2D(const vtkCellBased2D&); // Not implemented.
    void operator=(const vtkCellBased2D&); // Not implemented.
	string StateArray;
	bool onetime;
	vtkCellArray*triang;
	vispartdem::ReformCellStructure structura;
};
}
}
#endif /* vtkPolydataDelaunay2D_H_ */
