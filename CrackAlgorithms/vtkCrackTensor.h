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

#ifndef vtkCrackTensor_H_
#define vtkCrackTensor_H_
#include "vtkMultiProcessController.h"
#include "vtkPolyDataAlgorithm.h"
#include "vtkMultiProcessController.h"
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkDataArray.h"
#include "vtkDataSet.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkDoubleArray.h"
#include <vector>
#include <vector>
#include <deque>
#include <math.h>
#include <algorithm>
using namespace std;
namespace vispartdem {
namespace tensors {

//##################################################################################################

class PaliunisPoint {
public:
	float x;
	float y;
	float z;
	PaliunisPoint() :
			x(0), y(0), z(0) {
		;
	}
	;
};

struct PaliunisPointID {
	unsigned int id;
	bool IsInCircle;
};

struct PointWithNeighbours {
	unsigned int PointID;
	vector<PaliunisPointID> Neighbours;
	double TensorElements[4]; // vienos daleles tenzorius su visais jos kaimynais
};

struct CircleData // taplpyka
{
	double TotalTensorElements[4]; //sferos tenzorius
	vector<unsigned int> PointsAbroad; // saugom konteineri taskus kurie yra uz ribos
	vector<PointWithNeighbours> PointsInCircle; // ir kurie sferoje konteineriai
};

class VTK_GRAPHICS_EXPORT vtkCrackTensor: public vtkPolyDataAlgorithm {
public:
	static vtkCrackTensor *New();vtkTypeRevisionMacro(vtkCrackTensor, vtkPolyDataAlgorithm)
	;
	void PrintSelf(ostream& os, vtkIndent indent);
	void SetRadius(double radius) {
		Radius = radius;
	}
	;
	void setForceArray(std::string force) {
		this->force_ = force;
	}
	;
	const char *getForceArray() {
		return force_.c_str();
	}
	;

protected:
	vtkCrackTensor();
	~vtkCrackTensor() {
	}
	;
	// Usual data generation method
	virtual int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
	virtual int FillInputPortInformation(int port, vtkInformation *info);
private:
	vtkCrackTensor(const vtkCrackTensor&); // Not implemented.
	void operator=(const vtkCrackTensor&); // Not implemented.
	double Radius;
	std::string force_;
};
}
}
#endif /* vtkCrackTensor_H_ */
