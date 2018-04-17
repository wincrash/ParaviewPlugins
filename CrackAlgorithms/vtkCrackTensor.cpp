/*

 #     # ###  #####  ######     #    ######  ####### ######  ####### #     #
 #     #  #  #     # #     #   # #   #     #    #    #     # #       ##   ##
 #     #  #  #       #     #  #   #  #     #    #    #     # #       # # # #
 #     #  #   #####  ######  #     # ######     #    #     # #####   #  #  #
 #   #   #        # #       ####### #   #      #    #     # #       #     #
 # #    #  #     # #       #     # #    #     #    #     # #       #     #
 #    ###  #####  #       #     # #     #    #    ######  ####### #     #


 * Author       : Edgaras Paliunis
 * E-mail       : epaliunis@gmail.com
 * Vendor       : VGTU
 * Home page    : http://lsl.vgtu.lt/vispartdem
 */

#include "vtkCrackTensor.h"

namespace vispartdem {
namespace tensors {
vtkCxxRevisionMacro(vtkCrackTensor, "$Revision: 1.42 $");
vtkStandardNewMacro(vtkCrackTensor);

vtkCrackTensor::vtkCrackTensor() {

}

using namespace std;

//##################################################################################################
// ar dalelyte viduje sferos
bool IsPointInCircle(double circle[], double point[], double radius) {
	// double cx = circle[0];
	// double cy = circle[1];
	// double px = point[0];
	// double py = point[1];
	double first = pow((circle[0] - point[0]), 2);
	double second = pow((circle[1] - point[1]), 2);
	double distance = sqrt(first + second);
	return radius > distance;
}

//#################################################################################################

unsigned int GenerateCircles(double x_min, double x_max, double y_min, double y_max, double z_min, double z_max, unsigned int dimensions, double r, vtkPoints *points) {
	// float test = acos(-1.0/6.0);
	//  float test2 = cos(test);
	float h = cos(acos(-1.0) / 6.0) * 2 * r;
	unsigned int burbuluKiekisXAsyje = (unsigned int) (floor((float) (x_max - x_min) / (r + r * (sqrt(3.0) - 1.0)))); //(unsigned int)(floor((float)(x_max - x_min - 2.0*r)/h));
	unsigned int burbuluKiekisYAsyje = (unsigned int) (floor((float) (y_max - y_min) / (2 * r))); // (unsigned int)(floor((float)(y_max - y_min)/(2*r)));
	unsigned int burbuluKiekisZAsyje = (unsigned int) (floor((float) (z_max - z_min) / (2.0 * r)));
	unsigned long pointsTotal = 0;
	int imod;
	PaliunisPoint currentPoint;

	currentPoint.y = y_min;
	for (unsigned int j = 1; j <= burbuluKiekisYAsyje; j++) {
		currentPoint.x = x_min - r;
		for (unsigned int i = 1; i <= burbuluKiekisXAsyje; i++) {
			currentPoint.x = (float) (x_min + r + (i - 1) * h);
			currentPoint.y = (float) (y_min - r + r + (j - 1) * h);

			imod = i % 2;
			if (imod == 1) {
				currentPoint.y = (float) (y_min - r - r + r + j * 2 * r);
				//  currentPoint.x =  x_min + ((sqrt(3.0)-1.0) * r + r) * i;
				//  currentPoint.y = (float)(y_min + r + 2* r * j);
			} else {
				currentPoint.y = (float) (y_min + j * 2 * r);
				// currentPoint.x =  x_min + ((sqrt(3.0)-1.0) * r + r) * i;
				// currentPoint.y = (float)(y_min + 2*j*r);
			}
			if (dimensions == 3) {
				currentPoint.z = 0;
				for (unsigned int k = 0; k < burbuluKiekisZAsyje; k++) {
					currentPoint.z = currentPoint.z + 2.0 * r;
					points->InsertPoint(pointsTotal++, currentPoint.x, currentPoint.y, currentPoint.z);
				}
			} else
				points->InsertPoint(pointsTotal++, currentPoint.x, currentPoint.y, currentPoint.z);
		}
	}
	return points->GetNumberOfPoints();
}
;

//##################################################################################################

void GetForceBetwinTwoPoints(vtkDoubleArray * array, vtkPolyData *input, unsigned int id1, unsigned int id2, double *outForce, double *outLength)
// array-jegu masyvas,  id1-pirmo tasko id, antro veliau
		{
	//unsigned int cellNum = input->GetNumberOfCells();
	unsigned int temp_id1 = 0;
	unsigned int temp_id2 = 0;
	*outForce = 0;
	*(outForce + 1) = 0;
	*(outForce + 2) = 0;
	*outLength = 0;
	*(outLength + 1) = 0;
	*(outLength + 2) = 0;
//	double alpha = 0;
	//   input->Squeeze();
	vtkPoints * points = input->GetPoints();
	double point1[3];
	points->GetPoint(id1, point1);
	double point2[3];
	points->GetPoint(id2, point2);
	vtkIdList * list = vtkIdList::New();
	input->GetPointCells(id1, list);
	unsigned int numberOfIDs = list->GetNumberOfIds();
	for (unsigned int j = 0; j < numberOfIDs; j++) {
		temp_id1 = input->GetCell(list->GetId(j))->GetPointId(0);
		temp_id2 = input->GetCell(list->GetId(j))->GetPointId(1);
		if ((id1 == temp_id1 || id1 == temp_id2) && (id2 == temp_id1 || id2 == temp_id2)) {
			*outLength = (point1[0] - point2[0]) / 2;
			*(outLength + 1) = (point1[1] - point2[1]) / 2;
			double l_modul = sqrt(pow(outLength[0], 2.0) + pow(outLength[1], 2.0));
			double forceTemp = array->GetValue(list->GetId(j));
			*outForce = forceTemp * outLength[0] / l_modul;
			*(outForce + 1) = forceTemp * outLength[1] / l_modul;
			points->Squeeze();
			return;
		}
	}
}

//##################################################################################################

int vtkCrackTensor::RequestData(vtkInformation *vtkNotUsed(request), vtkInformationVector **inputVector, vtkInformationVector *outputVector) {
	// get the info objects
	vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
	vtkInformation *outInfo = outputVector->GetInformationObject(0);

	// get the input and ouptut
	vtkPolyData *input = vtkPolyData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
	vtkPolyData *output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

	int numPts = input->GetNumberOfPoints();
	vtkPoints * points = input->GetPoints();

	vtkDataArray*array = input->GetCellData()->GetArray(getForceArray());

	vtkPoints * circlePoints = vtkPoints::New();

	input->ComputeBounds();

	double* bounds = new double[6];
	input->GetBounds(bounds);
	double V = 3.14 * pow(Radius, 2);

	GenerateCircles(bounds[0], bounds[1], bounds[2], bounds[3], bounds[4], bounds[5], 2, Radius, circlePoints);

	unsigned int circleNum = circlePoints->GetNumberOfPoints();
	double currentCircle[3];
	double currentPoint[3];
	CircleData * circleData = new CircleData[circleNum];
	PointWithNeighbours current_point;
	vtkIdList * list = vtkIdList::New();
	//   vtkIdList * listOfNeighbours = vtkIdList::New();
	for (unsigned int i = 0; i < circleNum; i++) {
		circlePoints->GetPoint(i, currentCircle);
		for (int k = 0; k < numPts; k++) {
			points->GetPoint(k, currentPoint);
			if (IsPointInCircle(currentCircle, currentPoint, Radius)) {
				current_point.PointID = k;
				circleData[i].PointsInCircle.push_back(current_point);
			}
		}
		// apeinam viduje esancius taskus zinodami kokie ten yra ir ieskom ju kaimyniu
		for (unsigned int k = 0; k < circleData[i].PointsInCircle.size(); k++) {
			input->GetPointCells(circleData[i].PointsInCircle[k].PointID, list);
			for (unsigned int j = 0; j < list->GetNumberOfIds(); j++) {
				unsigned int id1 = input->GetCell(list->GetId(j))->GetPointId(0);
				unsigned int id2 = input->GetCell(list->GetId(j))->GetPointId(1);
				PaliunisPointID pointID;
				pointID.IsInCircle = false;
				if (id1 != circleData[i].PointsInCircle[k].PointID)
					pointID.id = id1;
				else
					pointID.id = id2;

				double PointTemp[3];
				points->GetPoint(pointID.id, PointTemp);
				if (IsPointInCircle(currentCircle, PointTemp, Radius))
					pointID.IsInCircle = true;
				else {
					// ar nebuvo irasytas i ta sfera dalelyte id
					vector<unsigned int>::iterator it = find(circleData[i].PointsAbroad.begin(), circleData[i].PointsAbroad.end(), pointID.id);
					if (it == circleData[i].PointsAbroad.end())
						circleData[i].PointsAbroad.push_back(pointID.id);
				}

				circleData[i].PointsInCircle[k].Neighbours.push_back(pointID);
			}
		}

	}
	// suminiu jegu skaiciavimas nuo daleles vienos
	double PointForces[3];
	double PointDistance[3];
	unsigned int temp_id1;
	unsigned int temp_id2;
	unsigned int another_id;
	unsigned int bendras_dalelyciu_skaicius = 0;
	unsigned int numberOfIDs;
	for (unsigned int i = 0; i < circleNum; i++) {
		//   // LOG_INFO("Sfera Nr " << i);
		unsigned int pointsInSphereNum = circleData[i].PointsInCircle.size();
		//  // LOG_INFO("Daleliu skaicius = " << pointsInSphereNum);
		bendras_dalelyciu_skaicius += pointsInSphereNum;
		// // LOG_INFO("Bendras skaicius = " << bendras_dalelyciu_skaicius);
		for (unsigned int k = 0; k < pointsInSphereNum; k++) {
			// vienos celes l ir f skaiciuojam tenzorius
			input->GetPointCells(circleData[i].PointsInCircle[k].PointID, list);
			numberOfIDs = list->GetNumberOfIds();
			memset(circleData[i].PointsInCircle[k].TensorElements, 0, sizeof(double) * 4);
			for (unsigned int j = 0; j < numberOfIDs; j++) {
				temp_id1 = input->GetCell(list->GetId(j))->GetPointId(0);
				temp_id2 = input->GetCell(list->GetId(j))->GetPointId(1);
				another_id = circleData[i].PointsInCircle[k].PointID == temp_id1 ? temp_id2 : temp_id1;

				GetForceBetwinTwoPoints((vtkDoubleArray *) array, input, circleData[i].PointsInCircle[k].PointID, another_id, PointForces, PointDistance);

				vector<unsigned int>::iterator it = find(circleData[i].PointsAbroad.begin(), circleData[i].PointsAbroad.end(), another_id);
				if (it == circleData[i].PointsAbroad.end()) {
					circleData[i].PointsInCircle[k].TensorElements[0] += PointForces[0] * PointDistance[0];
					circleData[i].PointsInCircle[k].TensorElements[1] += PointForces[1] * PointDistance[0];
					circleData[i].PointsInCircle[k].TensorElements[2] += PointForces[0] * PointDistance[1];
					circleData[i].PointsInCircle[k].TensorElements[3] += PointForces[1] * PointDistance[1];
				} else {
					circleData[i].PointsInCircle[k].TensorElements[0] += PointForces[0] * PointDistance[0];
					circleData[i].PointsInCircle[k].TensorElements[1] += PointForces[1] * PointDistance[0];
					circleData[i].PointsInCircle[k].TensorElements[2] += PointForces[0] * PointDistance[1];
					circleData[i].PointsInCircle[k].TensorElements[3] += PointForces[1] * PointDistance[1];
				}
			}
		}
	}

	for (unsigned int i = 0; i < circleNum; i++) {
		memset(circleData[i].TotalTensorElements, 0, sizeof(double) * 4);
		for (unsigned int k = 0; k < circleData[i].PointsInCircle.size(); k++) {
			circleData[i].TotalTensorElements[0] += circleData[i].PointsInCircle[k].TensorElements[0] / V;
			circleData[i].TotalTensorElements[1] += circleData[i].PointsInCircle[k].TensorElements[1] / V;
			circleData[i].TotalTensorElements[2] += circleData[i].PointsInCircle[k].TensorElements[2] / V;
			circleData[i].TotalTensorElements[3] += circleData[i].PointsInCircle[k].TensorElements[3] / V;
		}

	}

	vtkCellArray * pointscells = vtkCellArray::New();
	vtkDoubleArray * radiuss = vtkDoubleArray::New();
	radiuss->SetName("Radius");
	radiuss->SetNumberOfComponents(1);
	radiuss->SetNumberOfTuples(circleNum);

	for (int i = 0; i < circlePoints->GetNumberOfPoints(); i++) {
		pointscells->InsertNextCell(1);
		pointscells->InsertCellPoint(i);

	}

	vtkDoubleArray * resultArray = vtkDoubleArray::New();
	resultArray->SetName("Tensors");
	resultArray->SetNumberOfComponents(9);
	resultArray->SetNumberOfTuples(circleNum);

	vtkDoubleArray * Sigma11 = vtkDoubleArray::New();
	Sigma11->SetName("Sigma11");
	Sigma11->SetNumberOfComponents(1);
	Sigma11->SetNumberOfTuples(circleNum);

	vtkDoubleArray * Sigma12 = vtkDoubleArray::New();
	Sigma12->SetName("Sigma12");
	Sigma12->SetNumberOfComponents(1);
	Sigma12->SetNumberOfTuples(circleNum);

	vtkDoubleArray * Sigma21 = vtkDoubleArray::New();
	Sigma21->SetName("Sigma21");
	Sigma21->SetNumberOfComponents(1);
	Sigma21->SetNumberOfTuples(circleNum);

	vtkDoubleArray * Sigma22 = vtkDoubleArray::New();
	Sigma22->SetName("Sigma22");
	Sigma22->SetNumberOfComponents(1);
	Sigma22->SetNumberOfTuples(circleNum);

	for (unsigned int i = 0; i < circleNum; i++) {
		radiuss->SetValue(i, Radius);
		Sigma11->SetTuple1(i, circleData[i].TotalTensorElements[0]);
		Sigma12->SetTuple1(i, circleData[i].TotalTensorElements[1]);
		Sigma21->SetTuple1(i, circleData[i].TotalTensorElements[2]);
		Sigma22->SetTuple1(i, circleData[i].TotalTensorElements[3]);
		resultArray->SetTuple9(i, circleData[i].TotalTensorElements[0], circleData[i].TotalTensorElements[1], 0, circleData[i].TotalTensorElements[2], circleData[i].TotalTensorElements[3], 0, 0, 0, 0

		);

	}
	int kieka = 0;
	for (unsigned int i = 0; i < circleNum; i++) {
		kieka = kieka + circleData[i].PointsInCircle.size();
	}
	output->SetVerts(pointscells);
	output->SetPoints(circlePoints);
	output->GetPointData()->AddArray(radiuss);
	output->GetPointData()->SetTensors(resultArray);
	output->GetPointData()->AddArray(resultArray);
	output->GetPointData()->AddArray(Sigma11);
	output->GetPointData()->AddArray(Sigma12);
	output->GetPointData()->AddArray(Sigma21);
	output->GetPointData()->AddArray(Sigma22);
	output->GetPointData()->SetActiveTensors("Tensors");

	//------------------------------------------cia idomi vieta ar deapcopy ar shallowcopy
	//output->DeepCopy(input);
	output->Squeeze();
	return 1;
}

int vtkCrackTensor::FillInputPortInformation(int, vtkInformation *info) {
	info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
	return 1;
}

void vtkCrackTensor::PrintSelf(ostream& os, vtkIndent indent) {

}

}
}
