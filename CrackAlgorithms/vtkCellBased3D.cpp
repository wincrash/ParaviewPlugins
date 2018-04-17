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

#include "vtkCellBased3D.h"
namespace vispartdem {
namespace cellbased3D {
vtkCxxRevisionMacro(vtkCellBased3D, "$Revision: 1.42 $")
vtkStandardNewMacro(vtkCellBased3D)

vtkCellBased3D::vtkCellBased3D() {
	this->StateArray = "";
	onetime = true;
}
vtkCellBased3D::~vtkCellBased3D() {
}

int vtkCellBased3D::FillInputPortInformation(int, vtkInformation *info) {
	info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
	return 1;
}

void vtkCellBased3D::PrintSelf(ostream& os, vtkIndent indent) {
	this->Superclass::PrintSelf(os, indent);
}

void vtkCellBased3D::SetStateArray(string ResultArrayName) {
	this->StateArray = ResultArrayName;
}

const char *vtkCellBased3D::getStateArray() {
	return this->StateArray.c_str();
}

void vtkCellBased3D::CreateTetraUnstructuredGrid(vtkUnstructuredGrid *output, vtkPolyData* input, vtkDataArray*state) {
	double point1[3];
	double point2[3];
	double point3[3];

	double virsune1[3];
	double normal[3];
	double DD, d;
	vector<vispartdem::INT> tempas;
	tempas.resize(8, -1);
	vector<vispartdem::INT> tempas1;
	tempas1.resize(4, -1);
	vector<vispartdem::MyStruct> *strukturos = structura3D.getStruktura();
	for (vector<vispartdem::MyStruct>::iterator it = strukturos->begin(); it != strukturos->end(); it++) {
		if (it->ids.size() == 4) {
			input->GetPoint(it->ids[0], point1);
			input->GetPoint(it->ids[1], point2);
			input->GetPoint(it->ids[2], point3);
			input->GetPoint(it->ids[3], virsune1);
			vtkTriangle::ComputeNormal(point1, point2, point3, normal);

			d = point1[0] * normal[0] + point1[1] * normal[1] + point1[2] * normal[2];
			DD = (normal[0] * virsune1[0] + normal[1] * virsune1[1] + normal[2] * virsune1[2] + d) / vtkMath::Norm(normal);
			if (DD < 0) {
				swap(it->ids[0], it->ids[1]);
			}

			tempas1[0] = it->ids[0];
			tempas1[1] = it->ids[1];
			tempas1[2] = it->ids[2];
			tempas1[3] = it->ids[3];
			tempas[0] = it->cellIDs[0];
			tempas[1] = it->cellIDs[1];
			tempas[2] = it->cellIDs[2];
			tempas[3] = it->cellIDs[3];
			tempas[4] = it->cellIDs[4];
			tempas[5] = it->cellIDs[5];
			cellsIDTetra.push_back(tempas);
			IDTetra.push_back(tempas1);
		}
	}
}

void vtkCellBased3D::CreatePyramidUnstructuredGrid(vtkUnstructuredGrid *output, vtkPolyData* input, vtkDataArray*state) {
	double point1[3];
	double point2[3];
	double point3[3];
	double point4[3];
	double virsune1[3];
	double normal[3];
	double DD, d;
	vector<vispartdem::INT> tempas;
	tempas.resize(8, -1);
	vector<vispartdem::INT> tempas1;
	tempas1.resize(5, -1);
	vector<vispartdem::MyStruct> *strukturos = structura3D.getStruktura();
	for (vector<vispartdem::MyStruct>::iterator it = strukturos->begin(); it != strukturos->end(); it++) {

		if (it->ids.size() != 4) {
			input->GetPoint(it->ids[0], point1);
			input->GetPoint(it->ids[1], point2);
			input->GetPoint(it->ids[2], point3);
			input->GetPoint(it->ids[3], point4);
			input->GetPoint(it->ids[4], virsune1);
			vtkTriangle::ComputeNormal(point1, point2, point3, normal);
			d = point1[0] * normal[0] + point1[1] * normal[1] + point1[2] * normal[2];
			DD = (normal[0] * virsune1[0] + normal[1] * virsune1[1] + normal[2] * virsune1[2] + d) / vtkMath::Norm(normal);

			if (DD < 0) {
				tempas1[0] = it->ids[0];
				tempas1[1] = it->ids[3];
				tempas1[2] = it->ids[2];
				tempas1[3] = it->ids[1];
				tempas1[4] = it->ids[4];
			} else {
				tempas1[0] = it->ids[0];
				tempas1[1] = it->ids[1];
				tempas1[2] = it->ids[2];
				tempas1[3] = it->ids[3];
				tempas1[4] = it->ids[4];
			}

			tempas[0] = it->cellIDs[0];
			tempas[1] = it->cellIDs[1];
			tempas[2] = it->cellIDs[2];
			tempas[3] = it->cellIDs[3];
			tempas[4] = it->cellIDs[4];
			tempas[5] = it->cellIDs[5];
			tempas[6] = it->cellIDs[6];
			tempas[7] = it->cellIDs[7];
			cellsIDPyra.push_back(tempas);
			IDPyra.push_back(tempas1);

			//=========================================================
			if (it->ids[5] != -1) {
				input->GetPoint(it->ids[5], virsune1);
				DD = (normal[0] * virsune1[0] + normal[1] * virsune1[1] + normal[2] * virsune1[2] + d) / vtkMath::Norm(normal);

				if (DD < 0) {
					tempas1[0] = it->ids[0];
					tempas1[1] = it->ids[3];
					tempas1[2] = it->ids[2];
					tempas1[3] = it->ids[1];
					tempas1[4] = it->ids[5];
				} else {
					tempas1[0] = it->ids[0];
					tempas1[1] = it->ids[1];
					tempas1[2] = it->ids[2];
					tempas1[3] = it->ids[3];
					tempas1[4] = it->ids[5];
				}
				tempas[0] = it->cellIDs[0];
				tempas[1] = it->cellIDs[1];
				tempas[2] = it->cellIDs[2];
				tempas[3] = it->cellIDs[3];
				tempas[4] = it->cellIDs[8];
				tempas[5] = it->cellIDs[9];
				tempas[6] = it->cellIDs[10];
				tempas[7] = it->cellIDs[11];
				cellsIDPyra.push_back(tempas);
				IDPyra.push_back(tempas1);
			}
		}
	}

}

int vtkCellBased3D::RequestData(vtkInformation *vtkNotUsed(request), vtkInformationVector **inputVector, vtkInformationVector *outputVector) {
	// get the info objects
	vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
	vtkInformation *outInfo = outputVector->GetInformationObject(0);

	// get the input and ouptut
	vtkPolyData *input = vtkPolyData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
    vtkUnstructuredGrid *output1 = vtkUnstructuredGrid::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

    vtkUnstructuredGrid*output = vtkUnstructuredGrid::New();
	vtkDataArray*state = input->GetCellData()->GetArray(this->StateArray.c_str()); //gauna nurodyta state masyva
	if (state == NULL) {
		vtkErrorMacro( << "Nenurodytas state masyvas (turi buti idetas i celldata)!\n");
		return 0;
	}
	output->Allocate(1, 1);
	Timeris t_b;
	Timeris t;
	Timeris t1;
	t_b.start();
	if (onetime) {
		t.start();
        structura3D.Create3DStructure(input,StateArray,true);
		t.stop();
		onetime = false;
		t1.start();
		CreateTetraUnstructuredGrid(output, input, state);
		CreatePyramidUnstructuredGrid(output, input, state);
        cells = vtkUnstructuredGrid::New();
		cells->Allocate(1, 1);
		vtkIdList*list = vtkIdList::New();
		list->SetNumberOfIds(4);
        for (unsigned int k = 0; k < cellsIDTetra.size(); k++) {
			for (int i = 0; i < 4; i++) {
                list->SetId(i, IDTetra[k][i]);
			}            
            cells->InsertNextCell(VTK_TETRA, list);
		}
		list->SetNumberOfIds(5);
        for (unsigned int k = 0; k < cellsIDPyra.size(); k++) {
            for (int i = 0; i < 5; i++) {
				list->SetId(i, IDPyra[k][i]);
			}
            cells->InsertNextCell(VTK_PYRAMID, list);
		}

		t1.stop();
		list->Delete();        
	}

    t.print("vtkCellBased3D_Triangulation");
    t1.print("vtkCellBased3D_Structure_To_VTK");
	Timeris t2;
	t2.start();

	vtkPoints*points = vtkPoints::New(); //tie patys taskai
	vtkIntArray*newstate = vtkIntArray::New(); //bus laikoma naujas state masyvas
	newstate->SetName(state->GetName());
	newstate->SetNumberOfComponents(1);

	int suma;
	for (unsigned int k = 0; k < cellsIDTetra.size(); k++) {
		suma = 0;
		for (int i = 0; i < 6; i++) {
			suma = suma + state->GetTuple1(cellsIDTetra[k][i]);
		}
		suma = 6 - suma;
		newstate->InsertNextTuple1(suma);
	}

	for (unsigned int k = 0; k < cellsIDPyra.size(); k++) {
		suma = 0;
		for (int i = 0; i < 8; i++) {
			suma = suma + state->GetTuple1(cellsIDPyra[k][i]);
		}
		suma = 8 - suma;
		newstate->InsertNextTuple1(suma);
	}
    output->DeepCopy(cells);
	points->DeepCopy(input->GetPoints());
	output->SetPoints(points); //.sudedam mazgus
	output->GetPointData()->DeepCopy(input->GetPointData());
	output->GetCellData()->AddArray(newstate); //idedama nauja sate masyva
	output->Squeeze();    
	points->Delete();
	newstate->Delete();
	t2.stop();
    t2.print("vtkCellBased3D_Attribute_calculation");
    Timeris t3;
	t3.start();
    vtkExtractCells*extract = vtkExtractCells::New();
	extract->SetInput(output);
	vtkIdList *ptrIds = vtkIdList::New();

	for (int i = 0; i < output->GetNumberOfCells(); i++) {
		if (output->GetCellData()->GetArray("group/Scalar3")->GetTuple1(i) != 0) {
			ptrIds->InsertNextId(i);
		}
	}
	extract->SetCellList(ptrIds);
    extract->Update();
	t3.stop();
    t3.print("vtkCellBased3D_vtkExtractCells_veikimo_laikas");
    output1->DeepCopy(extract->GetOutput());
   //  output1->DeepCopy(output);
     output1->Squeeze();
    //extract->Delete();

	t_b.stop();
    t_b.print("vtkCellBased3D_bendras");
	return 1;
}
}
}
