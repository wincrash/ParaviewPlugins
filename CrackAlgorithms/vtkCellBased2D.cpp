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

#include "vtkCellBased2D.h"
namespace vispartdem {
namespace cellbased2D {

vtkCxxRevisionMacro(vtkCellBased2D, "$Revision: 1.42 $")
vtkStandardNewMacro(vtkCellBased2D)

vtkCellBased2D::vtkCellBased2D() {
	this->StateArray = "";
	onetime = true;
}
vtkCellBased2D::~vtkCellBased2D() {

}

int vtkCellBased2D::RequestData(vtkInformation *vtkNotUsed(request), vtkInformationVector **inputVector, vtkInformationVector *outputVector) {
	// get the info objects
	vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
	vtkInformation *outInfo = outputVector->GetInformationObject(0);

	// get the input and ouptut
	vtkPolyData *input = vtkPolyData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
	vtkPolyData *output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
	Timeris t_b;
	Timeris t;
	Timeris t1;
	t_b.start();
	if (onetime) {
		t.start();
        structura.Create2DStructure(input,StateArray,true);
		t.stop();

		t1.start();
		triang = vtkCellArray::New();
		vector<vispartdem::MyStruct>* strukturos = structura.getStruktura();
		for (vector<vispartdem::MyStruct>::iterator iterator = strukturos->begin(); iterator != strukturos->end(); iterator++) {
			triang->InsertNextCell(3);
			triang->InsertCellPoint((*iterator).ids[0]);
			triang->InsertCellPoint((*iterator).ids[1]);
			triang->InsertCellPoint((*iterator).ids[2]);
		}
		t1.stop();
		onetime = false;
	}
    t.print("vtkCellBased2D_Triangulation");
    t1.print("vtkCellBased2D_Structure_To_VTK");
	Timeris t2;
	t2.start();
	vtkPoints*points = vtkPoints::New(); //tie patys taskai
	vtkCellArray*faces = vtkCellArray::New(); /// trikiampiai bus dedami cia
	vtkIntArray*newstate = vtkIntArray::New(); //bus laikoma naujas state masyvas
	vtkDataArray*state = input->GetCellData()->GetArray(this->StateArray.c_str()); //gauna nurodyta state masyva
	if (state == NULL) {
		vtkErrorMacro( << "Nenurodytas state masyvas (turi buti idetas i celldata)!\n");
		return 0;
	}

	vector<vispartdem::MyStruct>* strukturos = structura.getStruktura();
	newstate->SetName(state->GetName());
	newstate->SetNumberOfComponents(1);
	newstate->SetNumberOfTuples(strukturos->size());

	int i = 0;
	for (vector<vispartdem::MyStruct>::iterator iterator = strukturos->begin(); iterator != strukturos->end(); iterator++) {
		newstate->SetTuple1(i, state->GetTuple1((*iterator).cellIDs[0]) + state->GetTuple1((*iterator).cellIDs[1]) + state->GetTuple1((*iterator).cellIDs[2]));
		i++;
	}

	faces->DeepCopy(triang);
	points->DeepCopy(input->GetPoints());
	output->SetPoints(points); //.sudedam mazgus
	points->Delete();
	output->GetCellData()->AddArray(newstate); //idedama nauja ste masyva
	newstate->Delete();
	output->SetPolys(faces); //sudedama trikampus
	faces->Delete();
	output->Squeeze();
	t2.stop();
    t2.print("vtkCellBased2D_Attribute_calculation");
	t_b.stop();
    t_b.print("vtkCellBased2D_bendras_laikas");
	return 1;
}

int vtkCellBased2D::FillInputPortInformation(int, vtkInformation *info) {
	info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
	return 1;
}

void vtkCellBased2D::PrintSelf(ostream& os, vtkIndent indent) {
	this->Superclass::PrintSelf(os, indent);
	os << "State array name " << StateArray << endl;
}

//-------------------------------------------

void vtkCellBased2D::SetStateArray(string ResultArrayName) {
	this->StateArray = ResultArrayName;
}

const char *vtkCellBased2D::getStateArray() {
	return this->StateArray.c_str();
}
}
}
