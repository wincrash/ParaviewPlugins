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

#include "vtkCellCut2D.h"
namespace vispartdem {
namespace CellCut3D {
vtkCxxRevisionMacro(vtkCellCut2D, "$Revision: 1.42 $")
vtkStandardNewMacro(vtkCellCut2D)

vtkCellCut2D::vtkCellCut2D() {
	this->StateArray = "";
	onetime = true;
}
vtkCellCut2D::~vtkCellCut2D() {

}

int vtkCellCut2D::RequestData(vtkInformation *vtkNotUsed(request), vtkInformationVector **inputVector, vtkInformationVector *outputVector) {
	// get the info objects
	vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
	vtkInformation *outInfo = outputVector->GetInformationObject(0);

	// get the input and ouptut
	vtkPolyData *input = vtkPolyData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
	vtkPolyData *output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
	Timeris t;
	Timeris t1;
	Timeris t_b;
	t_b.start();
    //if (onetime) {
		t.start();
        structura.Create2DStructure(input,StateArray,false);
		onetime = false;
		t.stop();
    //}
	t1.start();
	vtkPoints*points = vtkPoints::New(); //tie patys taskai
	vtkCellArray*lines = vtkCellArray::New(); /// linijos bus dedami cia
	vtkCellArray*vertex = vtkCellArray::New(); /// trikiampiai bus dedami cia
	vtkDataArray*state = input->GetCellData()->GetArray(this->StateArray.c_str()); //gauna nurodyta state masyva
	if (state == NULL) {
		vtkErrorMacro( << "Nenurodytas state masyvas (turi buti idetas i celldata)!\n");
		return 0;
	}
	vector<vispartdem::MyStruct>* strukturos = structura.getStruktura();
	int suma = 0;
	double p1[3];
	double p2[3];
	double p3[3];

	double cp1[3];
	double cp2[3];
	double cp3[3];
	double ccenter[3];
	vispartdem::INT centerID;
	vector<vispartdem::INT> pointsID;
	pointsID.resize(input->GetNumberOfCells(), -1);
    map<pair<int,int>, int> buvooo;
    pair<int,int> aaa;
    map<int,int> buvo1;
    map<pair<pair<int,int>,int>, int> buvo3;
    pair<pair<int,int>,int> trys;
	for (vector<vispartdem::MyStruct>::iterator iterator = strukturos->begin(); iterator != strukturos->end(); iterator++) {
		suma = state->GetTuple1((*iterator).cellIDs[0]) + state->GetTuple1((*iterator).cellIDs[1]) + state->GetTuple1((*iterator).cellIDs[2]);

		if (suma == 3) {
			continue;
		}
		input->GetPoint((*iterator).ids[0], p1);
		input->GetPoint((*iterator).ids[1], p2);
		input->GetPoint((*iterator).ids[2], p3);


        if (pointsID[(*iterator).cellIDs[0]] == -1 && state->GetTuple1((*iterator).cellIDs[0])==0) {
			vtkMath::Add(p1, p2, cp1); //1
			vtkMath::MultiplyScalar(cp1, 0.5);
			pointsID[(*iterator).cellIDs[0]] = points->InsertNextPoint(cp1);
		}
        if (pointsID[(*iterator).cellIDs[1]] == -1 && state->GetTuple1((*iterator).cellIDs[1])==0) {
			vtkMath::Add(p2, p3, cp2); //2
			vtkMath::MultiplyScalar(cp2, 0.5);
			pointsID[(*iterator).cellIDs[1]] = points->InsertNextPoint(cp2);
		}
        if (pointsID[(*iterator).cellIDs[2]] == -1 && state->GetTuple1((*iterator).cellIDs[2])==0) {
			vtkMath::Add(p1, p3, cp3); //3
			vtkMath::MultiplyScalar(cp3, 0.5);
			pointsID[(*iterator).cellIDs[2]] = points->InsertNextPoint(cp3);
		}

		if (suma == 0) //visi nutruke
				{


            trys.first.first=(*iterator).cellIDs[0];
            trys.first.second=(*iterator).cellIDs[1];
            trys.second=(*iterator).cellIDs[2];
            if(trys.first.first>trys.first.second)
            {
                swap(trys.first.first,trys.first.second);
            }
            if(trys.first.first>trys.second) {
                 swap(trys.first.first,trys.second);
            }
            if(trys.first.second>trys.second) {
                 swap(trys.first.second,trys.second);
            }
            if(buvo3[trys]==0)
            {
                buvo3[trys]=5;

            double a, b, c, abc;
            a = vtkMath::Distance2BetweenPoints(p1, p2);
            b = vtkMath::Distance2BetweenPoints(p2, p3);
            c = vtkMath::Distance2BetweenPoints(p1, p3);
			abc = a + b + c;
			ccenter[0] = (a * p3[0] + b * p1[0] + c * p2[0]) / abc;
			ccenter[1] = (a * p3[1] + b * p1[1] + c * p2[1]) / abc;
			ccenter[2] = (p1[2] + p2[2] + p3[2]) / 3;
			centerID = points->InsertNextPoint(ccenter);
			lines->InsertNextCell(2);
			lines->InsertCellPoint(centerID);
			lines->InsertCellPoint(pointsID[(*iterator).cellIDs[0]]);
			lines->InsertNextCell(2);
			lines->InsertCellPoint(centerID);
			lines->InsertCellPoint(pointsID[(*iterator).cellIDs[1]]);
			lines->InsertNextCell(2);
			lines->InsertCellPoint(centerID);
			lines->InsertCellPoint(pointsID[(*iterator).cellIDs[2]]);
            }

		} else if (state->GetTuple1((*iterator).cellIDs[0]) == 0 && state->GetTuple1((*iterator).cellIDs[1]) == 0) { //1,2
            aaa.first=(*iterator).cellIDs[0];
            aaa.second=(*iterator).cellIDs[1];
            if(aaa.first>aaa.second)
            {
                swap(aaa.first,aaa.second);
            }
            if(buvooo[aaa]==0)
            {
			lines->InsertNextCell(2);
			lines->InsertCellPoint(pointsID[(*iterator).cellIDs[0]]);
			lines->InsertCellPoint(pointsID[(*iterator).cellIDs[1]]);
            buvooo[aaa]=5;
            }

		} else if (state->GetTuple1((*iterator).cellIDs[1]) == 0 && state->GetTuple1((*iterator).cellIDs[2]) == 0) { //2,3
            aaa.first=(*iterator).cellIDs[1];
            aaa.second=(*iterator).cellIDs[2];
            if(aaa.first>aaa.second)
            {
                swap(aaa.first,aaa.second);
            }
            if(buvooo[aaa]==0)
            {
			lines->InsertNextCell(2);
			lines->InsertCellPoint(pointsID[(*iterator).cellIDs[1]]);
			lines->InsertCellPoint(pointsID[(*iterator).cellIDs[2]]);
            buvooo[aaa]=5;
            }

		} else if (state->GetTuple1((*iterator).cellIDs[2]) == 0 && state->GetTuple1((*iterator).cellIDs[0]) == 0) { //1,3
                aaa.first=(*iterator).cellIDs[0];
                aaa.second=(*iterator).cellIDs[2];
                if(aaa.first>aaa.second)
                {
                    swap(aaa.first,aaa.second);
                }
                if(buvooo[aaa]==0)
                {
			lines->InsertNextCell(2);
			lines->InsertCellPoint(pointsID[(*iterator).cellIDs[0]]);
			lines->InsertCellPoint(pointsID[(*iterator).cellIDs[2]]);
            buvooo[aaa]=5;
            }
		} else if (state->GetTuple1((*iterator).cellIDs[0]) == 0) {
            if(buvo1[(*iterator).cellIDs[0]]==0)
            {
                buvo1[(*iterator).cellIDs[0]]=5;
            vertex->InsertNextCell(1);
            vertex->InsertCellPoint(pointsID[(*iterator).cellIDs[0]]);

            }

		} else if (state->GetTuple1((*iterator).cellIDs[1]) == 0) {
            if(buvo1[(*iterator).cellIDs[1]]==0)
            {
                buvo1[(*iterator).cellIDs[1]]=5;
            vertex->InsertNextCell(1);
            vertex->InsertCellPoint(pointsID[(*iterator).cellIDs[1]]);
            }

		} else if (state->GetTuple1((*iterator).cellIDs[2]) == 0) {
                if(buvo1[(*iterator).cellIDs[2]]==0)
                {
                    buvo1[(*iterator).cellIDs[2]]=5;
            vertex->InsertNextCell(1);
            vertex->InsertCellPoint(pointsID[(*iterator).cellIDs[2]]);
                }
		}
	}
    t1.stop();
	output->SetPoints(points); //.sudedam mazgus
	points->Delete();
	output->SetLines(lines); //sudedama trikampus
	output->SetVerts(vertex);
	vertex->Delete();
	lines->Delete();
	output->Squeeze();

	t_b.stop();
    t.print("vtkCellCut2D_Triangulation_su_strukturomis");
    t1.print("vtkCellCut2D_Marching");
    t_b.print("vtkCellCut2D_bendras_laikas");


	return 1;
}

int vtkCellCut2D::FillInputPortInformation(int, vtkInformation *info) {
	info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
	return 1;
}

void vtkCellCut2D::PrintSelf(ostream& os, vtkIndent indent) {
	this->Superclass::PrintSelf(os, indent);
	os << "State array name " << StateArray << endl;
}
void vtkCellCut2D::SetStateArray(string ResultArrayName) {
	this->StateArray = ResultArrayName;
}
const char *vtkCellCut2D::getStateArray() {
	return this->StateArray.c_str();
}
}
}
