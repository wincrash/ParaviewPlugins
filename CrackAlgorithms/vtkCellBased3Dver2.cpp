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

#include "vtkCellBased3Dver2.h"
namespace vispartdem {
namespace cellbased3DRev2 {
vtkCxxRevisionMacro(vtkCellBased3Dver2, "$Revision: 1.42 $")
vtkStandardNewMacro(vtkCellBased3Dver2)

vtkCellBased3Dver2::vtkCellBased3Dver2() {
	this->StateArray = "";
	onetime = true;
}

int vtkCellBased3Dver2::RequestData(vtkInformation *vtkNotUsed(request), vtkInformationVector **inputVector, vtkInformationVector *outputVector) {
	// get the info objects
	vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
	vtkInformation *outInfo = outputVector->GetInformationObject(0);

	// get the input and ouptut
	vtkPolyData *input = vtkPolyData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
	vtkUnstructuredGrid *output1 = vtkUnstructuredGrid::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
	vtkUnstructuredGrid*output = vtkUnstructuredGrid::New();
	output->Allocate(1, 1);
	Timeris t_b;

	Timeris t;
	Timeris t1;
	t_b.start();
	if (onetime == true) {

		t.start();
		onetime = false;
		//-----get numbers of cells and points
		int numberofcells = input->GetNumberOfCells();
		vtkIdType* connections = (vtkIdType*) input->GetLines()->GetPointer();
		map<_line, int> linijos;
		map<_line, int>::iterator linijosIT;

		map<_trik, INT> trikampiai;
		map<_pyramid, INT> pyramid;
		map<_pyramid, INT>::iterator pyramidIT;
		map<_keturkampis, INT> keturkampiai;
		map<_keturkampis, INT>::iterator keturkampiaiIT;
		_keturkampis ket;

		_trik tria;
		_line l;
		_pyramid pyra;
		int id1, id2, id3, id4, id;
		vector<vector<int> > pointCells;
		vector<int> a;
		a.reserve(20);
		pointCells.resize(input->GetNumberOfPoints(), a);

		for (int i = 0; i < numberofcells; ++i) {
			l.first = connections[i * 3 + 1];
			l.second = connections[i * 3 + 2];
			if (l.first > l.second) {
				swap(l.first, l.second);
			}
			linijos[l] = i;
			pointCells[l.first].push_back(l.second);
			pointCells[l.second].push_back(l.first);
		}
		for (unsigned int i = 0; i < pointCells.size(); ++i) {
			sort(pointCells[i].begin(), pointCells[i].end());
		}

		map<int, int>::iterator tempIterator;
		for (id = 0; id < input->GetNumberOfPoints(); ++id) { //sakninis pointas

			for (unsigned int k = 0; k < pointCells[id].size(); ++k) {
				id1 = pointCells[id][k]; //pirmas kaimynas
				if (id == id1) {
					continue;
				}
				for (unsigned int v = k; v < pointCells[id].size(); ++v) {
					id2 = pointCells[id][v]; //antras kaimynas
					if (id2 == id || id2 == id1) {
						continue;
					}
					for (unsigned int z = v; z < pointCells[id].size(); ++z) {
						id3 = pointCells[id][z]; //trecias kaimynas
						if (id3 == id || id3 == id1 || id3 == id2) {
							continue;
						}
						//darom tetraedra
						id4 = -1;
						GetPyramidID(pyra, id, id1, id2, id3);
						if (pyramid.find(pyra) == pyramid.end()) {
							GetLineID(l, id1, id2);
							if (linijos.find(l) != linijos.end()) {
								MyStruct myStruct;
								myStruct.cellIDs.push_back(linijos[l]);
								GetLineID(l, id1, id3);
								if (linijos.find(l) != linijos.end()) {
									myStruct.cellIDs.push_back(linijos[l]);
									GetLineID(l, id3, id2);
									if (linijos.find(l) != linijos.end()) {
										myStruct.cellIDs.push_back(linijos[l]);
										pyramid[pyra] = 1;
										id4 = -2; // pazymime kad buvo sudarytas tetraedras ir neieskosime piramides
										myStruct.ids.push_back(id1); //pagrindas
										myStruct.ids.push_back(id2); //pagrindas
										myStruct.ids.push_back(id3); //pagrindas
										myStruct.ids.push_back(id); //virsune
										GetLineID(l, id, id1);
										myStruct.cellIDs.push_back(linijos[l]);
										GetLineID(l, id, id2);
										myStruct.cellIDs.push_back(linijos[l]);
										GetLineID(l, id, id3);
										myStruct.cellIDs.push_back(linijos[l]);
										strukturos.push_back(myStruct);
									}
								}
							}
						} else {
							id4 = -2;
						}
						int d1, d2, d3, d4;
						if (id4 == -1) {
							for (unsigned int x = z; x < pointCells[id].size(); ++x) {
								id4 = pointCells[id][x]; //ketvirtas kaimynas
								if (id4 == id3 || id4 == id2 || id4 == id1 || id4 == id) {
									continue;
								}
								MyStruct myStruct;
								vector<int> keturkapioIstrizaines; //ieskosim istrizainiu kuriu neturi buti
								keturkapioIstrizaines.reserve(4);
								GetLineID(l, id1, id2);
								if (linijos.find(l) != linijos.end()) {
									myStruct.cellIDs.push_back(linijos[l]);
								} else {
									keturkapioIstrizaines.push_back(id1);
									keturkapioIstrizaines.push_back(id2);
								}
								GetLineID(l, id2, id3);
								if (linijos.find(l) != linijos.end()) {
									myStruct.cellIDs.push_back(linijos[l]);
								} else {
									keturkapioIstrizaines.push_back(id2);
									keturkapioIstrizaines.push_back(id3);
								}
								GetLineID(l, id3, id4);
								if (linijos.find(l) != linijos.end()) {
									myStruct.cellIDs.push_back(linijos[l]);
								} else {
									keturkapioIstrizaines.push_back(id3);
									keturkapioIstrizaines.push_back(id4);
								}
								GetLineID(l, id4, id1);
								if (linijos.find(l) != linijos.end()) {
									myStruct.cellIDs.push_back(linijos[l]);
								} else {
									keturkapioIstrizaines.push_back(id4);
									keturkapioIstrizaines.push_back(id1);
								}
								GetLineID(l, id1, id3);
								if (linijos.find(l) != linijos.end()) {
									myStruct.cellIDs.push_back(linijos[l]);
								} else {
									keturkapioIstrizaines.push_back(id1);
									keturkapioIstrizaines.push_back(id3);
								}
								GetLineID(l, id2, id4);
								if (linijos.find(l) != linijos.end()) {
									myStruct.cellIDs.push_back(linijos[l]);
								} else {
									keturkapioIstrizaines.push_back(id2);
									keturkapioIstrizaines.push_back(id4);
								}
								if (keturkapioIstrizaines.size() != 4) { //tikrinsim     kad butu keturkampis ir neturetu skersai jungciu
									continue;
								}
								for (unsigned int a1 = 0; a1 < pointCells[id1].size(); ++a1) { //ieskom bendru kaimunu
									d1 = pointCells[id1][a1];
									if (d1 == id) { //svarbu kad nebutu jau pradinis mazgo id
										continue;
									}
									for (unsigned int a2 = 0; a2 < pointCells[id2].size(); ++a2) {
										d2 = pointCells[id2][a2];

										if (d1 != d2) {
											continue;
										}
										for (unsigned int a3 = 0; a3 < pointCells[id3].size(); ++a3) {
											d3 = pointCells[id3][a3];
											if (d1 != d3) {
												continue;
											}
											for (unsigned int a4 = 0; a4 < pointCells[id4].size(); ++a4) {
												d4 = pointCells[id4][a4];
												if (d1 != d4) {
													continue;
												}

												if (d1 == d2 && d1 == d3 && d1 == d4 && d1 != id && d1 != id1 && d1 != id2 && d1 != id3 && d1 != id4) {
													GetKeturkampis(ket, id, id1, id2, id3, id4, d1);
													if (keturkampiai.find(ket) == keturkampiai.end()) {
														keturkampiai[ket] = 1; //pazymime kad toki jau radome (tai yra 2 piramides su bendru pagrindu keturkampiu)
														//pagal istrizainiu nebuvimo indeksus sukonstrojame keturkampi kuri apiesimime tokia tvarka (0,2,1,3), bet mes dar nezinome ar apeiname desines rankos ar kaires rankos taisykle
														myStruct.ids.push_back(keturkapioIstrizaines[0]);
														myStruct.ids.push_back(keturkapioIstrizaines[2]);
														myStruct.ids.push_back(keturkapioIstrizaines[1]);
														myStruct.ids.push_back(keturkapioIstrizaines[3]);

														myStruct.ids.push_back(id);
														myStruct.ids.push_back(d1);
														GetLineID(l, id, id1);
														myStruct.cellIDs.push_back(linijos[l]);
														GetLineID(l, id, id2);
														myStruct.cellIDs.push_back(linijos[l]);
														GetLineID(l, id, id3);
														myStruct.cellIDs.push_back(linijos[l]);
														GetLineID(l, id, id4);
														myStruct.cellIDs.push_back(linijos[l]);
														GetLineID(l, d1, id1);
														myStruct.cellIDs.push_back(linijos[l]);
														GetLineID(l, d1, id2);
														myStruct.cellIDs.push_back(linijos[l]);
														GetLineID(l, d1, id3);
														myStruct.cellIDs.push_back(linijos[l]);
														GetLineID(l, d1, id4);
														myStruct.cellIDs.push_back(linijos[l]);
														strukturos.push_back(myStruct);
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		t.stop();
//////////////////
		t1.start();
		cells = vtkUnstructuredGrid::New();
		cells->Allocate(1, 1);

		vtkPyramid*vtkPyram = vtkPyramid::New();
		vtkIdList*vtkPyramList = vtkIdList::New();
		vtkTetra*vtkTetra = vtkTetra::New();
		vtkIdList*vtkTetraList = vtkIdList::New();
		vtkTetraList->SetNumberOfIds(4);
		vtkPyramList->SetNumberOfIds(5);

		double point1[3];
		double point2[3];
		double point3[3];
		double point4[3];
		double normal[3];
		double virsune1[3];
		double DD, d;
		for (unsigned int i = 0; i < strukturos.size(); i++) {
			if (strukturos[i].ids.size() == 4) {
				input->GetPoint(strukturos[i].ids[0], point1);
				input->GetPoint(strukturos[i].ids[1], point2);
				input->GetPoint(strukturos[i].ids[2], point3);
				input->GetPoint(strukturos[i].ids[3], virsune1);
				vtkTriangle::ComputeNormal(point1, point2, point3, normal);

				d = point1[0] * normal[0] + point1[1] * normal[1] + point1[2] * normal[2];
				DD = (normal[0] * virsune1[0] + normal[1] * virsune1[1] + normal[2] * virsune1[2] + d) / vtkMath::Norm(normal);
				if (DD < 0) {
					swap(strukturos[i].ids[0], strukturos[i].ids[1]);
				}
				vtkTetraList->SetId(0, strukturos[i].ids[0]);
				vtkTetraList->SetId(1, strukturos[i].ids[1]);
				vtkTetraList->SetId(2, strukturos[i].ids[2]);
				vtkTetraList->SetId(3, strukturos[i].ids[3]);
				cells->InsertNextCell(vtkTetra->GetCellType(), vtkTetraList);

			} else {
				input->GetPoint(strukturos[i].ids[0], point1);
				input->GetPoint(strukturos[i].ids[1], point2);
				input->GetPoint(strukturos[i].ids[2], point3);
				input->GetPoint(strukturos[i].ids[3], point4);
				input->GetPoint(strukturos[i].ids[4], virsune1);
				vtkTriangle::ComputeNormal(point1, point2, point3, normal);
				d = point1[0] * normal[0] + point1[1] * normal[1] + point1[2] * normal[2];
				DD = (normal[0] * virsune1[0] + normal[1] * virsune1[1] + normal[2] * virsune1[2] + d) / vtkMath::Norm(normal);
				if (DD < 0) {
					vtkPyramList->SetId(0, strukturos[i].ids[0]);
					vtkPyramList->SetId(1, strukturos[i].ids[3]);
					vtkPyramList->SetId(2, strukturos[i].ids[2]);
					vtkPyramList->SetId(3, strukturos[i].ids[1]);
					vtkPyramList->SetId(4, strukturos[i].ids[4]);
				} else {
					vtkPyramList->SetId(0, strukturos[i].ids[0]);
					vtkPyramList->SetId(1, strukturos[i].ids[1]);
					vtkPyramList->SetId(2, strukturos[i].ids[2]);
					vtkPyramList->SetId(3, strukturos[i].ids[3]);
					vtkPyramList->SetId(4, strukturos[i].ids[4]);
				}
				cells->InsertNextCell(vtkPyram->GetCellType(), vtkPyramList);

				//=========================================================
				if (strukturos[i].ids[5] != -1) {
					input->GetPoint(strukturos[i].ids[5], virsune1);
					DD = (normal[0] * virsune1[0] + normal[1] * virsune1[1] + normal[2] * virsune1[2] + d) / vtkMath::Norm(normal);
					if (DD < 0) {
						vtkPyramList->SetId(0, strukturos[i].ids[0]);
						vtkPyramList->SetId(1, strukturos[i].ids[3]);
						vtkPyramList->SetId(2, strukturos[i].ids[2]);
						vtkPyramList->SetId(3, strukturos[i].ids[1]);
						vtkPyramList->SetId(4, strukturos[i].ids[5]);
					} else {
						vtkPyramList->SetId(0, strukturos[i].ids[0]);
						vtkPyramList->SetId(1, strukturos[i].ids[1]);
						vtkPyramList->SetId(2, strukturos[i].ids[2]);
						vtkPyramList->SetId(3, strukturos[i].ids[3]);
						vtkPyramList->SetId(4, strukturos[i].ids[5]);
					}

					cells->InsertNextCell(vtkPyram->GetCellType(), vtkPyramList);

				}
			}
		}
		t1.stop();
/////////////////
	}
	t.print("vtkCracks3Drev2_Triangulation");
	t1.print("vtkCracks3Drev2_Structure_To_VTK");
	Timeris t2;
	t2.start();
	vtkPoints*points = vtkPoints::New(); //tie patys taskai
	vtkIntArray*newstate = vtkIntArray::New(); //bus laikoma naujas state masyvas
	vtkDataArray*state = input->GetCellData()->GetArray(this->StateArray.c_str()); //gauna nurodyta state masyva
	if (state == NULL) {
		vtkErrorMacro( << "Nenurodytas state masyvas (turi buti idetas i celldata)!\n");
		return 0;
	}
	newstate->SetName(state->GetName());
	newstate->SetNumberOfComponents(1);

	int suma;
	for (unsigned int i = 0; i < strukturos.size(); i++) {
		if (strukturos[i].ids.size() == 4) {
			suma = state->GetTuple1(strukturos[i].cellIDs[0]) + state->GetTuple1(strukturos[i].cellIDs[1]) + state->GetTuple1(strukturos[i].cellIDs[2]) + state->GetTuple1(strukturos[i].cellIDs[3]) + state->GetTuple1(strukturos[i].cellIDs[4]) + state->GetTuple1(strukturos[i].cellIDs[5]);
			suma = 6 - suma;
			newstate->InsertNextTuple1(suma);
		} else {
			suma = state->GetTuple1(strukturos[i].cellIDs[0]) + state->GetTuple1(strukturos[i].cellIDs[1]) + state->GetTuple1(strukturos[i].cellIDs[2]) + state->GetTuple1(strukturos[i].cellIDs[3]) + state->GetTuple1(strukturos[i].cellIDs[4]) + state->GetTuple1(strukturos[i].cellIDs[5]) + state->GetTuple1(strukturos[i].cellIDs[6]) + state->GetTuple1(strukturos[i].cellIDs[7]);
			suma = 8 - suma;
			newstate->InsertNextTuple1(suma);
			//=========================================================
			if (strukturos[i].ids[5] != -1) {
				suma = state->GetTuple1(strukturos[i].cellIDs[0]) + state->GetTuple1(strukturos[i].cellIDs[1]) + state->GetTuple1(strukturos[i].cellIDs[2]) + state->GetTuple1(strukturos[i].cellIDs[3]) + state->GetTuple1(strukturos[i].cellIDs[8]) + state->GetTuple1(strukturos[i].cellIDs[9]) + state->GetTuple1(strukturos[i].cellIDs[10]) + state->GetTuple1(strukturos[i].cellIDs[11]);
				suma = 8 - suma;
				newstate->InsertNextTuple1(suma);
			}
		}
	}
	output->DeepCopy(cells);
	points->DeepCopy(input->GetPoints());
	output->SetPoints(points); //.sudedam mazgus
	points->Delete();
	output->GetCellData()->AddArray(newstate); //idedama nauja ste masyva
	newstate->Delete();
	output->Squeeze();
	t2.stop();
    t2.print("vtkCellBased3Dver2_Attribute_calculation");
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
    t3.print("vtkCellBased3Dver2_vtkExtractCells_veikimo_laikas");
	output1->ShallowCopy(extract->GetOutput());
	extract->Delete();
	t_b.stop();
	t_b.print("vtkCracks3D_bendras");
	return 1;
}

int vtkCellBased3Dver2::FillInputPortInformation(int, vtkInformation *info) {
	info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
	return 1;
}

void vtkCellBased3Dver2::PrintSelf(ostream& os, vtkIndent indent) {
	this->Superclass::PrintSelf(os, indent);

}

void vtkCellBased3Dver2::SetStateArray(string ResultArrayName) {
	this->StateArray = ResultArrayName;
}

const char *vtkCellBased3Dver2::getStateArray() {
	return this->StateArray.c_str();
}
}
}
