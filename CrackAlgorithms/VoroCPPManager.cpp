/*
 * VoroCPPManager.cpp
 *
 *  Created on: Sep 24, 2012
 *      Author: ruslan
 */

#include "VoroCPPManager.h"

#include "voro/voro++.cpp"
VoroCPPManager::VoroCPPManager() {
	voronoicell_neighbor *v = new voronoicell_neighbor;
	this->voronoi_cell = v;

}
void VoroCPPManager::Voro_Cell_init(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax) {
	voronoicell_neighbor *v = static_cast<voronoicell_neighbor*>(voronoi_cell);
	v->init(xmin, xmax, ymin, ymax, zmin, zmax);
}
void VoroCPPManager::nplane(double x, double y, double z, long i) {
	voronoicell_neighbor *v = static_cast<voronoicell_neighbor*>(voronoi_cell);
	v->nplane(x, y, z, i);
}
vector<vector<int> > VoroCPPManager::output_vetice_relations() {
	voronoicell_neighbor *v = static_cast<voronoicell_neighbor*>(voronoi_cell);
	vector<vector<int> > realations;
	v->print_edges(realations);
	return realations;
}
void VoroCPPManager::output_neighbors(vector<int> &kaimynai) {
	voronoicell_neighbor *v = static_cast<voronoicell_neighbor*>(voronoi_cell);
	v->output_neighbors(kaimynai);
}
void VoroCPPManager::output_face_orders(vector<int> &faces_kiekis) {
	voronoicell_neighbor *v = static_cast<voronoicell_neighbor*>(voronoi_cell);
	v->output_face_orders(faces_kiekis);
}
void VoroCPPManager::output_face_vertices(int **indexes) {
	voronoicell_neighbor *v = static_cast<voronoicell_neighbor*>(voronoi_cell);
	v->output_face_vertices(indexes);
}
int VoroCPPManager::Voro_cell_p() {
	voronoicell_neighbor *v = static_cast<voronoicell_neighbor*>(voronoi_cell);
	return v->p;
}
double* VoroCPPManager::Voro_pts() {
	voronoicell_neighbor *v = static_cast<voronoicell_neighbor*>(voronoi_cell);
	return v->pts;

}
cele3d *VoroCPPManager::Voro_Cont_getInfo3D(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax,vtkPolyData*input) {
	container c(xmin, xmax, ymin, ymax, zmin, zmax, 8, 8, 8, false, false, false, 8);
	cele3d *faces = new cele3d[input->GetNumberOfPoints()];
	for (int i = 0; i < input->GetNumberOfPoints(); ++i) {
		c.put(i,input->GetPoint(i)[0],input->GetPoint(i)[1],input->GetPoint(i)[2]);
		faces[i].nrfaces = 0;
		faces[i].id = -1;
		faces[i].nrvertext = 0;
	}
	c.print_all_custom("%i %w %s %a %n %t %P", faces);
	return faces;
}
cele3d *VoroCPPManager::Voro_Cont_getInfo2D(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax,vtkPolyData*input) {
	container c(xmin, xmax, ymin, ymax, zmin, zmax, 8, 8, 8, false, false, false, 8);
	cele3d *faces = new cele3d[input->GetNumberOfPoints()];
	for (int i = 0; i < input->GetNumberOfPoints(); ++i) {
		c.put(i,input->GetPoint(i)[0],input->GetPoint(i)[1],0.0);
		faces[i].nrfaces = 0;
		faces[i].id = -1;
		faces[i].nrvertext = 0;
	}
	c.print_all_custom("%i %w %s %a %n %t %P", faces);
	return faces;
}

VoroCPPManager::~VoroCPPManager() {
	// TODO Auto-generated destructor stub
}

