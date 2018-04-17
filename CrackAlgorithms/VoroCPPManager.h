/*
 * VoroCPPManager.h
 *
 *  Created on: Sep 24, 2012
 *      Author: ruslan
 */

#ifndef VOROCPPMANAGER_H_
#define VOROCPPMANAGER_H_
#include <vector>
#include "vtkPolyData.h"
using namespace std;

struct cele3d {
	int id;
	int nrvertext;
	int nrfaces;
	int* facevertextnr;
	int * kaimynai;
	int **faceorder;
	float** vertexes;
	~cele3d() {

		int k = 0;

		for (k = 0; k < nrfaces; k++) {
			delete[] faceorder[k];

		}

		for (k = 0; k < nrvertext; k++) {

			delete[] vertexes[k];

		}
		if (nrfaces != 0 && nrvertext != 0) {
			delete[] vertexes;
			delete[] kaimynai;
			delete[] facevertextnr;
			delete[] faceorder;
		}

	}
	;
	cele3d() {
	}
	;

};
class VoroCPPManager {
public:
	VoroCPPManager();
	virtual ~VoroCPPManager();
	void Voro_Cell_init(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax);
	void nplane(double x, double y, double z, long i);
	void output_neighbors(vector<int> &kaimynai);
	void output_face_orders(vector<int> &faces_kiekis);
	void output_face_vertices(int **indexes);
	vector<vector<int> > output_vetice_relations();
	int Voro_cell_p();
	double* Voro_pts();
	cele3d * Voro_Cont_getInfo2D(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax,vtkPolyData*input);
	cele3d * Voro_Cont_getInfo3D(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax,vtkPolyData*input);

private:
	void* voronoi_cell;

};

#endif /* VOROCPPMANAGER_H_ */
