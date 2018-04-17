/*
 * CellSearch.h
 *
 *  Created on: Sep 19, 2011
 *      Author: ruslan
 */

#ifndef CELLSEARCH_H_
#define CELLSEARCH_H_
#include <vector>
#include "vtkType.h"
#include <iostream>
#include <algorithm>
#include "Timeris.h"
using namespace std;

namespace vispartdem {
typedef vtkIdType INT;
typedef vector< vector<INT> > cell;

class CellSearch {
public:
	CellSearch(INT pointCount,INT cellCount,INT* connections);
    CellSearch();
    void Fill(INT pointCount,INT cellCount,INT* connections,int times);
	virtual ~CellSearch();
	bool cellExist(INT &id1,INT &id2);
	INT cellID(INT &id1,INT &id2);
	vector<INT>* getPointNeighbours(INT & id1);
	vector<INT>* getPointNeighboursCells(INT & id1);
	long long getHash(long x,long y);
protected:
	cell data;
	cell dataID;
    int onetime;
};

}

#endif /* CELLSEARCH_H_ */
