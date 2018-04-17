/*
 * CellSearch.cpp
 *
 *  Created on: Sep 19, 2011
 *      Author: ruslan
 */

#include "CellSearch.h"

namespace vispartdem {


CellSearch::CellSearch()
{
    this->onetime=0;

}

void CellSearch::Fill(INT pointCount,INT cellCount,INT* connections,int times)
{
    if(this->onetime!=times)
    {
    Timeris t;
    t.start();
    vector<INT> temp;
    temp.reserve(20);
    INT tempID;
    INT id1, id2;
    data.resize(pointCount, temp);
    dataID.resize(pointCount, temp);
    for (int i = 0; i < cellCount; i++) {
        tempID = i * 3;
        id1 = connections[tempID + 1];
        id2 = connections[tempID + 2];
        data[id1].push_back(id2);
        data[id2].push_back(id1);
        dataID[id1].push_back(i);
        dataID[id2].push_back(i);
    }
    t.stop();
    t.printWithout("CellSearch_uzpildymas_info");
    this->onetime++;
    }
}




CellSearch::CellSearch(INT pointCount, INT cellCount, INT* connections) {
        Fill(pointCount,cellCount,connections,1);

}

CellSearch::~CellSearch() {
	data.clear();
	dataID.clear();
}

INT CellSearch::cellID(INT & id1, INT & id2) {
	for (unsigned int i = 0; i < data[id1].size(); i++) {
		if (data[id1][i] == id2) {
			return dataID[id1][i];
		}
	}
	return -1;
}

bool CellSearch::cellExist(INT & id1, INT & id2) {
	if (cellID(id1, id2) == -1) {
		return false;
	} else {
		return true;
	}
}
vector<INT> * CellSearch::getPointNeighboursCells(INT & id1) {
	return &dataID[id1];
}

vector<INT> * CellSearch::getPointNeighbours(INT & id1) {
	return &data[id1];
}


}
