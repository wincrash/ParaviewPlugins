/*
 * CombinationsUtils.cpp
 *
 *  Created on: Sep 27, 2012
 *      Author: ruslan
 */

#include "CombinationsUtils.h"

namespace vispartdem {

CombinationsUtils::CombinationsUtils() {
	N = 0;
	K = 0;
}

CombinationsUtils::~CombinationsUtils() {
}

template<typename Iterator>
bool CombinationsUtils::next_combination(const Iterator first, Iterator k, const Iterator last) {
	/* Credits: Mark Nelson http://marknelson.us */
	if ((first == last) || (first == k) || (last == k))
		return false;
	Iterator i1 = first;
	Iterator i2 = last;
	++i1;
	if (last == i1)
		return false;
	i1 = last;
	--i1;
	i1 = k;
	--i2;
	while (first != i1) {
		if (*--i1 < *i2) {
			Iterator j = k;
			while (!(*i1 < *j))
				++j;
			std::iter_swap(i1, j);
			++i1;
			++j;
			i2 = k;
			std::rotate(i1, j, last);
			while (last != j) {
				++j;
				++i2;
			}
			std::rotate(k, i2, last);
			return true;
		}
	}
	std::rotate(first, k, last);
	return false;
}

vector<int> CombinationsUtils::combinations(int n, int k) {
	vector<int> seed;
	for (int i = 0; i < n; i++) {
		seed.push_back(i);
	}
	vector<int> data;
	vector<int>::iterator it;
	do {
		for (it = seed.begin(); it != (seed.begin() + k); it++) {
			data.push_back(*it);
		}

	} while (next_combination(seed.begin(), seed.begin() + k, seed.end()));
	return data;
}
void CombinationsUtils::calcCombinations(int n, int k) {
	if(n==N && k==K)
	{
		return;
	}
	this->combination.clear();
    if(n<k)
    {
        return;
    }
	this->combination = combinations(n, k);
	this->N = n;
	this->K = k;
}
} /* namespace vispartdem */
