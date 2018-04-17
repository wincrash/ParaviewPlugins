/*
 * CombinationsUtils.h
 *
 *  Created on: Sep 27, 2012
 *      Author: ruslan
 */

#ifndef COMBINATIONSUTILS_H_
#define COMBINATIONSUTILS_H_
#include <vector>
#include <algorithm>

namespace vispartdem {

using namespace std;
class CombinationsUtils {
public:
	CombinationsUtils();
	virtual ~CombinationsUtils();
	template<typename Iterator>
	bool next_combination(const Iterator first, Iterator k, const Iterator last);
	vector<int> combinations(int n, int k);
	void calcCombinations(int n, int k);
	int N;
	int K;
	vector<int> combination;
};

} /* namespace vispartdem */
#endif /* COMBINATIONSUTILS_H_ */
