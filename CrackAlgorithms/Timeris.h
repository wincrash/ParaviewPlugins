/* 
 * File:   Timeris.h
 * Author: ruslan
 *
 * Created on June 21, 2010, 10:33 AM
 */

#ifndef _TIMERIS_H
#define	_TIMERIS_H
#include <iostream>
#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>

using namespace std;
class Timeris {
public:
	Timeris();
	Timeris(const Timeris& orig);
	void start();
	void stop();
	void print(string label);
    void printWithout(string label);
	virtual ~Timeris();
private:
	struct timeval stage1_start, stage1_stop, stage1_total, stage1_sum;
};

#endif	/* _TIMERIS_H */

