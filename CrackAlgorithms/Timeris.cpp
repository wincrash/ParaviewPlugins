/* 
 * File:   Timeris.cpp
 * Author: ruslan
 * 
 * Created on June 21, 2010, 10:33 AM
 */

#include "Timeris.h"

Timeris::Timeris() {
	timerclear(&stage1_sum);
}

Timeris::Timeris(const Timeris& orig) {
}

Timeris::~Timeris() {
}

void Timeris::start() {
	gettimeofday(&stage1_start, NULL);
}

void Timeris::stop() {
	gettimeofday(&stage1_stop, NULL);
	timersub(&stage1_stop, &stage1_start, &stage1_total);
	timeradd(&stage1_sum, &stage1_total, &stage1_sum);
	// printf("\n%li.%06li s", stage1_total.tv_sec, stage1_total.tv_usec);
}

void Timeris::print(string label) {
    printf("%s %li.%06li\n", label.c_str(), stage1_sum.tv_sec,
            stage1_sum.tv_usec);
}
void Timeris::printWithout(string label) {
    printf("%s %li.%06li\n", label.c_str(), stage1_sum.tv_sec,
            stage1_sum.tv_usec);
}
