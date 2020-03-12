#ifndef _SIMPLEMATH_H_
#define _SIMPLEMATH_H_

#include <cmath>
#include <iostream>
#include "object.h"
struct sol
{
    bool solvable;
    double solution[2];
};

sol QuadEq(double a, double b, double c);

double minpositive(double a, double b);

int gettreelevel(int i);

#endif