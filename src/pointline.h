
#ifndef __POINTLINE_H
#define __POINTLINE_H

//#include <R.h>
#include <fstream>
#include <algorithm>
#include <time.h>
#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <limits.h>

#define MAXINT INT_MAX

using namespace std;

class MyPoint
{
public:
    int comparator;
	double x,y,angle;
	int operator ==(MyPoint two)
	{
		if(two.x==x &&two.y==y)
			return 1;
		return 0;
	}
//	friend ostream& operator <<(ostream &os,const MyPoint &obj);
	friend bool operator<(const MyPoint &obj,const MyPoint &obj1);
};

//ostream& operator <<(ostream &os,const MyPoint &obj);


bool operator<(const MyPoint &obj,const MyPoint &obj1);


class MyLine
{
public:
	MyPoint a,b;
//	friend ostream& operator <<(ostream &os,const MyPoint &obj);
};

//ostream& operator <<(ostream &os,const MyLine &obj);

int IsUnique(MyPoint p);
MyPoint Intersects( MyLine first, MyLine second);
int InBoundingBox(MyLine a, MyPoint p);
double RSDepth(MyPoint q);
double RSDepth(double qx, double qy);
#endif
