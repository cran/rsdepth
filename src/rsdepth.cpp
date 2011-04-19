// RSsimulation.cpp : Defines the entry point for the console application.
//
#include <R.h>
#include <fstream>
#include <algorithm>
#include <time.h>
#include <math.h>
#define MAXPOINTS 100000

using namespace std;

class MyPoint
{
public:
	double x,y,angle;
	int operator ==(MyPoint two)
	{
		if(two.x==x &&two.y==y)
			return 1;
		return 0;
	}
	friend ostream& operator <<(ostream &os,const MyPoint &obj);
	friend bool operator<(const MyPoint &obj,const MyPoint &obj1);
};

ostream& operator <<(ostream &os,const MyPoint &obj)
{
      //os<<"("<<obj.x<<", "<<obj.x<<")";
	os<<obj.x<<endl<<obj.y;
      return os;
}


bool operator<(const MyPoint &obj,const MyPoint &obj1)
{
	return obj.angle<obj1.angle;
}


class MyLine
{
public:
	MyPoint a,b;
};

int numberOfPoints = 0;
MyPoint P[MAXPOINTS], q;
double Depth;
int algo;
int did2,did1;
int m=10,n=60;
char ch;



double RSDepth(MyPoint q)
{
	double angle[MAXPOINTS];

	int total[MAXPOINTS];

	//int endOf[MAXPOINTS];
	//ofstream oo("dump.txt");

	//plot points on unit cicle around q and calculate angles
	for(int i=0;i<numberOfPoints;i++)
	{
		double slope ;
		if(q.y-P[i].y==0&&(q.x-P[i].x)==0)
			slope=0;
		else if ((q.x-P[i].x)==0)
			slope = 99999999999999999999.0;
		else
		slope = (q.y-P[i].y)/(q.x-P[i].x);
		angle[i]=atan(slope)*57.35;

		if (angle[i] < 0)   {
            angle[i] += 180;
        }

        if (q.y > P[i].y) {
            angle[i] += 180;
        }
		P[i].angle=angle[i];
	}

	//sort
	sort(angle,angle+numberOfPoints);
	sort(P,P+numberOfPoints);


	//Fill array from number of points between a point and its antipodal

    int count = 0;//points in between
    int index=0;//of points for which we have found antipodals
	for (int i = 0; index < numberOfPoints  ; i++)
    {
		count++;
        if ( (i-index)>= numberOfPoints ||angle[i % numberOfPoints ] > (angle[index] + 180) ||

			((angle[i % numberOfPoints ] < (angle[index]) && (angle[i % numberOfPoints ] + 360 > (angle[index] + 180)))))

		{
			count--;

			total[index] = count - 1;
            //endOf[index] = i - 1;
            index++;
            count--;
            i--;

        }
    }

	//calculate depth for first ray between 1 and n points.

    double minDepth = 100000;

	double ray[MAXPOINTS];
    ray[0] = 0;
	did2=0;
	double current_angle=0,largest_angle = 0;
    for (int i = 0; i <= total[0]; i++)
		ray[0] += numberOfPoints - ( (total[i])+ i +1);
    minDepth = ray[0];

	if(ray[0]>=1.0/9.0*numberOfPoints*numberOfPoints)
	{
		//oo<<"1"<<endl<<P[numberOfPoints-1]<<endl<<P[0]<<endl;
		current_angle+=angle[0]+360-angle[numberOfPoints-1];
	}
	else;
		//oo<<"0"<<endl<<P[numberOfPoints-1]<<endl<<P[0]<<endl;


	for(int i=1;i<numberOfPoints+1;i++)
    {
        ray[i] = ray[i - 1] - numberOfPoints + 2*total[i-1]+1;

		if(ray[i]>=1.0/9.0*numberOfPoints*numberOfPoints)
		{
			/*
			if(i==numberOfPoints)
				oo<<"1"<<endl<<P[i-1]<<endl<<P[0]<<endl;
			else
				oo<<"1"<<endl<<P[i-1]<<endl<<P[i]<<endl;

				*/
			current_angle+=angle[i]-angle[i-1];
		}
		else
		{
			/*
			if(i==numberOfPoints)
				oo<<"0"<<endl<<P[i-1]<<endl<<P[0]<<endl;
			else
				oo<<"0"<<endl<<P[i-1]<<endl<<P[i]<<endl;
			*/
			if(current_angle > largest_angle)
			{
				largest_angle=current_angle;
				current_angle=0;

			}
		}


		if (minDepth > ray[i])
		{
            minDepth = ray[i];
			did2=i;
		}
    }
//	oo.close();
    return minDepth;

}

double RSDepth(double qx, double qy)
{
	MyPoint pt;
	pt.x=qx;
	pt.y=qy;
	return RSDepth(pt);
}
int InBoundingBox(MyLine a, MyPoint p)
{


    if (p.x <min(a.a.x, a.b.x))
        return 0;

    if (p.x > max(a.a.x, a.b.x))
        return 0;

    if (p.y <min(a.a.y, a.b.y))
        return 0;

    if (p.y > max(a.a.y, a.b.y))
        return 0;

    return 1;

}
MyPoint Intersects( MyLine first, MyLine second)
{

	MyPoint intersect;
	intersect.x=intersect.y=0;

	//this line

	double a1 = (first.b.y - first.a.y);

	double b1 = (first.a.x - first.b.x);

	double c1 = (first.b.x*first.a.y - first.a.x*first.b.y);


	//Line to test against
	double a2 = (second.b.y - second.a.y);
	double b2 = (second.a.x - second.b.x);
	double c2 = (second.b.x*second.a.y - second.a.x*second.b.y);

	double denom = a1*b2 - a2*b1;

	//Check for parallel lines
	if(denom == 0) { return intersect; }

	//Get the intersection point
	intersect.x = ( (b1*c2 - b2*c1)/denom);
	intersect.y = ( (a2*c1 - a1*c2)/denom);

	//InBoundingBox tests to see if the point is in the bounding box for
		//a line segment. The point must be in both lines' bounding boxes to
		//register and intersection.
	if( InBoundingBox( first,intersect )==0 || InBoundingBox( second,intersect )==0 )
		{
			intersect.x=0;
			intersect.y = 0;
		}

	return intersect;

}

int IsUnique(MyPoint p)
{
    for(int i=0;i<numberOfPoints;i++)    {
        if (p.x == P[i].x && p.y == P[i].y)
            return 0;
    }

    return 1;
}


 MyPoint RSMedian()
{
    MyPoint dum,ok;
    MyLine l1;
	MyLine l2;

	l1.a=P[0];
	l1.b=P[1];
	l1.a=P[2];
	l1.b=P[3];

    double max = 0;
    double t = 0;

    for (int i = 0; i < numberOfPoints; i++)
    {
        l1.a = P[i];
        for (int j = i + 1; j < numberOfPoints; j++)
        {
            l1.b = P[j];

            for (int k = 0; k < numberOfPoints; k++)
            {
                if (i == k || j == k)
                    continue;
                l2.a = P[k];

                for (int l = k + 1; l < numberOfPoints; l++)
                {

                    if (i == l || j == l)
                        continue;
                    l2.b = P[l];
                    dum = Intersects(l1, l2);

					if (dum.x != 0 && dum.y != 0)  {
                        t =RSDepth(dum.x , dum.y);
                        if (t > max && IsUnique(dum)) { max = t;ok = dum;}
                    }

                }
            }

        }

    }

	return ok;
}

 ifstream ii;

 void parse2()
{
	double one;

	for(int i=0;i<numberOfPoints;i++)
	{
		ch='0';
		while(ch!=']')
			ii>>ch;

		ii>>one;
		P[i].x=one;
		ii>>one;
		P[i].y=one;
	}
}

 void parse()
{
	ii.open("C:\\works\\60\\cauchy_60.txt");
	ofstream oo("C:\\works\\60\\cauchy.txt");
	numberOfPoints=n;
	MyPoint temp;

	for(int i=0;i<m;i++)
	{
		ch='0';
		while(ch!=']')
			ii>>ch;
		ch='0';
		while(ch!=']')
			ii>>ch;
		parse2();
		temp = RSMedian();
		oo<<temp.x<<"   ";
		oo<<temp.y<<endl;
		double th = (numberOfPoints*numberOfPoints)/RSDepth(temp);
		oo<<th<<endl;
	}


	ii.close();
	oo.close();
	return;
}

extern "C" {

void rs_depth(double *x,double *y,double *p, double *dp,int *nn)
{
	int count=0;
	numberOfPoints=nn[0];
	double dd;
	for(int i=0;i<numberOfPoints;i++)
	{
		P[i].x=x[count];
		P[i].y=y[count];
		count++;
	}
	dd=RSDepth(p[0],p[1]);
	dp[0]=dd;
}
void rs_med(double *x, double *y, double *md, int *nn)
{
	int count=0;
	numberOfPoints=nn[0];
	MyPoint dd;
	for(int i=0;i<numberOfPoints;i++)
	{
		P[i].x=x[count];
		P[i].y=y[count];
		count++;
	}
	dd=RSMedian();
	md[0]=dd.x;
	md[1]=dd.y;
}

} // extern "C"

/*
int _tmain(int argc, _TCHAR* argv[])
{
	parse();
	return 0;
}
*/

