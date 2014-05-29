// RSsimulation.cpp : Defines the entry point for the console application.
//
#include "pointline.h"

#define MAXPOINTS 10000

using namespace std;


class MyPoint;
class MyLine;


int numberOfPoints = 0;
MyPoint P[MAXPOINTS], q;
double Depth;
int algo;
int did2,did1=0;
int m=10,n=60;
char ch;

bool depthcompare(const MyPoint &obj,const MyPoint &obj1)
{
	return RSDepth(obj)<RSDepth(obj1);
}


int IsUnique(MyPoint p)
{
    for(int i=0;i<numberOfPoints;i++)    {
        if (p.x == P[i].x && p.y == P[i].y)
            return 0;
    }

    return 1;
}

double RSDepth(MyPoint q)
{
	double* angle;
	int *total;

	angle = (double*) malloc(sizeof(double)*numberOfPoints*2);
	total = (int*) malloc(sizeof(int)*numberOfPoints*2);;


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
	double *ray;
	ray = (double*) malloc(sizeof(double)*numberOfPoints*2);

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
    free(angle);
    free(ray);
    free(total);

    return minDepth;

}

double RSDepth(double qx, double qy)
{
	MyPoint pt;
	pt.x=qx;
	pt.y=qy;
	return RSDepth(pt);
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

					if (dum.x != MAXINT && dum.y != MAXINT)  {
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


MyPoint centroid;

MyPoint compute2DPolygonCentroid(double *x,double *y, int vertexCount)
{
    centroid.x=0;
    centroid.y=0;
    double signedArea = 0.0;
    double x0 = 0.0; // Current vertex X
    double y0 = 0.0; // Current vertex Y
    double x1 = 0.0; // Next vertex X
    double y1 = 0.0; // Next vertex Y
    double a = 0.0;  // Partial signed area

    // For all vertices except last
    int i=0;
    for (i=0; i<vertexCount-1; ++i)
    {
        x0 = x[i];
        y0 = y[i];
        x1 = x[i+1];
        y1 = y[i+1];
        a = x0*y1 - x1*y0;
        signedArea += a;
        centroid.x += (x0 + x1)*a;
        centroid.y += (y0 + y1)*a;
    }

    // Do last vertex
    x0 = x[i];
    y0 = y[i];
    x1 = x[0];
    y1 = y[0];
    a = x0*y1 - x1*y0;
    signedArea += a;
    centroid.x += (x0 + x1)*a;
    centroid.y += (y0 + y1)*a;

    signedArea *= 0.5;
    centroid.x /= (6.0*signedArea);
    centroid.y /= (6.0*signedArea);

    return centroid;
}

int FillBag(double *bag,double medDepth, double *center,int *cSize)
{
	int bagSize=0;
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

						if (dum.x != MAXINT && dum.y != MAXINT)  {
	                        t =RSDepth(dum.x , dum.y);

	                        if (t > (medDepth+numberOfPoints-1) && IsUnique(dum))
	                        {
	                        	bag[bagSize++]= dum.x;
	                        	bag[bagSize++]= dum.y;

	                        	if(t >= (numberOfPoints*numberOfPoints)/9.0)
	                        	{
	                        	    center[cSize[0]++]= dum.x;
                                    center[cSize[0]++]= dum.y;
	                        	}
	                        	max = t;ok = dum;
	                        }
	                    }

	                }
	            }

	        }

	    }


	return bagSize;

}


void GetCenter(double *center,int *cSize)
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

						if (dum.x != MAXINT && dum.y != MAXINT)  {
	                        t =RSDepth(dum.x , dum.y);

	                        	if(t >= (numberOfPoints*numberOfPoints)/9.0)
	                        	{
	                        	    center[cSize[0]++]= dum.x;
                                    center[cSize[0]++]= dum.y;
                                    //cout<<cSize[0]<<endl;
	                        	}

	                        	max = t;ok = dum;

	                    }

	                }
	            }

	        }

	    }


return ;
}


extern "C" {

void rs_centroid(double *x,double *y,int* count,double*cent)
{

    MyPoint centroid = compute2DPolygonCentroid(x,y, count[0]);
    cent[0]=centroid.x;
    cent[1]=centroid.y;
}
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

void rs_depthrings(double *inx,double *iny,double *outx, double *outy,int *size)
{
    int count=0;
	numberOfPoints=size[0];
	MyPoint*Q;
	Q=(MyPoint*)malloc(sizeof(MyPoint)*(numberOfPoints+1));

	for(int i=0;i<numberOfPoints;i++)
	{
		P[i].x=inx[count];
		P[i].y=iny[count];
		Q[i].x=inx[count];
		Q[i].y=iny[count];
		count++;
	}

    sort(Q,Q+numberOfPoints,depthcompare);
    for(int i=0;i<numberOfPoints;i++)
	{
		outx[i]=Q[i].x;
		outy[i]=Q[i].y;
	}

	free(Q);

    return;
}

void rs_getbag(double *ptX, double *ptY, double *bag, int *num, int *bagsz, double *center, int *cSize)
{

	int count=0;
	std::vector<double> depthList;
	numberOfPoints=num[0];

	for(int i=0;i<numberOfPoints;i++)
	{
		P[i].x=ptX[count];
		P[i].y=ptY[count];
		count++;
	}

	for(int i=0;i<numberOfPoints;i++)
	{
		depthList.push_back(RSDepth(P[i]));
		//printf("%d. %f\n",i,RSDepth(P[i]));
	}

	std::vector<double>::iterator first = depthList.begin();
	std::vector<double>::iterator last = depthList.end();
	std::vector<double>::iterator middle = first + (last - first) / 2;
	//middle = middle + (last - first) / 4;
	std::nth_element(first, middle, last); // can specify comparator as optional 4th arg
	double median = *middle;
	int bagSize = FillBag(bag,median,center,cSize);
	//printf("\ndepth of median is :%f\n",median);
	bagsz[0]=bagSize;
	return;

}


void rs_getcenter(double *ptX, double *ptY, int *num, double *center, int *cSize)
{
	int count=0;
	std::vector<double> depthList;
	numberOfPoints=num[0];

	for(int i=0;i<numberOfPoints;i++)
	{
		P[i].x=ptX[count];
		P[i].y=ptY[count];
		count++;
	}

	GetCenter(center,cSize);
	//printf("\ndepth of median is :%f\n",median);
	return;

}


} // extern "C"

