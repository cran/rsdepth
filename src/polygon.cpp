//============================================================================
// Name        : iotest.cpp
// Author      :
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include "pointline.h"

using namespace std;

//#define MAXPOINTS 5000


int orient(MyLine l, MyPoint p)
{

    double position = ( (l.b.x-l.a.x)*(p.y-l.a.y) - (l.b.y-l.a.y)*(p.x-l.a.x) );
    return ( position > 0 );
}

int isinsidetriangle(MyPoint *list,MyPoint p)
{
    MyLine l1,l2,l3;
    l1.a=list[0];
    l1.b=list[1];

    l2.a=list[1];
    l2.b=list[2];

    l3.a=list[2];
    l3.b=list[0];

    //cout<<l1<<endl<<l2<<endl<<l3<<endl;


    if(orient(l1,p)!=orient(l2,p)|| orient(l1,p)!=orient(l3,p)|| orient(l2,p)!=orient(l3,p) )
    return 0;

    return 1;
}

int isinsidetriangle(MyPoint a, MyPoint b,MyPoint c,MyPoint p)
{
	MyPoint list[3];
	list[0]=a;
	list[1]=b;
	list[2]=c;
	return isinsidetriangle(list,p);
}


int IsInsidePolygon(MyPoint *list, MyPoint pt, int sz)
{
    int start=0,mid=(sz-1)/2,last=sz-1;
    MyLine l;
    

    while( (last-start) > 2)
    {
		// There is some problem.
		if(mid < 0 || mid >= sz )
		return 1;
		
        l.a = list[0];
        l.b = list[mid];

        if( orient(l,pt) == orient(l,list[1]) )
        {
            last=mid;
            mid = mid-(last-start)/2;
        }
        else
        {
            start=mid;
            mid = mid+(last-start)/2;
        }



        //cout<<start<<", "<<mid<<", "<<last<<"\n";
    }


    //cout<<start<<", "<<mid<<", "<<last<<"\n";
    //cout<<list[start]<<"\n";

    if(last-start==1)
    	return isinsidetriangle(list[0],list[start],list[start+1],pt);

    return (isinsidetriangle(list[0],list[start],list[start+1],pt)||
    		isinsidetriangle(list[0],list[start+1],list[start+2],pt) );

}

int PolygonIntersection(MyPoint *p1, MyPoint *p2, double *p3x, double *p3y, int n, int m)
{
	int count = 0;

    
	for(int i=0;i<n;i++)
	{
		MyLine l1;

		if(IsInsidePolygon(p2,p1[i],m))
		{
		    //cout<<p1[i]<<endl;
			p3x[count]=p1[i].x;
			p3y[count]=p1[i].y;
			count++;
			
			// we don't have enough memory
			if(count >= 2*(n+m) )
			 return count;
		}
//		cout<<count<<" "<<i<<"  "<<n<<"   "<<m<<endl;

		l1.a=p1[i];
		if(i<n-1)
			l1.b=p1[i+1];
		else
			l1.b=p1[0];

		for(int j=0;j<m;j++)
		{


			if(i==0 && IsInsidePolygon(p1,p2[j],n))
			{
			    //cout<<p2[j]<<endl;
				p3x[count]=p2[j].x;
				p3y[count]=p2[j].y;
				count++;
			// we don't have enough memory
			if(count >= 2*(n+m) )
			 return count;
				
			}


		    if(p2[j]==l1.a||p2[j]==l1.b||p2[j+1]==l1.a||p2[j+1]==l1.b)
		    continue;

			MyLine l2;

			l2.a=p2[j];
			if(j<m-1)
				l2.b=p2[j+1];
			else
				l2.b=p2[0];

			MyPoint temp = Intersects(l1,l2);

			if(temp.x!=MAXINT||temp.y!=MAXINT)
			{
				//if(IsInsidePolygon(p2,temp,m)&&
				//	IsInsidePolygon(p1,temp,n)	)
				{
	//			    cout<<l1<<endl<<l2<<endl;
					p3x[count]=temp.x;
					p3y[count]=temp.y;
					count++;
			// we don't have enough memory
			if(count >= 2*(n+m) )
			 return count;
					
				}
			}

		}

	}
	return count;

}

extern "C" {

void polygonintersection(double *p1x, double *p1y,
		double *p2x, double *p2y,
		double *p3x, double *p3y,
		int *n, int *m, int *p)
{
	MyPoint *p1,*p2;
	MyPoint first_list[5001], second_list[5001];
	p1 = first_list;
	p2 = second_list;
	
	//p1=(MyPoint *) malloc(sizeof(MyPoint)*(n[0]+2) );
	//p2=(MyPoint *) malloc(sizeof(MyPoint)*(m[0])+2) ;
	
	// Improper use
	if(n[0]< 3 || n[0] > 5000 || m[0]< 3 || m[0] > 5000)
	return;

	for(int i=0;i<n[0];i++)
	{
		p1[i].x=p1x[i];
		p1[i].y=p1y[i];
		//cout<<p1[i]<<endl;
	}
	
	p1[n[0]].x=p1[0].x;
	p1[n[0]].y=p1[0].y;
	
	//cout<<endl<<endl;

	for(int i=0;i<m[0];i++)
	{
		p2[i].x=p2x[i];
		p2[i].y=p2y[i];
		//cout<<p2[i]<<endl;
	}
	p2[m[0]].x=p2[0].x;
	p2[m[0]].y=p2[0].y;
	
	p[0] = PolygonIntersection(p1, p2, p3x, p3y,  n[0], m[0]);
	
	//free(p1);
	//free(p2);
}

} // extern "C"

