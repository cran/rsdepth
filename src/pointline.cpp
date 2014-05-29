#include "pointline.h"

/*
ostream& operator <<(ostream &os,const MyPoint &obj)
{
      os<<"("<<obj.x<<", "<<obj.y<<")";
	//os<<obj.x<<endl<<obj.y;
      return os;
}

*/
bool operator<(const MyPoint &obj,const MyPoint &obj1)
{
	return obj.angle<obj1.angle;
}
/*
ostream& operator <<(ostream &os,const MyLine &obj)
{
      os<<obj.a<<"---"<<obj.b;
      return os;
}
*/

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
	intersect.x=intersect.y=MAXINT;

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
			intersect.x=MAXINT;
			intersect.y = MAXINT;
		}

	return intersect;

}

