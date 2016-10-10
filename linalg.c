#include "linalg.h"


ostream& operator<< (ostream& outwhere, const point& p)
{
    outwhere << p.x <<" "<<p.y<<" "<<p.z;
    return outwhere;
}

istream& operator>> (istream& infrom, point& p)
{
    infrom >> p.x;
    infrom >> p.y;
    infrom >> p.z;
    return infrom;
}

int operator== (point a, point b)
{
    int test=1;
    if (a.x != b.x) test=0;
    if (a.y != b.y) test=0;
    if (a.z != b.z) test=0;
    return test;
}

/* Inverts a point in a circle about the origin */

point invert(point testpt, double radius)
{
   
   point result;
   double scale;
   
   scale = radius*radius/(testpt.x*testpt.x+testpt.y*testpt.y);   
   result.x = scale * testpt.x;
   result.y = scale * testpt.y;
   result.z = testpt.z;
   return result;
}

/* Rotates point about z axis theta degrees counterclockwise */

point rotz(double theta, point coords)
{
    point newcoords;
    newcoords.x=cos(theta)*coords.x-sin(theta)*coords.y;
    newcoords.y=cos(theta)*coords.y+sin(theta)*coords.x;
    newcoords.z=coords.z;
    return newcoords;
}
