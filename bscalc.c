#include "ubiq.h"
#include "rkf.h"

point biot(point stseg, point endseg, point testpt)
{
      point l1, l2, l3, deldot;
      double dot11, dot12, dot22, dot33, dot23, factor1, factor2;

      l1 = testpt - stseg;
      l2 = endseg - stseg;
      l3 = testpt - endseg;      
      
      dot11 = dot(l1, l1);
      dot12 = dot(l1, l2);
      dot22 = dot(l2, l2);
      dot23 = dot12 - dot22;
      dot33 = dot(l3, l3);

      factor1 = dot11*dot22-dot12*dot12;

      if (fabs(factor1/dot11/dot22) > 1e-24)
      {
        factor2 = dot23/sqrt(dot33) - dot12/sqrt(dot11);
        deldot = (factor2 / factor1) * cross(l1, l2);
      }
      else
      {
        factor1 = dot22*dot33-dot23*dot23;
        if  (fabs(factor1/dot22/dot33) > 1e-24)
        {
            factor2 = dot23/sqrt(dot33) - dot12/sqrt(dot11);
            deldot = (factor2 / factor1) * cross(l3, l2);
        }
        else
        {
            if (dot11*dot33 > 1e-12) {
                factor2 = fabs(dot11-dot33)/dot11/dot33/sqrt(dot22);
                deldot = (factor2 / 2) * cross(l1, l2);
            }
            else {
                deldot.x = 0; deldot.y = 0; deldot.z = 0;
            }
        }
      }
      return deldot;
}
