      subroutine setkmt (alat1, slon1, elon1,
     &                   alat2, slon2, elon2, num)
c
c @(#) SCCS module: setkmt.F  version: 1.9
c     Creation date: 12/18/95
c
c------------------------------------------------------------------
c     set the topography mask "kmt(i,j)" = "num" within the area of 
c     the parallelogram bounded by vertices:
c     (alat1,slon1), (alat1,elon1), (alat2,slon1), & (alat2,elon2)
c------------------------------------------------------------------
c
#include "def_master.h"
#include "param.h"
#include "coord.h"
#include "levind.h"
c
c     convert the four vertices into model indices
c     (js,is1), (js,ie1), (je,is2), (je,ie2)
c
      j1 = indp (alat1, yt, JMT_M)
      j2 = indp (alat2, yt, JMT_M)
      js = min (j1,j2)
      je = max (j1,j2)
c
      i1  = indp (slon1, xt, IMT_M)
      i2  = indp (elon1, xt, IMT_M)
      is1 = min (i1,i2)
      ie1 = max (i1,i2)
c
      i1  = indp (slon2, xt, IMT_M)
      i2  = indp (elon2, xt, IMT_M)
      is2 = min (i1,i2)
      ie2 = max (i1,i2)
c
      is = is1
      ie = ie1
c     print *
c     print *, xt
c     print *, yt
c     print *,'alat1, slon1, elon1 = ',alat1, slon1, elon1
c     print *,'alat2, slon2, elon2 = ',alat2, slon2, elon2
c     print *,' num =',num
c     print *,'js, is1, is2 = ',js, is1, is2
c     print *,'je, ie1, ie2 = ',je, ie1, ie2
c
c     fill in the area bounded by (js,is1), (js,ie1), 
c                                 (je,is2), (je,ie2)
c
      if (js .eq. je) then
        rdj = c1
	iden=1
      else
	iden=je-js
        rdj = c1/(je-js)
      endif
      do 100 j=js,je
        do 90 i=is,ie
          kmt(i,j) = num
90      continue
#ifdef oldkmt
        is = nint(rdj*((j-js)*is2 + (je-j)*is1))
        ie = nint(rdj*((j-js)*ie2 + (je-j)*ie1))
#else
        inum1 = (j-js)*is2 + (je-j)*is1
        inum2 = (j-js)*ie2 + (je-j)*ie1
	is=int(inum1/iden)
	ie=int(inum2/iden)
	if (mod(inum1,iden) .ge. nint(1.0*iden/2)) is=is+1
	if (mod(inum2,iden) .ge. nint(1.0*iden/2)) ie=ie+1
#endif
100   continue
      return
      end
