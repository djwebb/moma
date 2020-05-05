c====================== include file "param.h" =========================
c
c     main parameter file which sets ocean characteristics:
c
c     imt    = number of grid points in longitudinal direction
c     jmt    = number of grid points in latitudinal direction
c     km     = number of grid points in the vertical direction
c     nt     = number of tracers
c     lseg   = maximum number of longitudinal stream function 
c              segments
c     nisle  = number of islands
c     maxipp = max number of island perimeter points
c
      parameter (imt=92, jmt=36, km=15, nt=2, lseg=5, nisle=2,
     &           maxipp = 1000, imu=imt)
      parameter  (imtp1=imt+1, imtm1=imt-1, imtm2=imt-2, 
     &            imum1=imu-1, imum2=imu-2, jmtp1=jmt+1, 
     &            jmtm1=jmt-1, jmtm2=jmt-2, 
     &            istart=2, iend=imtm1, jscan=jmtm2,
     &            kmp1=km+1, kmp2=km+2, kmm1=km-1)
c
c     add parameter constants
c
#include "pconst.h"
