c====================== include file "param.h" =========================
#ifdef hcomments
c
c     @(#) SCCS module: param.h, version 1.16
c
c     main parameter file which sets ocean characteristics:
c
c     imt    = number of grid points in longitudinal direction
c     jmt    = number of grid points in latitudinal direction
c     km     = number of grid points in the vertical direction
c     nt     = number of tracers
c     lseg   = maximum number of longitudinal stream function 
c              segments
c
#endif
      parameter  (imt=92, jmt=36, km=15, nt=2, imu=imt)
      parameter  (imtp1=imt+1, imtm1=imt-1, imtm2=imt-2, 
     &            imum1=imu-1, imum2=imu-2, jmtp1=jmt+1, 
     &            jmtm1=jmt-1, jmtm2=jmt-2, 
     &            istart=2, iend=imtm1, jscan=jmtm2,
     &            kmp1=km+1, kmp2=km+2, kmm1=km-1)
#ifdef hcomments
c
c     add parameter constants
c
#endif
#include "pconst.h"
