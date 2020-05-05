c====================== include file "param.h" =========================
#ifdef hcomments
c
c @(#) SCCS module: param.h  version: 1.9
c     Creation date: 12/18/95
c
c-----------------------------------------------------------------------
c     Main parameter file which sets ocean characteristics:
c
c     IMT_M    = number of grid points in longitudinal direction(master)
c     JMT_M    = number of grid points in latitudinal direction (master)
c     IMT_S    = maximum number of grid points in longitudinal direction
c                                                               (slave)
c     JMT_S    = maximum number of grid points in latitudinal direction
c                                                               (slave)
c     IMT_L    = maximum of IMT_M and IMT_S (combined code only)
c     JMT_L    = maximum of JMT_M and JMT_S (combined code only)
c     KM       = number of grid points in the vertical direction
c     NT       = number of tracers
c     MXSLAVE  = maximum number of processors to be used
c
c archive parameters:
c     NARCHV   = number of variable types to be requested from the
c                slaves for archiving.
c     NUM2D    = The number of the NARCHV variables which are 2-d fields
c     JSUB_M   = number of j-rows which can be temporarily held within
c                the master whilst archiving 3-d fields.        (master)
c     NZONE_M  = number of times the master will have to request each
c                3-d variable.
c
c snapshot parameters:
c     NSNAPS   = total number of horizontal snapshots
c
c  pvm parameters
c     LENBUFF  = size of buffer used with pvm
c
c  surface flux parameters:
c    Master:
c    NWX_M     = size of wind stress array in longitude direction
c    NWY_M     = size of wind stress array in latitude direction
c    NSX_M     = size of surface T/S array in longitude direction
c    NSY_M     = size of surface T/S array in latitude direction
c    Slave:
c    NWX_S     = size of wind stress array in longitude direction
c    NWY_S     = size of wind stress array in latitude direction
c    NSX_S     = size of surface T/S array in longitude direction
c    NSY_S     = size of surface T/S array in latitude direction
c-----------------------------------------------------------------------
c
#endif
      parameter  (IMT_M=90,  JMT_M=36,  KM=15,  NT=2,     NTP2=NT+2,
     &            IMTM1_M=IMT_M-1,      JMTM1_M=JMT_M-1,
     &            KMP1=KM+1,            KMM1=KM-1,
     &            NARCHV=7,             NUM2D=3,
     &            JSUB_M=10,            NZONE_M=(JMT_M-1)/JSUB_M + 1,
     &            NSNAPS=11                                         )
      parameter  (IMT_S=92,  JMT_S=36,              
     &            IMTM1_S=IMT_S-1,      JMTM1_S=JMT_S-1)
      parameter  (MXSLAVE=255, IBOUND=4*(IMT_S+JMT_S))
      parameter (NWX_M = 80, NWY_M = 41, NSX_M = 80, NSY_M = 40) 
      parameter (NWX_S = 80, NWY_S = 41, NSX_S = 80, NSY_S = 40) 
#ifdef pvm_buffer
      parameter (LENBUF1=3*IMT_S*JMT_S+10           )
      parameter (LENBUF2=IMT_S*(JSUB_M+2)*KM+10     )
      parameter (LENBUF3=LENBUF1/LENBUF2            )
      parameter (LENBUF4=LENBUF3/(LENBUF3-0.1)      )
      parameter (LENBUF =LENBUF2+LENBUF3*(LENBUF1-LENBUF2) )
#endif
#ifdef hcomments
c
c     add parameter constants
c
#endif
#include "pconst.h"
