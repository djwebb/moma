c====================== include file "levind.h" ========================
#ifdef hcomments
c
c     @(#) SCCS module: levind.h, version 1.8
c
c     vertical level indicators which define model geometry & 
c     bottom topography:
c
c     kmt = number of vertical boxes over "t" points
c     kmu = number of vertical boxes over "u,v" points
c
#endif
      common /levind/ kmt(imt,jmt), kmu(imt,jmt)

