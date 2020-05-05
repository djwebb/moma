c====================== include file "levind.h" ========================
c
c     vertical level indicators which define model geometry & 
c     bottom topography:
c
c     kmt = number of vertical boxes over "t" points
c     kmu = number of vertical boxes over "u,v" points
c
      common /levind/ kmt(imt,jmt), kmu(imt,jmt)
