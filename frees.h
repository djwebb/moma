c====================== include file "frees.h2" =========================
#ifdef hcomments
c
c     @(#) SCCS module: frees.h, version 1.13
c
c     variables for free surface model calculation of external mode
c
c     h0   =  free surface model sea surface height
c     u0   =  free surface model velocity east
c     v0   =  free surface model velocity north
c     h    = depth over "u,v" points
c     hr   = reciprocal depth over "u,v" points
c     zu   = vertically averaged zonal forcing 
c     zv   = vertically averaged meridional forcing 
c
#endif
      common /fields/
     & h0(imt,jmt,3), u0(imt,jmt,3), v0(imt,jmt,3)
     &,h(imt,jmt),    hr(imt,jmt),    zu(imt,jmt),  zv(imt,jmt)
