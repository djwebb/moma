c====================== include file "frees.h" =========================
#ifdef hcomments
c
c     @(#) SCCS module: frees.h, version 1.17
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
# ifndef free_eb
c     freeav = time average of free surface fields
# endif
c
#endif
      common /fields/
     & h0(imt,jmt,3), u0(imt,jmt,3), v0(imt,jmt,3)
     &,h(imt,jmt),    hr(imt,jmt),    zu(imt,jmt),  zv(imt,jmt)
#ifndef free_eb
     &,freeav(3,imt,jmt)
#endif
