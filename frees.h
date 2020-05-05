c====================== include file "frees.h" =========================
#ifdef hcomments
c
c @(#) SCCS module: frees.h  version: 1.8
c     Creation date: 08/23/95
c
c-----------------------------------------------------------------------
c     Variables for free surface model calculation of external mode
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
c-----------------------------------------------------------------------
c
#endif
      common /fields/
     & h0(IMT_S,JMT_S,3), u0(IMT_S,JMT_S,3), v0(IMT_S,JMT_S,3)
     &,h(IMT_S,JMT_S),    hr(IMT_S,JMT_S),   zu(IMT_S,JMT_S)
     &,zv(IMT_S,JMT_S)
#ifndef free_eb
     &,freeav(3,IMT_S,JMT_S)
#endif
