c====================== include file "grdvar.h" ========================
#ifdef hcomments
c
c @(#) SCCS module: grdvar.h  version: 1.7
c     Creation date: 06/30/95
c
c-----------------------------------------------------------------------
c  This version assumes all dx and dy boxes are of equal size when
c  measured in degrees.  Variables dxt, dxu etc replaced by dx, 
c  dyt,dyu etc by dy and dzt by dz.
c
c     dx     = longitudinal width of "t" grid box (in cm)
c     dxr    = reciprocal of "dxt"
c     dx2r   = reciprocal of "2*dxt"
c     dx4r   = reciprocal of "4*dxt"
c
c     dy     = latitudinal height of "t" grid box (in cm)
c     dyr    = reciprocal of "dyt"
c     dy2r   = reciprocal of "2*dyt"
c     dy4r   = reciprocal of "4*dyt"
c
c     csu     = cosine of "u,v" grid point latitude
c     csur    = reciprocal of "csu"
c     cst     = cosine of "t" grid point latitude
c     cstr    = reciprocal of "cst"
c     phi     = latitude of "u,v" grid point in radians
c     phit    = latitude of "t" grid point in radians
c     sine    = sine of "u,v" grid point latitude
c     tng     = tan of "u,v" grid point latitude
c     fcor    = 2*omega*sine(j)
c
c     dz(k)  = level thickness of "t" and "u,v" grid boxes (in cm)
c               (vertical separation between "w" velocity points)
c     c2dt(k)= "2*dzt"
c     dzr(k) = reciprocal of dzt
c     dz2r(k)= reciprocal of "2*dzt"
c     dzw(k)  = vertical separation between centers of levels k & k+1
c     dzwr(k) = reciprocal of dzw
c     dzw2r(k)= reciprocal of "2*dzw"
c     dzur(k)= upper diffusion grid factor = 1.0/(dzw(k-1)*dz(k))
c     dzlr(k)= lower diffusion grid factor = 1.0/(dzw(k)*dz(k))
c     dzwur(k)= upper diffusion grid factor = 1.0/(dz(k)*dzw(k))
c     dzwlr(k)= lower diffusion grid factor = 1.0/(dz(k+1)*dzw(k))
c     tanra(j)= tan(j)/radius on the "u,v" grid
c-----------------------------------------------------------------------
c
#endif
#ifdef Masterslave
      parameter (MS_DIF1 = 10*JMT_M - 10*JMT_S)
      parameter (ADIF1 = 0.5*MS_DIF1/((MS_DIF1+0.01)*(MS_DIF1+0.01)) )
      parameter (LOG_M1 = 1-ADIF1, LOG_S1 = 1+ADIF1)
      parameter (LBUF_S1 = 1 + MS_DIF1*LOG_S1)
      parameter (LBUF_M1 = 1 - MS_DIF1*LOG_M1)
#endif
#ifdef Master
      common /grdvar/ dx, dxr, dx2r, dx4r, dy, dyr, dy2r, dy4r,
     &                csu(JMT_M),  csur(JMT_M), cst(JMT_M), cstr(JMT_M),
     &                phi(JMT_M),  phit(JMT_M), sine(JMT_M), tng(JMT_M),
     &                fcor(JMT_M), tanra(JMT_M),
     &                dz(KM),    dzr(KM),    dz2r(KM),   c2dz(KM), 
     &                dzw(0:KM), dzwr(0:KM), dzw2r(0:KM)
#ifdef Masterslave
     &               ,buffgg(LBUF_M1)
#endif
#else
      common /grdvar/ dx, dxr, dx2r, dx4r, dy, dyr, dy2r, dy4r,
     &                csu(JMT_S),  csur(JMT_S), cst(JMT_S), cstr(JMT_S),
     &                phi(JMT_S),  phit(JMT_S), sine(JMT_S), tng(JMT_S),
     &                fcor(JMT_S), tanra(JMT_S),
     &                dz(KM),    dzr(KM),    dz2r(KM),   c2dz(KM), 
     &                dzw(0:KM), dzwr(0:KM), dzw2r(0:KM)
#ifdef Masterslave
     &               ,buffgg(LBUF_S1)
#endif
#endif


