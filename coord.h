c====================== include file "coord.h" =========================
#ifdef hcomments
c
c @(#) SCCS module: coord.h  version: 1.1
c     Creation date: 03/16/95
c
c-----------------------------------------------------------------------
c     Grid definition:
c
c     the model uses a staggard arakawa "b" grid. horizontally, 
c     tracer quantities are defined at the centers of "t" grid 
c     boxes and velocities are defined at the centers of "u,v" grid 
c     boxes.  the centers of "u,v" grid boxes are located at the 
c     northeast corners of "t" grid boxes. the 1st "t" grid box is 
c     located in the southwest corner of the "t" grid. this grid 
c     system is replicated and stacked vertically one on top of 
c     another from the surface of the ocean downward. vertically,
c     tracers & velocities are defined at the centers of their 
c     corresponding boxes and are at the same depths.
c
c     primary grid information:
c
c     stlon = starting longitude of the "t" grid in degrees. It 
c             defines the west face of the 1st "t" grid box.
c     stlat = starting latitude of the "t" grid in degrees. It 
c             defines the south face of the 1st "t" grid box
c     dxdeg = east-wast width of grid in degrees
c     dydeg = north-south width of grid in degrees
c
c     xt(i) = longitude of the ith "t" point in degrees.   i=1..imt
c     xu(i) = longitude of the ith "u,v" point in degrees. i=1..imt
c     yt(j) = latitude of the jth "t" point in degrees.    j=1..jmt
c     yu(j) = latitude of the jth "u,v" point in degrees.  j=1..jmt
c     zt(k) = distance (cm) from surface down to center of level k 
c             (depth of "t" and "u,v" grid points)
c     zw(k) = distance (cm) from surface down to bottom of level k 
c
c-----------------------------------------------------------------------
c  
#endif
#ifdef Masterslave
      parameter (MS_DIF2 = 2*IMT_M + 2*JMT_M - 2*IMT_S - 2*JMT_S)
      parameter (ADIF2 = 0.5*MS_DIF2/((MS_DIF2+0.01)*(MS_DIF2+0.01)) )
      parameter (LOG_M2 = 1-ADIF2, LOG_S2 = 1+ADIF2)
      parameter (LBUF_S2 = 1 + MS_DIF2*LOG_S2)
      parameter (LBUF_M2 = 1 - MS_DIF2*LOG_M2)
#endif
#ifdef Master
      common /coord/ stlon, stlat, dxdeg, dydeg, 
     &               xt(IMT_M), yt(JMT_M), xu(IMT_M), yu(JMT_M), 
     &               zw(KM), zt(KMP1)
#ifdef Masterslave
     &              ,buffcc(LBUF_M2)
#endif
#else
      common /coord/ stlon, stlat, dxdeg, dydeg, 
     &               xt(IMT_S), yt(JMT_S), xu(IMT_S), yu(JMT_S), 
     &               zw(KM), zt(KMP1)
#ifdef Masterslave
     &              ,buffcc(LBUF_S2)
#endif
#endif



