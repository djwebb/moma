c====================== include file "levind.h" ========================
#ifdef hcomments
c
c @(#) SCCS module: levind.h  version: 1.10
c     Creation date: 08/09/96
c
c-----------------------------------------------------------------------
c     Vertical level indicators which define model geometry & 
c     bottom topography:
c
c     kmt = number of vertical boxes over "t" points
c     kmu = number of vertical boxes over "u,v" points
c     kpn = processor responsible for this grid point
c     kmd = bottom level for vertical velocity calculations
c-----------------------------------------------------------------------
c
#endif
#ifdef Masterslave
      parameter (MS_DIF3 = 3*IMT_M*JMT_M - 4*IMT_S*JMT_S)
      parameter (ADIF3 = 0.5*MS_DIF3/((MS_DIF3+0.01)*(MS_DIF3+0.01)) )
      parameter (LOG_M3 = 1-ADIF3, LOG_S3 = 1+ADIF3)
      parameter (LBUF_S3 = 1 + MS_DIF3*LOG_S3)
      parameter (LBUF_M3 = 1 - MS_DIF3*LOG_M3)
#endif
#ifdef Master
      common /levind/ kmt(IMT_M,JMT_M), kmu(IMT_M,JMT_M)
     &              , kpn(IMT_M,JMT_M)
#ifdef Masterslave
     &              , buffll(LBUF_M3)
#endif
#else
      common /levind/ kmt(IMT_S,JMT_S), kmu(IMT_S,JMT_S)
     &              , kpn(IMT_S,JMT_S), kmd(IMT_S,JMT_S)
#ifdef Masterslave
     &              , buffll(LBUF_S3)
#endif
#endif
