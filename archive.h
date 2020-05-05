c-------------------------archive.h------------------------------
#ifdef hcomments
c
c @(#) SCCS module: archive.h  version: 1.8
c     Creation date: 08/23/95
c
c-----------------------------------------------------------------------
c Temporary storage for restart or archive data. Each slave
c copies the model state into its own copy of this storage and
c asynchronously sends it to the master in response to the
c master's requests. Progress is blocked should the slave reach
c the next archive step before the master has collected all the
c buffered data.
c
c The master has different storage requirements which are allowed
c for here by the use of preprocessor directives. At any one time
c the master will be processing either a single two-dimensional
c field or a sub-volume of one, three-dimensional field. Values of
c the following variables are already known to the master and are
c not requested from the slaves:
c
c   itt 			- integer variable
c   totsec, totday, years	- real variables
c   kmt				- IMT_M*JMT_M integer array
c
c Each of the following fields are requested in turn using the
c following identifiers:
c
c   identifier (variable type)      variable
c------------------------------+-----------------
c		1		free surface height (h0)
c		2		zonal velocity (u0)
c		3		meridional velocity (v0)
c		4		tracer type 1 (potential temperature)
c		5		tracer type 2 (salinity)
c		6		baroclinic zonal velocity (u)
c		7		baroclinic merid'al velocity (v)
c-----------------------------------------------------------------------
c
#endif
#ifdef Masterslave
      parameter (LEN_M4 = KM*IMT_M*JSUB_M)
      parameter (LEN_S4 = IMT_S*JMT_S*NUM2D+KM*IMT_S*JMT_S*NTP2)
      parameter (MS_DIF4 = LEN_M4 - LEN_S4)
      parameter (ADIF4 = 0.5*MS_DIF4/((MS_DIF4+0.01)*(MS_DIF4+0.01)) )
      parameter (LOG_M4 = 1-ADIF4, LOG_S4 = 1+ADIF4)
      parameter (LBUF_S4 = 1 + MS_DIF4*LOG_S4)
      parameter (LBUF_M4 = 1 - MS_DIF4*LOG_M4)
#endif
#ifdef Master
      common /resti/ ittsav
      common /restr/ totssav, totdsav, yrssav,
     &               rest3d(KM,IMT_M,JSUB_M)
#ifdef Masterslave
     &              ,buffaa(LBUF_M4)
#endif
      real rest2d(IMT_M,JMT_M)
      equivalence (rest3d(1,1,1),rest2d(1,1))
#else
      common /resti/ ittsav
      common /restr/ totssav, totdsav, yrssav,
     &               hsave(IMT_S,JMT_S,NUM2D)  
     &             , tsave(KM,IMT_S,JMT_S,NTP2)
#ifdef Masterslave
     &             ,buffaa(LBUF_S4)
#endif
#endif
