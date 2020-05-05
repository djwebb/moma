c====================== include file "snaps.h" =========================
#ifdef hcomments
c
c @(#) SCCS module: snaps.h  version: 1.8
c     Creation date: 08/23/95
c
c-----------------------------------------------------------------------
c Temporary storage for snapshot data. Each slave
c copies the model state into its own copy of this storage and
c asynchronously sends it to the master in response to the
c master's requests. Progress is blocked should the slave reach
c the next snapshot step before the master has collected all the
c buffered data.
c
c The master has different storage requirements which are allowed
c for here by the use of preprocessor directives. At any one time
c the Master will be processing a single two-dimensional field.
c
c ispvar   = integer array containing the variable identifiers
c            for each of the required horizontal snapshots.
c isplev   = integer array containing the vertical levels for
c            each of the corresponding snapshots
c
c The following identifiers are used:
c
c   identifier (variable type)      variable
c------------------------------+-----------------
c		1		free surface height (h0)
c		2		zonal velocity (u0)
c		3		meridional velocity (v0)
c		4		tracer type 1 (temperature)
c		5		tracer type 2 (salinity)
c		6		baroclinic zonal velocity (u)
c		7		baroclinic merid'al velocity (v)
c
c-----------------------------------------------------------------------
c
#endif
      common /snapi/ ispvar(NSNAPS),isplev(NSNAPS)
#ifdef Masterslave
      parameter (MS_DIF5 = IMT_M*JMT_M - IMT_S*JMT_S*NSNAPS)
      parameter (ADIF5 = 0.5*MS_DIF5/((MS_DIF5+0.01)*(MS_DIF5+0.01)) )
      parameter (LOG_M5 = 1-ADIF5, LOG_S5 = 1+ADIF5)
      parameter (LBUF_S5 = 1 + MS_DIF5*LOG_S5)
      parameter (LBUF_M5 = 1 - MS_DIF5*LOG_M5)
#endif
#ifdef Master
      common /snapr/ snap2d(IMT_M,JMT_M)
#ifdef Masterslave
     &         ,buffss(LBUF_M5)
#endif
     &         ,vmask(4)
#else
      common /snapr/ snap2d(IMT_S,JMT_S,NSNAPS)
#ifdef Masterslave
     &         ,buffss(LBUF_S5)
#endif
#endif
