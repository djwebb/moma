c====================== include file "moddata.h" =======================
#ifdef hcomments
c
c @(#) SCCS module: moddata.h  version: 1.7
c     Creation date: 06/30/95
c
c-----------------------------------------------------------------------
c This module contains various data statements for inline inclusion
c in master.F and slave.F. Data used only by the master process are
c bracketed by the "#ifdef Master" - "#endif" preprocesssor option.
c-----------------------------------------------------------------------
c
#endif
#ifdef Master
c
c Position descriptions and indexing
c
      data mesind/ 0, 8, 9, 9,10,11, 0, 0, 0, 0, 1,
     &                2, 0, 0, 3, 0, 0, 0, 0, 0, 4,
     &                5, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     &                0, 6, 0, 0, 7, 0, 0, 0, 0,12/
c
c position descriptions. If any descriptions are added to this
c list then the length of the messtr array will need to be
c changed in mesdta.h
c
      data messtr/ 'no description available'      ,
     &             'waiting before main loop'      ,
     &             'at top of main loop'           ,
     &             'entering clinic-tracer loop'   ,
     &             'entering free-surface loop'    ,
     &             'entering addv loop'            ,
     &             'preparing timestep info'       ,
     &             'processing master requests'    ,
     &             'aborting slave processes'      ,
     &             'receiving restart data'          ,
     &             'waiting for restart data'      ,
     &             'receiving free surface data'   ,
     &             'closing down'                  /
c
c snapshot variable identifiers and levels (see snaps.h)
c
      data ispvar/ 1, 2, 3, 4, 5, 6, 7,  4,  5,  6,  7/
      data isplev/ 1, 1, 1, 1, 1, 1, 1,  8,  8,  8,  8/
c
#endif
c
c variable descriptions, units and grid indicator (see varinfo.h)
c
      data (varname(n),varunit(n),vartgrd(n),n=1,NARCHV)/
     &  'free surface height      ', 'cm          ', .true.  ,
     &  'barotropic u-velocity    ', 'cm/s        ', .false. ,
     &  'barotropic v-velocity    ', 'cm/s        ', .false. ,
     &  'potential temperature    ', 'C           ', .true.  ,
     &  'salinity                 ', '.001(NSU-35)', .true.  ,
     &  'u-velocity               ', 'cm/s        ', .false. ,
     &  'v-velocity               ', 'cm/s        ', .false. /

