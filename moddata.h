c====================== include file "moddata.h" =======================
#ifdef hcomments
c
c @(#) SCCS module: moddata.h  version: 1.1
c     Creation date: 10/13/97
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
      data mesind/ 0, 8, 9, 9,10,11, 0,13, 0, 0, 1,
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
     &             'receiving restart data'        ,
     &             'waiting for restart data'      ,
     &             'receiving free surface data'   ,
     &             'closing down'                  ,
     &             'receiving meteorological data' /
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
c
c     time in days at which monthly met data is centred (see ctmngr.h):
c
      data eday / -15.5, 15.50, 45.125, 74.75,105.25,133.75,166.25,
     &            196.75,227.75,258.25,288.75,319.25,349.75, 380.75/
#ifdef Master
c
c  flux data
c
      data sstobs / -1.75, -1.75, -1.50, -1.50, -1.28,
     &              -0.55,  0.90,  2.92,  5.45,  8.62,
     &              12.27, 15.49, 18.30, 20.67, 22.64,
     &              24.14, 25.27, 26.37, 26.52, 26.16,
     &              26.85, 27.27, 26.82, 26.42, 25.53,
     &              24.03, 22.07, 19.73, 17.02, 12.77,
     &               8.93,  7.25,  6.22,  4.67,  4.57,
     &               3.03, -0.01, -1.05, -1.75, -1.75/
c
      data salobs / 34.30, 34.30, 34.30, 34.13, 33.98,
     &              33.97, 33.97, 33.98, 34.03, 34.24,
     &              34.61, 35.02, 35.37, 35.61, 35.72,
     &              35.68, 35.51, 35.22, 35.05, 35.12,
     &              34.80, 34.56, 34.71, 34.90, 35.27,
     &              35.67, 35.56, 35.49, 35.23, 34.28,
     &              33.57, 33.57, 33.60, 33.80, 34.04,
     &              34.05, 32.65, 32.30, 32.10, 32.00/
c
      data wsxobs /  0.00,
     &               0.00,  0.00, -0.02,  0.15,  0.31,
     &               0.50,  0.82,  1.08,  1.23,  1.16,
     &               0.84,  0.41,  0.02, -0.35, -0.55,
     &              -0.67, -0.64, -0.46, -0.29, -0.19,
     &              -0.16, -0.33, -0.52, -0.59, -0.55,
     &              -0.32,  0.09,  0.42,  0.56,  0.76,
     &               0.81,  0.65,  0.29,  0.06, -0.10,
     &              -0.05, -0.03,  0.05,  0.10,  0.01/
c
      data wsyobs /  .000,
     &               .000,  .009,  .032,  .005, -.023,
     &              -.075, -.155, -.202, -.230, -.179,
     &              -.049,  .093,  .214,  .294,  .344,
     &               .383,  .364,  .269,  .189,  .178,
     &               .125, -.122, -.213, -.251, -.259,
     &              -.202, -.189, -.179, -.183, -.009,
     &               .023,  .053, -.048, -.185, -.225,
     &              -.097, -.050, -.023, -.006,  .000/
#endif
