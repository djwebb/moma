c====================== include file "thick.h" ========================
#ifdef hcomments
c
c     @(#) SCCS module: thick.h, version 1.14
c
c     specifications for box thickness (cm) in the vertical
c
#endif
      data dzt /   30.00e2,   46.15e2,   68.93e2,   99.93e2,  140.63e2
     &      ,  192.11e2,  254.76e2,  327.95e2,  409.81e2,  497.11e2
     &      ,  585.36e2,  669.09e2,  742.41e2,  799.65e2,  836.10e2
     &   /
#ifdef hcomments
c
c     dzt = model level thicknesses computed as
c     dzt(km+1-k)=exp(-(k/sigma)**2) + offset) for k=1,km
c      km =  15
c      sigma=  8.204686
c      offset= 0.0000000
c      1st level thickness =  30.00000
c      # of levels with this thickness =  1
c      depth at bottom of last level =  5700.000
#endif

