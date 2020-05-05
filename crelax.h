c====================== include file "crelax.h" ========================
#ifdef hcomments
c
c     @(#) SCCS module: crelax.h, version 1.4
c
c     variables used in relaxation. see "stream.F" 
c
c     inputs:  (set through namelist)
c       mxscan  = max number of allowable scans
c       crit    = criterian used to stop relaxation
c       sor     = sucessive over-relaxation constant
c
c     output:  (from stream)
c       mscan   = actual number of scans taken
c
#endif
      common /crelax/  mxscan, crit, sor, mscan

