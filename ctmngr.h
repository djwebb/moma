c====================== include file "ctmngr.h" ========================
#ifdef hcomments
c
c @(#) SCCS module: ctmngr.h  version: 1.6
c     Creation date: 04/13/95
c
c-----------------------------------------------------------------------
c     Time manager clock parameters 
c
c     itt     = current baroclinic timestep
c     itbt    = current barotropic timestep
c     itbtp   = current pass of free surface model timestepping scheme
c     totsec  = accumulated time in seconds from jan 1, year 0
c     totsec  = accumulated time in days
c     years   = accumulated time in years
c     stamp   = date & time stamp corresponding to itt
c-----------------------------------------------------------------------
c
#endif
      character*32 stamp
      common /ctmngr/ totsec, totday, years
      common /ctmngi/ itt, itbt, itbtp
      common /ctmngc/ stamp

