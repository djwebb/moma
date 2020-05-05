c====================== include file "ctmngr.h" ========================
#ifdef hcomments
c
c @(#) SCCS module: ctmngr.h  version: 1.9
c     Creation date: 12/18/95
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
c
c     surface flux variables
c
c     mlast   = month number of previous (monthly) flux dataset
c     mnext   = month number of next flux dataset
c     mnextp  = month number of following flux dataset
c     eday    = central day for each of the (monthly) flux data sets
c
c-----------------------------------------------------------------------
c
#endif
      character*32 stamp
      real eday(0:13)
      common /ctmngr/ totsec, totday, years, twa, twb, eday
      common /ctmngi/ itt, itbt, itbtp, mlast, mnext, mnextp
      common /ctmngc/ stamp

