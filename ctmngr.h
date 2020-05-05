c====================== include file "ctmngr.h" ========================
#ifdef hcomments
c
c     @(#) SCCS module: ctmngr.h, version 1.7
c
c     time manager clock parameters 
c
c     itt     = current time step number
c     totsec  = accumulated time in seconds from jan 1, year 0
c     totsec  = accumulated time in days
c     years   = accumulated time in years
c     stamp   = date & time stamp corresponding to itt
c
#endif
      character*32 stamp
      common /ctmngr/ itt, totsec, totday, years
      common /ctmngc/ stamp

