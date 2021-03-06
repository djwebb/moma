c====================== include file "scalar.h" ========================
#ifdef hcomments
c
c @(#) SCCS module: scalar.h  version: 1.1
c     Creation date: 10/13/97
c
c-----------------------------------------------------------------------
c     Various scalar quantities:
c
c     dtts   = time step for density & tracers (in seconds)
c     dtuv   = time step for baroclinic velocity (in seconds)
c     dtbt   = time step for barotropic velocity (in seconds)
c              used by free surface model
c     c2dtts = 2*dtts
c     c2dtuv = 2*dtuv
c     area   = surface area of ocean (cm**2)
c     volume = volume of ocean (cm**3)
c     omega  = earth's rotation rate (radians/sec)
c     radius = earth's radius (cm)
c     grav   = earth's gravitational acceleration (cm/sec**2)
c     cdbot  = bottom drag coefficient
c
c     ncon   = number of  passes through convective code in tracer
c     ntbt   = number of free surface timesteps per baroclinic velocity
c              time steps.
c     ntbt2  = 2*ntbt
# ifdef de_checkbd
c   dchkbd   = weighting factor for del-plus-del-cross filter.
# endif
c-----------------------------------------------------------------------
c
#endif
      common /scalar/ dtts, dtuv, dtbt, c2dtts, c2dtuv, c2dtbt 
     &,               area, volume, omega, radius, grav, cdbot
#ifdef de_checkbd
     &,               dchkbd
#endif
      common /scalri_s/ it, jt, itm1, jtm1, iswm1, jswm1
      common /scalri/ ncon, ntbt, ntbt2
#ifdef hcomments
c
c     non dimensional quantities:
c
c     radian = degrees per radian
c     pi     = something good to eat
c
#endif
      common /ndcon/ radian, pi

