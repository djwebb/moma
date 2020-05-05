c====================== include file "switch.h" =========================
#ifdef hcomments
c
c     @(#) SCCS module: switch.h, version 1.12
c
c     all time dependent decisions are made in time manager
c     "tmngr.F" and communicated elsewhere in the model via 
c     logical switches .
c
c     inputs: (set through namelist)
c
c     days    = number of days to integrate
c     init    = true if this run is starts from scratch
c               false if restarting from an archived data set.
c     nmix    = number of time steps between time step mixing
c               to damp leap frog time splitting
c     eb      = true implies euler backward mixing, false
c               implies a forward timestep.
c     acor    = (>0, 0) = (implicit, explicit) treatment of coriolis
c               term
c     tsi     = number of days between printing of time step info
c     dgnstc  = number of days between diagnostic calculations:
c     snaps   = number of days between saving an archive dataset
c     restrt  = true if a restart data set is to be written
c               at the end of this run
c
c     outputs: (set in tmngr.F)
c
c     the following logical switches are set within "tmngr" every
c     time step based on the above requirements.
c
c     last    = true if this is the last timestep of a run
c     mixts   = true if this is a mixing time step
c     prntsi  = true if time step info to be printed
c     diagts  = true if diagnostics are to be printed
c     snapts  = true if this is an archive timestep
c
c     the following switches are set within the main program
c
c     first   = true if this is the first timestep of a run
c     eots    = end of a time step. always true except for first
c               pass of an euler backward time step
c     mxpas2  = second pass of mixing timestep
c     frpas1  = first pass of a free surface model timestep
c
#endif
      logical init, restrt, eb, first, last, mixts, eots,
     &        mxpas2, frpas1, prntsi, diagts, snapts 
      common /switch/ days, init, nmix, eb, acor, tsi, dgnstc,  
     &        snaps, restrt, first, last, mixts, eots, mxpas2,  
     &        frpas1, prntsi, diagts, snapts 

