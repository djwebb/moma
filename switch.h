c===================== include file "switch.h" =========================
#ifdef hcomments
c
c @(#) SCCS module: switch.h  version: 1.5
c     Creation date: 03/24/95
c
c-----------------------------------------------------------------------
c
c     All time dependent decisions are made in time manager
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
c     dgnstc  = number of days between extra idagnostic printing
c     snapd   = number of days between saving a snapshot dataset
c     archd   = number of days between saving an archive dataset
c     restrt  = true if a restart data set is to be written
c               at the end of this run
c     idebug  =  controls printing of debug information
c
c     outputs: (set in tmngr.F)
c
c     the following logical switches are set within "tmngr" every
c     time step based on the above requirements.
c
c     last    = true if this is the last timestep of a run
c     mixts   = true if this is a mixing time step
c     prntsi  = master: true if time step info to be printed 
c               slave: true if time step info needed but as yet
c               not sent to master
c     diagts  = true if this is an extra diagnostic timestep
c     snapts  = master :true if this is a snapshot timestep
c               slave: true if this is a snapshot timestep and
c               as yet not all data sent to master.
c     archts  = master: true if this is an archive timestep 
c               slaves: true if this is a snapshot timestep and
c               as yet not all data sent to master.
c
c     lastn   = true if the next timestep is the last timestep of a run
c     prntsn  = true if time step info to be printed next timestep
c     snaptn  = true if there is a snapshot next timestep
c     archtn  = true if there is an archive next timestep
c
c     the following switches are set within the main program
c
c     first   = true if this is the first timestep of a run
c     eots    = end of a time step. always true except for first
c               pass of an euler backward time step
c     mxpas2  = second pass of mixing timestep
c     frpas1  = first pass of a free surface model timestep
# ifdef de_checkbd
c     lchkbd  = logical variable used within the barotropic
c               time loop to decide whether or not to apply
c               the del-cross-del-plus filter
# endif
c-----------------------------------------------------------------------
c
#endif
      LOGICAL_N init, restrt, eb, first, last, mixts, eots
     &,         mxpas2, frpas1, prntsi, diagts, snapts, archts
     &,         lastn , prntsn, snaptn, archtn 
#ifdef de_checkbd
     &,       lchkbd
#endif
      common /switchl/init, eb, restrt, first, last, mixts, eots
     &,       mxpas2, frpas1, prntsi, diagts, snapts, archts
     &,       lastn , prntsn, snaptn, archtn 
#ifdef de_checkbd
     &,       lchkbd
#endif
      common /switchr/days,acor,tsi,dgnstc,snapd,archd
      common /switchi/nmix,idebug

