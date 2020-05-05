      subroutine tmngr (dt)
c
c @(#) SCCS module: tmngr.F  version: 1.1
c     Creation date: 03/16/95
c
c==================================================================
c
c                      t-i-m-e      m-a-n-a-g-e-r
c
c     inputs:
c
c     itt    = current model time step (from ctmngr.h).
c     dtts   = length of time step in seconds.
c     dgnstc = number of days between printing diagnostics
c     nmix   = number of timesteps between mixint timesteps
c     tsi    = number of days between printing timestep information
c     days   = day at which to end run.
c
c     outputs:
c
c     itt    = itt (on input)  + 1
c     totday = total number of days elapsed 
c              (0.5 = noon on first day).
c     years is number of years elapsed (365.5 days = 1 year)
c     diagts = true if time to print diagnostics
c     mixts  = true if a mixing timestep
c     prntsi = true if time to print timestep information
c     last   = true if last timestep of run
c
c==================================================================
c
#include "def_master.h"
#include "pconst.h"
#include "ctmngr.h"
#include "switch.h"
#include "timelv.h"
      logical*4 tmnset
c
c------------------------------------------------------------------
c     increment time step counters
c     and create a time stamp
c------------------------------------------------------------------
c
      itt = itt + 1
      totsec = totsec + dt
      totday = totsec*secday
      years  = totday/365.25
      write(stamp,11)totday,years
  11  format('totday =',f8.2,',years =',f8.3)
c
c------------------------------------------------------------------
c
c     set time dependent logical switches.
c     1. set last if it is the last timestep
c     2. set diagts if time to do diagnostics 
c     3. set mixts if time mixing to be done
c     4. set prntsi if time to print timestep info
c     5. set snapts if time to archive model
c
c------------------------------------------------------------------
c
      last   = (totday .ge. days) 
      mixts  = (mod(itt,nmix) .eq. 1) 
      prntsi = tmnset (tsi, totday, dt)
      diagts = tmnset (dgnstc, totday, dt)
      snapts = tmnset (snapd, totday, dt)
      archts = tmnset (archd, totday, dt) .or. (restrt .and. last)
c
c  set flags for next timestep (itt + 1)
c
      totdayp= (totsec+dt)*secday
      lastn  = (totdayp .ge. days) 
      prntsn = tmnset (tsi, totdayp, dt)
      snaptn = tmnset (snapd, totdayp, dt) 
      archtn = tmnset (archd, totdayp, dt) .or. (restrt .and. lastn)
c
      return
      end
