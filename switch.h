#ifdef hcomments
!
!======================================================================
!
!    File:  switch.h
!    Copyright 1992-1997, 2020  David J. Webb
!
!    This file is part of the Moma ocean model code.
!    Released under licence GPL-3.0-or-later
!
!    Moma is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    Moma is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with Moma.  If not, see <https://www.gnu.org/licenses/>.
!
!======================================================================
!
!     all time dependent decisions are made in time manager
!     "tmngr.F" and communicated elsewhere in the model via
!     logical switches .
!
!     inputs: (set through namelist)
!
!     days    = number of days to integrate
!     init    = true if this run is starts from scratch
!               false if restarting from an archived data set.
!     nmix    = number of time steps between time step mixing
!               to damp leap frog time splitting
!     eb      = true implies euler backward mixing, false
!               implies a forward timestep.
!     acor    = (>0, 0) = (implicit, explicit) treatment of coriolis
!               term
!     tsi     = number of days between printing of time step info
!     dgnstc  = number of days between diagnostic calculations:
!     snapd   = number of days between saving an archive dataset
!     archd   = number of days between saving an archive dataset
!     restrt  = true if a restart data set is to be written
!               at the end of this run
!
!     outputs: (set in tmngr.F)
!
!     the following logical switches are set within "tmngr" every
!     time step based on the above requirements.
!
!     last    = true if this is the last timestep of a run
!     mixts   = true if this is a mixing time step
!     prntsi  = true if time step info to be printed
!     diagts  = true if diagnostics are to be printed
!     snapts  = true if this is an archive timestep
!     archts  = true if this is an archive timestep
!
!     the following switches are set within the main program
!
!     first   = true if this is the first timestep of a run
!     eots    = end of a time step. always true except for first
!               pass of an euler backward time step
!     mxpas2  = second pass of mixing timestep
# ifdef de_checkbd
!     lchkbd  = logical variable used within the barotropic
!               time loop to decide whether or not to apply
!               the del-cross-del-plus filter
# endif
!
#endif
      LOGICAL_N init, restrt, eb, first, last, mixts, eots
     &,       mxpas2, prntsi, diagts, snapts, archts
#ifdef de_checkbd
     &,       lchkbd
#endif
      common /switchl/ init, eb, restrt, first, last, mixts, eots
     &,       mxpas2, prntsi, diagts, snapts, archts
#ifdef de_checkbd
     &,       lchkbd
#endif
      common /switchr/days,acor,tsi,dgnstc,snapd,archd
      common /switchi/nmix,idebug

