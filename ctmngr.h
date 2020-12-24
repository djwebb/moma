#ifdef hcomments
!
!======================================================================
!
!    File:  ctmngr.h
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
!     time manager clock parameters
!
!     itt     = current time step number
!     itbt    = current barotropic timestep
!     itbtp   = current pass of free surface model timestepping scheme
!     totsec  = accumulated time in seconds from jan 1, year 0
!     totsec  = accumulated time in days
!     years   = accumulated time in years
!     stamp   = date & time stamp corresponding to itt
!
#endif
      character*32 stamp
      common /ctmngr/ totsec, totday, years
      common /ctmngi/ itt, itbt, itbtp
      common /ctmngc/ stamp

