#ifdef hcomments
!
!======================================================================
!
!    File:  scalar.h
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
!     various scalar quantities:
!
!     dtts   = time step for density & tracers (in seconds)
!     dtuv   = time step for baroclinic velocity (in seconds)
!     dtbt   = time step for barotropic velocity (in seconds)
!              used by free surface model
!     c2dtts = 2*dtts
!     c2dtuv = 2*dtuv
!     c2dtbt = 2*dtbt
!     area   = surface area of ocean (cm**2)
!     volume = volume of ocean (cm**3)
!     omega  = earth's rotation rate (radians/sec)
!     radius = earth's radius (cm)
!     grav   = earth's gravitational acceleration (cm/sec**2)
!     cdbot  = bottom drag coefficient
!
!     ncon   = number of  passes through convective code in tracer
!     ntbt   = number of free surface timesteps per baroclinic velocity
!              time steps.
!     ntbt2  = 2*ntbt
# ifdef de_checkbd
!   dchkbd   = weighting factor for del-plus-del-cross filter.
# endif
!
#endif
      common /scalar/ dtts, dtuv, dtbt, c2dtts, c2dtuv, c2dtbt
     &,               area, volume, omega, radius, grav, cdbot
#ifdef de_checkbd
     &,               dchkbd
#endif
      common /scalri/ ncon, ntbt, ntbt2
#ifdef hcomments
!
!     non dimensional quantities:
!
!     radian = degrees per radian
!     pi     = something good to eat
!
#endif
      common /ndcon/ radian, pi

