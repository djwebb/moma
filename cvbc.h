#ifdef hcomments
!
!======================================================================
!
!    File:  cvbc.h
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
!     vertical boundary condition variables:
!
!     smf = surface momentum flux   (dynes/cm**2)
!          1 => zonal wind stress  2 => meridional wind stress
!     bmf = bottom momentum flux
!          1 => zonal bottom drag  2 => meridional bottom drag
!     stf = surface tracer flux
!          1 => surface heat flux (cal/cm**2/sec = cm*degC/sec
!               = ly/sec) (assuming rho*cp = 1 cal/degC/cm**3)
!          2 => surface water flux (ppt/cm**2/sec)
!     btf = bottom tracer flux (for consistency but normally zero!)
!          1 => bottom heat flux   2 => bottom water flux
!
#endif
       common /cvbc/ smf(2),  bmf(2), stf(nt), btf(nt)

