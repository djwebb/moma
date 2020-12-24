#ifdef hcomments
!
!======================================================================
!
!    File:  XXXXX
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
!     variables used for vertical diffusion
!
!     inputs: (set through namelist)
!
!     fkph   = vertical diffusion coefficient (cm**2/sec)
!     fkpm   = vertical viscosity coefficient (cm**2/sec)
!
!     derived quantities:
!
!     vvc  = vertical viscosity coeff
!     vdc  = vertical diffusion coeff
!
#endif
      common /cvmix/  fkph, fkpm, vvc(km), vdc(km)

