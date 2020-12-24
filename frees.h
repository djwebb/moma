#ifdef hcomments
!
!======================================================================
!
!    File:  frees.h
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
!     variables for free surface model calculation of external mode
!
!     h0   =  free surface model sea surface height
!     u0   =  free surface model velocity east
!     v0   =  free surface model velocity north
!     h    = depth over "u,v" points
!     hr   = reciprocal depth over "u,v" points
!     zu   = vertically averaged zonal forcing
!     zv   = vertically averaged meridional forcing
# ifndef free_eb
!     freeav = time average of free surface fields
# endif
!
#endif
      common /fields/
     & h0(imt,jmt,3), u0(imt,jmt,3), v0(imt,jmt,3)
     &,h(imt,jmt),    hr(imt,jmt),    zu(imt,jmt),  zv(imt,jmt)
#ifndef free_eb
     &,freeav(3,imt,jmt)
#endif
