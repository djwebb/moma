#ifdef hcomments
!
!======================================================================
!
!    File:  grdvar.h
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
!  This version assumes all dx and dy boxes are of equal size when
!  measured in degrees.  Variables dxt, dxu etc replaced by dx,
!  dyt,dyu etc by dy and dzt by dz.
!
!     dx     = longitudinal width of "t" grid box (in cm)
!     dxr    = reciprocal of "dxt"
!     dx2r   = reciprocal of "2*dxt"
!     dx4r   = reciprocal of "4*dxt"
!
!     dy     = latitudinal height of "t" grid box (in cm)
!     dyr    = reciprocal of "dyt"
!     dy2r   = reciprocal of "2*dyt"
!     dy4r   = reciprocal of "4*dyt"
!
!     csu     = cosine of "u,v" grid point latitude
!     csur    = reciprocal of "csu"
!     cst     = cosine of "t" grid point latitude
!     cstr    = reciprocal of "cst"
!     phi     = latitude of "u,v" grid point in radians
!     phit    = latitude of "t" grid point in radians
!     sine    = sine of "u,v" grid point latitude
!     tng     = tan of "u,v" grid point latitude
!     fcor    = 2*omega*sine(j)
!
!     dz(k)  = level thickness of "t" and "u,v" grid boxes (in cm)
!               (vertical separation between "w" velocity points)
!     c2dt(k)= "2*dzt"
!     dzr(k) = reciprocal of dzt
!     dz2r(k)= reciprocal of "2*dzt"
!     dzw(k)  = vertical separation between centers of levels k & k+1
!     dzwr(k) = reciprocal of dzw
!     dzw2r(k)= reciprocal of "2*dzw"
!     dzur(k)= upper diffusion grid factor = 1.0/(dzw(k-1)*dz(k))
!     dzlr(k)= lower diffusion grid factor = 1.0/(dzw(k)*dz(k))
!     dzwur(k)= upper diffusion grid factor = 1.0/(dz(k)*dzw(k))
!     dzwlr(k)= lower diffusion grid factor = 1.0/(dz(k+1)*dzw(k))
!     tanra(j)= tan(j)/radius on the "u,v" grid
!
#endif
      common /grdvar/ dx, dxr, dx2r, dx4r, dy, dyr, dy2r, dy4r,
     &                csu(jmt),  csur(jmt),  cst(jmt),   cstr(jmt),
     &                phi(jmt),  phit(jmt),  sine(jmt),  tng(jmt),
     &                fcor(jmt), tanra(jmt),
     &                dz(km),    dzr(km),    dz2r(km),   c2dz(km),
     &                dzw(0:km), dzwr(0:km), dzw2r(0:km)

