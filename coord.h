#ifdef hcomments
!
!======================================================================
!
!    File:  coord.h
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
!                    model grid point coordinates
!
!     grid definition:
!
!     The model uses a staggard Arakawa "b" grid. horizontally,
!     tracer quantities are defined at the centers of "t" grid
!     boxes and velocities are defined at the centers of "u,v" grid
!     boxes.  the centers of "u,v" grid boxes are located at the
!     northeast corners of "t" grid boxes. the 1st "t" grid box is
!     located in the southwest corner of the "t" grid. this grid
!     system is replicated and stacked vertically one on top of
!     another from the surface of the ocean downward. vertically,
!     tracers & velocities are defined at the centers of their
!     corresponding boxes and are at the same depths.
!
!     primary grid information:
!
!     stlon = starting longitude of the "t" grid in degrees. It
!             defines the west face of the 1st "t" grid box.
!     stlat = starting latitude of the "t" grid in degrees. It
!             defines the south face of the 1st "t" grid box
!     dxdeg = east-wast width of grid in degrees
!     dydeg = north-south width of grid in degrees
!
!     xt(i) = longitude of the ith "t" point in degrees.   i=1..imt
!     xu(i) = longitude of the ith "u,v" point in degrees. i=1..imt
!     yt(j) = latitude of the jth "t" point in degrees.    j=1..jmt
!     yu(j) = latitude of the jth "u,v" point in degrees.  j=1..jmt
!     zt(k) = distance (cm) from surface down to center of level k
!             (depth of "t" and "u,v" grid points)
!     zw(k) = distance (cm) from surface down to bottom of level k
!
#endif
      common /coord/ stlon, stlat, dxdeg, dydeg, xt(imt), yt(jmt),
     &               xu(imt), yu(jmt), zw(km), zt(kmp1)

