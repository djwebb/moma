#ifdef hcomments
!
!======================================================================
!
!    File:  chmix.h
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
!     horizontal mixing coefficients
!
!     am  = lateral eddy viscosity (cm**2/sec)
!     ah  = lateral eddy diffusivity (cm**2/sec)
!     aq  = lateral turbulent kinetic energy diffusivity
!           (cm**2/sec)
#endif
      common /chmix/ am, ah,
     &        bbu(jmt), ccu(jmt), ddu(jmt), ggu(jmt), hhu(jmt),
     &        bbt(jmt), cct(jmt), ddt(jmt)
