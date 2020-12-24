#ifdef hcomments
!
!======================================================================
!
!    File:  timelv.h
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
!     time level pointers (indices)
!
!     baroclinic pointers
!     nm = time level pointer (index) corresponding to "itt-2"
!     nc = time level pointer (index) corresponding to "itt-1"
!     np = time level pointer (index) corresponding to "itt"
!
!     free surface model
!     nbt = number of free surface model timesteps "dtbt"
!                               per baroclinic timestep "dtuv"
!     nm0 = time level pointer to "itt0-2"
!     nc0 = time level pointer to "itt0-1"
!     np0 = time level pointer to "itt0"
!       where "itt0" is the free surface model timestep.
!
#endif
      common /ctlevp/ nm, nc, np, nbt, nm0, nc0, np0

