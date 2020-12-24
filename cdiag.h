#ifdef hcomments
!
!======================================================================
!
!    File:  cdiag.h
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
!     variables used for computing diagnostics:
!
!     ektot    = total kinetic energy per unit volume
!     tddt     = rate of change of tracer per unit volume
!     dtabs    = average per unit volume of modulus of rate of change of tracer
!     tvar     = rate of change of tracer variance per unit volume
!======================================================================
!
#endif
      common /cdiag/ ektot, tddt(nt), dtabs(nt),tvar(nt)

