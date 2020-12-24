#ifdef hcomments
!
!======================================================================
!
!    File:  thick.h
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
!     specifications for box thickness (cm) in the vertical
!
#endif
      data dzt /   30.00e2,   46.15e2,   68.93e2,   99.93e2,  140.63e2
     &      ,  192.11e2,  254.76e2,  327.95e2,  409.81e2,  497.11e2
     &      ,  585.36e2,  669.09e2,  742.41e2,  799.65e2,  836.10e2
     &   /
#ifdef hcomments
!
!     dzt = model level thicknesses computed as
!     dzt(km+1-k)=exp(-(k/sigma)**2) + offset) for k=1,km
!      km =  15
!      sigma=  8.204686
!      offset= 0.0000000
!      1st level thickness =  30.00000
!      # of levels with this thickness =  1
!      depth at bottom of last level =  5700.000
#endif

