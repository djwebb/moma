#ifdef hcomments
!
!======================================================================
!
!    File:  param.h
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
!     main parameter file which sets ocean characteristics:
!
!     imt    = number of grid points in longitudinal direction
!     jmt    = number of grid points in latitudinal direction
!     km     = number of grid points in the vertical direction
!     nt     = number of tracers
!     lseg   = maximum number of longitudinal stream function
!              segments
!
#endif
      parameter  (imt=92, jmt=36, km=15, nt=2, imu=imt)
      parameter  (imtp1=imt+1, imtm1=imt-1, imtm2=imt-2,
     &            imum1=imu-1, imum2=imu-2, jmtp1=jmt+1,
     &            jmtm1=jmt-1, jmtm2=jmt-2,
     &            istart=2, iend=imtm1, jscan=jmtm2,
     &            kmp1=km+1, kmp2=km+2, kmm1=km-1)
#ifdef hcomments
!
!     add parameter constants
!
#endif
#include "pconst.h"
