#ifdef hcomments
!
!======================================================================
!
!    File:  slabs.h
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
!     This defines the basic array storage scheme.  In the MOM
!     model the data is stored on disk in slabs (hence the name).
!     For array processors it is stored in large arrays distributed
!     over the different processors
!
!     time levels defined in timelv.h
!
!     nm     = points to time level itt-2
!     nc     = points to current time level itt-1
!     np     = points to time level itt
!
!     the slab data is:
!
!     t = tracer quantites. the 1st is temperature &
!                           the 2nd is salinity
!         if nt > 2 then other tracers are allowed.
!     u = zonal component of velocity.
!     v = meridional component of velocity.
!
!     (note: only the baroclinic component of the velocity is
!      stored on the slab disks. in memory, the barotropic
!      component is added on to give the full velocity)
!
#endif
       common /slabs/ t(km,imt,jmt,nt,3),
     &                u(km,imt,jmt,3), v(km,imt,jmt,3)
#ifdef hcomments
!
!  the working arrays are
!
!     uclin  = internal mode u (used for advective fluxes &
!               diagnostics)
!     vclin  = internal mode v
!     rhoo   = density at (i  ,j  ) on the "t" grid
!     rhpo   = density at (i+1,j  ) on the "t" grid
!     rhop   = density at (i  ,j+1) on the "t" grid
!     rhpp   = density at (i+1,j+1) on the "t" grid
!     dpdx   = zonal gradient of pressure on "u,v" grid
!     dpdy   = meridional gradient of pressure on "u,v" grid
!     fue    = advective coeff for eastern face of "u,v" grid box
!              in "clinic" & "t" grid box in "tracer"
!     fuw    = advective coeff for western face of grid box
!     fvn    = advective coeff for northern face of grid box
!     fvs    = advective coeff for southern face of grid box
!     fm     = (0,1) over "t" grid (land,ocean) points
!     gm     = (0,1) over "u,v" grid (land,ocean) points
!     vmf    = array (1) used for vertical differences of u
!     vmf    = array (2) used for vertical differences of v
!     vtf    = array used for vertical tracer flux
!     fw     = vertical velocity defined at the bottom of "u,v"
!              boxes in "clinic" and "t" boxes in "tracer"
!     fwb1   = w * (quantity 1) defined at the bottoms of the boxes
!     fwb2   = w * (quantity 2) defined at the bottoms of the boxes
!     restr  = term for restoring surface tracers to prescribed
!              values via newtonain damping
!     rests  = time scale for restoring surface tracers (days)
!
#endif
       common /work/
     &     dpdx(km),  dpdy(km),
     &     fue(km), fuw(km),   fvn(km),  fvs(km),
     &     vmf(0:km,2),   vtf(0:km),
     &     fw(0:km),  fwb1(0:km), fwb2(0:km), temp (km),
     &     maskpo(km), maskmo(km), maskop(km), maskom(km),
     &     rhoo(km),   rhpo(km),   rhpp(km),   rhop(km)
#ifdef presetp
#ifdef hcomments
!
!  common press stores the pressure field precalculated for
!  subroutine clinic by subroutine setp
!
#endif
      common /press/  p(km,imt,jmt)
#endif


