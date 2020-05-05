c====================== include file "slabs.h" =========================
#ifdef hcomments
c
c @(#) SCCS module: slabs.h  version: 1.5
c     Creation date: 03/24/95
c
c-----------------------------------------------------------------------
c     This defines the basic array storage scheme.  In the MOM
c     model the data is stored on disk in slabs (hence the name).
c     For array processors it is stored in large arrays distributed
c     over the different processors
c
c     time levels defined in timelv.h
c
c     nm     = points to time level itt-2
c     nc     = points to current time level itt-1
c     np     = points to time level itt
c
c     the slab data is:
c
c     t = tracer quantites. the 1st is temperature & 
c                           the 2nd is salinity
c         if nt > 2 then other tracers are allowed.
c     u = zonal component of velocity.
c     v = meridional component of velocity.
c
c     (note: only the baroclinic component of the velocity is
c      stored on the slab disks. in memory, the barotropic
c      component is added on to give the full velocity)
c-----------------------------------------------------------------------
c
#endif
       common /slabs/ t(KM,IMT_S,JMT_S,NT,3),
     &                u(KM,IMT_S,JMT_S,3), v(KM,IMT_S,JMT_S,3)
#ifdef hcomments
c
c  the working arrays are
c 
c     uclin  = internal mode u (used for advective fluxes &
c               diagnostics)
c     vclin  = internal mode v 
c     rhoo   = density at (i  ,j  ) on the "t" grid
c     rhpo   = density at (i+1,j  ) on the "t" grid
c     rhop   = density at (i  ,j+1) on the "t" grid
c     rhpp   = density at (i+1,j+1) on the "t" grid
c     dpdx   = zonal gradient of pressure on "u,v" grid
c     dpdy   = meridional gradient of pressure on "u,v" grid
c     fue    = advective coeff for eastern face of "u,v" grid box
c              in "clinic" & "t" grid box in "tracer"
c     fuw    = advective coeff for western face of grid box
c     fvn    = advective coeff for northern face of grid box
c     fvs    = advective coeff for southern face of grid box
c     fm     = (0,1) over "t" grid (land,ocean) points
c     gm     = (0,1) over "u,v" grid (land,ocean) points
c     vmf    = array (1) used for vertical differences of u
c     vmf    = array (2) used for vertical differences of v
c     vtf    = array used for vertical tracer flux
c     fw     = vertical velocity defined at the bottom of "u,v" 
c              boxes in "clinic" and "t" boxes in "tracer"
c     fwb1   = w * (quantity 1) defined at the bottoms of the boxes
c     fwb2   = w * (quantity 2) defined at the bottoms of the boxes
c     restr  = term for restoring surface tracers to prescribed 
c              values via newtonain damping
c     rests  = time scale for restoring surface tracers (days)
c
#endif
       common /work/
     &     dpdx(KM),  dpdy(KM),
     &     fue(KM), fuw(KM),   fvn(KM),  fvs(KM),
     &     vmf(0:KM,2),   vtf(0:KM),
     &     fw(0:KM),  fwb1(0:KM), fwb2(0:KM), temp(KM),
     &     maskpo(KM), maskmo(KM), maskop(KM), maskom(KM),
     &     rhoo(KM),   rhpo(KM),   rhpp(KM),   rhop(KM)
#ifdef presetp
#ifdef hcomments
c
c  common press stores the pressure field precalculated for
c  subroutine clinic by subroutine setp
c
#endif
      common /press/  p(KM,IMT_S,JMT_S)
#endif
     

