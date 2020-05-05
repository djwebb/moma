c====================== include file "emode.h" =========================
c
c     variables for rigid lid calculation of external mode
c
c     p     = stream function 
c     ztd  = time change of vorticity
c     ptd  = time change of stream function
c     h    = depth over "u,v" points
c     hr   = reciprocal depth over "u,v" points
c     zu   = time change of vertically averaged zonal forcing 
c     zv   = time change of vertically averaged merid forcing 
c     sfu  = correction to u when calculating advection in clinic
c     sfv  = correction to v.  Calculated from the stream function.
c
c     cfn  = relaxation coefficient to the north of point (i,j)
c     cfs  = relaxation coefficient to the south of point (i,j)
c     cfe  = relaxation coefficient to the east of point (i,j)
c     cfw  = relaxation coefficient to the west of point (i,j)
c     cpf  = relaxation coefficient for point (i,j)
c     rmsk = relaxation mask of (0,1) on (land & perimeter, 
c            open ocean)
c     res   = relaxation residual
c     cof   = normalizing coefficient for 1
c     cofis = flow around 1
c
      common /emode/
     & cfn(imt,jmt), cfs(imt,jmt), cfe(imt,jmt), cfw(imt,jmt)
     &,ptd(imt,jmt,3), res(imt,jmt), cpf(imt,jmt)
     &,p  (imt,jmt,3), ztd(imt,jmt), hr(imt,jmt),   h(imt,jmt)
     &,zu(imt,jmt),    zv(imt,jmt),  sfu(imt,jmt),  sfv(imt,jmt)
     &,rmsk(imt,jmt)
     &,cof(imt,jmt), cofis(nisle)
