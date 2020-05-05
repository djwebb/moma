c====================== include file "cvbc.h" ==========================
#ifdef hcomments
c
c @(#) SCCS module: cvbc.h  version: 1.10
c     Creation date: 08/09/96
c
c-----------------------------------------------------------------------
c
c   Vertical boundary condition variables:
c
c     smf = surface momentum flux   (dynes/cm**2)
c           1 - zonal wind stress  
c           2 - meridional wind stress
c     bmf = bottom momentum flux
c           1 - zonal bottom drag  
c           2 - meridional bottom drag 
c     stf = surface tracer flux    
c           1 - surface heat flux (cal/cm**2/sec = cm*degC/sec
c               = ly/sec) (assuming rho*cp = 1 cal/degC/cm**3)
c           2 - surface water flux (ppt/cm**2/sec)
c     btf = bottom tracer flux (for consistency but normally zero!)
c           1 - bottom heat flux   
c           2 - bottom water flux
c
c-----------------------------------------------------------------------
c
#endif
       common /cvbc/ smf(2),  bmf(2), stf(nt), btf(nt)
#ifdef hcomments
c
c-----------------------------------------------------------------------
c  surface flux scheme variables
c-----------------------------------------------------------------------
c
#endif
      parameter (NOLAT = 40, NOLATP = NOLAT+1)
      logical        lcycw,lcycs
      common /metl/  lcycw,lcycs
      common /meti/  i_w,j_w,i_s,j_s
      common /metr/  stlon_w,stlat_w,dx_w,dy_w,
     &               stlon_s,stlat_s,dx_s,dy_s
#ifdef Master 
      logical        lcycm
      common /mmetl/ lcycm
      common /mmeti/ iw_s(MXSLAVE),iw_e(MXSLAVE),iw_l(MXSLAVE),
     &               jw_s(MXSLAVE),jw_e(MXSLAVE),jw_l(MXSLAVE),
     &               is_s(MXSLAVE),is_e(MXSLAVE),is_l(MXSLAVE),
     &               js_s(MXSLAVE),js_e(MXSLAVE),js_l(MXSLAVE)
      common /mmetr/ swlo_w(MXSLAVE),swla_w(MXSLAVE),
     &               swlo_s(MXSLAVE),swla_s(MXSLAVE),
     &               wstx(NWX_M,NWY_M),wsty(NWX_M,NWY_M),
     &               ssta(NSX_M,NSY_M),sssa(NSX_M,NSY_M),
     &               sstobs(NOLAT), salobs(NOLAT),
     &               wsxobs(NOLATP),wsyobs(NOLATP)
#endif
#ifdef Slave
      common /smetr/ wstx(NWX_S,NWY_S,3),wsty(NWX_S,NWY_S,3),
     &               ssta(NSX_S,NSY_S,3),sssa(NSX_S,NSY_S,3)
#endif
            
