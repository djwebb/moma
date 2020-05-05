c====================== include file "cvbc.h" ==========================
#ifdef hcomments
c
c     @(#) SCCS module: cvbc.h, version 1.8
c
c     vertical boundary condition variables:
c
c     smf = surface momentum flux   (dynes/cm**2)
c          1 => zonal wind stress  2 => meridional wind stress
c     bmf = bottom momentum flux
c          1 => zonal bottom drag  2 => meridional bottom drag 
c     stf = surface tracer flux    
c          1 => surface heat flux (cal/cm**2/sec = cm*degC/sec
c               = ly/sec) (assuming rho*cp = 1 cal/degC/cm**3)
c          2 => surface water flux (ppt/cm**2/sec)
c     btf = bottom tracer flux (for consistency but normally zero!)
c          1 => bottom heat flux   2 => bottom water flux
c
#endif
       common /cvbc/ smf(2),  bmf(2), stf(nt), btf(nt)

