c====================== include file "cvmix.h" =========================
#ifdef hcomments
c
c     @(#) SCCS module: cvmix.h, version 1.14
c
c     variables used for vertical diffusion
c
c     inputs: (set through namelist)
c
c     fkph   = vertical diffusion coefficient (cm**2/sec)
c     fkpm   = vertical viscosity coefficient (cm**2/sec)
c
c     derived quantities:
c
c     vvc  = vertical viscosity coeff
c     vdc  = vertical diffusion coeff
c
#endif
      common /cvmix/  fkph, fkpm, vvc(km), vdc(km)

