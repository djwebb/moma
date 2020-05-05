c====================== include file "cvmix.h" =========================
#ifdef hcomments
c
c @(#) SCCS module: cvmix.h  version: 1.6
c     Creation date: 04/13/95
c
c-----------------------------------------------------------------------
c     Variables used for vertical diffusion
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
c-----------------------------------------------------------------------
c
#endif
      common /cvmix/  fkph, fkpm, vvc(KM), vdc(KM)

