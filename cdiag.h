c====================== include file "cdiag.h" =========================
#ifdef hcomments
c
c @(#) SCCS module: cdiag.h  version: 1.1
c     Creation date: 03/16/95
c
c-----------------------------------------------------------------------
c     Variables used for computing diagnostics:
c
c     ektot    = total kinetic energy per unit volume
c     tddt     = rate of change of tracer per unit volume
c     dtabs    = average per unit volume of modulus of rate of change of tracer 
c     tvar     = rate of change of tracer variance per unit volume
c-----------------------------------------------------------------------
c
#endif
      common /cdiag/ ektot, tddt(nt), dtabs(nt),tvar(nt)

