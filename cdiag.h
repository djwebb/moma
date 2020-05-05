c====================== include file "cdiag.h" =========================
c
c     variables used for computing diagnostics:
c
c     ektot    = total kinetic energy per unit volume
c     tddt     = rate of change of tracer per unit volume
c     dtabs    = average per unit volume of modulus of rate of change of tracer 
c     tvar     = rate of change of tracer variance per unit volume
c
      common /cdiag/ ektot, tddt(nt), dtabs(nt),tvar(nt)
