c====================== include file "index.h" =========================
#ifdef hcomments
c
c     @(#) SCCS module: index.h, version 1.2
c
c    various starting & ending indices for controlling where quantities
c    are calculated:
c
c    isz   = starting longitudinal index for stream function calculation
c    iez   = ending longitudinal index for stream function calculation
c
c    alonis= longiude of arbitrary point within island.
c    alatis= latitude of arbitrary point within island.
c            coordinate "alonis", "alatis" must point to a land "t"
c            grid box within the island.
c    nippts= number of island perimeter points
c    iperm = "i" coordinate for the island perimeter point
c    jperm = "j" coordinate for the island perimeter point
c    iofs  = offset for indexing into the island perimeter 
c            coordinates
c
#endif
      common /index1/ isz(jmt,lseg), iez(jmt,lseg)
      common /index2/ alonis(nisle), alatis(nisle), nippts(nisle)
     &,               iofs(nisle), iperm(maxipp), jperm(maxipp)

