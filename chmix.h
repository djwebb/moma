c====================== include file "chmix.h" =========================
#ifdef hcomments
c
c     @(#) SCCS module: chmix.h, version 1.9
c
c     horizontal mixing coefficients
c
c     am  = lateral eddy viscosity (cm**2/sec)
c     ah  = lateral eddy diffusivity (cm**2/sec)
c     aq  = lateral turbulent kinetic energy diffusivity 
c           (cm**2/sec)
#endif
      common /chmix/ am, ah, 
     &        bbu(jmt), ccu(jmt), ddu(jmt), ggu(jmt), hhu(jmt),
     &        bbt(jmt), cct(jmt), ddt(jmt)

