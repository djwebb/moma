c====================== include file "chmix.h" =========================
#ifdef hcomments
c
c @(#) SCCS module: chmix.h  version: 1.9
c     Creation date: 12/18/95
c
c-----------------------------------------------------------------------
c     am  = lateral eddy viscosity (cm**2/sec)
c     ah  = lateral eddy diffusivity (cm**2/sec)
c     aq  = lateral turbulent kinetic energy diffusivity 
c           (cm**2/sec)
c-----------------------------------------------------------------------
c
#endif
      common /chmix/ am, ah, 
     &        bbu(JMT_S), ccu(JMT_S), ddu(JMT_S), ggu(JMT_S), 
     &        hhu(JMT_S), bbt(JMT_S), cct(JMT_S), ddt(JMT_S)

