c====================== include file "varinfo.h" =======================
#ifdef hcomments
c
c @(#) SCCS module: varinfo.h  version: 1.1
c     Creation date: 10/13/97
c
c common storage for variable descriptors
c
#endif
#ifdef cray-t3d
      logical*8 vartgrd
#else
      logical*4 vartgrd
#endif
      character varname*25,varunit*12
      common /varinfl/ vartgrd(NARCHV)
      common /varinfa/ varname(NARCHV),varunit(NARCHV)
