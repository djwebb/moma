c====================== include file "iounit.h" ========================
#ifdef hcomments
c
c     @(#) SCCS module: iounit.h, version 1.9
c
c     i/o units & related variables
c
c     iorest  = unit for history files (restarts)
c     fnrest  = name of history file used to start run
c   note: ansi needs character variables in a separate common block
c
#endif
      character*80    fnrest
      common /iounit/ iorest
      common /iounit/ fnrest

