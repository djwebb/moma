c====================== include file "iounit.h" ========================
#ifdef hcomments
c
c     @(#) SCCS module: iounit.h, version 1.15
c     Creation date: 03/14/95
c
c-----------------------------------------------------------------------
c     i/o units & related variables
c
c     iorest  = unit for archive files (restarts)
c     iosnap  = unit for snapshot files 
c     fnrest  = name of archive file used to start run
c     iokmt   = unit for kmt file (depths)
c     fnkmt   = name of kmt file used
c     ftrest  = file type of archive file used to start run
c               'std' := restart, 'hdf' := hdf file
c     ftarch  = file type used for archive files
c               'std' := restart, 'hdf' := hdf file
c     ftsnap  = file type for snapshot files
c               'std' := ascout, 'hdf' := hdf file
c   note: ansi needs character variables in a separate common block
c
#endif
      character    fnrest*80,fnkmt*80,ftrest*4,ftsnap*4,ftarch*4
      common /iouniti/ iorest,iokmt
      common /iounitc/ fnrest,fnkmt,ftrest,ftsnap,ftarch

