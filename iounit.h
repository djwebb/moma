c====================== include file "iounit.h" ========================
#ifdef hcomments
c
c @(#) SCCS module: iounit.h  version: 1.3
c     Creation date: 03/17/95
c
c-----------------------------------------------------------------------
c     i/o units & related variables
c
c     iorest  = unit for master archive files (restarts)
c     iosnap  = unit for master snapshot files 
c     fnrest  = name of master archive file used to start run
c     ftrest  = file type of archive file used to start run
c               'unf' := unformatted, 'hdf' := hdf file
c     ftarch  = file type used for archive files
c               'unf' := unformatted, 'hdf' := hdf file
c     ftsnap  = file type for snapshot files
c               'asc' := ascout, 'hdf' := hdf file
c     ioslave = unit used for standard output from slave
c     fnslave = name of file used by slave for output
c   note: ansi needs character variables in a separate common block
c-----------------------------------------------------------------------
c
#endif
      character    fnrest*80,ftrest*4,ftsnap*4,ftarch*4
      common /iouniti/ iorest,iosnap,ioslave
      common /iounitc/ fnrest,ftrest,ftsnap,ftarch
