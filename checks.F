      subroutine checks
c
c @(#) SCCS module: checks.F  version: 1.8
c     Creation date: 08/23/95
c
c==================================================================
c     do consistency checks
c==================================================================
c
#include "def_master.h"
#include "param.h"
#include "scalar.h"
#include "switch.h"
#include "iounit.h"
c
c   check timesteps
c
      if( dtuv.ne.dtts)then
        write (stdout,9000)
     $  '==> Warning: use of unequal time steps implies the transient '
     $, '             response is unimportant and multiple equilibria '
     $, '             do not exist.                                   '
      endif
c
c   check archive and restart file formats
c
      if(ftrest.ne.'ufm'.and.
     &   ftrest.ne.'hdf '.and.init)then
        print *,' ftrest has invalid value'
        print *,' Program stopping... '
        stop
      endif
      if(ftsnap.ne.'asc '.and.
     &   ftsnap.ne.'hdf '.and.snapd.ne.0.0)then
        print *,' ftsnap has invalid value'
        print *,' Program stopping... '
        stop
      endif
      if(ftarch.ne.'ufm'.and.
     &   ftarch.ne.'hdf '.and.archd.ne.0.0)then
        print *,' ftarch has invalid value'
        print *,' Program stopping... '
        stop
      endif
#ifndef ncsa_hdf
      if(ftrest.eq.'hdf '.and.init)then
        print *,' hdf restart file not supported.  "ncsa_hdf" not set'
        print *,' Program stopping... '
        stop
      endif
      if(ftsnap.eq.'hdf '.and.snapd.ne.0.0)then
        print *,' hdf snapshot file not supported.  "ncsa_hdf" not set'
        print *,' Program stopping... '
        stop
      endif
      if(ftarch.eq.'hdf '.and.(archd.ne.0.0.or.restrt))then
        print *,' hdf restart file not supported.  "ncsa_hdf" not set'
        print *,' Program stopping... '
        stop
      endif
#endif
      return
9000  format (/,(1x,a80))
      end

