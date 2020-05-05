            subroutine docmnt (nunit, wrdoc, wrc0, wrc1, wrc2)
c
c @(#) SCCS module: docmnt.F  version: 1.7
c      Creation date: 06/30/95
c
c=======================================================================
c
c     This code arranges data to be written out as a form of
c     documentation summarizing the model characteristics that
c     uniquely define this run.
c     In order to reap the full benefits of this feature, the user will
c     need to review this routine when configuring a model run or
c     adding new options to the code.
c     We envision this feature as being beneficial to those setting up
c     new model runs, comparing different model runs, and for analysis
c     purposes (much of the information used in setting up mapping
c     procedures and other analysis products are contained in the
c     data written by this routine).
c
c     Five namelist writes are available.  Which ones are written
c     is determined from the 4 logical variables passed into the
c     subroutine (wrdoc, wrc0, wrc1, & wrc2).  The unit to which
c     they are written is set by "nunit".
c
c     If "wrdoc" is true, then the "ndocs" & "docums" namelists
c     are written...
c       The 1st, "ndocs" simply contains one variable (ndoc); the
c     dimension of "docum(ndoc)" found in the 2nd namelist write.
c       The 2nd, "docums" consists of a series of 60 character-long
c     lines of information the user should use to list the variables
c     in the subsequent namelist writes, and to describe the model
c     run in general.
c
c     If "wrc0" is true, then the "cntrl0" namelist is written...
c       Namelist "cntrl0" contains the 2 variables needed to restart
c     the model (itt, years).
c
c       The "cntrl1" and "cntrl2" namelists contain information
c     describing the model's configuration and specifications.
c
c     If "wrc1" is true, then the "cntrl1" namelist is written...
c       "cntrl1" contains variables (non-dimensioned) that will be
c     common to most, if not all, incarnations of the model code.
c
c     If "wrc2" is true, then the "cntrl2" namelist is written...
c       "cntrl2" is intended to contain the larger (dimensioned)
c     data that will help identify the model specifics to users
c     and analysis programs alike.
c
c=======================================================================
c  This routine is based on the corresponding mom routine but
c  contains only a few of the standard features
c=======================================================================
c
#include "param.h"
#include "switch.h"
#include "versno.h"
      parameter (ndocum = 16, ifdmax = 50)
c
      character*16 cifdef(ifdmax)
      logical wrdoc, wrc0, wrc1, wrc2
c
c-----------------------------------------------------------------------
c     list options that are in effect on first timestep
c-----------------------------------------------------------------------
c
      if (first) then
c
        do 200 n=1,ifdmax
          cifdef(n) = '                '
 200    continue
        nifdef = 0
        cifdef(1) = '0 ifdefs set'
c
#ifdef Masterslave
        nifdef = nifdef + 1
        cifdef(nifdef) = 'Masterslave     '
#endif
#ifdef Master
        nifdef = nifdef + 1
        cifdef(nifdef) = 'Master          '
#endif
#ifdef Slave
        nifdef = nifdef + 1
        cifdef(nifdef) = 'Slave           '
#endif
#ifdef cyclic_master
        nifdef = nifdef + 1
        cifdef(nifdef) = 'cyclic_master   '
#endif
#ifdef oldad
        nifdef = nifdef + 1
        cifdef(nifdef) = 'oldav           '
#endif
#ifdef presetp
        nifdef = nifdef + 1
        cifdef(nifdef) = 'presetp         '
#endif
#ifdef inline
        nifdef = nifdef + 1
        cifdef(nifdef) = 'inline          '
#endif
#ifdef free_eb
        nifdef = nifdef + 1
        cifdef(nifdef) = 'free_eb         '
#endif
#ifdef de_checkbd
        nifdef = nifdef + 1
        cifdef(nifdef) = 'de_checkbd      '
#endif
#ifdef ncsa_hdf
        nifdef = nifdef + 1
        cifdef(nifdef) = 'ncsa_hdf        '
#endif
#ifdef DOUBLE
        nifdef = nifdef + 1
        cifdef(nifdef) = 'DOUBLE          '
#endif
#ifdef cray-t3d
        nifdef = nifdef + 1
        cifdef(nifdef) = 'cray-t3d        '
#endif
#ifdef no_namelist
        nifdef = nifdef + 1
        cifdef(nifdef) = 'no_namelist     '
#endif
#ifdef pvm_buffer
        nifdef = nifdef + 1
        cifdef(nifdef) = 'pvm_buffer      '
#endif
#ifdef SYNC_2
        nifdef = nifdef + 1
        cifdef(nifdef) = 'SYNC_2          '
#endif
#ifdef SYNC_3
        nifdef = nifdef + 1
        cifdef(nifdef) = 'SYNC_3      '
#endif
        write (stdout,'(/,a45)')
     $   ' The following "ifdef" options are in effect:'
        nn = (nifdef+4)/5
        do 300 i=1,nn
          ns = 1+(5*(i-1))
          ne = ns + 4
          if (ne .gt. nifdef) ne = nifdef
          write (stdout,9101) (cifdef(n),n=ns,ne)
300     continue
        write (stdout,'(a21,/)') ' end of options list.'
      endif
      return
9101  format(5(3x,a16))
      end

