c======================include file: mesdta.h==========================
#ifdef hcomments
c
c @(#) SCCS module: mesdta.h  version: 1.3
c     Creation date: 03/17/95
c
c-----------------------------------------------------------------------
c Common storage for miscellaneous message handling data
c
c  itida     -  array of process id's
c  nproc     -  number of processes in use
c  mytid     -  tid of this pricess
c  mtid      -  tid of master process
c  me        -  number of this process (0 = master)
c  icycl     -  (1 or 0) indicating (whether or not) slave is 
c               self-cyclic
c  infoa     -  array of progress information for each slave
c
c  tsi:
c
c  ltsip    - set true in master when waiting for tsi data 
c           - set true in slave when tsi data being processed
c  ltsia    - array set true by master when tsi information received 
c             from slave
c  nstsi    - number of slaves still to return tsi information
c  ltsiw    - set true when slave has tsi data waiting for master
c  itttsi   - timestep for which tsi information is being collected
c  timetsi  - timestamp for tsi information
c
c  snapshot & archive:
c
c  lregra   - logical register map. Set true in master for each
c             station filled in response to an archive or snapshot
c             request
c  nseat    - integer number of sea-stations (t-grid)
c  nseau    - integer number of sea-stations (u-grid)
c  nsubst   - integer array, number of sea-stations (t-grid) for 
c             each block
c  nsubsu   - integer array, number of sea-stations (u-grid) for
c             each block
c  nseachk  - integer containing the current reference number 
c             of sea points (set in archreq to one of nseat, nseau,
c             nsubst(n) or nsubsu(n) dependent on current archive 
c             variable)
c  idvar    - integer indicating which archive variable is being 
c             processed
c
c  snapshot:
c
c  lsnapp   - master: set true while snapshots being processed
c             slave: set true while snapshots in buffer waiting to
c             be sent to master
c  lsnapq   - second flag used when snapshots being processed
c  lsnapa   - array set true by master when snap information 
c             received from slave
c  nssnap   - number of slaves still to return snapshot information
c  ittsnap  - timestep for which snapshot information is being 
c             collected
c  idsnp    - integer indicating which snapshot variable is being
c             processed
c
c  archive:
c
c  larchp   - master: set true while archive being processed
c             slave: set true while archive in buffer waiting to
c             be sent to master
c  larchq   - second flag used when archive being processed. This
c             second flag is used recursively by the slave for each
c             partial request from the master.
c  lidvra   - array whose idvar'th element is set true by master when 
c             the idvar'th variable type is being collected.
c  nsarch   - number of slaves to broadcast archive request to
c  iarcha   - used to collected processor ids when broadcasting to
c	      a sub-set of processors
c  narcha   - used to collected processor indices when broadcasting to
c	      a sub-set of processors
c  ittarch  - timestep for which archive information is being collected
c  jarchl   - lower bound of the last sub-volume processed by master
c  jarchu   - upper bound of the last sub-volume processed by master
c
c  flux:
c
c  lfluxa   - array set true when flux data sent to slave
c
c  others:
c
c  utime     - array holding elasped execution time (master and slave)
c              1 = process, 2 = system, 3 = sum 
c  iotime    - time (in seconds since Jan 1st 1970) when status arrays 
c              last printed
c  intime    - time when master started
c  iposn     - in slave last reporting position 
c  lbrkpt    - in slave true when permission received to proceed 
c              from breakpoint (slave)
c  outstr    - in slave, character*80 used to send messages 
c
c  Slave parameters:
c
c The following are held in arrays on the master and elements are 
c passed to individual scalars (defined in scalar.h) on each slave.
c
c  it_s,jt_s - start position of each slave region in master indices
c  it_e,jt_e - end position of each slave region in master indices
c  it_l,jt_l - length of each slave region in master indices
c              Note: start and end indices should always be given as
c              positive integers. These will be converted to modulo 
c              IMT_M by the message handling routines.
c  swlo      - swlo and swla contain the long-lat co-ordinates of 
c  swla        the sw corner of the the most south-westerly halo point.
c              NOTE: For the "single" slave test case there are no top  
c              or bottom halo rows.
c  icycl_M   - (1 or 0) for each slave indicating (whether or not) it
c              is self-cyclic.
c  iswm1_M   - iswm1 and jswm1 are the master indices for the cell to 
c  jswm1_M     the south-west of the most south-westerly halo point 
c              (even it if doesn't exist). These values are used as 
c              offsets by the slaves to determine the local indices of
c              data passed through messages. I.e.:
c
c              slave(i_s,j_s) --> master(mod(iswm1+i_s,imt_m),jswm1+j_s)
c
c              Thus all messages passed from and to slaves which refer 
c              to specific cells should use the master indexing.
c
c The following are held and used only by the master:
c
c  messtr    - Character array containing descriptions of possible
c              reported positions within the model
c  mesind    - index array used to reference the messtr array. This
c              permits position numbers to be widely spaced without
c              requiring redundant storage within the character array.
c              Current maximum position number is: 40. The declarations
c              below and the data statements in caesar.F will need to
c              be changed if this number is exceeded.
c
c The following are held and used only by the slaves:
c
c  ia        - i values of interior and inner halo points
c  ja        - j values of interior and inner halo points
c  iota      - array set true when boundary data has been processed
c  ninner    - number of interior points
c  nhalo     - number of inner halo points
c  ncore     - ninner + nhalo
c
c  isd       - i, j values of inner halo points and processors to which
c              they must be sent
c  nbound    - number of boundary points to be sent
c  nsend     - array of number of points to be sent to each processor
c  nfalse    - counter for outstanding boundary point requests
c  nactive   - number of "active" neighbours (i.e. those neighbours
c              which receive data from "me"
c  nactme    - array containing the slave numbers of the active
c              neighbours
c  nactid    - array containing the process ids of the active
c              neighbours
c  itrpclr   - number of active neighbours waiting for barotropic data
c  iclnclr   - number of active neighbours waiting for baroclinic data
c
c  iout      - i values of outer halo points
c  jout      - j values of outer halo points
c  nouter    - number of outer halo points
c
c  bclnflg   - array of logical flags to indicate readiness to receive
c              baroclinic boundary information.
c  btrpflg   - array of logical flags to indicate readiness to receive
c              barotropic boundary information.
c  lbcln     - logical variable set true when all nactive elements of
c              bclnflg are set true
c  lbtrp     - logical variable set true when all nactive elements of 
c              btrpflg are set true
c  nmax      - max processor to which messages must be sent
c
c  itt3      = itt value when 3-D data sent, used as a check
c  ittp3     = pass of timestepping scheme when 3-D data sent
c  np3       = pointer to np value into which 3-D data to be received
c  itt2      = itt value when 2-D data sent, used as a check
c  itbt2     = itbt value when 2-D data sent, used as a check
c  itbtp2    = itbtp value when 2-D data sent, used as a check
c  np02      = pointer to np0 value into which 2-D data to be received
c
#endif
#include "fpvm3.h"
#ifdef DOUBLE 
# define REAL_PVM REAL8 
#else 
# define REAL_PVM REAL4 
#endif 
#ifdef cray-t3d
# define INTEGER_PVM INTEGER8
#else
# define INTEGER_PVM INTEGER4
#endif
      character*80    outstr
      character*32    timetsi
      logical*4       ltsia(MXSLAVE), lsnapa(MXSLAVE), lfluxa(MXSLAVE),
     &                ltsip, ltsiw, lsnapp, lsnapq, lbrkpt, 
     &                larchp, larchq, lidvra(NARCHV)
      common /mesdti/ itida(MXSLAVE),nproc,mtid,mytid,me,iposn,
     &                infoa(3,MXSLAVE),iotime,intime,itime_start,
     &                itttsi,ittsnap,ittarch,nstsi,nssnap,nsarch,
     &                nseat,nseau,idvar,idsnp,jarchl,jarchu
      common /mesdta/ outstr,timetsi
      common /mesdtl/ ltsia,  lsnapa, lfluxa, ltsip, ltsiw, 
     &                lsnapp, lsnapq, larchp, larchq, lidvra, lbrkpt
      common /mesdtr/ utime(3),time_start
#ifdef pvm_buffer
      common /pvmbuf/ bufout(LENBUF)
#endif
#ifdef Master
c
      character*30    messtr(0:12)
      logical*4       lregra(IMT_M,JMT_M)
      common /masdtal/lregra
      common /masdtai/it_s(MXSLAVE), jt_s(MXSLAVE),
     &                it_e(MXSLAVE), jt_e(MXSLAVE),
     &                it_l(MXSLAVE), jt_l(MXSLAVE), 
     &                icycl_M(MXSLAVE),
     &                iswm1_M(MXSLAVE), jswm1_M(MXSLAVE),
     &                nsubst(NZONE_M), nsubsu(NZONE_M),
     &                iarcha(MXSLAVE), narcha(MXSLAVE),
     &                mesind(0:40), nseapt, nseachk
      common /masdtar/swlo(MXSLAVE),swla(MXSLAVE)
      common /masdtaa/messtr
#else
c
      logical*4       iota(IMT_S,JMT_S,2), bclnflg(MXSLAVE),
     &                btrpflg(MXSLAVE), lbcln, lbtrp
      common /slvdtal/iota, bclnflg, btrpflg, lbcln, lbtrp
      common /slvdtai/ia(IMT_S*JMT_S), ja(IMT_S*JMT_S), isd(IBOUND,3), 
     &                nsend(MXSLAVE),  nfalse(2), nactme(MXSLAVE),
     &                nactid(MXSLAVE), iout(IBOUND), jout(IBOUND), 
     &                nactive, ninner, nhalo, ncore, nbound, nouter, 
     &                nmax, icycl, itrpclr, iclnclr,
     &                itt3,ittp3,np3,itt2,itbt2,itbtp2,np02
#endif
c 
c Message types defined as Fortran parameters. 
c
      parameter( MSG_ANY     = -1 )
      parameter( MSG_ABORT   =  1 )
      parameter( MSG_BYE     =  1 )
      parameter( MSG_POSN    =  2 )
      parameter( MSG_CONT    =  3 )
      parameter( MSG_OUTSTR  =  3 )
      parameter( MSG_IDS     = 10 )
      parameter( MSG_CONTROL = 11 )
      parameter( MSG_TOPOG   = 12 )
      parameter( MSG_RQ_TSI  = 20 )
      parameter( MSG_RY_TSI  = 20 )
      parameter( MSG_RQ_ARC  = 22 )
      parameter( MSG_RY_ARC  = 22 )
      parameter( MSG_ARC_CLR = 23 )
      parameter( MSG_TIMEVAR = 24 )
      parameter( MSG_RESTART = 25 )
      parameter( MSG_RQ_SNAP = 26 )
      parameter( MSG_RY_SNAP = 27 )
      parameter( MSG_BOUND1  = 100)
      parameter( MSG_BOUND2  = 101)
      parameter( MSG_CLR_BD1 = 102)
      parameter( MSG_CLR_BD2 = 103)
