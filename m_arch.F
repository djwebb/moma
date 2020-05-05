      subroutine archrd(ncntl)
c
c @(#) SCCS module: m_arch.F  version: 1.7
c     Creation date: 06/30/95
c
c=======================================================================
c subroutine to read restart data and send to slaves
c=======================================================================
c
#include "pconst.h"
#include "iounit.h"
c
      if(ftrest.eq.'ufm')then
        call archrd_ufm(ncntl)
#ifdef ncsa_hdf
      elseif(ftrest.eq.'hdf')then
        call archrd_hdf(ncntl)
#endif
      else
        write(stdout,*) 'Archive write: file type "',
     &                  ftrest,'" not supported'
      endif
      return
      end


      subroutine archrd_ufm(ncntl)
c
c=======================================================================
c  subroutine to read an unformatted file c=======================================================================
c
#include "def_master.h"
#include "param.h"
#include "archive.h"
#include "iounit.h"
#include "levind.h"
#include "ctmngr.h"
#include "mesdta.h"
c
c-----------------------------------------------------------------------
c  open archive file when ncntl zero (entry from routine init_kmt).
c  and read basic information
c-----------------------------------------------------------------------
c
      if(ncntl.eq.0)then
        if(fnrest(1:1).ne.' ')then
          open(iorest, file=fnrest, form='UNFORMATTED',
     &         access='SEQUENTIAL',status= 'OLD',err=666)
        else
          open(iorest, form='UNFORMATTED',
     &         access='SEQUENTIAL',status= 'OLD',err=666)
        endif
        read(iorest) itt, totsec, totday, years
        write(stdout,24) itt, iorest, fnrest
        read(iorest) kmt
        return
      endif
  24  format('  Reading restart file for itt = ',i7,/,
     &       '  unit = ',i2,'  filename = ',a)
c
c------------------------------------------------------------------
c   read 2-D archive data arrays and send to slaves
c------------------------------------------------------------------
c
      do 50 n=1,nproc
        iarcha(n)=itida(n)
        narcha(n)=n
  50  continue
      jarchl = 1
      jarchu = JMT_M
      nsarch = nproc
      do 100 idvar1=1,NUM2D
        idvar = idvar1
        read (iorest) ((rest2d(i,j),i=1,IMT_M),j=1,JMT_M)
        call m_send(MSG_RESTART)
        call mchkarch(2)
 100  continue  
c
c------------------------------------------------------------------
c   read 3-D archive data arrays and send to slaves
c------------------------------------------------------------------
c
      nsub = 1+(JMT_M-1)/JSUB_M
      do 400 idvar1=NUM2D+1,NARCHV
      idvar = idvar1
      do 400 jj = 1,nsub
        jarchl = 1+(jj-1)*JSUB_M
        jarchu = min(jarchl-1+JSUB_M, JMT_M)
        nsarch = 0
        do 200 n=1,nproc
          if(jarchu.gt.jt_s(n).and.jarchl.le.jt_e(n))then
            nsarch=nsarch+1
            iarcha(nsarch) = itida(n)
            narcha(nsarch) = n
          endif
 200    continue
        do 300 j=1,jarchu-jarchl+1
          read(iorest) ((rest3d(k,i,j),k=1,KM),i=1,IMT_M)
 300    continue
          call m_send(MSG_RESTART)
          call mchkarch(2)
 400  continue
c
c   tell slaves to continue past checkpoint 
c   and close restart file
c
      call m_send(MSG_CONT)
      close(unit=iorest)
      write(stdout,*) ' Restart complete'
      return
c
c  error
c
 666  write(stdout,*) ' Error: restart dataset not found'
      write(stdout,*) ' restart failed.....aborting'
      call m_send(MSG_ABORT)
      end



      subroutine archreq
c
c==================================================================
c
c subroutine to decide on the next partial archive request and
c prepare the required message information ready for packing
c
c==================================================================
c
#include "def_master.h"
#include "param.h"
#include "varinfo.h"
#include "mesdta.h"
c
c------------------------------------------------------------------
c  2-D arrays 
c  request from all processors
c------------------------------------------------------------------
c
      if(idvar.lt.NUM2D) then
        idvar = idvar + 1
        lidvra(idvar) = .true.
        if(idvar.ne.1) lidvra(idvar-1)=.false.
        call initregister(vartgrd(idvar))
        call initrest
        nseachk = nseat
        if(.not.vartgrd(idvar)) nseachk = nseau
        nseapt = 0
        do 50 n=1,nproc
         iarcha(n)=itida(n)
  50    continue
        nsarch=nproc
        jarchl = 1
        jarchu = JMT_M
c
c------------------------------------------------------------------
c  3-D arrays 
c  requestv a band at a time
c------------------------------------------------------------------
c
      else
c
c------------------------------------------------------------------
c check jarchu
c------------------------------------------------------------------
c
        if(jarchu.ge.JMT_M) then
          idvar = idvar+1
c
c------------------------------------------------------------------
c  if archive complete then close dataset, reset flags and inform 
c  the slaves
c------------------------------------------------------------------
c
          if(idvar.gt.NARCHV) then
            write(stdout,*) ' Archive complete'
            larchp = .false.
            call m_send(MSG_ARC_CLR)
            return
          endif
c
c------------------------------------------------------------------
c  new variable
c------------------------------------------------------------------
c

          lidvra(idvar) = .true.
          if(idvar.ne.1) lidvra(idvar-1)=.false.
          jarchu=0
          call initregister(vartgrd(idvar))
        endif
c
c------------------------------------------------------------------
c  new block
c------------------------------------------------------------------
c

        call initrest
        jarchl = jarchu+1
        jarchu = min(jarchu + JSUB_M,JMT_M)
c
        mzon   = (jarchu-1)/JSUB_M + 1
        nseachk = nsubst(mzon)
        nseapt = 0
        if(.not.vartgrd(idvar)) nseachk = nsubsu(mzon)
c
        nsarch = 0
        do 100 n=1,nproc
          if(jarchu.ge.jt_s(n).and.jarchl.le.jt_e(n)) then
            nsarch=nsarch+1
            iarcha(nsarch) = itida(n)
          endif
 100    continue
      endif
c
c------------------------------------------------------------------
c set "partial archive out-standing" flag and send request
c------------------------------------------------------------------
c
      larchq = .true.
      call m_send(MSG_RQ_ARC)
      return
      end
        
      subroutine initregister(tgrid)
c
c Subroutine to initialise register array which is used to keep tabs
c on slave responses to archive requests
c
c tgrid  - (true,false) if the register (should,shouldn't) be initialised
c          for a variable stored on the t-grid.
c
#include "param.h"
#include "levind.h"
#include "mesdta.h"
      logical tgrid
c
c Initialise register to false but set land points true, since these 
c won't be sent by the slaves
c
      if(tgrid) then
       do 20 j=1,JMT_M
       do 20 i=1,IMT_M
        lregra(i,j) = .false.
        if(kmt(i,j).eq.0) lregra(i,j) = .true.
  20   continue
      else
       do 30 j=1,JMT_M 
       do 30 i=1,IMT_M 
        lregra(i,j) = .false. 
        if(kmu(i,j).eq.0) lregra(i,j) = .true. 
  30   continue
      endif
c
      return
      end
      
      subroutine initrest
c
c subroutine to re-initialise buffer storage between uses
c
#include "param.h"
#include "archive.h"
#include "mesdta.h"
c
      if(idvar.le.NUM2D) then
c
c Initilaise 2-D storage to zero everywhere
c
      do 300 j=1,JMT_M
      do 300 i=1,IMT_M
        rest2d(i,j) = c0
300   continue
c
      elseif(idvar.eq.NUM2D+1) then
c
c Initialise temperature to -2. everywhere
c
      do 301 j=1,JSUB_M
      do 301 i=1,IMT_M
      do 301 k=1,KM
        rest3d(k,i,j)= -2.0
301   continue
c
      elseif(idvar.eq.NUM2D+2) then
c
c Initialise salinity to 0.01 everywhere
c
      do 302 j=1,JSUB_M
      do 302 i=1,IMT_M
      do 302 k=1,KM  
        rest3d(k,i,j)= 0.01 
302   continue
c
      else
c 
c Initialise to zero everywhere
c 
      do 303 j=1,JSUB_M 
      do 303 i=1,IMT_M 
      do 303 k=1,KM   
        rest3d(k,i,j)= c0
303   continue
c
      endif
      return
      end

      subroutine archwr
c
c=======================================================================
c subroutine to write out the next partial archive
c=======================================================================
c
#include "pconst.h"
#include "iounit.h"
c
      if(ftarch.eq.'ufm')then
        call archwr_ufm
#ifdef ncsa_hdf
      elseif(ftarch.eq.'hdf')then
        call archwr_hdf
#endif
      else
        write(stdout,*) 'Archive write: file type "',
     &                  ftarch,'" not supported'
      endif
      return
      end
 
      subroutine archwr_ufm
c
c=======================================================================
c  subroutine to write an unformatted file 
c=======================================================================
c
#include "def_master.h"
#include "param.h"
#include "archive.h"
#include "iounit.h"
#include "levind.h"
#include "ctmngr.h"
#include "mesdta.h"
      character*40 rsfile
      logical around
c
c------------------------------------------------------------------
c  create archive file and write base data
c------------------------------------------------------------------
c
      if(idvar.eq.0) then
        write (rsfile,'(i10.10)')nint(totday)
        rsfile = 'd'//rsfile(7:10)//'.ufm'
        inquire(file=rsfile,exist=around)
        if(around) then
          write (rsfile,'(i10.10)') itt
          rsfile = 'd'//rsfile(7:10)//'.restart2'
        endif
c
        open(iorest, file=rsfile, form='UNFORMATTED',
     &       access='SEQUENTIAL')
c
        write(stdout,24) itt, iorest, rsfile
        write(iorest) itt, totsec, totday, years
        write(iorest) kmt
        return
      endif
c 
c------------------------------------------------------------------
c 2-D field 
c------------------------------------------------------------------
c
      if(idvar.le.NUM2D) then
        write (iorest) ((rest2d(i,j),i=1,IMT_M),j=1,JMT_M)
        larchq = .false.
      else
c
c------------------------------------------------------------------
c 3-D field 
c------------------------------------------------------------------
c
       mzon   = (jarchu-1)/JSUB_M + 1
       jupper = JSUB_M
       if(mzon*JSUB_M.gt.JMT_M) jupper = JMT_M - (mzon-1)*JSUB_M
c
       do 10 j=1,jupper
        write(iorest) ((rest3d(k,i,j),k=1,KM),i=1,IMT_M)
  10   continue
        larchq = .false.
c
c------------------------------------------------------------------
c at end close file 
c------------------------------------------------------------------
c
        if(jarchu.ge.JMT_M)then
          close(unit=iorest)
        endif 
      endif
      return
  24  format('  Writing restart file for itt = ',i7,/,
     &       '  unit = ',i2,'  filename = ',a)
      end
        
