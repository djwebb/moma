      subroutine snapreq
c
c @(#) SCCS module: m_snap.F  version: 1.5
c     Creation date: 03/24/95
c
c==================================================================
c subroutine to prepare for the next snapshot request
c==================================================================
c
#include "def_master.h"
#include "param.h"
#include "snaps.h"
#include "mesdta.h"
#include "iounit.h"
#include "varinfo.h"
c
      real vmask(4)
      data vmask/1.e7,2.e7,3.e7,4.e7/
c
      idsnp = idsnp + 1
c
      if(idsnp.gt.NSNAPS) then
        write(stdout,*) ' Snapshot complete'
        lsnapp=.false.
        return
      endif
c
c  initialise snap2d array
c
      if(ftsnap.eq.'asc')then
        call initsnp(vmask)
      else
c
c Initialise temperature to -2. everywhere
c
        if(ispvar(idsnp).eq.4) then
          do 300 jc=1,JMT_M
          do 300 ic=1,IMT_M
            snap2d(ic,jc)= -2.0
 300      continue
c
c Initialise salinity to 0.01 everywhere
c
        elseif(ispvar(idsnp).eq.5) then
          do 400 jc=1,JMT_M
          do 400 ic=1,IMT_M
            snap2d(ic,jc)= 0.01 
 400      continue
c 
c Initialise other fields to zero everywhere
c 
        else
          do 500 jc=1,JMT_M 
          do 500 ic=1,IMT_M 
            snap2d(ic,jc)= c0
500       continue
        endif
      endif
c
      lsnapq = .true.
      nssnap = nproc
      call m_send(MSG_RQ_SNAP)
      return
      end

      subroutine initsnp(vmsk)
c
c==================================================================
c subroutine to initialise the snapshot array to produce
c the correct land masking for routine ascout0.
c==================================================================
c
c
#include "def_master.h"
#include "param.h"
#include "snaps.h"
#include "levind.h"
#include "varinfo.h"
#include "mesdta.h"
      real vmsk(4)
c
      idsvar=ispvar(idsnp)
      do 10 j=1,JMT_M
      do 10 i=1,IMT_M
       snap2d(i,j)=vmsk(2)
       if(vartgrd(idsvar)) then
        if(kmt(i,j).eq.0) snap2d(i,j)=vmsk(1)
       else
        if(kmu(i,j).eq.0) snap2d(i,j)=vmsk(1)
       endif
 10   continue
      return
      end

      subroutine snapwr
c
c @(#) SCCS module: m_snap.F  version: 1.5
c     Creation date: 03/24/95
c
c------------------------------------------------------------------
c  subroutine to write out the next snapshot field
c------------------------------------------------------------------
c
#include "iounit.h"
c
      if(ftsnap.eq.'asc')then
        call snap_asc
#ifdef ncsa_hdf
      elseif(ftsnap.eq.'hdf')then
        call snap_hdf
#endif
      else
        write(stdout,*) 'Snapshot write: file type "',
     &                  ftsnap,'" not supported'
      endif
      return
      end
 
      subroutine snap_asc
#include "def_master.h"
#include "param.h"
#include "snaps.h"
#include "iounit.h"
#include "levind.h"
#include "ctmngr.h"
#include "varinfo.h"
#include "mesdta.h"
      character*40 spfile
      logical*4 around
      character*2 opform,namrun*50
c Load some strings used by the header routine
      data namrun/'FF Occam pvm model'/opform/'CD'/
c
c+->
c|
      if(idsnp.eq.0) then
c
c------------------------------------------------------------------
c  create snapshot file name
c------------------------------------------------------------------
c
        write (spfile,'(i10.10)')nint(totday)
        spfile = 'd'//spfile(7:10)//'.snap'
        inquire(file=spfile,exist=around)
        if(around) then
         write (spfile,'(i10.10)') itt
         spfile = 'd'//spfile(7:10)//'.snap2'
        endif
c
        write(stdout,24) itt, iosnap, spfile
        open(iosnap, file=spfile)
c
c Write out the header
c
        call header(iosnap,'HUVTSUVTSUV','DEPTH',1,IMT_M,1,
     +                  JMT_M,1,1,opform,namrun)
        return
c|
c+->
      endif
c
c write out 2-D field 
c
      call ascout0(snap2d,IMT_M,IMT_M,JMT_M,vmask,2,iosnap)
c
      if(idsnp.eq.NSNAPS) then
        close(unit=iosnap)
      endif
c
c clear the  "partial snapshot out-standing" flag
c
      lsnapq = .false.
c
      return
  24  format('  Writing snapshot file for itt = ',i7,/,
     &       '  unit = ',i2,'  filename = ',a)
      end
        

      subroutine header(op,trac,depvar,ibase1,itop1,ibase2,itop2,tpts,
     &                  intlv,opform,namrun)
c
c Subroutine to produce a header for ascout files.
c tpts - 1 if t-grid 0 if u-grid
c
#include "param.h"
#include "scalar.h"
#include "coord.h"
#include "ctmngr.h"
      character trac*(*),depvar*9,opform*(*),namrun*(*)
      character*9 quan(3),from(3),incr(3),to(3)
      integer nop(3),op,tpts
c
c convert input strings to upper-case
      call conv2up(trac)
      call conv2up(depvar)
      call conv2up(opform)
      call conv2up(namrun)
      rintlv=float(intlv)
      xincr=dxdeg
      yincr=dydeg
c
          quan(1) = 'LONGITUDE'
          quan(2) = ' LATITUDE'
          write(from(1),'(f9.3)')
     +                   (ibase1-1)*xincr+(2-tpts)*0.5*xincr+stlon
          write(from(2),'(f9.3)')
     +                   (ibase2-1)*yincr+(2-tpts)*0.5*yincr+stlat
          write(incr(1),'(f9.3)')   xincr
          write(incr(2),'(f9.3)')   yincr
          write(to(1),'(f9.3)')
     +                    (itop1-1)*xincr+(2-tpts)*0.5*xincr+stlon
          write(to(2),'(f9.3)')
     +                    (itop2-1)*yincr+(2-tpts)*0.5*yincr+stlat
          nop(1) = itop1-ibase1+1
          nop(2) = itop2-ibase2+1
c
c
      quan(3)=' TIMESTEP'
      ndinc=0
      ndfir1 = itt
      ndlas = itt
      write(from(3),'(i9)')ndfir1
      write(incr(3),'(i9)')ndinc
      write(to(3),'(i9)')ndlas
      nop(3) = 1
c
        write(op,5100) trac,depvar,opform
c
      write(op,5102) namrun(1:2),namrun
      write(op,5103) (i,i=1,3)
      write(op,5104) (quan(i),i=1,3)
      write(op,5105) (from(i),i=1,3)
      write(op,5106) (incr(i),i=1,3)
      write(op,5107) (to(i),i=1,3)
      write(op,5108) (nop(i),i=1,3)
      write(op,5110) rintlv
      write(op,5111) ttsec ,dtts
 5100 format('VARIABLE :',a15,2x,a9,t41,'FORMAT   :',a2)
 5101 format('VARIABLE :',a15,t41,'FORMAT   :',a2)
 5102 format('MODEL : ',a2,'      COMMENTS:',a55)
 5103 format('INDEX ',9x,':',3('   ',i1,'     :'))
 5104 format('QUANTITY ',6x,':',a9,':',a9,':',a9,':')
 5105 format('FROM     ',6x,':',a9,':',a9,':',a9,':')
 5106 format('INCREMENT',6x,':',a9,':',a9,':',a9,':')
 5107 format('TO       ',6x,':',a9,':',a9,':',a9,':')
 5108 format('NO.OF POINTS ',2x,':',i9,':',i9,':',i9,':')
 5109 format(a50)
 5110 format(f7.3)
 5111 format('FIRST TTSEC ',f12.0,' dtts ',f5.0)
      return
      end
      subroutine conv2up(mixcase)
      character mixcase*(*)
      integer upa
c
#ifdef cray-t3d
      return
#else
      nchar=lnblnk(mixcase)
c
      lowa=ichar('a')
      lowz=ichar('z')
      upa =ichar('A')
c
      do 10 n=1,nchar
      nc=ichar(mixcase(n:n))
      if(nc.ge.lowa.and.nc.le.lowz) then
      nc=upa+(nc-lowa)
      mixcase(n:n) = char(nc)
      endif
 10   continue
c
#endif
      end

      subroutine ascout0 (array,idim,id,jd,vmask,nchar,nout)            
c
c     version with single precision real arguments for 64-bit machine   
c                                                                       
c     subroutine to encode a section of an array as sets of 'nchar'     
c     printable characters, and write as a formatted card-image dump.   
c     (uses ascii characters 0-9 , a-z , lower case a-z and brackets)   
c                                                                       
c     array - 2-d array of values to be converted                       
c     idim  - declared i-dimension of array in calling program          
c     id,jd - specify section of array to be converted                  
c     vmask - 4-element array whose values indicate 'masked' points.    
c             such points are denoted by one of the four possible       
c             combinations of full stop and comma, padded out to nchar  
c             characters by repetition of the last character of the pair
c             these values are ignored in finding max and mins for scali
c             the vmask values are normally much larger than other value
c     nchar - number of characters to be used to represent an array valu
c     nout  - fortran channel number of output dataset.                 
c                                                                       
c     m. a. rowe  sept. 1987 ( rewritten j. r. blundell  07/07/1988 )   
c     this version (internally declared character array) 14/12/1988     
c     modified to allow for four types of masked point   07/02/1989     
c     n.b. internal arithmetic always done in 64-bit mode               
c                                                                       
c     internal parameters:                                              
c                                                                       
c     lrecl - max. length of data record to be output                   
c     nascc - number of different ascii characters used in              
c             representation of numbers (at unmasked points)            
c     ncmax - max. number of characters which can be used               
c             to represent an array element                             
c                                                                       
      integer lrecl,nascc,ncmax                                         
c                                                                       
      parameter ( lrecl=80, nascc=64, ncmax=5 )                         
c                                                                       
c     local variables                                                   
c                                                                       
      integer icode(ncmax),idim,id,jd,nchar,nout,                       
     *        i,j,nnum,ic,integ,ncbuff,linlen,mtype                     
      real array(idim,jd),vmask(4)                                      
      real fmin,fmax,range,arang,scale                                  
      character*1 ascarr(lrecl),lkup(nascc),cmask(2),mask(ncmax,4)      
      character*(nascc) charep                                          
c                                                                       
      equivalence (charep(1:1),lkup(1))                                 
c                                                                       
c     specify the nascc characters to be used in the number             
c     representation, and the characters denoting masked points         
c                                                                       
      charep( 1:10) = '0123456789'                                      
      CHAREP(11:36) = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'                      
      CHAREP(37:62) = 'abcdefghijklmnopqrstuvwxyz'                      
      CHAREP(63:64) = '()'                                              
      CMASK(1)='.'                                                      
      CMASK(2)=','                                                      
c                                                                       
c     write out coding info in first data record                        
c     (write warning to unit 6 if it won't fit)                         
c                                                                       
      if ( nascc.gt.72 ) write(6,50) nascc                              
   50 format(/,2x,'**ascout warning: overlength coding record,',        
     *       ' nascc =',i3)                                             
      write(nout,'(i4,1x,2a1,1x,72a1)') nascc,cmask,(lkup(i),i=1,nascc) 
c                                                                       
c     check input value of nchar                                        
c                                                                       
      if ( nchar.lt.2 .or. nchar.gt.ncmax ) then                        
        write(6,100) nchar                                              
  100   format(/,2x,'**ascout warning: routine called',                 
     *         ' with incorrect nchar =',i4)                            
        return                                                          
      end if                                                            
c                                                                       
c     check input values of vmask are all different,                    
c     otherwise masking will be ambiguous                               
c                                                                       
      do 110 j=1,3                                                      
        do 110 i=j+1,4                                                  
          if ( vmask(i).eq.vmask(j) ) write(6,120) i,j                  
  110 continue                                                          
  120 format(/,2x,'**ascout warning: vmask(',i1,') = vmask(',i1,')',    
     *       /,2x,'**masking produced will be ambiguous')               
c                                                                       
c     create the 4 types of mask, including padding characters          
c                                                                       
      do 130 ic=1,nchar                                                 
        mask(ic,1) = cmask(1)                                           
        mask(ic,2) = cmask(1)                                           
        mask(ic,3) = cmask(2)                                           
        mask(ic,4) = cmask(2)                                           
  130 continue                                                          
      mask(2,2) = cmask(2)                                              
      mask(2,3) = cmask(1)                                              
c                                                                       
c     establish range of data and scaling for conversion                
c     ( typical size of values assumed o(10**5) )                       
c                                                                       
      fmax = -9999999.9                                                 
      fmin =  9999999.9                                                 
      do 150 i=1,id                                                     
        do 150 j=1,jd                                                   
          do 140 mtype=1,4                                              
            if ( array(i,j).eq.vmask(mtype) ) goto 145                  
  140     continue                                                      
            fmin = min( fmin,array(i,j) )                               
            fmax = max( fmax,array(i,j) )                               
  145     continue                                                      
  150 continue                                                          
      if ( fmax.lt.-99999.9 .or. fmin.gt.99999.9 )                      
     *  write(6,200) fmin,fmax                                          
  200 format(/,2x,'**ascout warning: large +ve minimum or large',       
     *       ' -ve maximum value',/,2x,'fmax, fmin = ',1p,2e16.5)       
c                                                                       
      nnum = id*jd                                                      
      write(nout,'(1p,2e20.12,4i10)') fmin,fmax,id,jd,nnum,nchar        
c                                                                       
c     calculate scaling factors                                         
c                                                                       
      arang = real( nascc**nchar - 1 )                                  
      range = fmax - fmin                                               
      size  = max(abs(fmax),abs(fmin))                                  
      if (range .gt. 1e-10)then                                         
       scale = arang/range                                              
      else                                                              
        scale = 0.0                                                     
      endif                                                             
c                                                                       
c     exit if range too small                                           
c                                                                       
      if(scale.eq.0.0)then                                              
        temp = 0.0                                                      
      else                                                              
        temp = abs(range/scale)                                         
      endif                                                             
      if(temp.lt.1e-10)then                                             
        write(6,250)fmax,fmin                                           
  250   format('**ascout warning: field approx. constant,',             
     *       ' not character coded',/,'fmax, fmin = ',1p,2e16.5)        
c                                                                       
c     scale array and encode as nchar printable characters              
c                                                                       
      else                                                              
      ncbuff = 0                                                        
      if ( nchar.eq.3 ) linlen=78                                       
      if ( nchar.ne.3 ) linlen=80                                       
      do 500 j=1,jd                                                     
        do 500 i=1,id                                                   
c                                                                       
          do 350 mtype = 1,4                                            
            if ( array(i,j).eq.vmask(mtype) ) then                      
c             type mtype masked point; copy from mask(ncmax,mtype)      
              do 300 ic = 1,nchar                                       
                ascarr(ncbuff+ic) = mask(ic,mtype)                      
  300         continue                                                  
              goto 450                                                  
            endif                                                       
  350     continue                                                      
c         normal point; encode as nchar characters                      
          integ = nint( (array(i,j)-fmin)*scale )                       
          do 400 ic=nchar,1,-1                                          
            icode(ic) = 1 + mod( integ, nascc )                         
            ascarr(ncbuff+ic) = lkup( icode(ic) )                       
            integ = integ/nascc                                         
  400     continue                                                      
  450     continue                                                      
          ncbuff = ncbuff + nchar                                       
c                                                                       
          if ( ncbuff.eq.linlen ) then                                  
c           buffer ascarr full; write to                                
c           channel nout (card-image format)                            
            if ( nchar.ne.3 ) then                                      
              write(nout,'(80a1)') (ascarr(ic),ic=1,ncbuff)             
             else                                                       
              write(nout,'(1x,78a1,1x)') (ascarr(ic),ic=1,ncbuff)       
            end if                                                      
            ncbuff = 0                                                  
          end if                                                        
c                                                                       
  500 continue                                                          
c                                                                       
c     flush character buffer if not empty                               
c                                                                       
      if ( ncbuff.ne.0 ) then                                           
        if ( nchar.ne.3 ) then                                          
          write(nout,'(80a1)') (ascarr(ic),ic=1,ncbuff)                 
         else                                                           
          write(nout,'(1x,78a1,1x)') (ascarr(ic),ic=1,ncbuff)           
        end if                                                          
        ncbuff = 0                                                      
      end if                                                            
c                                                                       
        end if                                                          
c                                                                       
      return                                                            
      end                                                               

