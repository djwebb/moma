#ifdef ncsa_hdf
      subroutine archrd_hdf(ncntl)
c
c @(#) SCCS module: hdf.F  version: 1.4
c     Creation date: 03/23/95
c
c=======================================================================
c  subroutine to read hdf version of the restart file
c  written for HDF3.3r3 routines
c=======================================================================
c
#include "def_master.h"
#include "param.h"
#include "archive.h"
#include "coord.h"
#include "ctmngr.h"
#include "iounit.h"
#include "levind.h"
#include "mesdta.h"
c
c  hdf '.h' files
c
#include "constant.i"
#include "hdf.h"
c
c  local variables
c
      logical*4 around
      integer*4 file_id,sds_id,ret,
     &          start(MAXDIMS),stride(MAXDIMS),edge(MAXDIMS)
c
c=======================================================================
c  open archive file when ncntl zero (entry from routine init_kmt).
c  and read basic information
c=======================================================================
c
      if(ncntl.eq.0)then
        if(fnrest(1:1).eq.' ')then
          write(stdout,*)' archrd_hdf: restart filename not given'
          write(stdout,*)' Program stopping ...'
          stop
        endif
        inquire(file=fnrest,exist=around)
        if(.not.around)then
          write(stdout,*)' archrd_hdf: restart filename not present'
          write(stdout,*)' Filename = ',fnrest
          write(stdout,*)' Program stopping ...'
          stop
        endif
        file_id = sfstart(fnrest, DFACC_RDONLY)
          call hdf_error( file_id,'cannot open restart dataset', 0)
c
c  read dataset information
c
        start(1) =0
        start(2) =0
        start(3) =0
        stride(1)=1
        stride(2)=1
        stride(3)=1
        edge(1)  =1
        edge(2)  =1
        edge(3)  =1
c
        index = sfn2index(file_id,tlabel(1))
        call hdf_error(index,'cannot find id for TIMESTEP', 0)
        sds_id = sfselect(file_id,index)
        call hdf_error(sds_id,'error getting index for TIMESTEP', 0)
        ret = sfrdata(sds_id,start, stride,edge,itt)
        call hdf_error(ret,'error reading TIMESTEP', 0)
        ret = sfendacc(sds_id)
        call hdf_error( ret, 'cannot close TIMESTEP dataset', 0)
c
        index = sfn2index(file_id,tlabel(2))
        call hdf_error( index,'cannot find id for totsec', 0)
        sds_id = sfselect(file_id,index)
        call hdf_error(sds_id,'error getting index for totsec', 0)
        ret = sfrdata(sds_id,start, stride,edge,totsec)
        call hdf_error(ret,'error reading totsec', 0)
        ret = sfendacc(sds_id)
        call hdf_error( ret, 'cannot close totsec dataset', 0)
c
        index = sfn2index(file_id,tlabel(3))
        call hdf_error( index,'cannot find id for totday', 0)
        sds_id = sfselect(file_id,index)
        call hdf_error(sds_id,'error getting index for totday', 0)
        ret = sfrdata(sds_id,start, stride,edge,totday)
        call hdf_error(ret,'error reading totday', 0)
        ret = sfendacc(sds_id)
        call hdf_error( ret, 'cannot close totday dataset', 0)
c
        index = sfn2index(file_id,tlabel(4))
        call hdf_error( index,'cannot find id for years', 0)
        sds_id = sfselect(file_id,index)
        call hdf_error(sds_id,'error getting index for years', 0)
        ret = sfrdata(sds_id,start, stride,edge,years)
        call hdf_error(ret,'error reading years', 0)
        ret = sfendacc(sds_id)
        call hdf_error( ret, 'cannot close years dataset', 0)
c
c  read kmt
c
        edge(1) = IMT_M
        edge(2) = JMT_M
        index = sfn2index(file_id,'KMT')
        call hdf_error( index,'cannot find id for KMT', 0)
        sds_id = sfselect(file_id,index)
        call hdf_error(sds_id,'error getting index for KMT', 0)
        ret = sfrdata(sds_id,start, stride,edge,kmt)
        call hdf_error(ret,'error reading KMT', 0)
        ret = sfendacc(sds_id)
        call hdf_error( ret, 'cannot close KMT dataset', 0)
        return
      endif
c
c=======================================================================
c  read the 2-D datasets
c=======================================================================
c
      do 50 n=1,nproc
        iarcha(n)=itida(n)
        narcha(n)=n
  50  continue
      jarchl = 1
      jarchu = JMT_M
      nsarch = nproc
c
      start(1) =0
      start(2) =0
      start(3) =0
      stride(1)=1
      stride(2)=1
      stride(3)=1
      edge(1)  =IMT_M
      edge(2)  =JMT_M
c
      do 100 idvar = 1, NUM2D
c
        index = sfn2index(file_id,varid(idvar))
        call hdf_error( index,'cannot find id for 2-D dataset', 0)
        sds_id = sfselect(file_id,index)
        call hdf_error(sds_id,'error getting index for 2-D data', 0)
        ret = sfrdata(sds_id,start,stride,edge,rest2d)
        call hdf_error(ret,'error reading 2-D data', 0)
        ret = sfendacc(sds_id)
        call hdf_error( ret, 'cannot close 2-D dataset', 0)
c
        call m_send(MSG_RESTART)
        call mchkarch(2)
 100  continue
c
c=======================================================================
c  read the 3-D datasets
c=======================================================================
c
      edge(1) = KM
      edge(2) = IMT_M
c
      nsub = 1+(JMT_M-1)/JSUB_M
      do 400 idvar = NUM2D+1,NARCHV
        do 300 jj = 1,nsub
          jarchl = 1+(jj-1)*JSUB_M
          jarchu = min(jarchl-1+JSUB_M, JMT_M)
          nsarch = 0
          do 200 n=1,nproc
            if(jarchu.gt.jt_s(n).and.jarchl.le.jt_e(n))then
              nsarch=nsarch+1
              iarcha(nsarch) = itida(n)
              narcha(nsarch) = n
            endif
 200      continue
c
          start(3) = jarchl-1
          edge(3)  = jarchu-jarchl+1
c      
          index = sfn2index(file_id,varid(idvar))
          call hdf_error( index,'cannot find id for 3-D dataset', 0)
          sds_id = sfselect(file_id,index)
          call hdf_error(sds_id,'error getting index for 3-D data', 0)
          ret = sfrdata(sds_id,start,stride,edge,rest3d)
          call hdf_error(ret,'error reading 3-D data', 0)
          ret = sfendacc(sds_id)
          call hdf_error( ret, 'cannot close 3-D dataset', 0)
c
          call m_send(MSG_RESTART)
          call mchkarch(2)
 300    continue
 400  continue
c
c   tell slaves to continue past checkpoint 
c   and close restart file
c
      call m_send(MSG_CONT)
      ret = sfend(file_id)
      call hdf_error( ret,' error closing SD interface', 0)
      write(stdout,*) ' Restart complete'
      return
      end


      subroutine archwr_hdf
c
c=======================================================================c  subroutine to save hdf version of the restart file
c  The generated file name starts with 'd' and ends with '.hdf'.
c=======================================================================c
#include "def_master.h"
#include "param.h"
#include "archive.h"
#include "coord.h"
#include "ctmngr.h"
#include "iounit.h"
#include "levind.h"
#include "mesdta.h"
#include "switch.h"
#include "timelv.h"
#include "varinfo.h"
c
c  hdf '.h' files
c
# include "constant.i"
# include "hdf.h"
c
c  local variables
c
      logical*4 around
      integer*4 write_id,sdsid,dimid1,dimid2,dimid3,ret,ret2,
     &          srcdims(MAXDIMS),stride(MAXDIMS),edge(MAXDIMS),
     &          irank
c
c  initialise the hdf common blocks
c
      data lonlabel(1), lonlabel(2) /'LONGITUDE_T', 'LONGITUDE_U' /
      data latlabel(1), latlabel(2) /'LATITUDE_T',  'LATITUDE_U'  /
      data labels2(1), units2(1), fmats2(1) 
     &     /'Longitude', 'degrees', 'F10.4'/
      data labels2(2), units2(2), fmats2(2) 
     &     /'Latitude' , 'degrees', 'F10.4'/
      data labels3(1), units3(1), fmats3(1) 
     &     /'Depth'    , 'cm'    , 'F5.2' /
      data labels3(2), units3(2), fmats3(2) 
     &     /'Longitude', 'degrees', 'F10.4'/
      data labels3(3), units3(3), fmats3(3) 
     &     /'Latitude' , 'degrees', 'F10.4'/
      data (dims2(k),k=1,2) / IMT_M,JMT_M /
      data (dims3(k),k=1,3) / KM,IMT_M,JMT_M /
      data (cfext(k),k=1,2) /'.hdf','.hdfa'/
      data (tlabel(n),n=1,4) /'TIMESTEP', 'SECONDS', 'DAY', 'YEARS'/
      data itsize,itdims(1) /1,1/
c      data (varname(n),varunit(n),vartgrd(n),n=1,NARCHV)/
c     &  'free surface height      ', 'cm          ', .true.  ,
c     &  'barotropic u-velocity    ', 'cm/s        ', .false. ,
c     &  'barotropic v-velocity    ', 'cm/s        ', .false. ,
c     &  'potential temperature    ', 'C           ', .true.  ,
c     &  'salinity                 ', '.001(NSU-35)', .true.  ,
c     &  'u-velocity               ', 'cm/s        ', .false. ,
c     &  'v-velocity               ', 'cm/s        ', .false. /
        data (varid(n), varcordsys(n), n=1, NARCHV)/
     &  'FREE SURFACE HEIGHT',          'T_GRID'        ,
     &  'BAROTROPIC U VELOCITY',        'U_GRID'        ,
     &  'BAROTROPIC V VELOCITY',        'U_GRID'        ,
     &  'POTENTIAL TEMPERATURE',        'T_GRID'        ,
     &  'SALINITY',                     'T_GRID'        ,
     &  'U VELOCITY',                   'U_GRID'        ,
     &  'V VELOCITY',                   'U_GRID'        /
c
c=======================================================================c open file and write fixed data
c=======================================================================c
      if(idvar.eq.0)then
c
        write(hdffile, '(i10.10)') nint(totday)
        hdffile = 'd'//hdffile(7:10)//cfext(1)
        inquire(file = hdffile, exist = around )
        if(around) then
          write( hdffile,'(i10.10)') itt
          hdffile = 'd'// hdffile(7:10)//cfext(2)
        endif
        
        write_id = sfstart(hdffile, DFACC_CREATE)
        if(write_id.eq.-1) then
          write(stdout,*),' error opening the file'
          return
        endif
        write (stdout,24) itt, hdffile
  24  format('  Writing hdf archive file for itt = ',i7,/,
     &       '  hdf filename = ',a) 
c
c write the dataset containing 'itt'
c
        irank = 1
        sdsid = sfcreate( write_id,tlabel(1),
     &                  INTEGER_HDF,irank,itdims )
        call hdf_error( sdsid,'cannot create itt dataset', 0)
        ret = sfsdtstr( sdsid, 'TIMESTEP', '', '','')
        call hdf_error( ret, 'cannot set datastrings', 0)
        ret = sfwdata( sdsid, 0,1,itdims, itt) 
        call hdf_error( ret, 'cannot write dataset', 0)
        ret = sfendacc(sdsid)
        call hdf_error( ret, 'cannot close dataset', 0)
c
c write the dataset containing 'totsec'
c
        sdsid = sfcreate(write_id,tlabel(2),
     &                 REAL_HDF,itsize,itdims)
        call hdf_error( sdsid,'cannot create totsec dataset', 0)
        ret = sfsdtstr( sdsid, 'Model Seconds', '', '')
        call hdf_error( ret, 'cannot set datastrings', 0)
        ret = sfwdata( sdsid, 0, 1, itdims, totsec)
        call hdf_error( ret, 'cannot write dataset', 0)
        ret = sfendacc(sdsid)
        call hdf_error( ret, 'cannot close dataset', 0)
c
c write the dataset containing 'totday'
c
        sdsid = sfcreate(write_id,tlabel(3),
     &                 REAL_HDF,itsize,itdims)
        call hdf_error( sdsid,'cannot create totday dataset', 0)
        ret = sfsdtstr( sdsid, 'Model Day', '', '')
        call hdf_error( ret, 'cannot set datastrings', 0)
        ret = sfwdata(  sdsid, 0, 1,itdims, totday)
        call hdf_error( ret, 'cannot write dataset', 0)
        ret = sfendacc(sdsid)
        call hdf_error( ret, 'cannot close dataset', 0)
c
c write the dataset containing 'years'
c
        sdsid = sfcreate(write_id,tlabel(4),
     &                 REAL_HDF,itsize,itdims)
        call hdf_error( sdsid,'cannot create years dataset', 0)
        ret = sfsdtstr( sdsid, 'Model Year', '', '')
        call hdf_error( ret, 'cannot set datastrings', 0)
        ret = sfwdata(  sdsid,0,1,itdims,years)
        call hdf_error( ret, 'cannot write dataset', 0)
        ret = sfendacc(sdsid)
        call hdf_error( ret, 'cannot close dataset', 0)
c
c initialise the KMT array field.
c
        irank = 2
        srcdims(1)= 0
        srcdims(2)= 0
        stride(1) = 1
        stride(2) = 1
c
c call the HDF function to create the dataset.
c
        sdsid = sfcreate( write_id,'KMT',
     &                  INTEGER_HDF, irank, dims2) 
        call hdf_error( sdsid,'cannot create kmt dataset', 0)
        ret = sfsdtstr( sdsid, 'KMT', '', '')
        call hdf_error( ret, ' cannot set datastrings', 0)
        ifill_value = 999
        ret = sfsfill(sdsid,ifill_value)
        call hdf_error( ret, ' cannot set fill value', 0)
c
        dimid1 = sfdimid(sdsid,0)
        call hdf_error( dimid1, ' error getting dimid', 0)
        dimid2 = sfdimid(sdsid,1)
        call hdf_error( dimid2, ' error getting dimid', 0)
        ret = sfsdmname(dimid1,lonlabel(1))
        call hdf_error( ret, ' cannot set dimname', 0)
        ret = sfsdmname(dimid2,latlabel(1))
        call hdf_error( ret, ' cannot set dimname', 0)
        ret = sfsdmstr(dimid1,labels2(1),units2(1),fmats2(1))
        call hdf_error( ret, ' setting dimension string', 0)
        ret = sfsdmstr(dimid2,labels2(2),units2(2),fmats2(2))
        call hdf_error( ret, ' setting dimension string', 0)
        ret = sfsdscale(dimid1, dims2(1),REAL_HDF, xt)
        call hdf_error( ret, 'cannot set scales', dimid1)
        ret = sfsdscale(dimid2, dims2(2),REAL_HDF, yt)
        call hdf_error( ret, 'cannot set scales', dimid2)
        ret = sfwdata( sdsid,srcdims(1),stride(1),dims2(1),kmt)
        call hdf_error( ret,'cannot write kmt dataset',0)
        ret = sfendacc(sdsid)
        call hdf_error( ret, 'cannot close dataset', 0)
c
c initialise the KMU array field.
c
        sdsid = sfcreate( write_id,'KMU',
     &                  INTEGER_HDF,irank,dims2)
        call hdf_error( sdsid,'cannot create kmu dataset', 0)
        ret = sfsdtstr( sdsid, 'KMU', '', '')
        call hdf_error( ret, 'cannot set datastrings', 0)
        ifill_value = 999
        ret = sfsfill(sdsid,ifill_value)
        call hdf_error( ret, 'cannot set fill value', 0)
c
        dimid1 = sfdimid(sdsid,0)
        call hdf_error( dimid1, 'error getting dimid', 0)
        dimid2 = sfdimid(sdsid,1)
        call hdf_error( dimid2, 'error getting dimid', 0)
        ret = sfsdmname(dimid1,lonlabel(2))
        call hdf_error( ret, 'cannot set dimname', 0)
        ret = sfsdmname(dimid2,latlabel(2))
        call hdf_error( ret, 'cannot set dimname', 0)
c
        ret = sfsdmstr(dimid1,labels2(1),units2(1),fmats2(1))
        call hdf_error( ret, 'setting dimension string', 0)
        ret = sfsdmstr(dimid2,labels2(2),units2(2),fmats2(2))
        call hdf_error( ret, 'setting dimension string', 0)
        ret = sfsdscale(dimid1, dims2(1),REAL_HDF,xu)
        call hdf_error( ret,'cannot set scales', dimid1)
        ret = sfsdscale(dimid2, dims2(2),REAL_HDF,yu)
        call hdf_error( ret,'cannot set scales', dimid2)
        ret = sfwdata( sdsid,srcdims(1),stride(1),dims2(1),kmu)
        call hdf_error( ret,'cannot write kmu dataset',0)
        ret = sfendacc(sdsid)
        call hdf_error( ret, ' cannot close dataset', 0)
        return
c
c=======================================================================
c write two dimensional datasets to file
c=======================================================================c
      elseif(idvar.le.NUM2D)then
        irank = 2
c
        sdsid = sfcreate( write_id,varid(idvar),
     &                    REAL_HDF,irank,dims2)
        call hdf_error( sdsid, 'cannot create dataset', 0)
        ret = sfsdtstr(sdsid,varname(idvar),varunit(idvar),'','')
        call hdf_error( ret, 'cannot set data strings',0)
        fill_value = 0.0
        ret = sfsfill(sdsid,fill_value)
         call hdf_error( ret, 'cannot set fill value', 0)
c
c set the scales
c
        dimid1 = sfdimid( sdsid,0)
        call hdf_error( dimid1, 'error getting dimid', 0)
        dimid2 = sfdimid( sdsid,1)
        call hdf_error( dimid2, 'error getting dimid', 0)
        if(vartgrd(idvar))then
          ret = sfsdmname(dimid1, lonlabel(1))
          call hdf_error( ret, 'cannot set dimname', 0)
          ret = sfsdmname(dimid2, latlabel(1))
          call hdf_error( ret, 'cannot set dimname', 0)
          ret = sfsdscale(dimid1,IMT_M,REAL_HDF,xt)
          ret2= sfsdscale(dimid2,JMT_M,REAL_HDF,yt)
        else
          ret = sfsdmname(dimid1, lonlabel(2))
          call hdf_error( ret, 'cannot set dimname', 0)
          ret = sfsdmname(dimid2, latlabel(2))
          call hdf_error( ret, 'cannot set dimname', 0)
          ret = sfsdscale(dimid1,IMT_M,REAL_HDF,xu)
          ret2= sfsdscale(dimid2,JMT_M,REAL_HDF,yu) 
        endif
        call hdf_error( ret, 'cannot set x-scales', dimid1)
        call hdf_error( ret2, 'cannot set y-scales',dimid2)
c
        ret = sfsdmstr(dimid1,labels2(1),units2(1),fmats2(1))
        call hdf_error( ret, 'setting dimension string', 0)
        ret = sfsdmstr(dimid2,labels2(2),units2(2),fmats2(2))
        call hdf_error( ret, 'setting dimension string', 0)
c
c  write out the dataset
c
        srcdims(1)= 0
        srcdims(2)= 0
        stride(1) = 1
        stride(2) = 1
        ret = sfwdata(sdsid,srcdims(1),stride(1),dims2(1),rest2d)
        call hdf_error( ret, 'writing dataset ', 0)
        ret = sfsattr(sdsid, LEVEL_ATTR, INTEGER_HDF, 1, 0)
        call hdf_error( ret, 'setting attribute', 0)
        ret = sfsattr(sdsid, TGRID_ATTR, INTEGER_HDF, 1,
     &                                      vartgrd(idvar))
        call hdf_error( ret, 'setting attribute', 0)
        ret = sfendacc(sdsid)
        call hdf_error( ret, 'cannot close dataset', 0)
c
        larchq = .false.
        return
      endif
c
c=======================================================================
c archive the 3-D datasets
c=======================================================================c
c
      if(idvar.gt.NARCHV)then
        write(stdout,*)' archwr_hdf: routine called with idvar .gt.',
     &                 ' NARCHV '
        write(stdout,*)' program stopping ...'
        call m_send(MSG_ABORT)
      endif
      irank = 3
c
c  initialise dataset (if jarchl = 1)
c
      if(jarchl.eq.1)then
        sdsid = sfcreate(write_id,varid(idvar),
     &                   REAL_HDF,irank,dims3)
        call hdf_error( sdsid, 'cannot create dataset', 0)
        ret = sfsdtstr(sdsid,varname(idvar),varunit(idvar),'','')
        call hdf_error( ret, 'cannot set data strings', 0)
        fill_value = 0.0
        ret = sfsfill(sdsid,fill_value)
        call hdf_error( ret, 'cannot set fill value', 0)
c
        dimid1 = sfdimid(sdsid,0)
        call hdf_error( dimid1, 'Cannot get dimid:',dimid1) 
        dimid2 = sfdimid(sdsid,1)
        call hdf_error( dimid2, 'Cannot get dimid:',dimid2) 
        dimid3 = sfdimid(sdsid,2)
        call hdf_error( dimid3, 'Cannot get dimid:',dimid3) 
c
c set the scales for each dimensions
c
        ret = sfsdmname(dimid1, 'DEPTH')
        call hdf_error(ret, 'cannot set DEPTH dimname', 0)
        ret = sfsdscale(dimid1,KM,REAL_HDF,zt)
        call hdf_error(ret, 'cannot set DEPTH dimension scales' , 0)
        if( vartgrd(idvar) ) then
          ret = sfsdmname(dimid2,lonlabel(1))
          call hdf_error( ret, 'cannot set dimname', 0)
          ret = sfsdmname(dimid3, latlabel(1))
          call hdf_error( ret, 'cannot set dimname', 0)
          ret = sfsdscale(dimid2,IMT_M,REAL_HDF, xt)
          ret2= sfsdscale(dimid3,JMT_M,REAL_HDF, yt)
        else
          ret = sfsdmname(dimid2,lonlabel(2))
          call hdf_error( ret, 'cannot set dimname', 0)
          ret = sfsdmname(dimid3,latlabel(2))
          call hdf_error( ret, 'cannot set dimname', 0)
          ret = sfsdscale(dimid2,IMT_M,REAL_HDF,xu)
          ret2= sfsdscale(dimid3,JMT_M,REAL_HDF,yu)
        endif
        call hdf_error( ret, 'cannot set dimension scales' , 0)
        call hdf_error( ret2, 'cannot set dimension scales' ,0)
        ret = sfsdmstr(dimid1,labels3(1),units3(1),fmats3(1))
        call hdf_error( ret, 'setting dim strings',dimid1, 0)
        ret = sfsdmstr(dimid2,labels3(2),units3(2),fmats3(2))
        call hdf_error( ret, 'setting dim strings',dimid2, 0)
        ret = sfsdmstr(dimid3,labels3(3),units3(3),fmats3(3))
        call hdf_error( ret, 'setting dim strings',dimid3, 0)
      endif
c
c set the start and stride of the dataset to written
c
      srcdims(1) = 0
      srcdims(2) = 0
      srcdims(3) = jarchl-1
      stride(1) = 1
      stride(2) = 1
      stride(3) = 1
      edge(1) = KM
      edge(2) = IMT_M
      edge(3) = 1+jarchu-jarchl
      ret = sfwdata(sdsid,srcdims,stride,edge,rest3d)
      call hdf_error(ret,'error writing 3-D dataset ',idvar)
c
c on the last call, write the dataset attributes and close dataset
c
      if(jarchu.ge.JMT_M)then
        ret = sfsattr(sdsid, LEVEL_ATTR, 
     &                INTEGER_HDF, 1, 0)
        call hdf_error(ret,' error writing attribute', 0)
        ret = sfsattr(sdsid, TGRID_ATTR, 
     &                INTEGER_HDF, 1,vartgrd(idvar))
        call hdf_error(ret,' error writing attribute', 0)
        ret = sfendacc(sdsid)
        call hdf_error( ret, 'cannot close dataset', 0)
      endif
c
c------------------------------------------------------------------
c End of file.  Write the global attributes ... 
c------------------------------------------------------------------
c
      if(idvar.eq.NARCHV.and.jarchu.ge.JMT_M)then
        ret = sfsattr(write_id, FMODE_ATTR, 
     &                INTEGER_HDF, 1, MODELLED)
        call hdf_error(ret,' error writing global attribute', 0)
        ret = sfsattr(write_id, FTYPE_ATTR, 
     &              INTEGER_HDF, 1, RESTART)
        call hdf_error(ret,' error writing global attribute', 0)
c
c set rotation attribute.  1= non-rotated 2=rotated
c
        ret = sfsattr(write_id, ROTATION_ATTR, 
     &              INTEGER_HDF, 1, 1)
        call hdf_error(ret,' error writing global attribute', 0)
c
c------------------------------------------------------------------
c ... and close file 
c------------------------------------------------------------------
c
        ret = sfend(write_id)
        call hdf_error( ret,' error closing SD interface', 0)
        print *,' archwr_hdf:  hdf archive file written'
      endif
c 
      larchq = .false.
      return
      end

      subroutine snap_hdf
c
c=======================================================================c  This routine produces hdf snapshots.  The present version of the
c  routine saves the fields h0, u0, v0.
c  The generated file name starts with 's' and ends with '.hdf'.
c  At present only supports 2-D data
c=======================================================================c
#include "def_master.h"
#include "param.h"
#include "coord.h"
#include "ctmngr.h"
#include "iounit.h"
#include "levind.h"
#include "mesdta.h"
#include "snaps.h"
#include "switch.h"
#include "timelv.h"
#include "varinfo.h"
c
c using multi-file interface HDF3.3r3
c
# include "constant.i"
# include "hdf.h"
c
c  local variables
c
      logical*4 around
      integer*4 write_id,sdsid,dimid1,dimid2,ret,ret2,
     &          srcdims(MAXDIMS),stride(MAXDIMS),irank
c
c=======================================================================c open file and write fixed data
c=======================================================================c
      if(idsnp.eq.0)then
c
c  zero dataset
c
        write(hdffile, '(i10.10)') nint(totday)
        hdffile = 's'//hdffile(7:10)//cfext(1)
        inquire(file = hdffile, exist = around )
        if(around) then
          write( hdffile,'(i10.10)') itt
          hdffile = 's'// hdffile(7:10)//cfext(2)
        endif
c        
        write_id = sfstart(hdffile, DFACC_CREATE)
        if(write_id.eq.-1) then
          write(stdout,*),' error opening the file'
          return
        endif
        write (stdout,24) itt, hdffile
  24  format('  Writing hdf file for itt = ',i7,/,
     &       '  hdf filename = ',a) 
c
c write the dataset containing 'itt'
c
      irank = 1
      sdsid = sfcreate( write_id,tlabel(1),
     &                  INTEGER_HDF,irank,itdims )
      call hdf_error( sdsid,'cannot create itt dataset', 0)
      ret = sfsdtstr( sdsid, 'TIMESTEP', '', '','')
      call hdf_error( ret, 'cannot set datastrings', 0)
      ret = sfwdata( sdsid, 0,1,itdims, itt) 
      call hdf_error( ret, 'cannot write dataset', 0)
      ret = sfendacc(sdsid)
      call hdf_error( ret, 'cannot close dataset', 0)
c
c write the dataset containing 'totsec'
c
      sdsid = sfcreate(write_id,tlabel(2),
     &                 REAL_HDF,itsize,itdims)
      call hdf_error( sdsid,'cannot create totsec dataset', 0)
      ret = sfsdtstr( sdsid, 'Model Seconds', '', '')
      call hdf_error( ret, 'cannot set datastrings', 0)
      ret = sfwdata( sdsid, 0, 1, itdims, totsec)
      call hdf_error( ret, 'cannot write dataset', 0)
      ret = sfendacc(sdsid)
      call hdf_error( ret, 'cannot close dataset', 0)
c
c write the dataset containing 'totday'
c
      sdsid = sfcreate(write_id,tlabel(3),
     &                 REAL_HDF,itsize,itdims)
      call hdf_error( sdsid,'cannot create totday dataset', 0)
      ret = sfsdtstr( sdsid, 'Model Day', '', '')
      call hdf_error( ret, 'cannot set datastrings', 0)
      ret = sfwdata(  sdsid, 0, 1,itdims, totday)
      call hdf_error( ret, 'cannot write dataset', 0)
      ret = sfendacc(sdsid)
      call hdf_error( ret, 'cannot close dataset', 0)
c
c write the dataset containing 'years'
c
      sdsid = sfcreate(write_id,tlabel(4),
     &                 REAL_HDF,itsize,itdims)
      call hdf_error( sdsid,'cannot create years dataset', 0)
      ret = sfsdtstr( sdsid, 'Model Year', '', '')
      call hdf_error( ret, 'cannot set datastrings', 0)
      ret = sfwdata(  sdsid,0,1,itdims,years)
      call hdf_error( ret, 'cannot write dataset', 0)
      ret = sfendacc(sdsid)
      call hdf_error( ret, 'cannot close dataset', 0)
c
c initialise the KMT array field.
c
      irank = 2
      srcdims(1)= 0
      srcdims(2)= 0
      stride(1) = 1
      stride(2) = 1
c
c call the HDF function to create the dataset.
c
      sdsid = sfcreate( write_id,'KMT',
     &                  INTEGER_HDF, irank, dims2) 
      call hdf_error( sdsid,'cannot create kmt dataset', 0)
      ret = sfsdtstr( sdsid, 'KMT', '', '')
      call hdf_error( ret, ' cannot set datastrings', 0)
      ifill_value = 999
      ret = sfsfill(sdsid,ifill_value)
      call hdf_error( ret, ' cannot set fill value', 0)
c
      dimid1 = sfdimid(sdsid,0)
      call hdf_error( dimid1, ' error getting dimid', 0)
      dimid2 = sfdimid(sdsid,1)
      call hdf_error( dimid2, ' error getting dimid', 0)
      ret = sfsdmname(dimid1,lonlabel(1))
      call hdf_error( ret, ' cannot set dimname', 0)
      ret = sfsdmname(dimid2,latlabel(1))
      call hdf_error( ret, ' cannot set dimname', 0)
      ret = sfsdmstr(dimid1,labels2(1),units2(1),fmats2(1))
      call hdf_error( ret, ' setting dimension string', 0)
      ret = sfsdmstr(dimid2,labels2(2),units2(2),fmats2(2))
      call hdf_error( ret, ' setting dimension string', 0)
      ret = sfsdscale(dimid1, dims2(1),REAL_HDF, xt)
      call hdf_error( ret, 'cannot set scales', dimid1)
      ret = sfsdscale(dimid2, dims2(2),REAL_HDF, yt)
      call hdf_error( ret, 'cannot set scales', dimid2)
      ret = sfwdata( sdsid,srcdims(1),stride(1),dims2(1),kmt)
      call hdf_error( ret,'cannot write kmt dataset',0)
      ret = sfendacc(sdsid)
      call hdf_error( ret, 'cannot close dataset', 0)
c
c initialise the KMU array field.
c
      sdsid = sfcreate( write_id,'KMU',
     &                  INTEGER_HDF,irank,dims2)
      call hdf_error( sdsid,'cannot create kmu dataset', 0)
      ret = sfsdtstr( sdsid, 'KMU', '', '')
      call hdf_error( ret, 'cannot set datastrings', 0)
      ifill_value = 999
      ret = sfsfill(sdsid,ifill_value)
      call hdf_error( ret, 'cannot set fill value', 0)
c
      dimid1 = sfdimid(sdsid,0)
      call hdf_error( dimid1, 'error getting dimid', 0)
      dimid2 = sfdimid(sdsid,1)
      call hdf_error( dimid2, 'error getting dimid', 0)
      ret = sfsdmname(dimid1,lonlabel(2))
      call hdf_error( ret, 'cannot set dimname', 0)
      ret = sfsdmname(dimid2,latlabel(2))
      call hdf_error( ret, 'cannot set dimname', 0)
c
      ret = sfsdmstr(dimid1,labels2(1),units2(1),fmats2(1))
      call hdf_error( ret, 'setting dimension string', 0)
      ret = sfsdmstr(dimid2,labels2(2),units2(2),fmats2(2))
      call hdf_error( ret, 'setting dimension string', 0)
      ret = sfsdscale(dimid1, dims2(1),REAL_HDF,xu)
      call hdf_error( ret,'cannot set scales', dimid1)
      ret = sfsdscale(dimid2, dims2(2),REAL_HDF,yu)
      call hdf_error( ret,'cannot set scales', dimid2)
      ret = sfwdata( sdsid,srcdims(1),stride(1),dims2(1),kmu)
      call hdf_error( ret,'cannot write kmu dataset',0)
      ret = sfendacc(sdsid)
      call hdf_error( ret, ' cannot close dataset', 0)
c
c  set snap2d array to fill value
c
        do 110 j=1,JMT_M
        do 110 i=1,IMT_M
          snap2d(i,j) = 0.0
 110    continue
        return
      endif
c
c
c write two dimensional datasets to file
c
      if(idsnp.le.NSNAPS)then
        irank = 2
        iddsnp = ispvar(idsnp)
c
        sdsid = sfcreate( write_id,varid(iddsnp),
     &                    REAL_HDF,irank,dims2)
        call hdf_error( sdsid, 'cannot create dataset', 0)
        ret = sfsdtstr(sdsid,varname(iddsnp),varunit(iddsnp),'','')
        call hdf_error( ret, 'cannot set data strings',0)
        fill_value = 0.0
        ret = sfsfill(sdsid,fill_value)
         call hdf_error( ret, 'cannot set fill value', 0)
c
c set the scales for the first dimension
c
        dimid1 = sfdimid( sdsid,0)
        call hdf_error( dimid1, 'error getting dimid', 0)
        dimid2 = sfdimid( sdsid,1)
        call hdf_error( dimid2, 'error getting dimid', 0)
c
        if(vartgrd(iddsnp))then
          ret = sfsdmname(dimid1, lonlabel(1))
          call hdf_error( ret, 'cannot set dimname', 0)
          ret = sfsdmname(dimid2, latlabel(1))
          call hdf_error( ret, 'cannot set dimname', 0)
          ret = sfsdscale(dimid1,IMT_M,REAL_HDF,xt)
          ret2= sfsdscale(dimid2,JMT_M,REAL_HDF,yt)
        else
          ret = sfsdmname(dimid1, lonlabel(2))
          call hdf_error( ret, 'cannot set dimname', 0)
          ret = sfsdmname(dimid2, latlabel(2))
          call hdf_error( ret, 'cannot set dimname', 0)
          ret = sfsdscale(dimid1,IMT_M,REAL_HDF,xu)
          ret2= sfsdscale(dimid2,JMT_M,REAL_HDF,yu) 
        endif
        call hdf_error( ret, 'cannot set x-scales', dimid1)
        call hdf_error( ret2, 'cannot set y-scales',dimid2)
c
        ret = sfsdmstr(dimid1,labels2(1),units2(1),fmats2(1))
        call hdf_error( ret, 'setting dimension string', 0)
        ret = sfsdmstr(dimid2,labels2(2),units2(2),fmats2(2))
        call hdf_error( ret, 'setting dimension string', 0)
c
c  write out the dataset
c
      srcdims(1)= 0
      srcdims(2)= 0
      stride(1) = 1
      stride(2) = 1
        ret = sfwdata(sdsid,srcdims(1),stride(1),dims2(1),snap2d)
        call hdf_error( ret, 'writing dataset ', 0)
        ret = sfsattr(sdsid, LEVEL_ATTR, INTEGER_HDF, 1, 0)
        call hdf_error( ret, 'setting attribute', 0)
        ret = sfsattr(sdsid, TGRID_ATTR, INTEGER_HDF, 1,
     &                                     vartgrd(iddsnp))
        call hdf_error( ret, 'setting attribute', 0)
        ret = sfendacc(sdsid)
        call hdf_error( ret, 'cannot close dataset', 0)
c
c  set snap2d array to fill value before using it again
c
        do 410 j=1,JMT_M
        do 410 i=1,IMT_M
          snap2d(i,j) = 0.0
 410    continue
        lsnapq = .false.
c
c  dataset not supported yet !
c
      else
        lsnapq = .false.
      endif
c
c------------------------------------------------------------------
c End of file.  Write the global attributes ... 
c------------------------------------------------------------------
c
      if(idsnp.eq.NSNAPS)then
        ret = sfsattr(write_id, FMODE_ATTR, 
     &                INTEGER_HDF, 1, MODELLED)
        call hdf_error(ret,' error writing global attribute', 0)
        ret = sfsattr(write_id, FTYPE_ATTR, 
     &              INTEGER_HDF, 1, RESTART)
        call hdf_error(ret,' error writing global attribute', 0)
c
c set rotation attribute.  1= non-rotated 2=rotated
c
        ret = sfsattr(write_id, ROTATION_ATTR, 
     &              INTEGER_HDF, 1, 1)
        call hdf_error(ret,' error writing global attribute', 0)
c
c------------------------------------------------------------------
c ... and close file
c------------------------------------------------------------------
c
        ret = sfend(write_id)
        call hdf_error( ret,' error closing SD interface', 0)
        lsnapq = .false.
        print *,' snap_hdf:  hdf snapshot written'
      endif
      return
      end


        subroutine hdf_error(istat,cstr,info)
c
c subroutine to check return codes from the hdf calls
c   istat - return code
c   cstr  - character string explaining the error
c   info  - optional integer information (printed if non-zero)
c
#include "def_master.h"
#include "param.h"
#include "mesdta.h"
        character*(*) cstr
        integer istat, info
c
        if(istat .eq. -1) then
          if(info.ne.0)then
            write(stdout,*)'HDF Error: ',cstr, info
          else
            write(stdout,*)'HDF Error: ',cstr
          endif
          call m_send(MSG_ABORT)
        endif
c
        return
        end
#endif
