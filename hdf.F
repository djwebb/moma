      subroutine archwr_hdf
!
!======================================================================
!
!    File:  hdf.F
!    Copyright 1992-1997, 2020  David J. Webb
!
!    This file is part of the Moma ocean model code.
!    Released under licence GPL-3.0-or-later
!
!    Moma is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    Moma is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with Moma.  If not, see <https://www.gnu.org/licenses/>.
!
!======================================================================
!
!
!=======================================================================
!  Subroutine to save hdf version of the restart file
!  using multi-file HDF3.3r3 routines.
!=======================================================================
!
#include "param.h"
#include "ctmngr.h"
#include "frees.h"
#include "iounit.h"
#include "levind.h"
#include "slabs.h"
#include "switch.h"
#include "timelv.h"
!
!  hdf '.h' files
!
#include "constant.i"
#include "coord.h"
#include "hdf.h"
!
!  local variables
!
      logical*4 around
      integer*4 write_id,sds_id,dimid1,dimid2,dimid3,ret,ret2,
     &          irank,start(MAXDIMS),stride(MAXDIMS),edge(MAXDIMS)
!
!  initialise the hdf common blocks
!
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
      data (cfext(k),k=1,2) /'.hdf','.hdfa'/
      data (tlabel(n),n=1,4) /'TIMESTEP', 'SECONDS', 'DAY', 'YEARS'/
      data itsize,itdims(1) /1,1/
      data (varname(n),varunit(n),vartgrd(n),n=1,NARCHV)/
     &  'free surface height      ', 'cm          ', .true.  ,
     &  'barotropic u-velocity    ', 'cm/s        ', .false. ,
     &  'barotropic v-velocity    ', 'cm/s        ', .false. ,
     &  'potential temperature    ', 'C           ', .true.  ,
     &  'salinity                 ', '.001(NSU-35)', .true.  ,
     &  'u-velocity               ', 'cm/s        ', .false. ,
     &  'v-velocity               ', 'cm/s        ', .false. /
        data (varid(n), varcordsys(n), n=1, NARCHV)/
     &  'FREE SURFACE HEIGHT',          'T_GRID'        ,
     &  'BAROTROPIC U VELOCITY',        'U_GRID'        ,
     &  'BAROTROPIC V VELOCITY',        'U_GRID'        ,
     &  'POTENTIAL TEMPERATURE',        'T_GRID'        ,
     &  'SALINITY',                     'T_GRID'        ,
     &  'U VELOCITY',                   'U_GRID'        ,
     &  'V VELOCITY',                   'U_GRID'        /
!
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
!
! write the dataset containing 'itt'
!
      irank = 1
      sds_id = sfcreate( write_id,tlabel(1),
     &                  INTEGER_HDF,irank,itdims )
      call hdf_error(sds_id,'cannot create itt dataset', 0)
      ret = sfsdtstr(sds_id, 'TIMESTEP', '', '','')
      call hdf_error( ret, 'cannot set datastrings', 0)
      ret = sfwdata(sds_id, 0,1,itdims, itt)
      call hdf_error( ret, 'cannot write dataset', 0)
      ret = sfendacc(sds_id)
      call hdf_error( ret, 'cannot close dataset', 0)
!
! write the dataset containing 'totsec'
!
      sds_id = sfcreate(write_id,tlabel(2),
     &                 REAL_HDF,itsize,itdims)
      call hdf_error(sds_id,'cannot create totsec dataset', 0)
      ret = sfsdtstr(sds_id, 'Model Seconds', '', '')
      call hdf_error( ret, 'cannot set datastrings', 0)
      ret = sfwdata(sds_id, 0, 1, itdims, totsec)
      call hdf_error( ret, 'cannot write dataset', 0)
      ret = sfendacc(sds_id)
      call hdf_error( ret, 'cannot close dataset', 0)
!
! write the dataset containing 'totday'
!
      sds_id = sfcreate(write_id,tlabel(3),
     &                 REAL_HDF,itsize,itdims)
      call hdf_error(sds_id,'cannot create totday dataset', 0)
      ret = sfsdtstr(sds_id, 'Model Day', '', '')
      call hdf_error( ret, 'cannot set datastrings', 0)
      ret = sfwdata( sds_id, 0, 1,itdims, totday)
      call hdf_error( ret, 'cannot write dataset', 0)
      ret = sfendacc(sds_id)
      call hdf_error( ret, 'cannot close dataset', 0)
!
! write the dataset containing 'years'
!
      sds_id = sfcreate(write_id,tlabel(4),
     &                 REAL_HDF,itsize,itdims)
      call hdf_error(sds_id,'cannot create years dataset', 0)
      ret = sfsdtstr(sds_id, 'Model Year', '', '')
      call hdf_error( ret, 'cannot set datastrings', 0)
      ret = sfwdata( sds_id,0,1,itdims,years)
      call hdf_error( ret, 'cannot write dataset', 0)
      ret = sfendacc(sds_id)
      call hdf_error( ret, 'cannot close dataset', 0)
!
! initialise the KMT array field.
!
      irank = 2
      start(1)= 0
      start(2)= 0
      stride(1) = 1
      stride(2) = 1
#ifdef cyclic
      i1 = 2
      dims2(1) = imtm2
      dims2(2) = jmt
      edge(1)  = imtm2
      edge(2)  = 1
#else
      i1 = 1
      dims2(1) = imt
      dims2(2) = jmt
      edge(1)  = imt
      edge(2)  = jmt
#endif
!
! call the HDF function to create the dataset.
!
      sds_id = sfcreate( write_id,'KMT',
     &                  INTEGER_HDF, irank, dims2)
      call hdf_error(sds_id,'cannot create kmt dataset', 0)
      ret = sfsdtstr(sds_id, 'KMT', '', '')
      call hdf_error( ret, ' cannot set kmt datastrings', 0)
      ifill_value = 999
      ret = sfsfill(sds_id,ifill_value)
      call hdf_error( ret, ' cannot set fill value', 0)
!
      dimid1 = sfdimid(sds_id,0)
      call hdf_error( dimid1, ' error getting dimid 1', 0)
      dimid2 = sfdimid(sds_id,1)
      call hdf_error( dimid2, ' error getting dimid 2', 0)
      ret = sfsdmname(dimid1,lonlabel(1))
      call hdf_error( ret, ' cannot set dimname 1', 0)
      ret = sfsdmname(dimid2,latlabel(1))
      call hdf_error( ret, ' cannot set dimname 2', 0)
      ret = sfsdmstr(dimid1,labels2(1),units2(1),fmats2(1))
      call hdf_error( ret, ' setting dimension string 1', 0)
      ret = sfsdmstr(dimid2,labels2(2),units2(2),fmats2(2))
      call hdf_error( ret, ' setting dimension string 2', 0)
#ifndef REAL_8
!
!---------------------------------------------------------------------
!  skip writing scales when REAL_8 defined.  This is a fix to overcome
!  an error in pvm when lengths of scale and main variable (kmt) differ
!---------------------------------------------------------------------
!
      ret = sfsdscale(dimid1, dims2(1),REAL_HDF, xt(i1))
      call hdf_error( ret, 'cannot set scales xt', dimid1)
      ret = sfsdscale(dimid2, dims2(2),REAL_HDF, yt)
      call hdf_error( ret, 'cannot set scales yt', dimid2)
#endif
#ifdef cyclic
      do 110 j=1,jmt
        start(2) = j-1
        ret = sfwdata(sds_id,start,stride,edge,kmt(2,j))
        call hdf_error( ret,'cannot write kmt dataset',0)
 110  continue
#else
      ret = sfwdata(sds_id,start,stride,edge,kmt)
      call hdf_error( ret,'cannot write kmt dataset',0)
#endif
      ret = sfendacc(sds_id)
      call hdf_error( ret, 'cannot close kmt dataset', 0)
!
! initialise the KMU array field.
!
      sds_id = sfcreate( write_id,'KMU',
     &                  INTEGER_HDF,irank,dims2)
      call hdf_error(sds_id,'cannot create kmu dataset', 0)
      ret = sfsdtstr(sds_id, 'KMU', '', '')
      call hdf_error( ret, 'cannot set kmu datastrings', 0)
      ifill_value = 999
      ret = sfsfill(sds_id,ifill_value)
      call hdf_error( ret, 'cannot set fill value', 0)
!
      dimid1 = sfdimid(sds_id,0)
      call hdf_error( dimid1, 'error getting dimid 1', 0)
      dimid2 = sfdimid(sds_id,1)
      call hdf_error( dimid2, 'error getting dimid 2', 0)
      ret = sfsdmname(dimid1,lonlabel(2))
      call hdf_error( ret, 'cannot set dimname 1', 0)
      ret = sfsdmname(dimid2,latlabel(2))
      call hdf_error( ret, 'cannot set dimname 2', 0)
!
      ret = sfsdmstr(dimid1,labels2(1),units2(1),fmats2(1))
      call hdf_error( ret, 'setting dimension string 1', 0)
      ret = sfsdmstr(dimid2,labels2(2),units2(2),fmats2(2))
      call hdf_error( ret, 'setting dimension string 2', 0)
#ifndef REAL_8
!
!---------------------------------------------------------------------
!  skip writing scales when REAL_8 defined.  This is a fix to overcome
!  an error in pvm when lengths of scale and main variable (kmu) differ
!---------------------------------------------------------------------
!
      ret = sfsdscale(dimid1, dims2(1),REAL_HDF,xu(i1))
      call hdf_error( ret,'cannot set scales xu', dimid1)
      ret = sfsdscale(dimid2, dims2(2),REAL_HDF,yu)
      call hdf_error( ret,'cannot set scales yu', dimid2)
#endif
#ifdef cyclic
      do 120 j=1,jmt
        start(2) = j-1
        ret = sfwdata(sds_id,start,stride,edge,kmu(2,j))
        call hdf_error( ret,'cannot write kmu dataset',0)
 120  continue
#else
      ret = sfwdata(sds_id,start,stride,edge,kmu)
      call hdf_error( ret,'cannot write kmu dataset',0)
#endif
      ret = sfendacc(sds_id)
      call hdf_error( ret, ' cannot close kmu dataset', 0)
!
! write two dimensional datasets to file
!
      do 10 idvar = 1, NUM2D
!
        sds_id = sfcreate( write_id,varid(idvar),
     &                    REAL_HDF,irank,dims2)
        call hdf_error(sds_id, 'cannot create 2-D dataset', 0)
        ret = sfsdtstr(sds_id,varname(idvar),varunit(idvar),'','')
        call hdf_error( ret, 'cannot set 2-D data strings',0)
        fill_value = 0.0
        ret = sfsfill(sds_id,fill_value)
         call hdf_error( ret, 'cannot set 2-D fill value', 0)
!
! set the scales
!
        dimid1 = sfdimid(sds_id,0)
        call hdf_error( dimid1, 'error getting dimid 1', 0)
        dimid2 = sfdimid(sds_id,1)
        call hdf_error( dimid2, 'error getting dimid 2', 0)
        if(vartgrd(idvar))then
          ret = sfsdmname(dimid1, lonlabel(1))
          call hdf_error( ret, 'cannot set dimname 1', 0)
          ret = sfsdmname(dimid2, latlabel(1))
          call hdf_error( ret, 'cannot set dimname 2', 0)
          ret = sfsdscale(dimid1,dims2(1),REAL_HDF,xt(i1))
          ret2= sfsdscale(dimid2,dims2(2),REAL_HDF,yt)
        else
          ret = sfsdmname(dimid1, lonlabel(2))
          call hdf_error( ret, 'cannot set dimname 1', 0)
          ret = sfsdmname(dimid2, latlabel(2))
          call hdf_error( ret, 'cannot set dimname 2', 0)
          ret = sfsdscale(dimid1,dims2(1),REAL_HDF,xu(i1))
          ret2= sfsdscale(dimid2,dims2(2),REAL_HDF,yu)
        endif
        call hdf_error( ret, 'cannot set x-scales', dimid1)
        call hdf_error( ret2, 'cannot set y-scales',dimid2)
!
        ret = sfsdmstr(dimid1,labels2(1),units2(1),fmats2(1))
        call hdf_error( ret, 'setting dimension string 1', 0)
        ret = sfsdmstr(dimid2,labels2(2),units2(2),fmats2(2))
        call hdf_error( ret, 'setting dimension string 2', 0)
!
!  write out the dataset
!
#ifdef cyclic
      do 130 j=1,jmt
        start(2) = j-1
        if(idvar.eq.1) then
          ret = sfwdata(sds_id,start,stride,edge,h0(2,j,np0))
        elseif(idvar.eq.2) then
          ret = sfwdata(sds_id,start,stride,edge,u0(2,j,np0))
        elseif(idvar.eq.3) then
          ret = sfwdata(sds_id,start,stride,edge,v0(2,j,np0))
        endif
 130  continue
#else
        if(idvar.eq.1) then
          ret = sfwdata(sds_id,start,stride,edge,h0(1,1,np0))
        elseif(idvar.eq.2) then
          ret = sfwdata(sds_id,start,stride,edge,u0(1,1,np0))
        elseif(idvar.eq.3) then
          ret = sfwdata(sds_id,start,stride,edge,v0(1,1,np0))
        endif
#endif
        call hdf_error( ret, 'writing 2-D dataset ', 0)
        ret = sfsattr(sds_id, LEVEL_ATTR, INTEGER_HDF, 1, 0)
        call hdf_error( ret, 'setting 2-D attribute', 0)
        ret = sfsattr(sds_id, TGRID_ATTR, INTEGER_HDF, 1,
     &                                    vartgrd(idvar))
        call hdf_error( ret, 'setting 2-D attribute', 0)
        ret = sfendacc(sds_id)
        call hdf_error( ret, 'cannot close 2-D dataset', 0)
!
 10   continue
!
! archive the 3-D datasets
!
      irank = 3
      start(1) = 0
      start(2) = 0
      start(3) = 0
      stride(1) = 1
      stride(2) = 1
      stride(3) = 1
#ifdef cyclic
      i1 = 2
      dims3(1) = km
      dims3(2) = imtm2
      dims3(3) = jmt
      edge(1) = km
      edge(2) = imtm2
      edge(3) = 1
#else
      i1 = 1
      dims3(1) = km
      dims3(2) = imt
      dims3(3) = jmt
      edge(1) = km
      edge(2) = imt
      edge(3) = jmt
#endif
      do 20 idvar = NUM2D+1,NARCHV
!
        sds_id = sfcreate(write_id,varid(idvar),
     &                   REAL_HDF,irank,dims3)
        call hdf_error(sds_id, 'cannot create 3-D dataset', 0)
        ret = sfsdtstr(sds_id,varname(idvar),varunit(idvar),'','')
        call hdf_error( ret, 'cannot set 3-D data strings', 0)
        fill_value = 0.0
        ret = sfsfill(sds_id,fill_value)
        call hdf_error( ret, 'cannot set 3-D fill value', 0)
!
        dimid1 = sfdimid(sds_id,0)
        call hdf_error( dimid1, 'Cannot get dimid 1:',dimid1)
        dimid2 = sfdimid(sds_id,1)
        call hdf_error( dimid2, 'Cannot get dimid 2:',dimid2)
        dimid3 = sfdimid(sds_id,2)
        call hdf_error( dimid3, 'Cannot get dimid 3:',dimid3)
!
! set the scales for each dimensions
!
        ret = sfsdmname(dimid1,'DEPTH')
        ret = sfsdscale(dimid1,KM,REAL_HDF,zt)
        if( vartgrd(idvar) ) then
          ret = sfsdmname(dimid2,lonlabel(1))
          call hdf_error( ret, 'cannot set dimname 2', 0)
          ret = sfsdmname(dimid3, latlabel(1))
          call hdf_error( ret, 'cannot set dimname 3', 0)
          ret = sfsdscale(dimid2,dims3(2),REAL_HDF,xt(i1))
          ret2= sfsdscale(dimid3,dims3(3),REAL_HDF,yt)
        else
          ret = sfsdmname(dimid2,lonlabel(2))
          call hdf_error( ret, 'cannot set dimname 2', 0)
          ret = sfsdmname(dimid3,latlabel(2))
          call hdf_error( ret, 'cannot set dimname 3', 0)
          ret = sfsdscale(dimid2,dims3(2),REAL_HDF,xu(i1))
          ret2= sfsdscale(dimid3,dims3(3),REAL_HDF,yu)
        endif
        call hdf_error( ret, 'cannot set dimension scale 2' , 0)
        call hdf_error( ret2, 'cannot set dimension scale 3' ,0)
        ret = sfsdmstr(dimid1,labels3(1),units3(1),fmats3(1))
        call hdf_error( ret, 'setting dim strings 1',dimid1)
        ret = sfsdmstr(dimid2,labels3(2),units3(2),fmats3(2))
        call hdf_error( ret, 'setting dim strings 2',dimid2)
        ret = sfsdmstr(dimid3,labels3(3),units3(3),fmats3(3))
        call hdf_error( ret, 'setting dim strings 3',dimid3)
!
#ifdef cyclic
        do 220 j=1,jmt
          start(3) = j-1
          if(idvar.eq.4) then
            ret = sfwdata(sds_id,start,stride,edge,t(1,2,j,1,np))
          elseif(idvar.eq.5) then
            ret = sfwdata(sds_id,start,stride,edge,t(1,2,j,2,np))
          elseif(idvar.eq.6) then
            ret = sfwdata(sds_id,start,stride,edge,u(1,2,j,np))
          elseif(idvar.eq.7) then
            ret = sfwdata(sds_id,start,stride,edge,v(1,2,j,np))
          endif
          if(ret.ne.0)then
            write(stdout,*)' error writing 3-D dataset at j = ',j
            write(stdout,*)' start = ',start
            write(stdout,*)' stride = ',stride
            write(stdout,*)' edge = ',edge
          endif
          call hdf_error(ret,'error writing 3-D dataset ',idvar)
 220    continue
#else
        if(idvar.eq.4) then
          ret = sfwdata(sds_id,start,stride,edge,t(1,1,1,1,np))
        elseif(idvar.eq.5) then
          ret = sfwdata(sds_id,start,stride,edge,t(1,1,1,2,np))
        elseif(idvar.eq.6) then
          ret = sfwdata(sds_id,start,stride,edge,u(1,1,1,np))
        elseif(idvar.eq.7) then
          ret = sfwdata(sds_id,start,stride,edge,v(1,1,1,np))
        endif
        call hdf_error(ret, 'error writing 3-D dataset:', idvar)
#endif
!
! write the dataset attributes
!
        ret = sfsattr(sds_id, LEVEL_ATTR,
     &                INTEGER_HDF, 1, 0)
        call hdf_error(ret,' error writing 3-D attribute', 0)
        ret = sfsattr(sds_id, TGRID_ATTR,
     &                INTEGER_HDF, 1,vartgrd(idvar))
        call hdf_error(ret,' error writing 3-D attribute', 0)
!
! End access to this dataset
!
        ret = sfendacc(sds_id)
        call hdf_error( ret, 'cannot close 3-D dataset', 0)
!
!
 20   continue
!
! Now write the global attributes for the whole file
!
      ret = sfsattr(write_id, FMODE_ATTR,
     &              INTEGER_HDF, 1, MODELLED)
      call hdf_error(ret,' error writing global attribute', 0)
      ret = sfsattr(write_id, FTYPE_ATTR,
     &              INTEGER_HDF, 1, RESTART)
      call hdf_error(ret,' error writing global attribute', 0)
!
! set rotation attribute.  1= non-rotated 2=rotated
!
      ret = sfsattr(write_id, ROTATION_ATTR,
     &              INTEGER_HDF, 1, 1)
      call hdf_error(ret,' error writing global attribute', 0)
      ret = sfend(write_id)
      call hdf_error( ret,' error closing SD interface', 0)
!
        end


      subroutine snap_hdf
!
!=======================================================================
!  This routine produces hdf snapshots.  The present version of the
!  routine saves the fields h0, u0, v0.
!  The generated file name starts with 's' and ends with '.hdf'.
!=======================================================================
!
#include "param.h"
#include "ctmngr.h"
#include "frees.h"
#include "iounit.h"
#include "levind.h"
#include "slabs.h"
#include "switch.h"
#include "timelv.h"
!
! using multi-file interface HDF3.3r3 routines for the data
!
# include "constant.i"
# include "coord.h"
# include "hdf.h"
!
!  local variables
!
      logical*4 around
      integer*4 write_id,sds_id,dimid1,dimid2,ret,ret2,
     &          irank,start(MAXDIMS),stride(MAXDIMS),edge(MAXDIMS)
!
      write(hdffile, '(i10.10)') nint(totday)
      hdffile = 's'//hdffile(7:10)//cfext(1)
      inquire(file = hdffile, exist = around )
      if(around) then
         write( hdffile,'(i10.10)') itt
         hdffile = 's'// hdffile(7:10)//cfext(2)
      endif
!
      write_id = sfstart(hdffile, DFACC_CREATE)
      if(write_id.eq.-1) then
        write(stdout,*),' error opening the file'
        return
      endif
      write (stdout,24) itt, hdffile
  24  format('  Writing hdf snapshot file for itt = ',i7,/,
     &       '  hdf filename = ',a)
!
! write the dataset containing 'itt'
!
      irank = 1
      sds_id = sfcreate( write_id,tlabel(1),
     &                  INTEGER_HDF,irank,itdims )
      call hdf_error(sds_id,'cannot create itt dataset', 0)
      ret = sfsdtstr(sds_id, 'TIMESTEP', '', '','')
      call hdf_error( ret, 'cannot set itt datastrings', 0)
      ret = sfwdata(sds_id, 0,1,itdims, itt)
      call hdf_error( ret, 'cannot write itt dataset', 0)
      ret = sfendacc(sds_id)
      call hdf_error( ret, 'cannot close itt dataset', 0)
!
! write the dataset containing 'totsec'
!
      sds_id = sfcreate(write_id,tlabel(2),
     &                 REAL_HDF,itsize,itdims)
      call hdf_error(sds_id,'cannot create totsec dataset', 0)
      ret = sfsdtstr(sds_id, 'Model Seconds', '', '')
      call hdf_error( ret, 'cannot set totsec datastrings', 0)
      ret = sfwdata(sds_id, 0, 1, itdims, totsec)
      call hdf_error( ret, 'cannot write totsec dataset', 0)
      ret = sfendacc(sds_id)
      call hdf_error( ret, 'cannot close totsec dataset', 0)
!
! write the dataset containing 'totday'
!
      sds_id = sfcreate(write_id,tlabel(3),
     &                 REAL_HDF,itsize,itdims)
      call hdf_error(sds_id,'cannot create totday dataset', 0)
      ret = sfsdtstr(sds_id, 'Model Day', '', '')
      call hdf_error( ret, 'cannot set totday datastrings', 0)
      ret = sfwdata( sds_id, 0, 1,itdims, totday)
      call hdf_error( ret, 'cannot write totday dataset', 0)
      ret = sfendacc(sds_id)
      call hdf_error( ret, 'cannot close totday dataset', 0)
!
! write the dataset containing 'years'
!
      sds_id = sfcreate(write_id,tlabel(4),
     &                 REAL_HDF,itsize,itdims)
      call hdf_error(sds_id,'cannot create years dataset', 0)
      ret = sfsdtstr(sds_id, 'Model Year', '', '')
      call hdf_error( ret, 'cannot set years datastrings', 0)
      ret = sfwdata( sds_id,0,1,itdims,years)
      call hdf_error( ret, 'cannot write years dataset', 0)
      ret = sfendacc(sds_id)
      call hdf_error( ret, 'cannot close years dataset', 0)
!
! initialise the KMT array field.
!
      irank = 2
      start(1)  = 0
      start(2)  = 0
      stride(1) = 1
      stride(2) = 1
#ifdef cyclic
      i1 = 2
      dims2(1) = imtm2
      dims2(2) = jmt
      edge(1)  = imtm2
      edge(2)  = 1
#else
      i1 = 1
      dims2(1) = imt
      dims2(2) = jmt
      edge(1)  = imt
      edge(2)  = jmt
#endif
!
! call the HDF function to create the dataset.
!
      sds_id = sfcreate( write_id,'KMT',
     &                  INTEGER_HDF, irank, dims2)
      call hdf_error(sds_id,'cannot create kmt dataset', 0)
      ret = sfsdtstr(sds_id, 'KMT', '', '')
      call hdf_error( ret, ' cannot set kmt datastrings', 0)
      ifill_value = 999
      ret = sfsfill(sds_id,ifill_value)
      call hdf_error( ret, ' cannot set kmt fill value', 0)
!
      dimid1 = sfdimid(sds_id,0)
      call hdf_error( dimid1, ' error getting dimid', 1)
      dimid2 = sfdimid(sds_id,1)
      call hdf_error( dimid2, ' error getting dimid', 2)
      ret = sfsdmname(dimid1,lonlabel(1))
      call hdf_error( ret, ' cannot set dimname', 1)
      ret = sfsdmname(dimid2,latlabel(1))
      call hdf_error( ret, ' cannot set dimname', 2)
      ret = sfsdmstr(dimid1,labels2(1),units2(1),fmats2(1))
      call hdf_error( ret, ' setting dimension string', 1)
      ret = sfsdmstr(dimid2,labels2(2),units2(2),fmats2(2))
      call hdf_error( ret, ' setting dimension string', 2)
#ifndef REAL_8
!
!---------------------------------------------------------------------
!  skip writing scales when REAL_8 defined.  This is a fix to overcome
!  an error in pvm when lengths of scale and main variable (kmt) differ
!---------------------------------------------------------------------
!
      ret = sfsdscale(dimid1, dims2(1),REAL_HDF, xt(i1))
      call hdf_error( ret, 'cannot set scales 1', dimid1)
      ret = sfsdscale(dimid2, dims2(2),REAL_HDF, yt)
      call hdf_error( ret, 'cannot set scales 2', dimid2)
#endif
#ifdef cyclic
      do 110 j=1,jmt
        start(2) = j-1
        ret = sfwdata(sds_id,start,stride,edge,kmt(2,j))
        call hdf_error( ret,'cannot write kmt dataset',0)
 110  continue
#else
      ret = sfwdata(sds_id,start,stride,edge,kmt)
      call hdf_error( ret,'cannot write kmt dataset',0)
#endif
      ret = sfendacc(sds_id)
      call hdf_error( ret, 'cannot close dataset', 0)
!
! initialise the KMU array field.
!
      sds_id = sfcreate( write_id,'KMU',
     &                  INTEGER_HDF,irank,dims2)
      call hdf_error(sds_id,'cannot create kmu dataset', 0)
      ret = sfsdtstr(sds_id, 'KMU', '', '')
      call hdf_error( ret, 'cannot set kmu datastrings', 0)
      ifill_value = 999
      ret = sfsfill(sds_id,ifill_value)
      call hdf_error( ret, 'cannot set kmu fill value', 0)
!
      dimid1 = sfdimid(sds_id,0)
      call hdf_error( dimid1, 'error getting dimid', 1)
      dimid2 = sfdimid(sds_id,1)
      call hdf_error( dimid2, 'error getting dimid', 2)
      ret = sfsdmname(dimid1,lonlabel(2))
      call hdf_error( ret, 'cannot set dimname', 1)
      ret = sfsdmname(dimid2,latlabel(2))
      call hdf_error( ret, 'cannot set dimname', 2)
!
      ret = sfsdmstr(dimid1,labels2(1),units2(1),fmats2(1))
      call hdf_error( ret, 'setting dimension string', 1)
      ret = sfsdmstr(dimid2,labels2(2),units2(2),fmats2(2))
      call hdf_error( ret, 'setting dimension string', 2)
#ifndef REAL_8
!
!---------------------------------------------------------------------
!  skip writing scales when REAL_8 defined.  This is a fix to overcome
!  an error in pvm when lengths of scale and main variable (kmt) differ
!---------------------------------------------------------------------
!
      ret = sfsdscale(dimid1, dims2(1),REAL_HDF,xu(i1))
      call hdf_error( ret,'cannot set scales 1', dimid1)
      ret = sfsdscale(dimid2, dims2(2),REAL_HDF,yu)
      call hdf_error( ret,'cannot set scales 2', dimid2)
#endif
#ifdef cyclic
      do 120 j=1,jmt
        start(2) = j-1
        ret = sfwdata(sds_id,start,stride,edge,kmu(2,j))
        call hdf_error( ret,'cannot write kmu dataset',0)
 120  continue
#else
      ret = sfwdata(sds_id,start,stride,edge,kmu)
      call hdf_error( ret,'cannot write kmu dataset',0)
#endif
      ret = sfendacc(sds_id)
      call hdf_error( ret, ' cannot close kmu dataset', 0)
!
!
! write two dimensional datasets to file
!
      do 10 idvar = 1, NUM2D
!
        sds_id = sfcreate( write_id,varid(idvar),
     &                    REAL_HDF,irank,dims2)
        call hdf_error(sds_id, 'cannot create 2-D dataset', 0)
        ret = sfsdtstr(sds_id,varname(idvar),varunit(idvar),'','')
        call hdf_error( ret, 'cannot set 2-D data strings',0)
        fill_value = 0.0
        ret = sfsfill(sds_id,fill_value)
        call hdf_error( ret, 'cannot set 2-D fill value', 0)
!
! set the scales for the first dimension
!
        dimid1 = sfdimid(sds_id,0)
        call hdf_error( dimid1, 'error getting dimid', 1)
        dimid2 = sfdimid(sds_id,1)
        call hdf_error( dimid2, 'error getting dimid', 2)
!
        if(vartgrd(idvar))then
          ret = sfsdmname(dimid1, lonlabel(1))
          call hdf_error( ret, 'cannot set dimname', 1)
          ret = sfsdmname(dimid2, latlabel(1))
          call hdf_error( ret, 'cannot set dimname', 2)
          ret = sfsdscale(dimid1,dims2(1),REAL_HDF,xt(i1))
          ret2= sfsdscale(dimid2,dims2(2),REAL_HDF,yt)
        else
          ret = sfsdmname(dimid1, lonlabel(2))
          call hdf_error( ret, 'cannot set dimname', 1)
          ret = sfsdmname(dimid2, latlabel(2))
          call hdf_error( ret, 'cannot set dimname', 2)
          ret = sfsdscale(dimid1,dims2(1),REAL_HDF,xu(i1))
          ret2= sfsdscale(dimid2,dims2(2),REAL_HDF,yu)
        endif
        call hdf_error( ret, 'cannot set x-scales', dimid1)
        call hdf_error( ret2, 'cannot set y-scales',dimid2)
!
        ret = sfsdmstr(dimid1,labels2(1),units2(1),fmats2(1))
        call hdf_error( ret, 'setting dimension string', 1)
        ret = sfsdmstr(dimid2,labels2(2),units2(2),fmats2(2))
        call hdf_error( ret, 'setting dimension string', 2)
!
!  write out the dataset
!
#ifdef cyclic
      do 130 j=1,jmt
        start(2) = j-1
        if(idvar.eq.1) then
          ret = sfwdata(sds_id,start,stride,edge,h0(2,j,np0))
        elseif(idvar.eq.2) then
          ret = sfwdata(sds_id,start,stride,edge,u0(2,j,np0))
        elseif(idvar.eq.3) then
          ret = sfwdata(sds_id,start,stride,edge,v0(2,j,np0))
        endif
 130  continue
#else
        if(idvar.eq.1) then
          ret = sfwdata(sds_id,start,stride,edge,h0(1,1,np0))
        elseif(idvar.eq.2) then
          ret = sfwdata(sds_id,start,stride,edge,u0(1,1,np0))
        elseif(idvar.eq.3) then
          ret = sfwdata(sds_id,start,stride,edge,v0(1,1,np0))
        endif
#endif
        call hdf_error( ret, 'writing 2-D dataset ', 0)
        ret = sfsattr(sds_id, LEVEL_ATTR, INTEGER_HDF, 1, 0)
        call hdf_error( ret, 'setting attribute', 0)
        ret = sfsattr(sds_id, TGRID_ATTR, INTEGER_HDF, 1,vartgrd(idvar))
        call hdf_error( ret, 'setting attribute', 0)
        ret = sfendacc(sds_id)
        call hdf_error( ret, 'cannot close dataset', 0)
!
 10   continue
!
! Now write the global attributes for the whole file
!
      ret = sfsattr(write_id, FMODE_ATTR,
     &              INTEGER_HDF, 1, MODELLED)
      call hdf_error(ret,' error writing global attribute', 0)
      ret = sfsattr(write_id, FTYPE_ATTR,
     &              INTEGER_HDF, 1, RESTART)
      call hdf_error(ret,' error writing global attribute', 0)
!
! set rotation attribute.  1= non-rotated 2=rotated
!
      ret = sfsattr(write_id, ROTATION_ATTR,
     &              INTEGER_HDF, 1, 1)
      call hdf_error(ret,' error writing global attribute', 0)
      ret = sfend(write_id)
      call hdf_error( ret,' error closing SD interface', 0)
!
        end

        subroutine hdf_error(istat,cstr,info)
!
!=======================================================================
! subroutine to check return codes from the hdf calls
!   istat - return code
!   cstr  - character string explaining the error
!   info  - optional integer information (printed if non-zero)
!=======================================================================
!
        character*(*) cstr
        integer istat, info
!
        if(istat .eq. -1) then
          if(info.ne.0)then
            write(stdout,*)'HDF Error: ',cstr, info
          else
            write(stdout,*)'HDF Error: ',cstr
          endif
          stop
        endif
!
        return
        end

      subroutine archrd_hdf
!
!=======================================================================
!  subroutine to read hdf version of the restart file
!=======================================================================
!
#include "param.h"
#include "ctmngr.h"
#include "frees.h"
#include "iounit.h"
#include "levind.h"
#include "slabs.h"
#include "switch.h"
#include "timelv.h"
!
!  hdf '.h' files
!
#include "constant.i"
#include "coord.h"
#include "hdf.h"
!
!  local variables
!
      logical*4 around
      integer*4 file_id,sds_id,ret,
     &          start(MAXDIMS),stride(MAXDIMS),edge(MAXDIMS)
!
!  open file
!
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
      file_id = sfstart(fnrest, DFACC_READ)
      call hdf_error( file_id,'cannot open restart dataset', 0)
!
!  read dataset information
!
      start(1) =0
      start(2) =0
      start(3) =0
      stride(1)=1
      stride(2)=1
      stride(3)=1
      edge(1)  =1
      edge(2)  =1
      edge(3)  =1
!
      index = sfn2index(file_id,tlabel(1))
      call hdf_error(index,'cannot find id for TIMESTEP', 0)
      sds_id = sfselect(file_id,index)
      call hdf_error(sds_id,'error getting index for TIMESTEP', 0)
      ret = sfrdata(sds_id,start, stride,edge,itt)
      call hdf_error(ret,'error reading TIMESTEP', 0)
      ret = sfendacc(sds_id)
      call hdf_error( ret, 'cannot close dataset', 0)
!
      index = sfn2index(file_id,tlabel(2))
      call hdf_error( index,'cannot find id for totsec', 0)
      sds_id = sfselect(file_id,index)
      call hdf_error(sds_id,'error getting index for totsec', 0)
      ret = sfrdata(sds_id,start, stride,edge,totsec)
      call hdf_error(ret,'error reading totsec', 0)
      ret = sfendacc(sds_id)
      call hdf_error( ret, 'cannot close totsec dataset', 0)
!
      index = sfn2index(file_id,tlabel(3))
      call hdf_error( index,'cannot find id for totday', 0)
      sds_id = sfselect(file_id,index)
      call hdf_error(sds_id,'error getting index for totday', 0)
      ret = sfrdata(sds_id,start, stride,edge,totday)
      call hdf_error(ret,'error reading totday', 0)
      ret = sfendacc(sds_id)
      call hdf_error( ret, 'cannot close totday dataset', 0)
!
      index = sfn2index(file_id,tlabel(4))
      call hdf_error( index,'cannot find id for years', 0)
      sds_id = sfselect(file_id,index)
      call hdf_error(sds_id,'error getting index for years', 0)
      ret = sfrdata(sds_id,start, stride,edge,years)
      call hdf_error(ret,'error reading years', 0)
      ret = sfendacc(sds_id)
      call hdf_error( ret, 'cannot close years dataset', 0)
!
!  read kmt
!
#ifdef cyclic
      i1 = 2
      start(1) = 0
      start(2) = 0
      edge(1) = imtm2
      edge(2) = 1
#else
      i1 = 1
      start(1) = 0
      start(2) = 0
      edge(1) = imt
      edge(2) = jmt
#endif
      index = sfn2index(file_id,'KMT')
      call hdf_error( index,'cannot find id for KMT', 0)
      sds_id = sfselect(file_id,index)
      call hdf_error(sds_id,'error getting index for KMT', 0)
#ifdef cyclic
      do 210 j=1,jmt
        start(2) = j-1
        ret = sfrdata(sds_id,start,stride,edge,kmt(i1,j))
        call hdf_error(ret,'error reading KMT', 0)
 210  continue
#else
      ret = sfrdata(sds_id,start,stride,edge,kmt)
      call hdf_error(ret,'error reading KMT', 0)
#endif
      ret = sfendacc(sds_id)
      call hdf_error( ret, 'cannot close KMT dataset', 0)
!
!  read the 2-D datasets
!
      do 10 idvar = 1, NUM2D
      index = sfn2index(file_id,varid(idvar))
      call hdf_error( index,'cannot find id for 2-D dataset', 0)
      sds_id = sfselect(file_id,index)
      call hdf_error(sds_id,'error getting index for 2-D data', 0)
#ifdef cyclic
      do 220 j=1,jmt
      start(2) = j-1
      if(idvar.eq.1) then
        ret = sfrdata(sds_id,start,stride,edge,h0(2,j,np0))
      elseif(idvar.eq.2) then
        ret = sfrdata(sds_id,start,stride,edge,u0(2,j,np0))
      elseif(idvar.eq.3) then
        ret = sfrdata(sds_id,start,stride,edge,v0(2,j,np0))
      endif
      call hdf_error(ret,'error reading 2-D data', 0)
 220  continue
#else
      if(idvar.eq.1) then
        ret = sfrdata(sds_id,start,stride,edge,h0(1,1,np0))
      elseif(idvar.eq.2) then
        ret = sfrdata(sds_id,start,stride,edge,u0(1,1,np0))
      elseif(idvar.eq.3) then
        ret = sfrdata(sds_id,start,stride,edge,v0(1,1,np0))
      endif
      call hdf_error(ret,'error reading 2-D data', 0)
#endif
      ret = sfendacc(sds_id)
      call hdf_error( ret, 'cannot close 2-D dataset', 0)
  10  continue
!
!  read the 3-D dataset
!
#ifdef cyclic
      i1 = 2
      start(1) = 0
      start(2) = 0
      start(3) = 0
      edge(1) = km
      edge(2) = imtm2
      edge(3) = 1
#else
      i1 = 1
      edge(1) = km
      edge(2) = imt
      edge(3) = jmt
#endif
      do 20 idvar = NUM2D+1,NARCHV
        index = sfn2index(file_id,varid(idvar))
        call hdf_error( index,'cannot find id for 3-D dataset', 0)
        sds_id = sfselect(file_id,index)
        call hdf_error(sds_id,'error getting index for 3-D data', 0)
#ifdef cyclic
        do 230 j=1,jmt
        start(3) = j-1
        if(idvar.eq.4) then
          ret = sfrdata(sds_id,start,stride,edge,t(1,i1,j,1,np))
        elseif(idvar.eq.5) then
          ret = sfrdata(sds_id,start,stride,edge,t(1,i1,j,2,np))
        elseif(idvar.eq.6) then
          ret = sfrdata(sds_id,start,stride,edge,u(1,i1,j,np))
        elseif(idvar.eq.7) then
          ret = sfrdata(sds_id,start,stride,edge,v(1,i1,j,np))
        endif
 230    continue
#else
        if(idvar.eq.4) then
          ret = sfrdata(sds_id,start,stride,edge,t(1,1,1,1,np))
        elseif(idvar.eq.5) then
          ret = sfrdata(sds_id,start,stride,edge,t(1,1,1,2,np))
        elseif(idvar.eq.6) then
          ret = sfrdata(sds_id,start,stride,edge,u(1,1,1,np))
        elseif(idvar.eq.7) then
          ret = sfrdata(sds_id,start,stride,edge,v(1,1,1,np))
        endif
#endif
        call hdf_error(ret,'error reading 2-D data', 0)
        ret = sfendacc(sds_id)
        call hdf_error( ret, 'cannot close 2-D dataset', 0)
  20  continue
!
      ret = sfend(file_id)
      call hdf_error( ret,' error closing SD interface', 0)
#ifdef cyclic
!
!  cyclic boundaries
!
      do 410 j=1,jmt
        kmt(1,j)   = kmt(imtm1,j)
        kmt(imt,j) = kmt(2,j)
        h0(1,j,np0)   = h0(imtm1,j,np0)
        h0(imt,j,np0) = h0(2,j,np0)
        u0(1,j,np0)   = u0(imtm1,j,np0)
        u0(imt,j,np0) = u0(2,j,np0)
        v0(1,j,np0)   = v0(imtm1,j,np0)
        v0(imt,j,np0) = v0(2,j,np0)
        do 410 k=1,km
          t(k,1,j,1,np)   = t(k,imtm1,j,1,np)
          t(k,imt,j,1,np) = t(k,2,j,1,np)
          t(k,1,j,2,np)   = t(k,imtm1,j,2,np)
          t(k,imt,j,2,np) = t(k,2,j,2,np)
          u(k,1,j,np)   = u(k,imtm1,j,np)
          u(k,imt,j,np) = u(k,2,j,np)
          v(k,1,j,np)   = v(k,imtm1,j,np)
          v(k,imt,j,np) = v(k,2,j,np)
 410    continue
#endif
      return
      end
