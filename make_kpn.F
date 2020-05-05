      program make_kpn
c
c @(#) SCCS module: make_kpn.F  version: 1.8
c     Creation date: 08/23/95
c
c==================================================================
c
c  moma.pvm program to generate knp fields
c  D.J. Webb,  July 1994.
c
c==================================================================
c
c  Program generates a knp array for use with the moma.pvm program
c
c  Input parameters:
c    npx = number of partitions in x direction
c    npy = number of partitions in y direction
c    name_kmt = name of kmt file
c
#include "def_master.h"
#include "param.h"
      dimension  kmt(IMT_M,JMT_M),kpn(IMT_M,JMT_M)
      character*80 name,name_kmt
      character*8 string
      logical lexist
c
      print *,' Program make_kpn'
10    print *,' Maximum dimensions if kpn array with present version'
      print *,' of make_kpn is: ',IMT_M,' by ',JMT_M
      print *
15    print *,' Input name of kmt file'
      read  *,name_kmt
      inquire(file=name_kmt,exist=lexist,iostat=ios)
      if(.not.lexist.or.iostat.ne.0)then
        print *
        if(.not.lexist)then
          print *,' File does not exist'
        else
          print *,' Error opening file'
          print *,' iostat = ',iostat
        endif
        print *,' Try again:'
        print *
        goto 15
      endif
c
      open(11,file=name_kmt,status='old',iostat = iostat)
      if(iostat.ne.0)then
        print *
        print *,' Error opening kmt file'
        print *,' iostat = ',iostat
        print *,' Try again:'
        print *
        goto 15
      endif
      read(11,62)string,imt,jmt
      print *,' kmt file id = ',string
      print *,' Values of imt and jmt read from kmt file are:',imt,jmt

c
      if(imt.lt.2.or.imt.gt.IMT_M.or.jmt.lt.2.or.jmt.gt.JMT_M
     &            .or.string.ne.'kmt     ')then
        print *
        if(imt.lt.2)print *,' imt too small'
        if(jmt.lt.2)print *,' jmt too small'
        if(imt.gt.IMT_M)print *,' imt too large'
        if(jmt.gt.JMT_M)print *,' jmt too large'
        if(string.ne.'kmt     ')then
          print *,'file identifier not equal to "kmt     ".'
          print *,'file identifier equals       "',string,'".'
        endif
        print *,' Try again:'
        print *
        goto 10
      endif
c
      do 20 j=jmt,1,-1
        read(11,63)(kmt(i,j),i=1,imt)
  20  continue
      close(11)
c
  30  print *,' Input number of partitions in x and y directions:'
      read *,npx,npy
      if(npx.le.0.or.npy.le.0)then
        print *
        print *,' Out of range. Try again:'
        print *
        goto 30
      endif
c
c  allocate equally in i and j directions 
c   - but because the top and bottom rows need no work do not
c     include them when balancing workload
c   - for a rectangular ocean the same should be done in the x 
c     direction
c
      dx = float(imt)/float(npx)
      dy = float(jmt-2)/float(npy)
c
      do 40 j=jmt,1,-1
      do 40 i=1,imt
        if(kmt(i,j).eq.0)then
          kpn(i,j)=0
        else
          ii = 1+(i-0.5)/dx
          jj = max(1,min(npy,1+(j-1.5)/dy))
          kpn(i,j) = (jj-1)*npx+ii
        endif
  40  continue
c
  60  print *,' Input name of kpn output file:'
      read (*,*,err=60)name
      open(10,file=name,status='unknown',iostat=ios)
      if(ios.ne.0)then
        print *
        print *,' Error opening file ',name
        print *,' iostat = ',ios
        print *,' Try again:'
        print *
        goto 60
      endif
c
      string = 'kpn     '
      write(10,62)string,imt,jmt,npx*npy
      print 62,string,imt,jmt,npx*npy
      do 61 j=jmt,1,-1
      write(10,63)(kpn(i,j),i=1,imt)
  61  continue
  62  format(1x,a8,6i10)
  63  format(1x,30i3)
      stop
c
      end

        
