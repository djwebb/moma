      program make_kmt
c
c @(#) SCCS module: make_kmt.F  version: 1.10
c     Creation date: 08/09/96
c
c==================================================================
c
c  moma.pvm program to generate kmt fields
c  D.J. Webb,  July 1994.
c
c==================================================================
c
c  Program generates a kmt array for use with the moma.pvm program
c
c  Input parameters:
c    imt = first dimension of generated kmt array
c    imt = second dimension of generated kmt array
c    km  = maximum level number of generated array
c    iopt= type of model ocean to generate
c          1 = box type ocean of depth equal to km
c          2 = channel type ocean of depth equal to km
c          3 = ocean region based on mom test global ocean
c
#define Master
#include "param.h"
#include "coord.h"
#include "grdvar.h"
#include "levind.h"
      dimension  kmt(IMT_M,JMT_M)
      character*80 name
      character*8  string
c
      print *,' Program make_kmt'
10    print *,' Maximum dimensions if kmt array with present version'
      print *,' of make_kmt is: ',IMT_M,' by ',JMT_M
      print *,' Number of levels km when make_kmt compiled was ',km
      print *
      print *,' Input imt, jmt, km:'
      read  *,imt,jmt,kmm
      if(imt.lt.2.or.imt.gt.IMT_M.or.jmt.lt.2.or.jmt.gt.JMT_M)then
        print *
        print *,' Out of range. Try again:'
        print *
        goto 10
      endif
c
20    print *,' Input shape option. (1=box ocean, 2= channel ocean,',
     &        ' 3=mom test global ocean):'
      read *,iopt
      if(iopt.lt.1.or.iopt.gt.3)then
        print *
        print *,' Out of range. Try again:'
        print *
        goto 20
      endif
c
      do 30 j=1,jmt
      do 30 i=1,imt
        kmt(i,j)=kmm
30    continue
      if(iopt.eq.1)then
        do 41 i=1,imt
          kmt(i,1)=0
          kmt(i,jmt)=0
  41    continue
        do 42 j=1,jmt
          kmt(1,j)=0
          kmt(imt,j)=0
  42    continue
      elseif(iopt.eq.2)then
        do 51 i=1,imt
          kmt(i,1)=0
          kmt(i,jmt)=0
  51    continue
      elseif(iopt.eq.3)then
        print *
        print *,'  stlon and stlat refer to edge of model region.'
        print *,' Input stlon, stlat, dxdeg, dydeg:'
        read *,stlon, stlat, dxdeg, dydeg
c      do 100 i=1,IMT_M
c        xt(i) = stlon+(i-p5)*dxdeg
c100   continue
c      do 200 j=1,JMT_M
c        yt(j) = stlat+(j-p5)*dydeg
c200   continue
      pi     = c4*atan(c1)
      radian = c360/(c2*pi)
      omega  = pi/43082.0
      grav =980.6
      radius =6370.e5

      dz( 1) = 30.00e2
      dz( 2) = 46.15e2
      dz( 3) = 68.93e2
      dz( 4) = 99.93e2
      dz( 5) = 140.63e2
      dz( 6) = 192.11e2
      dz( 7) = 254.76e2
      dz( 8) = 327.95e2
      dz( 9) = 409.81e2
      dz(10) = 497.11e2
      dz(11) = 585.36e2
      dz(12) = 669.09e2
      dz(13) = 742.41e2
      dz(14) = 799.65e2
      dz(15) = 836.10e2
      call m_grids
      call topog
      endif
c
  60  print *,' Input name of kmt output file:'
      read (*,16,err=60)name
  16  format(a)
      open(10,file=name,status='unknown',err=80)
      string = 'kmt     '
      write(10,62)string,imt,jmt,km
      print 62,string,imt,jmt,km
      do 61 j=jmt,1,-1
      write(10,63)(kmt(i,j),i=1,imt)
  61  continue
  62  format(1x,a8,6i10)
  63  format(1x,30i3)
      stop
c
  80  print *
      print *,' Error opening file ',name
      print *,' Try again:'
      print *
      goto 60
      end

#include "topog.F"
#include "m_grids.F"
#include "setkmt.F"
#include "indp.F"
