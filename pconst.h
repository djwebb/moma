c====================== include file "pconst.h" ========================
#ifdef hcomments
c
c @(#) SCCS module: pconst.h  version: 1.1
c     Creation date: 10/13/97
c
c-----------------------------------------------------------------------
c     Rules for parameter constants
c
c     use prefix of "c" for whole real numbers (ie: c57 for 57.0)
c     use "m" after prefix for negative values (ie: cm7 for -7.0)
c     use prefix "p" for non repeating fractions (ie: p5 for 0.5)
c     use prefix "r" for reciprocals (ie: r3 for 1/3.0)
c     combine with "e" for scientific notation, 
c       (ie: c5e4 for 5.0e4, c1em10 for 1.0e-10)
c-----------------------------------------------------------------------
c
#endif
      parameter (c0=0.0, c1=1.0, c2=2.0, c4=4.0, c5=5.0, c8=8.0)
      parameter (c16=16.0, c360=360.0, c100=100.0)
      parameter (p25=0.25, p5=0.5, p75=0.75, p125=0.125)
      parameter (epsln=1.0e-20)
      parameter (c24=24.0, c60=60.0, c1440=1440.0)
      parameter (r24=c1/c24, r60=c1/c60, r1440=c1/c1440)
      parameter (secday=c1/(c60*c1440))
#ifdef hcomments
c
c     stdin, stdout, and stderr set the Fortran units for
c     standard input, standard output, and standard error messages
c     (used in read and write statements)
c
#endif
      integer stdin, stdout, stderr
      parameter (stdin = 5, stdout = 6, stderr = 6)
