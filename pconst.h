#ifdef hcomments
!
!======================================================================
!
!    File:  pconst.h
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
!     rules for parameter constants
!
!     use prefix of "c" for whole real numbers (ie: c57 for 57.0)
!     use "m" after prefix for negative values (ie: cm7 for -7.0)
!     use prefix "p" for non repeating fractions (ie: p5 for 0.5)
!     use prefix "r" for reciprocals (ie: r3 for 1/3.0)
!     combine with "e" for scientific notation,
!       (ie: c5e4 for 5.0e4, c1em10 for 1.0e-10)
!
#endif
      parameter (c0=0.0, c1=1.0, c2=2.0, c4=4.0, c5=5.0, c8=8.0)
      parameter (c16=16.0, c360=360.0, c100=100.0)
      parameter (p25=0.25, p5=0.5, p75=0.75, p125=0.125)
      parameter (epsln=1.0e-20)
      parameter (c24=24.0, c60=60.0, c1440=1440.0)
      parameter (r24=c1/c24, r60=c1/c60, r1440=c1/c1440)
      parameter (secday=c1/(c60*c1440))
#ifdef hcomments
!
!     stdin, stdout, and stderr set the Fortran units for
!     standard input, standard output, and standard error messages
!     (used in read and write statements)
!
#endif
      integer stdin, stdout, stderr
      parameter (stdin = 5, stdout = 6, stderr = 6)
#ifdef hcomments
!
!     define length of words on different machines
!
#endif
#ifdef cray-t3d
# define LOGICAL_N logical*8
#else
# define LOGICAL_N logical*4
#endif

