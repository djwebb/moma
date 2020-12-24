#ifdef hcomments
!
!======================================================================
!
!    File:  iounit.h
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
!     i/o units & related variables
!
!     iorest  = unit for archive files (restarts)
!     iosnap  = unit for snapshot files
!     fnrest  = name of archive file used to start run
!     iokmt   = unit for kmt file (depths)
!     fnkmt   = name of kmt file used
!     ftrest  = file type of archive file used to start run
!               'std' := restart, 'hdf' := hdf file
!     ftarch  = file type used for archive files
!               'std' := restart, 'hdf' := hdf file
!     ftsnap  = file type for snapshot files
!               'std' := ascout, 'hdf' := hdf file
!   note: ansi needs character variables in a separate common block
!
#endif
      character    fnrest*80,fnkmt*80,ftrest*4,ftsnap*4,ftarch*4
      common /iouniti/ iorest,iokmt
      common /iounitc/ fnrest,fnkmt,ftrest,ftsnap,ftarch

