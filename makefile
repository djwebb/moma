#
#=======================================================================
# @(#) SCCS module: makefile  version: 1.1
#     Creation date: 10/13/97
#=======================================================================
#
#     Moma.pvm make file for 64-bit SGI workstation 
#              running IRIX64 and mpich version of MPI
#
#   Needs the following environment variables to be set and
#   the PVM and HDF include files and libraries to be installed
#   in the corresponding locations:
#
#     ARCH        - architecture type
#     INSTALL_DIR - Root directory for MPI
#     HDF_ROOT    - Root directory for HDF

ARCH        = IRIX64
INSTALL_DIR = /working/jrd/occam/MPI
HDF_ROOT    = /working/jrd/occam/HDF4

#

F77         = mpif77 -64 -mips4 -O2
FLINKER     = mpif77 -64 -mips4 -O2
GPP         = cpp -P  
OPTFLAGS    = 

#
#   With HDF
#
INCLUDE_DIR = -I$(INSTALL_DIR)/include -I$(HDF_ROOT)/include/$(ARCH)
LIB_PATH    = -L/working/jrd/occam/MPI/lib/IRIX64/ch_p4 -L$(HDF_ROOT)/lib/$(ARCH)
LIB_LIST    = -lmfhdf -ldf -lz -lmpi -lm
#
#  Without HDF
#
#INCLUDE_DIR = -I$(INSTALL_DIR)/include 
#LIB_PATH    = -L$(INSTALL_DIR)/lib/$(ARCH)/ch_p4 
#LIB_LIST    = -lmpi -lm

#

FFLAGS =  -cckr -64 -mips4 -woff 1552,1174 $(INCLUDE_DIR) $(OPTFLAGS)

EXECS  =  master slave

# define '.F' suffix for the SGI system

.SUFFIXES:
.SUFFIXES: .o .F .h 

MOFILES=master.o   checks.o   docmnt.o   hdf.o      \
        init_kmt.o init_kpn.o m_arch.o   m_flx.o    \
        m_grids.o  m_msg.o    m_snap.o   setkmt.o   \
        topog.o     

SOFILES=slave.o    addv.o     bcest.o    bound.o    \
        clinic.o   diag.o     frees.o    matrix.o   \
        ocn1st.o   pdiag.o    s_arch.o   s_grids.o  \
        s_msg.o    s_snap.o   setvbc.o   setp.o     \
        state.o    step.o     theta0.o   tracer.o   

BOFILES=sysmsg.o   tmngr.o    tmnset.o   indp.o

HFILES=archive.h cdiag.h  chmix.h   coord.h   ctmngr.h   \
       cvbc.h    cvmix.h  dncoef.h  frees.h   grdvar.h   \
       iounit.h  levind.h mesdta.h  moddata.h param.h    \
       pconst.h  scalar.h slabs.h   snaps.h   switch.h   \
       thick.h   timelv.h varinfo.h versno.h  hdf.h

# Add -DMasterslave if necessary

DOPTNS= -Dde_checkbd -Dcyclic_master -Dncsa_hdf      \
        -DREAL_8     -Dinline        -Dpresetp 

#suppress SCCS retrieval

.SCCS_GET:

#

default:	

all:	$(EXECS)

masterslave: masterslave.o $(MOFILES) $(SOFILES) $(BOFILES) $(INSTALL_DIR)/include/mpi.h
	$(FLINKER) $(OPTFLAGS) -o masterslave masterslave.o $(OFILES) $(LIB_PATH) $(LIB_LIST)

master:	$(MOFILES) $(BOFILES)	
	$(FLINKER) $(OPTFLAGS) $(MOFILES) $(BOFILES) -o master $(LIB_PATH) $(LIB_LIST)

slave : $(SOFILES) $(BOFILES)	
	$(FLINKER) $(OPTFLAGS) $(SOFILES) $(BOFILES) -o slave  $(LIB_PATH) $(LIB_LIST)

make_kmt :
	$(F77) -o make_kmt make_kmt.F

make_kpn :
	$(F77) -o make_kpn make_kpn.F

# SGI system needs to call preprocessor separately

.F.o:   $(HFILES)
	$(GPP) $(DOPTNS) $(INCLUDE_DIR) $< > $*.f
	$(F77) $(FFLAGS)                  -c $*.f

clean:
	/bin/rm -f *.f *% *~ PI* 

cleanest:
	/bin/rm -f *.o *.f *% *~ PI* $(EXECS)

