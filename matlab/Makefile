SRVF_LIBDIR = ../build/src
SRVF_LIBS   = -lsrvf
SRVF_INCDIR = ../include

MEXSUFFIX  = mex
MATLABHOME = /usr/local/MATLAB/R2016b/
MEX        = /usr/local/MATLAB/R2016b/bin/mex
CXX        = g++

LIBS      = -L$(SRVF_LIBDIR) $(SRVF_LIBS) -lm
INCLUDE   = -I$(SRVF_INCDIR) -I$(MATLABHOME)/extern/include
CFLAGS    = -fPIC
MEXFLAGS  = -cxx CC='$(CXX)' CXX='$(CXX)' LD='$(CXX)'

MEXFILES= \
	srvf_fa_groupwise_reparam.mex \
	srvf_fa_karcher_mean.mex \
	srvf_fa_optimal_reparam.mex \
	srvf_karcher_mean.mex \
	srvf_optimal_reparam_mex.mex \
	srvf_pmatch_find_matches.mex

%.$(MEXSUFFIX): %.o
	$(MEX) $(MEXFLAGS) $(INCLUDE) $(LIBS) $< -lgsl -lgslcblas

%.o: %.cc
	$(CXX) $(INCLUDE) $(CFLAGS) -c $<


all: $(MEXFILES)

clean:
	rm -f $(MEXFILES)
	rm -f *.o
