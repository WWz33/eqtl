# eqtl — cis/trans eQTL (C++17)
# Vendored: third_party/htslib @ 1.24, third_party/gffsub, third_party/eigen,
#           third_party/OpenBLAS @ 0.3.28 (static, serial BLAS/LAPACKE for Eigen)
CXX ?= g++
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -Wno-unused-parameter -fopenmp
CPPFLAGS += -Iinclude -Ithird_party/eigen -Ithird_party/gffsub/src

# ---- OpenBLAS (vendored static; Eigen BLAS/LAPACKE backend) ----
OPENBLAS_SRC := third_party/OpenBLAS
OPENBLAS_LIB := $(OPENBLAS_SRC)/libopenblas.a
# serial BLAS: app uses OpenMP at SNP/gene level
OPENBLAS_MAKE_FLAGS := NO_SHARED=1 USE_THREAD=0 USE_OPENMP=0 BINARY=64 DYNAMIC_ARCH=1
CPPFLAGS += -DEIGEN_USE_BLAS -DEIGEN_USE_LAPACKE -I$(OPENBLAS_SRC)
# gfortran runtime required by static OpenBLAS LAPACK objects
LDFLAGS += $(OPENBLAS_LIB) -lgfortran -lpthread

# ---- htslib (default: vendored static) ----
USE_SYSTEM_HTS ?= 0
HTS_SRC := third_party/htslib
HTS_LIB := $(HTS_SRC)/libhts.a

ifeq ($(USE_SYSTEM_HTS),1)
  HTS_CFLAGS := $(shell pkg-config --cflags htslib 2>/dev/null)
  HTS_LIBS := $(shell pkg-config --libs htslib 2>/dev/null)
  ifeq ($(HTS_LIBS),)
    HTS_CFLAGS := -I/usr/local/include
    HTS_LIBS := -L/usr/local/lib -lhts -lz -lbz2 -llzma -lcurl -lcrypto -lpthread -lm
  endif
  HTS_REQ :=
else
  HTS_CFLAGS := -I$(HTS_SRC)
  -include $(HTS_SRC)/htslib_static.mk
  HTSLIB_static_LIBS ?= -lz -lbz2 -llzma -lcurl -lcrypto -ldeflate -lpthread -lm
  HTS_LIBS := $(HTS_LIB) $(HTSLIB_static_LIBS)
  HTS_REQ := $(HTS_LIB)
endif

# ---- gffsub (library objects, no CLI) ----
GFFSUB_SRC := third_party/gffsub
GFFSUB_CPP := \
  $(GFFSUB_SRC)/src/gff3_parser.cpp \
  $(GFFSUB_SRC)/src/gtf_parser.cpp \
  $(GFFSUB_SRC)/src/attributes.cpp
GFFSUB_OBJ := $(GFFSUB_CPP:.cpp=.o)

CPPFLAGS += $(HTS_CFLAGS)
LDFLAGS += $(HTS_LIBS) -fopenmp -lm

SRC := \
  src/main.cpp \
  src/options.cpp \
  src/util.cpp \
  src/pheno.cpp \
  src/annot.cpp \
  src/vcf_session.cpp \
  src/plink_bed.cpp \
  src/grm.cpp \
  src/model_lm.cpp \
  src/model_lmm.cpp \
  src/model_glm.cpp \
  src/model_glmm.cpp \
  src/stats_extra.cpp \
  src/scan.cpp \
  src/output.cpp

OBJ := $(SRC:.cpp=.o)
BIN := eqtl

.PHONY: all clean smoke htslib openblas

all: $(BIN)

htslib: $(HTS_LIB)
openblas: $(OPENBLAS_LIB)

$(OPENBLAS_LIB):
	@if [ ! -f $(OPENBLAS_SRC)/Makefile ]; then \
	  echo "[E] $(OPENBLAS_SRC) missing. Run: git submodule update --init --recursive"; \
	  exit 1; \
	fi
	@echo "[I] building vendored OpenBLAS (static, serial) ..."
	$(MAKE) -C $(OPENBLAS_SRC) $(OPENBLAS_MAKE_FLAGS) -j$$(nproc 2>/dev/null || echo 4)

$(HTS_LIB) $(HTS_SRC)/htslib_static.mk:
	@if [ ! -e $(HTS_SRC)/htslib/vcf.h ] && [ ! -e $(HTS_SRC)/vcf.h ]; then \
	  echo "[E] $(HTS_SRC) incomplete. Run: git submodule update --init --recursive"; \
	  exit 1; \
	fi
	@if [ ! -f $(HTS_SRC)/config.mk ]; then \
	  echo "[I] configuring vendored htslib ..."; \
	  cd $(HTS_SRC) && \
	    git submodule update --init --recursive htscodecs 2>/dev/null || true; \
	    if [ ! -f config.guess ]; then autoreconf -i || true; fi && \
	    if [ ! -x configure ]; then autoheader && autoconf; fi && \
	    ./configure; \
	fi
	$(MAKE) -C $(HTS_SRC) -j$$(nproc 2>/dev/null || echo 4) lib-static htslib_static.mk

$(BIN): $(OBJ) $(GFFSUB_OBJ) $(HTS_REQ) $(OPENBLAS_LIB)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJ) $(GFFSUB_OBJ) $(LDFLAGS)

src/%.o: src/%.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $<

$(GFFSUB_SRC)/src/%.o: $(GFFSUB_SRC)/src/%.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $<

clean:
	rm -f $(OBJ) $(GFFSUB_OBJ) $(BIN)

smoke: $(BIN)
	./scripts/run_smoke.sh
