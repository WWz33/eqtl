# eqtl — cis/trans eQTL (C++17)
CXX ?= g++
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -Wno-unused-parameter -fopenmp
CPPFLAGS += -Iinclude -Ithird_party/eigen

USE_SYSTEM_HTS ?= 1

ifeq ($(USE_SYSTEM_HTS),1)
  HTS_CFLAGS := $(shell pkg-config --cflags htslib 2>/dev/null)
  HTS_LIBS := $(shell pkg-config --libs htslib 2>/dev/null)
  ifeq ($(HTS_LIBS),)
    HTS_CFLAGS := -I/usr/local/include
    HTS_LIBS := -L/usr/local/lib -lhts -lz -lbz2 -llzma -lcurl -lcrypto -lpthread -lm
  endif
  HTS_REQ :=
else
  HTS_SRC := third_party/htslib
  HTS_LIB := $(HTS_SRC)/libhts.a
  HTS_CFLAGS := -I$(HTS_SRC)
  -include $(HTS_SRC)/htslib_static.mk
  HTSLIB_static_LIBS ?= -lz -lbz2 -llzma -lcurl -lcrypto -ldeflate -lpthread -lm
  HTS_LIBS := $(HTS_LIB) $(HTSLIB_static_LIBS)
  HTS_REQ := $(HTS_LIB)
endif

CPPFLAGS += $(HTS_CFLAGS)
LDFLAGS += $(HTS_LIBS) -fopenmp -lm

SRC := \
  src/main.cpp \
  src/options.cpp \
  src/util.cpp \
  src/pheno.cpp \
  src/annot.cpp \
  src/vcf_session.cpp \
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

.PHONY: all clean smoke

all: $(BIN)

$(BIN): $(OBJ) $(HTS_REQ)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJ) $(LDFLAGS)

src/%.o: src/%.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $<

clean:
	rm -f $(OBJ) $(BIN)

smoke: $(BIN)
	./scripts/run_smoke.sh
