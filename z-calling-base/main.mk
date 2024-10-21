ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../$(OSTYPE)-$(MACHINETYPE)/bin
endif

TARGET   := z-calling-base
SOURCES  := build_mod_bam.cpp main.cpp make_read_features.cpp necat_info.cpp

SRC_INCDIRS  := . ${LIBTORCH}/include ${LIBTORCH}/include/torch/csrc/api/include

TGT_LDFLAGS := -L${TARGET_DIR}
TGT_LDLIBS  := -lnecat
TGT_PREREQS := libnecat.a

SUBMAKEFILES :=
