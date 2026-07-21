# SpeedCHEM Library Makefile (Intel ifx only)
#
# Usage:
#   make                       # ifx build -> ifx/libSpeedCHEM64.a
#   make FC=mpiifx             # explicit (default)
#
# Output layout:
#   ifx/
#     libSpeedCHEM64.a
#     mod/*.mod
#     build/*.o

# ifx-only. Default to the Intel MPI wrapper. A command-line FC (e.g. `make FC=ifx`)
# wins; make's built-in default (FC=f77) or an inherited environment FC is ignored.
ifneq ($(origin FC),command line)
  FC := mpiifx
endif
AR ?= ar

SRCDIR    = src
OUTDIR    = ifx
BUILDDIR  = $(OUTDIR)/build
MODULEDIR = $(OUTDIR)/mod
BUILDINFO = $(BUILDDIR)/build_flags.env
LIBNAME   = libSpeedCHEM64.a
TARGET    = $(OUTDIR)/$(LIBNAME)

FFLAGS = -c \
  -extend-source 132 \
  -module $(MODULEDIR) \
  -O2

# Canonical compile order (mirrors scripts/ifx.sh).
# working_precision must come first so modules that `use working_precision`
# find the .mod. A strict serial sequence keeps intra-module dependencies
# satisfied — the Makefile should be invoked with -j1.
SRCS = \
  $(SRCDIR)/working_precision.f90 \
  $(SRCDIR)/SCutilities.f90 \
  $(SRCDIR)/SCsparse_definitions.f90 \
  $(SRCDIR)/SCsparse.f90 \
  $(SRCDIR)/sparse_MPI.f90 \
  $(SRCDIR)/dvode_f90_m.f90 \
  $(SRCDIR)/SCmodule.f90 \
  $(SRCDIR)/chemkin_module.f90 \
  $(SRCDIR)/SCconV.f90 \
  $(SRCDIR)/SCsetup.f90 \
  $(SRCDIR)/SCallocate.f90 \
  $(SRCDIR)/SCcklink.f90 \
  $(SRCDIR)/gam.f90 \
  $(SRCDIR)/gamsub.f90 \
  $(SRCDIR)/opkdmain.f90 \
  $(SRCDIR)/opkda1.f90 \
  $(SRCDIR)/opkda2.f90 \
  $(SRCDIR)/ddaspk.f90 \
  $(SRCDIR)/rodas.f90 \
  $(SRCDIR)/vode.f90 \
  $(SRCDIR)/MEBDFSO.f90 \
  $(SRCDIR)/rowmap.f90 \
  $(SRCDIR)/radau5.f90 \
  $(SRCDIR)/radaua.f90 \
  $(SRCDIR)/radau_sparse.f90 \
  $(SRCDIR)/chemistry_input.f90 \
  $(SRCDIR)/SCbroadcast.f90

OBJS = $(patsubst $(SRCDIR)/%.f90,$(BUILDDIR)/%.o,$(SRCS))

# デフォルトターゲット（$(eval ...) より先に宣言して .DEFAULT_GOAL を確保）
.DEFAULT_GOAL := all
all: $(BUILDINFO) $(TARGET)

# Force strict serial ordering via a chain of prerequisites so -j >1 still works.
PREV :=
define ORDER_rule
$1: $(PREV)
PREV := $1
endef
$(foreach o,$(OBJS),$(eval $(call ORDER_rule,$o)))

$(BUILDINFO): | $(BUILDDIR) $(MODULEDIR)
	@printf 'FC=%s\nAR=%s\nFFLAGS=%s\nTARGET=%s\n' \
		"$(FC)" "$(AR)" "$(strip $(FFLAGS))" "$(TARGET)" > $@

# ライブラリ作成
$(TARGET): $(BUILDINFO) $(BUILDDIR) $(MODULEDIR) $(OBJS)
	$(AR) rcs $@ $(OBJS)
	@echo "Library created: $@"

# ビルドディレクトリ作成
$(BUILDDIR):
	mkdir -p $(BUILDDIR)

$(MODULEDIR):
	mkdir -p $(MODULEDIR)

# コンパイル規則（.f90用）
$(BUILDDIR)/%.o: $(SRCDIR)/%.f90 | $(BUILDDIR) $(MODULEDIR)
	$(FC) $(FFLAGS) -o $@ $<

# 掃除用
clean:
	rm -rf $(OUTDIR)

.PHONY: all clean
