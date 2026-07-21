# SpeedCHEM Library Makefile
#
# Usage:
#   make                       # gfortran build -> gfortran/libSpeedCHEM64.a
#   make FC=mpif90             # MPI gfortran wrapper (treated as gfortran)
#   make FC=mpiifx             # Intel ifx via MPI wrapper -> ifx/libSpeedCHEM64.a
#   make COMPILER_TAG=ifx FC=mpiifx
#
# Output layout:
#   $(COMPILER_TAG)/
#     libSpeedCHEM64.a
#     mod/*.mod
#     build/*.o

FC ?= gfortran
AR ?= ar

# Decide COMPILER_TAG from FC when not given explicitly.
ifeq ($(origin COMPILER_TAG),undefined)
  ifneq (,$(findstring ifx,$(FC)))
    COMPILER_TAG := ifx
  else ifneq (,$(findstring ifort,$(FC)))
    COMPILER_TAG := ifx
  else
    COMPILER_TAG := gfortran
  endif
endif

SRCDIR    = src
OUTDIR    = $(COMPILER_TAG)
BUILDDIR  = $(OUTDIR)/build
MODULEDIR = $(OUTDIR)/mod
BUILDINFO = $(BUILDDIR)/build_flags.env
LIBNAME   = libSpeedCHEM64.a
TARGET    = $(OUTDIR)/$(LIBNAME)

# Compiler-specific flags.
ifeq ($(COMPILER_TAG),ifx)
  FFLAGS = -c \
    -extend-source 132 \
    -module $(MODULEDIR) \
    -O2
else
  FFLAGS = -c \
    -ffixed-line-length-none \
    -ffree-line-length-none \
    -fallow-argument-mismatch \
    -J$(MODULEDIR) \
    -O2
endif

# Canonical compile order (mirrors scripts/compile_gfort64opt.sh).
# working_precision must come first so .f files that `use working_precision`
# find the .mod. A strict serial sequence keeps intra-module dependencies
# satisfied — the Makefile should be invoked with -j1.
SRCS = \
  $(SRCDIR)/working_precision.f90 \
  $(SRCDIR)/SCutilities.f \
  $(SRCDIR)/SCsparse_definitions.f90 \
  $(SRCDIR)/SCsparse.f \
  $(SRCDIR)/sparse_MPI.f \
  $(SRCDIR)/dvode_f90_m.f90 \
  $(SRCDIR)/SCmodule.f \
  $(SRCDIR)/chemkin_module.f90 \
  $(SRCDIR)/SCconV.f90 \
  $(SRCDIR)/SCsetup.f \
  $(SRCDIR)/SCallocate.f \
  $(SRCDIR)/SCcklink.f90 \
  $(SRCDIR)/gam.f \
  $(SRCDIR)/gamsub.f \
  $(SRCDIR)/opkdmain.f \
  $(SRCDIR)/opkda1.f \
  $(SRCDIR)/opkda2.f \
  $(SRCDIR)/ddaspk.f \
  $(SRCDIR)/rodas.f \
  $(SRCDIR)/vode.f \
  $(SRCDIR)/MEBDFSO.f \
  $(SRCDIR)/rowmap.f \
  $(SRCDIR)/radau5.f \
  $(SRCDIR)/radaua.f \
  $(SRCDIR)/radau_sparse.f90 \
  $(SRCDIR)/chemistry_input.f90 \
  $(SRCDIR)/SCbroadcast.f90

OBJS = $(patsubst $(SRCDIR)/%.f90,$(BUILDDIR)/%.o, \
         $(patsubst $(SRCDIR)/%.f,$(BUILDDIR)/%.o,$(SRCS)))

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
	@printf 'FC=%s\nAR=%s\nCOMPILER_TAG=%s\nFFLAGS=%s\nTARGET=%s\n' \
		"$(FC)" "$(AR)" "$(COMPILER_TAG)" "$(strip $(FFLAGS))" "$(TARGET)" > $@

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

# コンパイル規則（.f用）
$(BUILDDIR)/%.o: $(SRCDIR)/%.f | $(BUILDDIR) $(MODULEDIR)
	$(FC) $(FFLAGS) -o $@ $<

# 掃除用（現在の COMPILER_TAG のみ）
clean:
	rm -rf $(OUTDIR)

# 全コンパイラ成果物を掃除
distclean:
	rm -rf gfortran ifx build

.PHONY: all clean distclean
