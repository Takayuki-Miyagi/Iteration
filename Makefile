#--------------------------------------------------
# Make file for Optimizer code
#
# When we use ifort,
# -heap-arrays option is needed for built in transpose function with
#  large dimension matrix
#--------------------------------------------------
TARGET=test
MODDIR = mod
Host= $(shell if hostname|grep -q apt1; then echo apt; \
  elif hostname|grep -q oak; then echo oak; \
  elif hostname|grep -q cedar; then echo cedar; \
  else echo other; fi)
HOST=$(strip $(Host))
DEBUG_MODE=on

OS = Linux
ifneq (,$(findstring arwin,$(shell uname)))
  OS = OSX
endif
$(info HOST is $(HOST).)
$(info OS is $(OS).)
$(info Debug mode is $(DEBUG_MODE).)

FDEP=
FC=
LFLAGS=  # library
FFLAGS=  # option
DFLAGS=  # debug
LINT=    # 8-byte integer

#--------------------------------------------------
# Default Parameters
#--------------------------------------------------
ifeq ($(strip $(HOST)),other)
  FDEP=makedepf90
  FC=gfortran
  ifeq ($(OS), OSX)
    LFLAGS+= -I/usr/local/include -L/usr/local/lib
  endif
  LFLAGS+= -lblas -llapack -lgsl -lz
  FFLAGS=-O3
  CFLAGS=-O3
  FFLAGS+= -fopenmp
  FFLAGS+= -Dsingle_precision_three_body_file
  #FFLAGS+= -ff2c # for dot product (LinAlgf90)
  ifeq ($(DEBUG_MODE),on)
    DFLAGS+=-Wall -pedantic -fbounds-check -O -Wuninitialized -fbacktrace
    #FDFLAGS+=-ffpe-trap=invalid,zero,overflow # Note: gsl larguerre signal
    ifneq ($(OS), OSX)
      DFLAGS+= -pg -g
    endif
  endif
endif

ifeq ($(DEBUG_MODE),on)
  DFLAGS+=-DIterationDebug
endif



#--------------------------------------------------
# Source Files
#--------------------------------------------------

SRCDIR = src
SRCDIR_LinAlg = submodule/LinAlgf90/src
SRCDIR_Opt = main
DEPDIR = .
OBJDIR = obj

SRCS=
OBJS=
MODS=

SRCC:=$(wildcard $(SRCDIR)/*.c)
SRCF77:=$(wildcard $(SRCDIR)/*.f)
SRCF90:=$(wildcard $(SRCDIR)/*.f90)
SRCF95:=$(wildcard $(SRCDIR)/*.F90)
OBJC:=$(addprefix $(OBJDIR)/, $(patsubst %c, %o, $(notdir $(SRCC))))
OBJF77:=$(addprefix $(OBJDIR)/, $(patsubst %f, %o, $(notdir $(SRCF77))))
OBJF90:=$(addprefix $(OBJDIR)/, $(patsubst %f90, %o, $(notdir $(SRCF90))))
OBJF95:=$(addprefix $(OBJDIR)/, $(patsubst %F90, %o, $(notdir $(SRCF95))))
SRCS= $(SRCC) $(SRCF77) $(SRCF90) $(SRCF95)
OBJS= $(OBJC) $(OBJF77) $(OBJF90) $(OBJF95)

SRCC_LinAlg:=$(wildcard $(SRCDIR_LinAlg)/*.c)
SRCF77_LinAlg:=$(wildcard $(SRCDIR_LinAlg)/*.f)
SRCF90_LinAlg:=$(wildcard $(SRCDIR_LinAlg)/*.f90)
SRCF95_LinAlg:=$(wildcard $(SRCDIR_LinAlg)/*.F90)
OBJC_LinAlg:=$(addprefix $(OBJDIR)/, $(patsubst %c, %o, $(notdir $(SRCC_LinAlg))))
OBJF77_LinAlg:=$(addprefix $(OBJDIR)/, $(patsubst %f, %o, $(notdir $(SRCF77_LinAlg))))
OBJF90_LinAlg:=$(addprefix $(OBJDIR)/, $(patsubst %f90, %o, $(notdir $(SRCF90_LinAlg))))
OBJF95_LinAlg:=$(addprefix $(OBJDIR)/, $(patsubst %F90, %o, $(notdir $(SRCF95_LinAlg))))
SRCS_LinAlg= $(SRCC_LinAlg) $(SRCF77_LinAlg) $(SRCF90_LinAlg) $(SRCF95_LinAlg)
OBJS_LinAlg= $(OBJC_LinAlg) $(OBJF77_LinAlg) $(OBJF90_LinAlg) $(OBJF95_LinAlg)

SRCC_Opt:=$(wildcard $(SRCDIR_Opt)/*.c)
SRCF77_Opt:=$(wildcard $(SRCDIR_Opt)/*.f)
SRCF90_Opt:=$(wildcard $(SRCDIR_Opt)/*.f90)
SRCF95_Opt:=$(wildcard $(SRCDIR_Opt)/*.F90)
OBJC_Opt:=$(addprefix $(OBJDIR)/, $(patsubst %c, %o, $(notdir $(SRCC_Opt))))
OBJF77_Opt:=$(addprefix $(OBJDIR)/, $(patsubst %f, %o, $(notdir $(SRCF77_Opt))))
OBJF90_Opt:=$(addprefix $(OBJDIR)/, $(patsubst %f90, %o, $(notdir $(SRCF90_Opt))))
OBJF95_Opt:=$(addprefix $(OBJDIR)/, $(patsubst %F90, %o, $(notdir $(SRCF95_Opt))))
SRCS_Opt= $(SRCC_Opt) $(SRCF77_Opt) $(SRCF90_Opt) $(SRCF95_Opt)
OBJS_Opt= $(OBJC_Opt) $(OBJF77_Opt) $(OBJF90_Opt) $(OBJF95_Opt)


SRCS_ALL = $(SRCS) $(SRCS_LinAlg) $(SRCS_Opt)
OBJS_ALL = $(OBJS) $(OBJS_LinAlg) $(OBJS_Opt)

MODOUT=
ifeq ($(strip $(HOST)),other)
  MODOUT=-J$(MODDIR)
endif

#--------------------------------------------------
# Rules
#--------------------------------------------------
all: dirs $(TARGET)
$(TARGET): $(OBJS_ALL)
	$(FC) $(FFLAGS) $(DFLAGS) -o $(TARGET).exe $^ $(LFLAGS)
	if test -d $(EXEDIR); then \
		: ; \
	else \
		mkdir -p $(EXEDIR); \
	fi

$(OBJDIR)/%.o:$(SRCDIR)/%.c
	$(CC) $(CFLAGS) -o $@ -c $<
$(OBJDIR)/%.o:$(SRCDIR)/%.f
	$(FC) $(FFLAGS) $(MODOUT) -o $@ -c $<
$(OBJDIR)/%.o:$(SRCDIR)/%.f90
	$(FC) $(FFLAGS) $(DFLAGS) $(MODOUT) -o $@ -c $<
$(OBJDIR)/%.o:$(SRCDIR)/%.F90
	$(FC) $(FFLAGS) $(DFLAGS) $(MODOUT) -o $@ -c $<

$(OBJDIR)/%.o:$(SRCDIR_LinAlg)/%.c
	$(CC) $(CFLAGS) -o $@ -c $<
$(OBJDIR)/%.o:$(SRCDIR_LinAlg)/%.f
	$(FC) $(FFLAGS) $(MODOUT) -o $@ -c $<
$(OBJDIR)/%.o:$(SRCDIR_LinAlg)/%.f90
	$(FC) $(FFLAGS) $(DFLAGS) $(MODOUT) -o $@ -c $<
$(OBJDIR)/%.o:$(SRCDIR_LinAlg)/%.F90
	$(FC) $(FFLAGS) $(DFLAGS) $(MODOUT) -o $@ -c $<

$(OBJDIR)/%.o:$(SRCDIR_Opt)/%.c
	$(CC) $(CFLAGS) -o $@ -c $<
$(OBJDIR)/%.o:$(SRCDIR_Opt)/%.f
	$(FC) $(FFLAGS) $(MODOUT) -o $@ -c $<
$(OBJDIR)/%.o:$(SRCDIR_Opt)/%.f90
	$(FC) $(FFLAGS) $(DFLAGS) $(MODOUT) -o $@ -c $<
$(OBJDIR)/%.o:$(SRCDIR_Opt)/%.F90
	$(FC) $(FFLAGS) $(DFLAGS) $(MODOUT) -o $@ -c $<

dep:
	$(FDEP) $(SRCS_ALL) -b $(OBJDIR)/ > $(DEPDIR)/makefile.d

dirs:
	if test -d $(OBJDIR); then \
		: ; \
	else \
		mkdir $(OBJDIR); \
	fi
	if test -d $(MODDIR); then \
		: ; \
	else \
		mkdir $(MODDIR); \
	fi

clean:
	rm -f $(TARGET).exe
	rm -f $(MODDIR)/*.mod
	rm -f $(OBJS_ALL)

#--------------------------------------------------
-include $(wildcard $(DEPDIR)/*.d)
