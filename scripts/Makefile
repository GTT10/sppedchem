# *********************************************************************
# **                                                                 **
# **   Subdirectory containing SpeedCHEM library                     **
# **                                                                 **
# **   Last modified by: Federico Perini                             **
# **   Last update     : 05/05/2012                                  **
# **                                                                 **
# *********************************************************************

include ../Makefile.inc

# Choose compilation flags
SCLIB     = ../libSpeedCHEM.a
OBJLIBS	  = $(SCLIB)

# NB order of the objects is important! do not change it
OBJS	= working_precision.o\
	      SCutilities.o\
	      SCsparse_definitions.o\
	      SCsparse.o\
	      sparse_MPI.o\
	      dvode_f90_m.o\
           SCmodule.o\
           chemkin_module.o\
           SCconV.o\
           SCsetup.o\
           SCallocate.o\
           SCcklink.o\
	       gam.o\
           gamsub.o\
           opkdmain.o\
           opkda1.o\
           opkda2.o\
	   ddaspk.o\
           rodas.o\
           vode.o\
           MEBDFSO.o\
           rowmap.o\
           radau5.o\
           radaua.o\
           chemistry_input.o\
	   SCbroadcast.o

# Decide on compilation flags	   
all    : $(OBJLIBS)
all    : COMPILESC = $(COMPILELIB) -c $(PRJCFLAGS) $(MPILIBS)

debug : $(OBJLIBS)
debug : COMPILESC = $(COMPILELIBDBG) -c $(PRJCFLAGS) $(MPILIBS)

COMPILESC = $(COMPILEFAST) $(MPILIBS) 
CFLAGS    = $(PRJCFLAGS) -I.. $(MPILIBS)

../libSpeedCHEM.a : $(OBJS) 
	$(ECHO) $(AR) $(ARFLAGS) rv $(SCLIB) $?
	$(AR) $(ARFLAGS) cr $(SCLIB) $?
	$(ECHO) $(RANLIB) $(SCLIB)
	$(RANLIB) $(SCLIB)

clean :
	$(ECHO) cleaning up in SpeedCHEM
	-$(RM) -f $(OBJS) *.mod


#force_look :
#	true
	
# COMPILE OBJECTS

%.o: %.f90 
	$(COMPILESC) $<

%.o: %.F90
	$(COMPILESC) $<

%.o: %.f
	$(COMPILESC) $<

%.o: %.f95
	$(COMPILESC) $<

%.o: %.f03
	$(COMPILESC) $<

%.o: %.f08
	$(COMPILESC) $<

%.o: %.for
	$(COMPILESC) $<

%.o: %.f77
	$(COMPILESC) $<	
