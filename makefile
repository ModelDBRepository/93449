CC = mpiCC
LD = mpiCC
CFLAGS = -Wall -O3 -ansi 
C_OBJ_BUILD = $(CC) $(CFLAGS) -c
PROGRAM = test.exe
OBJS = control.o integrator.o convlv.o realft.o twofft.o four1.o nrutil.o

$(PROGRAM) : $(OBJS)
	$(LD) $(OBJS) /home/xfeng/lib/odepack/liblsoda.a /usr/lib/gcc-lib/i386-redhat-linux/3.3.3/libg2c.a -o $@

control.o : control.C control.h integrator.h
	$(C_OBJ_BUILD) $<
integrator.o : integrator.C control.h integrator.h 
	$(C_OBJ_BUILD) $< 
convlv.o : convlv.c nrutil.h 
	$(C_OBJ_BUILD) $<
realft.o : realft.c nrutil.h 
	$(C_OBJ_BUILD) $<
twofft.o : twofft.c nrutil.h 
	$(C_OBJ_BUILD) $<
four1.o : four1.c 
	$(C_OBJ_BUILD) $<
nrutil.o : nrutil.c nrutil.h
	$(C_OBJ_BUILD) $<
clean :
	rm *.o
