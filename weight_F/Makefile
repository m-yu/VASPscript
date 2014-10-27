.SUFFIXES: .f90

FC = ifort

FFLAGS = -r8 -O

LINK = -static

OBJS = prec.o       \
       options.o    \
       matrix.o     \
       chgcar.o     \
       qsort.o      \
       wsflux.o     \
       weight.o     \
       main.o

TARGET = bader	 

%.o %.mod : %.f90
	$(FC) $(FFLAGS) -c $*.f90

$(TARGET): $(OBJS) 
	rm -f $(TARGET)
	$(FC) -o $(TARGET) $(LINK) $(OBJS) 

clean:
	rm -f *.o *.mod $(TARGET)
