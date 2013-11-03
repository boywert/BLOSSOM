THIS = $(CURDIR)
SRC = $(THIS)/src
LIB = $(THIS)/lib
MICLIB = $(THIS)/miclib
BIN = $(THIS)/bin
CF = mpif90
MICCF = mpif90

CFFLAGS = -fpp -O2 -g -openmp 
MICCFFLAGS = -fpp -O3  -openmp  -mmic

CFFLAGS += -L/$(LIB) -I/$(LIB)
MICCFFLAGS += -L/$(MICLIB) -I/$(MICLIB)
OPT  = 
OPT += -DCLEARTEMP 
OPT += -DUSEMAXSOURCESIZE 
#OPT += -DINCLUDEPROTOGALACTIC 
#OPT += -DDEBUG 
OPT += -DUSERHO178
#OPT = -DCALCORR_RR -DCALCORR_RG -DGENLOS_RR -DGENLOS_RG -DAPSORBTION -DSUBCELL  -DTEST -DDEBUG 
#OPT += -DCALCORR_RR -DCALCORR_RG 
#OPT += -DSUBCELL
#OPT += -DGENLOS_RG #-DDEBUG
#OPT += -DGENLOS_RR #-DDEBUG
#OPT += -DABSORPTION_RG #-DDEBUG
#OPT += -DABSORPTION_RR
#OPT += -DGENRANDOM
CFFLAGS += $(OPT)
MICCFFLAGS += $(OPT)

all: $(BIN)/genrandom #$(BIN)/subcell $(BIN)/genlos_rg $(BIN)/genlos_rr $(BIN)/genrandom $(BIN)/absorption_rg  $(BIN)/absorption_rr 

$(BIN)/genrandom: $(LIB)/read_parameters.a $(LIB)/libcommon_vars.a $(LIB)/librunprocs.a
	mkdir -p $(BIN)
	$(CF) $(CFFLAGS) $(SRC)/main.f90 -o $(BIN)/genrandom -lread_parameters -lcommon_vars -lrunprocs -lmpitools -lio_tools -lvectortools -larraytools -ldatatools -labsorptiontools -lconversiontools -lutilities_serial

$(LIB)/read_parameters.a: $(LIB)/libcommon_vars.a
	$(CF) $(CFFLAGS) -c $(SRC)/read_parameters.f90 -o $(LIB)/libread_parameters.a  -lcommon_vars

$(LIB)/libcommon_vars.a: $(LIB)/libmpitools.a
	$(CF) $(CFFLAGS) -c $(SRC)/common_vars.f90 -o $(LIB)/libcommon_vars.a -lmpitools

$(LIB)/libutilities_serial.a: $(LIB)/libcommon_vars.a
	$(CF) $(CFFLAGS) -c $(SRC)/utilities_serial.f90 -o $(LIB)/libutilities_serial.a -lcommon_vars

$(LIB)/libmpitools.a: 
	$(CF) $(CFFLAGS) -c $(SRC)/mpitools.f90 -o $(LIB)/libmpitools.a 

$(LIB)/libio_tools.a: $(LIB)/libmpitools.a
	$(CF) $(CFFLAGS) -c $(SRC)/io_tools.f90 -o $(LIB)/libio_tools.a -lmpitools

$(LIB)/libvectortools.a:
	$(CF) $(CFFLAGS) -c $(SRC)/vectortools.f90 -o $(LIB)/libvectortools.a

$(LIB)/libarraytools.a:
	$(CF) $(CFFLAGS) -c $(SRC)/arraytools.f90 -o $(LIB)/libarraytools.a

$(LIB)/libdatatools.a:
	$(CF) $(CFFLAGS) -c $(SRC)/datatools.f90 -o $(LIB)/libdatatools.a

$(LIB)/libabsorptiontools.a: $(LIB)/libcommon_vars.a
	$(CF) $(CFFLAGS) -c $(SRC)/absorptiontools.f90 -o $(LIB)/libabsorptiontools.a -lcommon_vars

$(LIB)/libconversiontools.a: $(LIB)/libcommon_vars.a
	$(CF) $(CFFLAGS) -c $(SRC)/conversiontools.f90 -o $(LIB)/libconversiontools.a -lcommon_vars

$(LIB)/librunprocs.a: $(LIB)/libconversiontools.a $(LIB)/libabsorptiontools.a $(LIB)/libdatatools.a $(LIB)/libarraytools.a $(LIB)/libvectortools.a  $(LIB)/libio_tools.a $(LIB)/libmpitools.a $(LIB)/libcommon_vars.a $(LIB)/libutilities_serial.a
	$(CF) $(CFFLAGS) -c $(SRC)/runprocs.f90 -o $(LIB)/librunprocs.a -lconversiontools -ldatatools -larraytools -lvectortools -lio_tools -lmpitools -lcommon_vars -lutilities_serial 

clean:
	rm -f $(THIS)/*.mod
	rm -f $(SRC)/*.o
	rm -f $(LIB)/*
	rm -f $(MICLIB)/*
	rm -f *.o
	rm -f *.exe
git:
	git add .
	git commit -am "make `date +%F-%T`"
	git push 
