
MAKEFILE_OPT = $(PWD)/Makefile_opt.in
include $(MAKEFILE_OPT)

MAKEFILE_H5  = $(PWD)/Makefile_dir.in
include $(MAKEFILE_H5)

CC = mpicc
FLAGS = -O3 -Wall -g

INC = -I$(H55)/include
LIB = -L$(H55)/lib -lm -lhdf5

OBJ = main.o readpar.o timestep.o onestep.o riemann.o mpisetup.o gridsetup.o domain.o misc.o geometry.o faces_alt.o exchange.o plm.o report.o profiler.o planet.o omega.o analysis.o bfields.o hlld.o rotframe.o $(INITIAL).o $(OUTPUT).o $(HYDRO).o $(BOUNDARY).o $(RESTART).o $(PLANETS).o $(METRIC).o $(FRAME).o #snapshot.o

default: disco

%.o: %.c andrea.h
	$(CC) $(FLAGS) $(INC) -c $<

$(TIMESTEP).o: Timestep/$(TIMESTEP).c andrea.h
	$(CC) $(FLAGS) $(INC) -c Timestep/$(TIMESTEP).c

$(INITIAL).o : Initial/$(INITIAL).c andrea.h
	$(CC) $(FLAGS) $(INC) -c Initial/$(INITIAL).c

$(HYDRO).o : Hydro/$(HYDRO).c andrea.h
	$(CC) $(FLAGS) $(INC) -c Hydro/$(HYDRO).c

$(PLANETS).o : Planets/$(PLANETS).c andrea.h
	$(CC) $(FLAGS) $(INC) -c Planets/$(PLANETS).c

$(BOUNDARY).o : Boundary/$(BOUNDARY).c andrea.h
	$(CC) $(FLAGS) $(INC) -c Boundary/$(BOUNDARY).c

$(OUTPUT).o : Output/$(OUTPUT).c andrea.h
	$(CC) $(FLAGS) $(INC) -c Output/$(OUTPUT).c

$(RESTART).o : Restart/$(RESTART).c andrea.h
	$(CC) $(FLAGS) $(INC) -c Restart/$(RESTART).c

$(METRIC).o : Hydro/Metric/$(METRIC).c andrea.h Hydro/metric.h
	$(CC) $(FLAGS) $(INC) -c Hydro/Metric/$(METRIC).c

$(FRAME).o : Hydro/Frame/$(FRAME).c andrea.h Hydro/frame.h
	$(CC) $(FLAGS) $(INC) -c Hydro/Frame/$(FRAME).c

disco: $(OBJ) andrea.h
	$(CC) $(FLAGS) $(LIB) -o disco $(OBJ)

clean:
	rm -f *.o disco
