################################################
# Makefile for compiling csf2det               #
#                                              #
################################################

FC = gfortran 
FFLAGS = -fopenmp
#FFLAGS = -fcheck=all -fopenmp
#FC = ifort
#FFLAGS = -check bounds -openmp
SRCS = csf2det.F90
OBJS = csf2det.o
EXEC = csf2det.x

$(EXEC): $(OBJS)
	$(FC) $(FFLAGS) -o $(EXEC) $(OBJS)

%.o: %.F90
	$(FC) $(FFLAGS) -c $<

mcdensity.o: csf2det.F90
