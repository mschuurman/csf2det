################################################
# Makefile for compiling NADVIBS               #
#                                              #
################################################

FC = gfortran 
FFLAGS = -fbounds-check 
SRCS = csf2det.F90
OBJS = csf2det.o
EXEC = csf2det.x

$(EXEC): $(OBJS)
	$(FC) $(FFLAGS) -o $(EXEC) $(OBJS)

%.o: %.F90
	$(FC) $(FFLAGS) -c $<

mcdensity.o: csf2det.F90
