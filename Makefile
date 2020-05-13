# Makefile 
# Nom du compilateur
FC = mpiifort

# Options de compilation: optimisation, debug etc...
OPT =  -O3 -heap-arrays  -fopenmp -mcmodel=large 



# -ipo -xCORE-AVX2 -align array32byte -qopenmp -Ofast
# -mcmodel=large
# -heap-arrays 64



EXE = sol0

LINKOPT = 

DIR_OBJ = ./obj
DIR_BIN = ./bin

# Defining the objects (OBJS) variables
OBJS =  \
    variables_module.o \
    main.o \
    reading_data.o \
    gridder.o \
    initial_conditions.o \
    boundary_conditions.o \
    discretisation_QUICK_centre.o \
    calcul_new_velocity.o \
    check_steady.o \
    filer.o \
    virtualForceIntegrator.o \
    Mpi_division.o \
    GS.o \
    airfoil.o \
    smagorinsky.o \
    plasma.o

# Linking object files
exe :   $(OBJS)
	$(FC) $(OPT)   -o $(EXE) \
    $(OBJS) \
    $(OPT)



	@echo "   ***************** successful *****************   "                                                                                      
	#    |  Author:  Zi-Hsuan Wei                                                 
	#    |  Version: 1.8                                                         
	#    |  Web:     http://smetana.me.ntust.edu.tw/     
	#    |  Editor : Hsuan
    #    




%.o:%.f90
	$(FC) $(OPT) -c $<

main.o : variables_module.o

reading_data.o : variables_module.o

gridder.o : variables_module.o

initial_conditions.o : variables_module.o

boundary_conditions.o : variables_module.o

discretisation_upwind.o : variables_module.o

discretisation_QUICK_centre.o : variables_module.o

discretisation_QUICK.o : variables_module.o

calcul_new_velocity.o : variables_module.o

check_steady.o : variables_module.o

filer.o : variables_module.o

virtualForceIntegrator.o : variables_module.o

Mpi_division.o : variables_module.o

airfoil.o : variables_module.o

smagorinsky.o : variables_module.o


# Removing object files
clean :
	/bin/rm -f *.dat
	/bin/rm -f *.x
	/bin/rm -f *.q

cleanall : 
	/bin/rm -f $(OBJS) $(EXE)  *.mod
	/bin/rm -f *.dat
	/bin/rm -f *.x
	/bin/rm -f *.q
    
config :
	if [ ! -d obj ] ; then mkdir obj ; fi
	if [ ! -d run ] ; then mkdir bin ; fi


