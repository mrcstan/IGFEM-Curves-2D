# The generated "Makefile" from "Makefile.cmake" is only usable after
# the Armadillo library has been configured and installed by CMake.

CXX=mex
#CXX=g++-4.2
## Under MacOS you may have an old compiler as default (e.g. GCC 4.0).
## However, GCC 4.2 or later is available and preferable due to better
## handling of template code.

#CXX=CC
## When using the Sun Studio compiler


# flags configured by CMake
ifeq (unix,macos)
  EXTRA_LIB_FLAGS = -framework Accelerate
endif

#EXTRA_LIB_FLAGS = -library=sunperf
## When using the Sun Studio compiler


LIB_FLAGS =  -larmadillo  -lgfortran $(EXTRA_LIB_FLAGS) 
## NOTE: on Ubuntu and Debian based systems you may need to add 
## -lgfortran to LIB_FLAGS


## As the Armadillo library uses recursive templates,
## compilation times depend on the level of optimisation:
##
## -O0: quick compilation, but the resulting program will be slow
## -O1: good trade-off between compilation time and execution speed
## -O2: produces programs which have almost all possible speedups,
##      but compilation takes longer


#OPT = -xO4 -xannotate=no
## When using the Sun Studio compiler


#EXTRA_OPT = -fwhole-program
## Uncomment the above line if you're compiling 
## all source files into one program in a single hit.


#DEBUG = -DARMA_EXTRA_DEBUG 
## Uncomment the above line to enable low-level
## debugging.  Lots of debugging information will
## be printed when a compiled program is run.
## Please enable this option when reporting bugs.


#FINAL = -DARMA_NO_DEBUG
## Uncomment the above line to disable Armadillo's checks.
## DANGEROUS!  Not recommended unless your code has been
## thoroughly tested.

CXXFLAGS  = -c -cxx -largeArrayDims -I/usr/local/include

EXTRA_CXXFLAGS = CXXFLAGS="$$CXXFLAGS -Wall -fopenmp"

#LDFLAGS = LDFLAGS="$$LDFLAGS -fopenmp" 
LDFLAGS = 

#MEXFUNC = mx_assemble_triplets
MEXFUNC = mx_assemble_sparse

SOURCES = $(MEXFUNC).cpp\
          assemble.cpp\
          compute_regular_element.cpp\
          shape_function.cpp\
          partition_stiffness_mat_n_load_vec.cpp\
          local_coord_2D.cpp\
          conductance.cpp\
          compute_IGFEM_element.cpp\
          child_element_NBJ.cpp\
          impose_Neumann_BC.cpp\
          impose_Dirichlet_BC.cpp\
          supg.cpp\
          body_source_functions.cpp
         # solve_mat_eqn.cpp

HEADERS = assemble.h


OBJECTS = $(SOURCES:.cpp=.o)

OBJDIR = 

EXES = $(MEXFUNC)

DIROBJECTS = $(addprefix $(OBJDIR), $(OBJECTS))

#$(EXES): $(DIROBJECTS) 
#	$(CXX) $(LDFLAGS)  $(DIROBJECTS)  -o $@ $(LIB_FLAGS) 
$(EXES): $(DIROBJECTS) 
	$(CXX) $(LDFLAGS)  $(DIROBJECTS)  $(LIB_FLAGS) 
	
$(OBJDIR)$(MEXFUNC).o: $(MEXFUNC).cpp $(HEADERS) 
	$(CXX) $(CXXFLAGS) $(EXTRA_CXXFLAGS) $<  

$(OBJDIR)assemble.o: assemble.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(EXTRA_CXXFLAGS) $< 

$(OBJDIR)compute_regular_element.o: compute_regular_element.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(EXTRA_CXXFLAGS) $< 

$(OBJDIR)shape_function.o: shape_function.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(EXTRA_CXXFLAGS) $<  

$(OBJDIR)partition_stiffness_mat_n_load_vec.o: partition_stiffness_mat_n_load_vec.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(EXTRA_CXXFLAGS) $< 

$(OBJDIR)local_coord_2D.o: local_coord_2D.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(EXTRA_CXXFLAGS) $< 

$(OBJDIR)conductance.o: conductance.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(EXTRA_CXXFLAGS) $<  

$(OBJDIR)compute_IGFEM_element.o: compute_IGFEM_element.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(EXTRA_CXXFLAGS) $< 
 
$(OBJDIR)child_element_NBJ.o: child_element_NBJ.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(EXTRA_CXXFLAGS) $<  

$(OBJDIR)impose_Neumann_BC.o: impose_Neumann_BC.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(EXTRA_CXXFLAGS) $<  

$(OBJDIR)impose_Dirichlet_BC.o: impose_Dirichlet_BC.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(EXTRA_CXXFLAGS) $<  

$(OBJDIR)supg.o: supg.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(EXTRA_CXXFLAGS) $<  

$(OBJDIR)body_source_functions.o: body_source_functions.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(EXTRA_CXXFLAGS) $<  



#$(OBJDIR)solve_mat_eqn.o: solve_mat_eqn.cpp $(HEADERS)
#	$(CXX) $(CXXFLAGS) $<  

.PHONY: clean

clean:
	rm -f $(EXES) $(OBJDIR)*.o

