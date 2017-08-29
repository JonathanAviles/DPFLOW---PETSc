all: main
CFLAGS	         =
FFLAGS	         =
CPPFLAGS         =
FPPFLAGS         =

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

OBJECTS = main.o pf.o complexMatOperations.o

main: $(OBJECTS) chkopts
	-${CLINKER} -o dpflow $(OBJECTS) ${PETSC_SNES_LIB}
	${RM} $(OBJECTS)
