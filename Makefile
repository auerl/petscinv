ALL: petscinv

CFLAGS = -I${PETSC_DIR}/include
FFLAGS = -I${PETSC_DIR}/include/finclude
SOURCESC = 
SOURCESF = petscinv.F90
OBJ = ${SOURCESF:.F90=.o}
CLEANFILES = ${OBJ} test

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

petscinv: ${OBJ} chkopts
	-${FLINKER} -o petscinv ${OBJ} ${PETSC_SYS_LIB} ${PETSC_KSP_LIB}
	${RM} -f *.mod *.o

