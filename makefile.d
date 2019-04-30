obj/Iteration.o : src/Iteration.F90 obj/LinAlgLib.o 
obj/LinAlgLib.o : submodule/LinAlgf90/src/LinAlgLib.f90 obj/MatVecComplex.o obj/MatVecDouble.o obj/MatVecSingle.o obj/MatrixComplex.o obj/MatrixDouble.o obj/MatrixSingle.o obj/VectorComplex.o obj/VectorDouble.o obj/VectorSingle.o obj/SingleDoubleComplex.o obj/LinAlgParameters.o 
obj/MatVecDouble.o : submodule/LinAlgf90/src/MatVecDouble.f90 obj/MatrixDouble.o obj/VectorDouble.o obj/LinAlgParameters.o 
obj/VectorSingle.o : submodule/LinAlgf90/src/VectorSingle.f90 obj/LinAlgParameters.o 
obj/LinAlgParameters.o : submodule/LinAlgf90/src/LinAlgParameters.f90 
obj/VectorComplex.o : submodule/LinAlgf90/src/VectorComplex.f90 obj/LinAlgParameters.o 
obj/MatrixSingle.o : submodule/LinAlgf90/src/MatrixSingle.f90 obj/VectorSingle.o obj/LinAlgParameters.o 
obj/MatrixDouble.o : submodule/LinAlgf90/src/MatrixDouble.f90 obj/VectorDouble.o obj/LinAlgParameters.o 
obj/MatVecSingle.o : submodule/LinAlgf90/src/MatVecSingle.f90 obj/MatrixSingle.o obj/VectorSingle.o obj/LinAlgParameters.o 
obj/VectorDouble.o : submodule/LinAlgf90/src/VectorDouble.f90 obj/LinAlgParameters.o 
obj/MatrixComplex.o : submodule/LinAlgf90/src/MatrixComplex.f90 obj/VectorComplex.o obj/LinAlgParameters.o 
obj/SingleDoubleComplex.o : submodule/LinAlgf90/src/SingleDoubleComplex.f90 obj/MatrixComplex.o obj/VectorComplex.o obj/MatrixDouble.o obj/VectorDouble.o obj/MatrixSingle.o obj/VectorSingle.o 
obj/MatVecComplex.o : submodule/LinAlgf90/src/MatVecComplex.f90 obj/MatrixComplex.o obj/VectorComplex.o obj/LinAlgParameters.o 
obj/test.o : main/test.f90 
