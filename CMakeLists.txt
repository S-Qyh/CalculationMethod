cmake_minimum_required(VERSION 3.21)
project(CalculationMethod)

set(CMAKE_CXX_STANDARD 14 Debug)

#set(CMAKE_BUILD_TYPE Debug)
#add_executable(CalculationMethod main_1.cpp Vector.h Matrix.h LinearEquationsSolve.cpp LinearEquationsSolve.h)
#add_executable(CalculationMethod EquationSolve.cpp EquationSolve.h Equation/main.cpp)
add_executable(CalculationMethod Interpolation.cpp Interpolation.h Matrix.h debug.cpp)
#add_executable(CalculationMethod NumericalIntegration.cpp NumericalIntegration.h debug_6.cpp)
