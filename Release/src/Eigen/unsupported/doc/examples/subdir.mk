################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/Eigen/unsupported/doc/examples/BVH_Example.cpp \
../src/Eigen/unsupported/doc/examples/FFT.cpp \
../src/Eigen/unsupported/doc/examples/MatrixExponential.cpp \
../src/Eigen/unsupported/doc/examples/MatrixFunction.cpp \
../src/Eigen/unsupported/doc/examples/MatrixLogarithm.cpp \
../src/Eigen/unsupported/doc/examples/MatrixPower.cpp \
../src/Eigen/unsupported/doc/examples/MatrixPower_optimal.cpp \
../src/Eigen/unsupported/doc/examples/MatrixSine.cpp \
../src/Eigen/unsupported/doc/examples/MatrixSinh.cpp \
../src/Eigen/unsupported/doc/examples/MatrixSquareRoot.cpp \
../src/Eigen/unsupported/doc/examples/PolynomialSolver1.cpp \
../src/Eigen/unsupported/doc/examples/PolynomialUtils1.cpp 

OBJS += \
./src/Eigen/unsupported/doc/examples/BVH_Example.o \
./src/Eigen/unsupported/doc/examples/FFT.o \
./src/Eigen/unsupported/doc/examples/MatrixExponential.o \
./src/Eigen/unsupported/doc/examples/MatrixFunction.o \
./src/Eigen/unsupported/doc/examples/MatrixLogarithm.o \
./src/Eigen/unsupported/doc/examples/MatrixPower.o \
./src/Eigen/unsupported/doc/examples/MatrixPower_optimal.o \
./src/Eigen/unsupported/doc/examples/MatrixSine.o \
./src/Eigen/unsupported/doc/examples/MatrixSinh.o \
./src/Eigen/unsupported/doc/examples/MatrixSquareRoot.o \
./src/Eigen/unsupported/doc/examples/PolynomialSolver1.o \
./src/Eigen/unsupported/doc/examples/PolynomialUtils1.o 

CPP_DEPS += \
./src/Eigen/unsupported/doc/examples/BVH_Example.d \
./src/Eigen/unsupported/doc/examples/FFT.d \
./src/Eigen/unsupported/doc/examples/MatrixExponential.d \
./src/Eigen/unsupported/doc/examples/MatrixFunction.d \
./src/Eigen/unsupported/doc/examples/MatrixLogarithm.d \
./src/Eigen/unsupported/doc/examples/MatrixPower.d \
./src/Eigen/unsupported/doc/examples/MatrixPower_optimal.d \
./src/Eigen/unsupported/doc/examples/MatrixSine.d \
./src/Eigen/unsupported/doc/examples/MatrixSinh.d \
./src/Eigen/unsupported/doc/examples/MatrixSquareRoot.d \
./src/Eigen/unsupported/doc/examples/PolynomialSolver1.d \
./src/Eigen/unsupported/doc/examples/PolynomialUtils1.d 


# Each subdirectory must supply rules for building sources it contributes
src/Eigen/unsupported/doc/examples/%.o: ../src/Eigen/unsupported/doc/examples/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


