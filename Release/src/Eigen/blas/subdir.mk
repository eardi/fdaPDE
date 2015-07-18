################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/Eigen/blas/complex_double.cpp \
../src/Eigen/blas/complex_single.cpp \
../src/Eigen/blas/double.cpp \
../src/Eigen/blas/single.cpp \
../src/Eigen/blas/xerbla.cpp 

OBJS += \
./src/Eigen/blas/complex_double.o \
./src/Eigen/blas/complex_single.o \
./src/Eigen/blas/double.o \
./src/Eigen/blas/single.o \
./src/Eigen/blas/xerbla.o 

CPP_DEPS += \
./src/Eigen/blas/complex_double.d \
./src/Eigen/blas/complex_single.d \
./src/Eigen/blas/double.d \
./src/Eigen/blas/single.d \
./src/Eigen/blas/xerbla.d 


# Each subdirectory must supply rules for building sources it contributes
src/Eigen/blas/%.o: ../src/Eigen/blas/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


