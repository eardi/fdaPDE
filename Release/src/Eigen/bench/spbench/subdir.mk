################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/Eigen/bench/spbench/sp_solver.cpp \
../src/Eigen/bench/spbench/spbenchsolver.cpp \
../src/Eigen/bench/spbench/test_sparseLU.cpp 

OBJS += \
./src/Eigen/bench/spbench/sp_solver.o \
./src/Eigen/bench/spbench/spbenchsolver.o \
./src/Eigen/bench/spbench/test_sparseLU.o 

CPP_DEPS += \
./src/Eigen/bench/spbench/sp_solver.d \
./src/Eigen/bench/spbench/spbenchsolver.d \
./src/Eigen/bench/spbench/test_sparseLU.d 


# Each subdirectory must supply rules for building sources it contributes
src/Eigen/bench/spbench/%.o: ../src/Eigen/bench/spbench/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


