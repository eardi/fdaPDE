################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CXX_SRCS += \
../src/Eigen/bench/btl/data/mean.cxx \
../src/Eigen/bench/btl/data/regularize.cxx \
../src/Eigen/bench/btl/data/smooth.cxx 

OBJS += \
./src/Eigen/bench/btl/data/mean.o \
./src/Eigen/bench/btl/data/regularize.o \
./src/Eigen/bench/btl/data/smooth.o 

CXX_DEPS += \
./src/Eigen/bench/btl/data/mean.d \
./src/Eigen/bench/btl/data/regularize.d \
./src/Eigen/bench/btl/data/smooth.d 


# Each subdirectory must supply rules for building sources it contributes
src/Eigen/bench/btl/data/%.o: ../src/Eigen/bench/btl/data/%.cxx
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


