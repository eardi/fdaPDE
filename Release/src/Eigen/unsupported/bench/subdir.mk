################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/Eigen/unsupported/bench/bench_svd.cpp 

OBJS += \
./src/Eigen/unsupported/bench/bench_svd.o 

CPP_DEPS += \
./src/Eigen/unsupported/bench/bench_svd.d 


# Each subdirectory must supply rules for building sources it contributes
src/Eigen/unsupported/bench/%.o: ../src/Eigen/unsupported/bench/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


