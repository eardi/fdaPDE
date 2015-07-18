################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/Eigen/bench/btl/libs/gmm/main.cpp 

OBJS += \
./src/Eigen/bench/btl/libs/gmm/main.o 

CPP_DEPS += \
./src/Eigen/bench/btl/libs/gmm/main.d 


# Each subdirectory must supply rules for building sources it contributes
src/Eigen/bench/btl/libs/gmm/%.o: ../src/Eigen/bench/btl/libs/gmm/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


