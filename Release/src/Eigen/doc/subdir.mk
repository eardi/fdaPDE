################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/Eigen/doc/tutorial.cpp 

OBJS += \
./src/Eigen/doc/tutorial.o 

CPP_DEPS += \
./src/Eigen/doc/tutorial.d 


# Each subdirectory must supply rules for building sources it contributes
src/Eigen/doc/%.o: ../src/Eigen/doc/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


