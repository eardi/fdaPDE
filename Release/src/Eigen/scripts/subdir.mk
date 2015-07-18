################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/Eigen/scripts/eigen_gen_credits.cpp 

OBJS += \
./src/Eigen/scripts/eigen_gen_credits.o 

CPP_DEPS += \
./src/Eigen/scripts/eigen_gen_credits.d 


# Each subdirectory must supply rules for building sources it contributes
src/Eigen/scripts/%.o: ../src/Eigen/scripts/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


