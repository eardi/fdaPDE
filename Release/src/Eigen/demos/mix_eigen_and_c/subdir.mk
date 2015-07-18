################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../src/Eigen/demos/mix_eigen_and_c/example.c 

CPP_SRCS += \
../src/Eigen/demos/mix_eigen_and_c/binary_library.cpp 

OBJS += \
./src/Eigen/demos/mix_eigen_and_c/binary_library.o \
./src/Eigen/demos/mix_eigen_and_c/example.o 

C_DEPS += \
./src/Eigen/demos/mix_eigen_and_c/example.d 

CPP_DEPS += \
./src/Eigen/demos/mix_eigen_and_c/binary_library.d 


# Each subdirectory must supply rules for building sources it contributes
src/Eigen/demos/mix_eigen_and_c/%.o: ../src/Eigen/demos/mix_eigen_and_c/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

src/Eigen/demos/mix_eigen_and_c/%.o: ../src/Eigen/demos/mix_eigen_and_c/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: Cross GCC Compiler'
	gcc -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


