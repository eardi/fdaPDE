################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/Eigen/demos/mandelbrot/mandelbrot.cpp 

OBJS += \
./src/Eigen/demos/mandelbrot/mandelbrot.o 

CPP_DEPS += \
./src/Eigen/demos/mandelbrot/mandelbrot.d 


# Each subdirectory must supply rules for building sources it contributes
src/Eigen/demos/mandelbrot/%.o: ../src/Eigen/demos/mandelbrot/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


