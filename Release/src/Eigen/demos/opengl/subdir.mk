################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/Eigen/demos/opengl/camera.cpp \
../src/Eigen/demos/opengl/gpuhelper.cpp \
../src/Eigen/demos/opengl/icosphere.cpp \
../src/Eigen/demos/opengl/quaternion_demo.cpp \
../src/Eigen/demos/opengl/trackball.cpp 

OBJS += \
./src/Eigen/demos/opengl/camera.o \
./src/Eigen/demos/opengl/gpuhelper.o \
./src/Eigen/demos/opengl/icosphere.o \
./src/Eigen/demos/opengl/quaternion_demo.o \
./src/Eigen/demos/opengl/trackball.o 

CPP_DEPS += \
./src/Eigen/demos/opengl/camera.d \
./src/Eigen/demos/opengl/gpuhelper.d \
./src/Eigen/demos/opengl/icosphere.d \
./src/Eigen/demos/opengl/quaternion_demo.d \
./src/Eigen/demos/opengl/trackball.d 


# Each subdirectory must supply rules for building sources it contributes
src/Eigen/demos/opengl/%.o: ../src/Eigen/demos/opengl/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


