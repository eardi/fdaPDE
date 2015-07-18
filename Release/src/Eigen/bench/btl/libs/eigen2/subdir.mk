################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/Eigen/bench/btl/libs/eigen2/btl_tiny_eigen2.cpp \
../src/Eigen/bench/btl/libs/eigen2/main_adv.cpp \
../src/Eigen/bench/btl/libs/eigen2/main_linear.cpp \
../src/Eigen/bench/btl/libs/eigen2/main_matmat.cpp \
../src/Eigen/bench/btl/libs/eigen2/main_vecmat.cpp 

OBJS += \
./src/Eigen/bench/btl/libs/eigen2/btl_tiny_eigen2.o \
./src/Eigen/bench/btl/libs/eigen2/main_adv.o \
./src/Eigen/bench/btl/libs/eigen2/main_linear.o \
./src/Eigen/bench/btl/libs/eigen2/main_matmat.o \
./src/Eigen/bench/btl/libs/eigen2/main_vecmat.o 

CPP_DEPS += \
./src/Eigen/bench/btl/libs/eigen2/btl_tiny_eigen2.d \
./src/Eigen/bench/btl/libs/eigen2/main_adv.d \
./src/Eigen/bench/btl/libs/eigen2/main_linear.d \
./src/Eigen/bench/btl/libs/eigen2/main_matmat.d \
./src/Eigen/bench/btl/libs/eigen2/main_vecmat.d 


# Each subdirectory must supply rules for building sources it contributes
src/Eigen/bench/btl/libs/eigen2/%.o: ../src/Eigen/bench/btl/libs/eigen2/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


