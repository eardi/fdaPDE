################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/Eigen/failtest/block_nonconst_ctor_on_const_xpr_0.cpp \
../src/Eigen/failtest/block_nonconst_ctor_on_const_xpr_1.cpp \
../src/Eigen/failtest/block_nonconst_ctor_on_const_xpr_2.cpp \
../src/Eigen/failtest/block_on_const_type_actually_const_0.cpp \
../src/Eigen/failtest/block_on_const_type_actually_const_1.cpp \
../src/Eigen/failtest/const_qualified_block_method_retval_0.cpp \
../src/Eigen/failtest/const_qualified_block_method_retval_1.cpp \
../src/Eigen/failtest/const_qualified_diagonal_method_retval.cpp \
../src/Eigen/failtest/const_qualified_transpose_method_retval.cpp \
../src/Eigen/failtest/diagonal_nonconst_ctor_on_const_xpr.cpp \
../src/Eigen/failtest/diagonal_on_const_type_actually_const.cpp \
../src/Eigen/failtest/failtest_sanity_check.cpp \
../src/Eigen/failtest/map_nonconst_ctor_on_const_ptr_0.cpp \
../src/Eigen/failtest/map_nonconst_ctor_on_const_ptr_1.cpp \
../src/Eigen/failtest/map_nonconst_ctor_on_const_ptr_2.cpp \
../src/Eigen/failtest/map_nonconst_ctor_on_const_ptr_3.cpp \
../src/Eigen/failtest/map_nonconst_ctor_on_const_ptr_4.cpp \
../src/Eigen/failtest/map_on_const_type_actually_const_0.cpp \
../src/Eigen/failtest/map_on_const_type_actually_const_1.cpp \
../src/Eigen/failtest/ref_1.cpp \
../src/Eigen/failtest/ref_2.cpp \
../src/Eigen/failtest/ref_3.cpp \
../src/Eigen/failtest/ref_4.cpp \
../src/Eigen/failtest/ref_5.cpp \
../src/Eigen/failtest/transpose_nonconst_ctor_on_const_xpr.cpp \
../src/Eigen/failtest/transpose_on_const_type_actually_const.cpp 

OBJS += \
./src/Eigen/failtest/block_nonconst_ctor_on_const_xpr_0.o \
./src/Eigen/failtest/block_nonconst_ctor_on_const_xpr_1.o \
./src/Eigen/failtest/block_nonconst_ctor_on_const_xpr_2.o \
./src/Eigen/failtest/block_on_const_type_actually_const_0.o \
./src/Eigen/failtest/block_on_const_type_actually_const_1.o \
./src/Eigen/failtest/const_qualified_block_method_retval_0.o \
./src/Eigen/failtest/const_qualified_block_method_retval_1.o \
./src/Eigen/failtest/const_qualified_diagonal_method_retval.o \
./src/Eigen/failtest/const_qualified_transpose_method_retval.o \
./src/Eigen/failtest/diagonal_nonconst_ctor_on_const_xpr.o \
./src/Eigen/failtest/diagonal_on_const_type_actually_const.o \
./src/Eigen/failtest/failtest_sanity_check.o \
./src/Eigen/failtest/map_nonconst_ctor_on_const_ptr_0.o \
./src/Eigen/failtest/map_nonconst_ctor_on_const_ptr_1.o \
./src/Eigen/failtest/map_nonconst_ctor_on_const_ptr_2.o \
./src/Eigen/failtest/map_nonconst_ctor_on_const_ptr_3.o \
./src/Eigen/failtest/map_nonconst_ctor_on_const_ptr_4.o \
./src/Eigen/failtest/map_on_const_type_actually_const_0.o \
./src/Eigen/failtest/map_on_const_type_actually_const_1.o \
./src/Eigen/failtest/ref_1.o \
./src/Eigen/failtest/ref_2.o \
./src/Eigen/failtest/ref_3.o \
./src/Eigen/failtest/ref_4.o \
./src/Eigen/failtest/ref_5.o \
./src/Eigen/failtest/transpose_nonconst_ctor_on_const_xpr.o \
./src/Eigen/failtest/transpose_on_const_type_actually_const.o 

CPP_DEPS += \
./src/Eigen/failtest/block_nonconst_ctor_on_const_xpr_0.d \
./src/Eigen/failtest/block_nonconst_ctor_on_const_xpr_1.d \
./src/Eigen/failtest/block_nonconst_ctor_on_const_xpr_2.d \
./src/Eigen/failtest/block_on_const_type_actually_const_0.d \
./src/Eigen/failtest/block_on_const_type_actually_const_1.d \
./src/Eigen/failtest/const_qualified_block_method_retval_0.d \
./src/Eigen/failtest/const_qualified_block_method_retval_1.d \
./src/Eigen/failtest/const_qualified_diagonal_method_retval.d \
./src/Eigen/failtest/const_qualified_transpose_method_retval.d \
./src/Eigen/failtest/diagonal_nonconst_ctor_on_const_xpr.d \
./src/Eigen/failtest/diagonal_on_const_type_actually_const.d \
./src/Eigen/failtest/failtest_sanity_check.d \
./src/Eigen/failtest/map_nonconst_ctor_on_const_ptr_0.d \
./src/Eigen/failtest/map_nonconst_ctor_on_const_ptr_1.d \
./src/Eigen/failtest/map_nonconst_ctor_on_const_ptr_2.d \
./src/Eigen/failtest/map_nonconst_ctor_on_const_ptr_3.d \
./src/Eigen/failtest/map_nonconst_ctor_on_const_ptr_4.d \
./src/Eigen/failtest/map_on_const_type_actually_const_0.d \
./src/Eigen/failtest/map_on_const_type_actually_const_1.d \
./src/Eigen/failtest/ref_1.d \
./src/Eigen/failtest/ref_2.d \
./src/Eigen/failtest/ref_3.d \
./src/Eigen/failtest/ref_4.d \
./src/Eigen/failtest/ref_5.d \
./src/Eigen/failtest/transpose_nonconst_ctor_on_const_xpr.d \
./src/Eigen/failtest/transpose_on_const_type_actually_const.d 


# Each subdirectory must supply rules for building sources it contributes
src/Eigen/failtest/%.o: ../src/Eigen/failtest/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


