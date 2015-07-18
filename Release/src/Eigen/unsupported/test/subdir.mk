################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/Eigen/unsupported/test/BVH.cpp \
../src/Eigen/unsupported/test/FFT.cpp \
../src/Eigen/unsupported/test/FFTW.cpp \
../src/Eigen/unsupported/test/NonLinearOptimization.cpp \
../src/Eigen/unsupported/test/NumericalDiff.cpp \
../src/Eigen/unsupported/test/alignedvector3.cpp \
../src/Eigen/unsupported/test/autodiff.cpp \
../src/Eigen/unsupported/test/bdcsvd.cpp \
../src/Eigen/unsupported/test/dgmres.cpp \
../src/Eigen/unsupported/test/forward_adolc.cpp \
../src/Eigen/unsupported/test/gmres.cpp \
../src/Eigen/unsupported/test/jacobisvd.cpp \
../src/Eigen/unsupported/test/kronecker_product.cpp \
../src/Eigen/unsupported/test/levenberg_marquardt.cpp \
../src/Eigen/unsupported/test/matrix_exponential.cpp \
../src/Eigen/unsupported/test/matrix_function.cpp \
../src/Eigen/unsupported/test/matrix_power.cpp \
../src/Eigen/unsupported/test/matrix_square_root.cpp \
../src/Eigen/unsupported/test/minres.cpp \
../src/Eigen/unsupported/test/mpreal_support.cpp \
../src/Eigen/unsupported/test/openglsupport.cpp \
../src/Eigen/unsupported/test/polynomialsolver.cpp \
../src/Eigen/unsupported/test/polynomialutils.cpp \
../src/Eigen/unsupported/test/sparse_extra.cpp \
../src/Eigen/unsupported/test/splines.cpp 

OBJS += \
./src/Eigen/unsupported/test/BVH.o \
./src/Eigen/unsupported/test/FFT.o \
./src/Eigen/unsupported/test/FFTW.o \
./src/Eigen/unsupported/test/NonLinearOptimization.o \
./src/Eigen/unsupported/test/NumericalDiff.o \
./src/Eigen/unsupported/test/alignedvector3.o \
./src/Eigen/unsupported/test/autodiff.o \
./src/Eigen/unsupported/test/bdcsvd.o \
./src/Eigen/unsupported/test/dgmres.o \
./src/Eigen/unsupported/test/forward_adolc.o \
./src/Eigen/unsupported/test/gmres.o \
./src/Eigen/unsupported/test/jacobisvd.o \
./src/Eigen/unsupported/test/kronecker_product.o \
./src/Eigen/unsupported/test/levenberg_marquardt.o \
./src/Eigen/unsupported/test/matrix_exponential.o \
./src/Eigen/unsupported/test/matrix_function.o \
./src/Eigen/unsupported/test/matrix_power.o \
./src/Eigen/unsupported/test/matrix_square_root.o \
./src/Eigen/unsupported/test/minres.o \
./src/Eigen/unsupported/test/mpreal_support.o \
./src/Eigen/unsupported/test/openglsupport.o \
./src/Eigen/unsupported/test/polynomialsolver.o \
./src/Eigen/unsupported/test/polynomialutils.o \
./src/Eigen/unsupported/test/sparse_extra.o \
./src/Eigen/unsupported/test/splines.o 

CPP_DEPS += \
./src/Eigen/unsupported/test/BVH.d \
./src/Eigen/unsupported/test/FFT.d \
./src/Eigen/unsupported/test/FFTW.d \
./src/Eigen/unsupported/test/NonLinearOptimization.d \
./src/Eigen/unsupported/test/NumericalDiff.d \
./src/Eigen/unsupported/test/alignedvector3.d \
./src/Eigen/unsupported/test/autodiff.d \
./src/Eigen/unsupported/test/bdcsvd.d \
./src/Eigen/unsupported/test/dgmres.d \
./src/Eigen/unsupported/test/forward_adolc.d \
./src/Eigen/unsupported/test/gmres.d \
./src/Eigen/unsupported/test/jacobisvd.d \
./src/Eigen/unsupported/test/kronecker_product.d \
./src/Eigen/unsupported/test/levenberg_marquardt.d \
./src/Eigen/unsupported/test/matrix_exponential.d \
./src/Eigen/unsupported/test/matrix_function.d \
./src/Eigen/unsupported/test/matrix_power.d \
./src/Eigen/unsupported/test/matrix_square_root.d \
./src/Eigen/unsupported/test/minres.d \
./src/Eigen/unsupported/test/mpreal_support.d \
./src/Eigen/unsupported/test/openglsupport.d \
./src/Eigen/unsupported/test/polynomialsolver.d \
./src/Eigen/unsupported/test/polynomialutils.d \
./src/Eigen/unsupported/test/sparse_extra.d \
./src/Eigen/unsupported/test/splines.d 


# Each subdirectory must supply rules for building sources it contributes
src/Eigen/unsupported/test/%.o: ../src/Eigen/unsupported/test/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


