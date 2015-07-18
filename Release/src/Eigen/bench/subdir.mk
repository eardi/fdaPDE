################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/Eigen/bench/basicbenchmark.cpp \
../src/Eigen/bench/benchBlasGemm.cpp \
../src/Eigen/bench/benchCholesky.cpp \
../src/Eigen/bench/benchEigenSolver.cpp \
../src/Eigen/bench/benchFFT.cpp \
../src/Eigen/bench/benchGeometry.cpp \
../src/Eigen/bench/benchVecAdd.cpp \
../src/Eigen/bench/bench_gemm.cpp \
../src/Eigen/bench/bench_norm.cpp \
../src/Eigen/bench/bench_reverse.cpp \
../src/Eigen/bench/bench_sum.cpp \
../src/Eigen/bench/benchmark.cpp \
../src/Eigen/bench/benchmarkSlice.cpp \
../src/Eigen/bench/benchmarkX.cpp \
../src/Eigen/bench/benchmarkXcwise.cpp \
../src/Eigen/bench/check_cache_queries.cpp \
../src/Eigen/bench/eig33.cpp \
../src/Eigen/bench/geometry.cpp \
../src/Eigen/bench/product_threshold.cpp \
../src/Eigen/bench/quat_slerp.cpp \
../src/Eigen/bench/quatmul.cpp \
../src/Eigen/bench/sparse_cholesky.cpp \
../src/Eigen/bench/sparse_dense_product.cpp \
../src/Eigen/bench/sparse_lu.cpp \
../src/Eigen/bench/sparse_product.cpp \
../src/Eigen/bench/sparse_randomsetter.cpp \
../src/Eigen/bench/sparse_setter.cpp \
../src/Eigen/bench/sparse_transpose.cpp \
../src/Eigen/bench/sparse_trisolver.cpp \
../src/Eigen/bench/spmv.cpp \
../src/Eigen/bench/vdw_new.cpp 

OBJS += \
./src/Eigen/bench/basicbenchmark.o \
./src/Eigen/bench/benchBlasGemm.o \
./src/Eigen/bench/benchCholesky.o \
./src/Eigen/bench/benchEigenSolver.o \
./src/Eigen/bench/benchFFT.o \
./src/Eigen/bench/benchGeometry.o \
./src/Eigen/bench/benchVecAdd.o \
./src/Eigen/bench/bench_gemm.o \
./src/Eigen/bench/bench_norm.o \
./src/Eigen/bench/bench_reverse.o \
./src/Eigen/bench/bench_sum.o \
./src/Eigen/bench/benchmark.o \
./src/Eigen/bench/benchmarkSlice.o \
./src/Eigen/bench/benchmarkX.o \
./src/Eigen/bench/benchmarkXcwise.o \
./src/Eigen/bench/check_cache_queries.o \
./src/Eigen/bench/eig33.o \
./src/Eigen/bench/geometry.o \
./src/Eigen/bench/product_threshold.o \
./src/Eigen/bench/quat_slerp.o \
./src/Eigen/bench/quatmul.o \
./src/Eigen/bench/sparse_cholesky.o \
./src/Eigen/bench/sparse_dense_product.o \
./src/Eigen/bench/sparse_lu.o \
./src/Eigen/bench/sparse_product.o \
./src/Eigen/bench/sparse_randomsetter.o \
./src/Eigen/bench/sparse_setter.o \
./src/Eigen/bench/sparse_transpose.o \
./src/Eigen/bench/sparse_trisolver.o \
./src/Eigen/bench/spmv.o \
./src/Eigen/bench/vdw_new.o 

CPP_DEPS += \
./src/Eigen/bench/basicbenchmark.d \
./src/Eigen/bench/benchBlasGemm.d \
./src/Eigen/bench/benchCholesky.d \
./src/Eigen/bench/benchEigenSolver.d \
./src/Eigen/bench/benchFFT.d \
./src/Eigen/bench/benchGeometry.d \
./src/Eigen/bench/benchVecAdd.d \
./src/Eigen/bench/bench_gemm.d \
./src/Eigen/bench/bench_norm.d \
./src/Eigen/bench/bench_reverse.d \
./src/Eigen/bench/bench_sum.d \
./src/Eigen/bench/benchmark.d \
./src/Eigen/bench/benchmarkSlice.d \
./src/Eigen/bench/benchmarkX.d \
./src/Eigen/bench/benchmarkXcwise.d \
./src/Eigen/bench/check_cache_queries.d \
./src/Eigen/bench/eig33.d \
./src/Eigen/bench/geometry.d \
./src/Eigen/bench/product_threshold.d \
./src/Eigen/bench/quat_slerp.d \
./src/Eigen/bench/quatmul.d \
./src/Eigen/bench/sparse_cholesky.d \
./src/Eigen/bench/sparse_dense_product.d \
./src/Eigen/bench/sparse_lu.d \
./src/Eigen/bench/sparse_product.d \
./src/Eigen/bench/sparse_randomsetter.d \
./src/Eigen/bench/sparse_setter.d \
./src/Eigen/bench/sparse_transpose.d \
./src/Eigen/bench/sparse_trisolver.d \
./src/Eigen/bench/spmv.d \
./src/Eigen/bench/vdw_new.d 


# Each subdirectory must supply rules for building sources it contributes
src/Eigen/bench/%.o: ../src/Eigen/bench/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


