#include "cuda.h"

#include <vector>
#include <iostream>

#include <cuda.h>

#define checkCudaErrors(val) check( (val), #val, __FILE__, __LINE__)

template<typename T>
void check(T err, const char* const func, const char* const file, const int line) {
  if (err != cudaSuccess) {
    std::cerr << "CUDA error at: " << file << ":" << line << std::endl;
    std::cerr << cudaGetErrorString(err) << " " << func << std::endl;
    exit(1);
  }
}

__global__ 
void cuda_test(float* left, float* right, float* result)
{
  int idx = threadIdx.x;
  
  result[idx] = left[idx] + right[idx];
}

void test() {
    int size = 128;
    std::vector<float> vec1 (size, 1.0);
    std::vector<float> vec2 (size, 3.0);
    std::vector<float> vec3 (size, 0.0);

    //Указатели на память видеокарте
    float* devVec1;
    float* devVec2;
    float* devVec3;

    checkCudaErrors(cudaMalloc((void**)&devVec1, sizeof(float) * size));
    checkCudaErrors(cudaMalloc((void**)&devVec2, sizeof(float) * size));
    checkCudaErrors(cudaMalloc((void**)&devVec3, sizeof(float) * size));

    cudaMemcpy(devVec1, vec1.data(), sizeof(float) * size, cudaMemcpyHostToDevice);
    cudaMemcpy(devVec2, vec2.data(), sizeof(float) * size, cudaMemcpyHostToDevice);

    dim3 gridSize = dim3(1, 1, 1);
    dim3 blockSize = dim3(size, 1, 1);

    //Выполняем вызов функции ядра
    cuda_test<<<gridSize, blockSize>>>(devVec1, devVec2, devVec3);

    cudaEvent_t syncEvent;

    cudaEventCreate(&syncEvent);    //Создаем event
    cudaEventRecord(syncEvent, 0);  //Записываем event
    cudaEventSynchronize(syncEvent);  //Синхронизируем event
  
    //Только теперь получаем результат расчета
    cudaMemcpy(vec3.data(), devVec3, sizeof(float) * size, cudaMemcpyDeviceToHost);

    std::cout << vec3[10] << std::endl;
}