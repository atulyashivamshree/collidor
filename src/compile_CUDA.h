/*
 *  compile_CUDA.h
    Description: main header file containing information regarding the CUDA
 library files

    @author Atulya Shivam Shree
    Created on: Dec 12, 2017
    Copyright (c) 2017 Atulya Shivam Shree
 */

#ifndef SRC_COMPILE_CUDA_H_
#define SRC_COMPILE_CUDA_H_

#include <sys/time.h>
#include <time.h>
#include <cstdlib>

#define CUDA_PREFIX __device__
#define INLINE_PREFIX __inline__
#define HOST_PREFIX __host__

HOST_PREFIX double get_wall_time() {
  struct timeval time;
  if (gettimeofday(&time, NULL)) {
    //  Handle error
    return 0;
  }
  return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

#endif  // SRC_COMPILE_CUDA_H_
