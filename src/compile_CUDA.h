#ifndef COLLIDOR_SRC_COMPILE_CUDA_H_
#define COLLIDOR_SRC_COMPILE_CUDA_H_

#include <cstdlib>
#include <time.h>
#include <sys/time.h>

#define CUDA_PREFIX  __device__
#define INLINE_PREFIX __inline__
#define HOST_PREFIX __host__


HOST_PREFIX double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

#endif
