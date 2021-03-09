#pragma once

#define CHUNKSIZE 15

#ifdef USE_ONEAPI

#include "dpc_common.hpp"
#include <CL/sycl.hpp>
#include <oneapi/mkl/vm.hpp>

#define SYCL_LINK SYCL_EXTERNAL
#define __sqrt sycl::sqrt
#define __nearbyint std::nearbyint

#else

#include <cmath>
#define SYCL_LINK
#define __sqrt sqrt
#define __nearbyint nearbyint

#endif