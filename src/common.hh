#pragma once

#define CHUNKSIZE 15

#ifdef USE_ONEAPI

#include "dpc_common.hpp"
#include <CL/sycl.hpp>

#include "Exceptions.hh"

#define SYCL_LINK SYCL_EXTERNAL
#define kernel_printf sycl::ONEAPI::experimental::printf
#else

#define SYCL_LINK
#define kernel_printf std::printf
#endif