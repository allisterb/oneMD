sycl::stream out(1024, 256, h);
#define __kernel_logger(K) \
    auto printd = [out] (auto m) {out << "[kernel] [" << K << "] " << m << sycl::endl;}; \
    auto printd1 = [out] (auto m, sycl::id<1> i) {out << "[kernel] [" << K <<"] [index:" << i << "] " << m << sycl::endl;}; \
    auto printd2 = [out] (auto m, sycl::id<2> i) {out << "[kernel] [" << K << "] [index:" << i << "] " << m << sycl::endl;};