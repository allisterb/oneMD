bool debug_log = Simulator::debug_log;
int debug_log_level = Simulator::debug_log_level;
sycl::stream out(1024, 256, h);
#define __kernel_logger(K) \
    auto print = [out] (auto m) {out << "[kernel] [" << K << "] " << m << sycl::endl;}; \
    auto print2 = [out] (auto m, auto n) {out << "[kernel] [" << K << "] " << m << n << sycl::endl;}; \
    auto print3 = [out] (auto m, auto n, auto o) {out << "[kernel] [" << K << "] " << m << n << o << sycl::endl;}; \
    auto print4 = [out] (auto m, auto n, auto o, auto p) {out << "[kernel] [" << K << "] " << m << n << o << p << sycl::endl;}; \
    auto printd = [out, debug_log] (auto m) {if (debug_log) out << "[kernel] [" << K << "] " << m << sycl::endl;}; \
    auto printdd = [out, debug_log, debug_log_level] (auto m) {if (debug_log && debug_log_level >=2) out << "[kernel] [" << K << "] " << m << sycl::endl;}; \
    auto printi = [out] (sycl::id<1> i, auto m) {out << "[kernel] [" << K <<"] [index:" << i << "] " << m << sycl::endl;}; \
    auto printid = [out, debug_log] (sycl::id<1> i, auto m) {if (debug_log) out << "[kernel] [" << K <<"] [index:" << i << "] " << m << sycl::endl;}; \
    auto printidd = [out, debug_log, debug_log_level] (sycl::id<1> i, auto m) {if (debug_log && debug_log_level >= 2) out << "[kernel] [" << K <<"] [index:" << i << "] " << m << sycl::endl;}; \
    auto printij = [out] (sycl::id<2> i, auto m) {out << "[kernel] [" << K << "] [index:" << i << "] " << m << sycl::endl;}; \
    auto printijd = [out, debug_log] (sycl::id<2> i, auto m) {if (debug_log) out << "[kernel] [" << K << "] [index:" << i << "] " << m << sycl::endl;}; \
    auto printijdd = [out, debug_log, debug_log_level] (sycl::id<2> i, auto m) {if (debug_log && debug_log_level >= 2) out << "[kernel] [" << K << "] [index:" << i << "] " << m << sycl::endl;}; \
    auto printijdd2 = [out, debug_log, debug_log_level] (sycl::id<2> i, auto m, auto n) {if (debug_log && debug_log_level >= 2) out << "[kernel] [" << K << "] [index:" << i << "] " << m << n << sycl::endl;}; 