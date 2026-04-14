// cuda_utils.cpp — Rcpp exports for CUDA device query functions.
//
// These wrappers are always compiled regardless of CUDA availability.
// When the package is built without CUDA (HAVE_CUDA not defined) they
// return 0 / empty list so that R-level code can check safely.

#ifdef HAVE_CUDA
#include <cuda_runtime.h>
#endif

#include <Rcpp.h>
#include "BioclimEngine.hpp"

//' Count available CUDA GPU devices
//'
//' Returns the number of CUDA-capable GPUs available on this machine.
//' Returns \code{0} when the package was built without CUDA support or when
//' no CUDA-capable device is found.
//'
//' @return Non-negative integer: number of CUDA devices detected.
//' @seealso \code{\link{has_cuda}}, \code{\link{cuda_device_info}}
//' @export
// [[Rcpp::export]]
int cuda_device_count() {
#ifdef HAVE_CUDA
    int count = 0;
    cudaError_t err = cudaGetDeviceCount(&count);
    if (err != cudaSuccess) return 0;
    return count;
#else
    return 0;
#endif
}

//' Query properties of the first CUDA GPU device
//'
//' Returns a named list with hardware information about the first
//' CUDA-capable GPU.  Returns an empty list when no CUDA device is
//' available or when the package was built without CUDA.
//'
//' @return Named list with:
//'   \describe{
//'     \item{name}{Character: GPU model name.}
//'     \item{memory_gb}{Numeric: total global memory in gigabytes.}
//'     \item{compute_capability}{Character: e.g. \code{"8.0"} for A100.}
//'   }
//'   An empty list when no CUDA device is detected.
//' @seealso \code{\link{has_cuda}}, \code{\link{cuda_device_count}}
//' @export
// [[Rcpp::export]]
Rcpp::List cuda_device_info() {
#ifdef HAVE_CUDA
    int count = 0;
    cudaError_t err = cudaGetDeviceCount(&count);
    if (err != cudaSuccess || count == 0) return Rcpp::List();
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    return Rcpp::List::create(
        Rcpp::Named("name")               = std::string(prop.name),
        Rcpp::Named("memory_gb")          =
            prop.totalGlobalMem / (1024.0 * 1024.0 * 1024.0),
        Rcpp::Named("compute_capability") =
            std::to_string(prop.major) + "." + std::to_string(prop.minor)
    );
#else
    return Rcpp::List();
#endif
}

//' Set the compute device for a BioclimEngine instance
//'
//' Controls whether the computation runs on a CUDA GPU or the CPU.
//' When \code{"auto"} is selected the engine uses the GPU if at least one
//' CUDA device is present, otherwise it falls back to the CPU.  GPU
//' requests on systems without a CUDA device silently fall back to the CPU.
//'
//' @param xptr External pointer returned by \code{\link{engine_create}}.
//' @param device Character scalar: one of \code{"auto"}, \code{"cpu"},
//'   or \code{"gpu"}.
//' @return \code{NULL} invisibly.
//' @seealso \code{\link{engine_create}}, \code{\link{has_cuda}},
//'   \code{\link{bioclim_engine}}
//' @export
// [[Rcpp::export]]
void engine_set_device(SEXP xptr, std::string device) {
    Rcpp::XPtr<rxbioclim::BioclimEngine>(xptr)->set_device(device);
}
