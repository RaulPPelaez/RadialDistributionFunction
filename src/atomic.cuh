/*Raul P.Pelaez 2018, atomic for doubles
 */
#ifndef ATOMIC_CUH
#define ATOMIC_CUH

// namespace gdr{
//   namespace detail{
//     template<class T>
//     __device__ T atomicAdd(T* address, T val){
//       return ::atomicAdd(address, val);
//     }
//     template<>
#if !defined(__CUDA_ARCH__) || __CUDA_ARCH__ >= 600
#else
__device__ double atomicAdd(double* address, double val)
    {
      unsigned long long int* address_as_ull =
	(unsigned long long int*)address;
      unsigned long long int old = *address_as_ull, assumed;
      do {
	assumed = old;
	old = atomicCAS(address_as_ull, assumed,
			__double_as_longlong(val +
					     __longlong_as_double(assumed)));
      } while (assumed != old);
      return __longlong_as_double(old);
    }
#endif
//   }
// }
#endif