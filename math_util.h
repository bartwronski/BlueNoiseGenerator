#include "debug_opt.h"
#if !defined(OLD_VER)


#ifndef MATH_UTIL_H
#define MATH_UTIL_H
#define MATH_UTIL_INCLUDE_START

#include "config.h"

inline uint32_t	NextPow2(uint32_t value);
inline size_t	IntPow(size_t base, size_t exp);
inline int		FastModulo(int dividand, int divisor);
inline int		FastDiv(int dividand, int divisor);
inline float	FastPowScalar(float value, float exponent);
inline double	FastPowScalar(double value, double exponent);
inline double	FastExp(double x);
inline float	Clamp(float value, float v0, float v1) { return std::max(v0, std::min(v1, value)); }
inline float    Saturate(float value) { return Clamp(value, 0.f, 1.f);  }

inline double Sqr(double val) { return val * val;  }

#include "math_util.inl"
#undef MATH_UTIL_INCLUDE_START
#endif // SSE_UTIL_H


#endif
