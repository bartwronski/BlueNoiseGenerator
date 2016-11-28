#include "debug_opt.h"
#if !defined(OLD_VER)


#ifndef MATH_UTIL_INCLUDE_START
	#error Please do not include this file directly, include "math_util.h" instead.
#endif

inline uint32_t NextPow2(uint32_t value)
{
	uint32_t power = 1;
	while (power < value) power *= 2u;
	return power;
}

size_t IntPow(size_t base, size_t exp)
{
	size_t ret = 1;
	while (exp--)
	{
		ret *= base;
	}
	return ret;
}

inline int FastModulo(int dividand, int divisor)
{
#ifdef USE_FAST_MODULO
	if ((divisor & (divisor - 1)) == 0)
	{
		return dividand & (divisor - 1);
	}
	else
	{
		return dividand % divisor;
	}
#else
	return dividand % divisor;
#endif
}

inline int FastDiv(int dividand, int divisor)
{
#ifdef USE_FAST_DIV
	if ((divisor & (divisor - 1)) == 0)
	{
		unsigned long fbs;
		_BitScanReverse(&fbs, divisor);
		return dividand >> fbs;
	}
	else
	{
		return dividand / divisor;
	}
#else
	return dividand / divisor;
#endif
}

inline float FastPowScalar(float value, float exponent)
{
#if defined(USE_SSE) && defined(USE_FAST_POW_SSE)
	if (exponent == 0.5f)
	{
		return _mm_cvtss_f32(_mm_rcp_ss(_mm_rsqrt_ss(_mm_set_ss(value))));
	}
	else if (exponent == 1.f)
	{
		return value;
	}
	else
	{
		return powf(value, exponent);
	}
#else
	return powf(value, exponent);
#endif
}

inline double	FastPowScalar(double value, double exponent)
{
	return pow(value, exponent);
}

// from https://codingforspeed.com/using-faster-exponential-approximation/
inline float FastExp(float x)
{
#ifdef USE_FAST_EXP
	x = x / 1024.f + 1.f;
	x *= x;
	x *= x;
	x *= x;
	x *= x;
	x *= x;
	x *= x;
	x *= x;
	x *= x;
	x *= x;
	x *= x;
	return float(x);
#else
	return float(expf(x));
#endif
}

inline double FastExp(double x)
{
	return exp(x);
}

#endif
