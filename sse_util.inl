#include "debug_opt.h"
#if !defined(OLD_VER)

#ifndef SSE_UTIL_INCLUDE_START
	#error Please do not include this file directly, include "sse_util.h" instead.
#endif

#ifdef USE_SSE

inline __m128 Negate(__m128 value)
{
	return _mm_xor_ps(value, _mm_set1_ps(-0.0));
}

inline __m128 Mad(__m128 v0, __m128 v1, __m128 v2)
{
	return _mm_add_ps(_mm_mul_ps(v0, v1), v2);
}

inline __m128 FastPowSSEVector(__m128 value, float exponent)
{
	if (exponent == 0.5f)
	{
		return _mm_rcp_ps(_mm_rsqrt_ps(value));
	}
	else if (exponent == 1.f)
	{
		return value;
	}
	const SSERegister &sseReg = (const SSERegister &)value;
	return _mm_set_ps(powf(sseReg.s[0], exponent), powf(sseReg.s[1], exponent), powf(sseReg.s[2], exponent), powf(sseReg.s[3], exponent));
}

inline __m128 FastExpSSEVector(__m128 x)
{
	#ifdef USE_FAST_EXP
		x = Mad(x, _mm_set_ps1(1.f / 1024.f), _mm_set_ps1(1.0));
		x = _mm_mul_ps(x, x);
		x = _mm_mul_ps(x, x);
		x = _mm_mul_ps(x, x);
		x = _mm_mul_ps(x, x);
		x = _mm_mul_ps(x, x);
		x = _mm_mul_ps(x, x);
		x = _mm_mul_ps(x, x);
		x = _mm_mul_ps(x, x);
		x = _mm_mul_ps(x, x);
		x = _mm_mul_ps(x, x);
		return x;
	#else
		const SSERegister &src = (const SSERegister &)x;
		SSERegister result;
		result.s[0] = expf(src.s[0]);
		result.s[1] = expf(src.s[1]);
		result.s[2] = expf(src.s[2]);
		result.s[3] = expf(src.s[3]);
		return result.v;
	#endif
}

#endif // USe_SSE

#endif
