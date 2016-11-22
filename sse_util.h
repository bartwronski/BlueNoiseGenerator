#include "debug_opt.h"
#if !defined(OLD_VER)


#ifndef SSE_UTIL_H
#define SSE_UTIL_H
#define SSE_UTIL_INCLUDE_START

#include "config.h"

#ifdef USE_SSE
	union SSERegister
	{
		__m128 v;
		float  s[4];
	};

	inline __m128 Negate(__m128 value);
	inline __m128 Mad(__m128 v0, __m128 v1, __m128 v2);
	inline __m128 FastPowSSEVector(__m128 value, float exponent);
	inline __m128 FastExpSSEVector(__m128 x);
	inline __m128 Abs(__m128 x);

#endif // USE_SSE

#include "sse_util.inl"
#undef SSE_UTIL_INCLUDE_START
#endif // SSE_UTIL_H

#endif
