#include <random>
#include <string>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <chrono>
#include <assert.h>
#include <intrin.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

union SSERegister
{
	__m128 v;
	float  s[4];
};

inline __m128 Negate(__m128 value)
{
	return _mm_xor_ps(value, _mm_set1_ps(-0.0));
}

inline __m128 Mad(__m128 v0, __m128 v1, __m128 v2)
{
	return _mm_add_ps(_mm_mul_ps(v0, v1), v2);
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
	if ((divisor & (divisor - 1)) == 0)
	{
		return dividand & (divisor - 1);
	}
	else
	{
		return dividand % divisor;
	}
}

inline int FastDiv(int dividand, int divisor)
{
	if ((divisor & (divisor - 1)) == 0)
	{
		unsigned long fbs;
		_BitScanReverse(&fbs, divisor);
		return dividand >> fbs;
	}
	else
	{
		return dividand % divisor;
	}
}

inline float FastPowSSEScalar(float value, float exponent)
{
	if (exponent == 0.5f)
	{
		return _mm_cvtss_f32(_mm_rcp_ss(_mm_rsqrt_ss(_mm_set_ss(value))));
	}
	else if (exponent == 1.f)
	{
		return value;
	}
	return powf(value, exponent);
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
	const SSERegister &sseReg = (const SSERegister &) value;
	return _mm_set_ps(powf(sseReg.s[0], exponent), powf(sseReg.s[1], exponent), powf(sseReg.s[2], exponent), powf(sseReg.s[3], exponent));
}

// from https://codingforspeed.com/using-faster-exponential-approximation/
inline float FastExp(double x) 
{
	x = x / 1024 + 1.0;
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
}

inline __m128 FastExpSSEVector(__m128 x)
{
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
}

size_t ComputeElementCount(size_t dimCount, const size_t sizePerDim[])
{
	size_t elemCount = 1;
	for (size_t currDim = 0; currDim < dimCount; ++currDim)
	{
		elemCount *= sizePerDim[currDim];
	}
	return elemCount;
}

enum EMethod
{
	Method_SolidAngle,
	Method_HighPass,
};

const bool   useIncrementalUpdate = true;
const size_t N_dimensions = 2;
const size_t dimensionSize[N_dimensions] = { 32, 32 };
const size_t N_valuesPerItem = 1;
const size_t totalElements = ComputeElementCount(N_dimensions, dimensionSize);
size_t dimensionElementCount[N_dimensions + 1] = { 0 };

const EMethod chosenMethod = Method_SolidAngle;

// Solid Angle method parameters
const size_t distanceToCheck = 3; // original paper mentioned looking at whole array; however this is N^2 and super expensive, while exp(-4*4/2.2) ~= 0.000694216
const size_t distanceToCheckBoth = distanceToCheck * 2 + 1; // in both directions
const size_t elementsToCheck = IntPow(distanceToCheckBoth, N_dimensions);
const size_t numIterationsToFindDistribution = 256 * 1024;

// Highpass filter parameters
const size_t convSize				 = 3;
const size_t convSizeTotal			 = IntPow(convSize, N_dimensions);

const float convWeights1D[3]		 = { -1, 2, -1 };

const float convWeights2D[3*3]		 = { -1, -2, -1, 
									 	 -2, 12, -2, 
									 	 -1, -2, -1 };

const float convWeights3D[3 * 3 * 3] = { -1, -2, -1, 
										 -2, -4, -2, 
									     -1, -2, -1, 

										 -2, -4, -2,
										 -4, 56, -4,
										 -2, -4, -2,

										 -1, -2, -1,
										 -2, -4, -2,
										 -1, -2, -1,
};

void PrintCodeOutput(const std::string& fileName, const std::vector<float>& arr, const std::string& arrName, bool mathematica, const size_t dimensionSize[N_dimensions], size_t N_dimensions, size_t N_valuesPerItem)
{
	std::ofstream outFile;
	outFile.open(fileName, std::ios::out | std::ios::trunc);

	const size_t arrSize = arr.size();

	if (!mathematica)
	{
		outFile << "static const float " << arrName;
		for (size_t d = 0; d < N_dimensions; ++d)
		{
			outFile << "[" << dimensionSize[d] << "]";
		}
		if (N_valuesPerItem > 1)
		{
			outFile << "[" << N_valuesPerItem << "]";
		}
		outFile << " = " << std::endl;
	}

	for (size_t i = 0; i < arrSize / N_valuesPerItem; ++i)
	{
		for (size_t d = 0; d < N_dimensions; ++d)
		{
			size_t dim = (i / dimensionElementCount[d]) % dimensionSize[d];

			if (dim == 0)
			{
				outFile << "{";
			}
			else
			{
				break;
			}
		}

		if (N_valuesPerItem == 1)
		{
			outFile << std::setprecision(8) << std::fixed << arr[i];
		}
		else
		{
			outFile << "{";
			for (size_t v = 0; v < N_valuesPerItem; ++v)
			{
				outFile << std::setprecision(8) << std::fixed << arr[i * N_valuesPerItem + v];
				if (v < N_valuesPerItem - 1)
				{
					outFile << ", ";
				}
			}
			outFile << "}";
		}

		for (size_t d = 0; d < N_dimensions; ++d)
		{
			size_t dim = (i / dimensionElementCount[d]) % dimensionSize[d];

			if (dim == dimensionSize[d] - 1)
			{
				outFile << "}";
			}
			else
			{
				outFile << ",";
				if (d > 0)
				{
					outFile << std::endl;
				}
				break;
			}
		}
	}
	if (!mathematica)
	{
		outFile << ";";
	}
	outFile << std::endl << std::endl;
}

std::string dimNames[4] = { "x", "y", "z", "w" };

void PrintWebGLOutputRecursive(std::ofstream& outFile, const std::vector<float>& arr, const std::string& arrName, size_t N_dimensions, size_t N_valuesPerItem, size_t lo, size_t high)
{
	if (high - lo == 1)
	{
		outFile << "return ";

		if (N_valuesPerItem > 1)
		{
			outFile << "vec" << N_valuesPerItem << "(";
			for (size_t v = 0; v < N_valuesPerItem; ++v)
			{
				outFile << arr[lo * N_valuesPerItem + v];
				if (v < N_valuesPerItem - 1)
				{
					outFile << ", ";
				}
			}

			outFile << ");";
		}
		else
		{
			outFile << arr[lo] << ";";
		}
	}
	else
	{
		size_t mid = (lo + high) / 2;

		outFile << "if(" << arrName << " < " << mid << ") " << std::endl << "{" << std::endl;
		PrintWebGLOutputRecursive(outFile, arr, arrName, N_dimensions, N_valuesPerItem, lo, mid);
		outFile << "} else {" << std::endl;
		PrintWebGLOutputRecursive(outFile, arr, arrName, N_dimensions, N_valuesPerItem, mid, high);
		outFile << std::endl << "}";
	}
};

void PrintWebGLOutput(const std::string& fileName, const std::vector<float>& arr, const std::string& arrName, size_t N_dimensions, size_t N_valuesPerItem, size_t lo, size_t high)
{
	std::ofstream outFile;
	outFile.open(fileName, std::ios::out | std::ios::trunc);

	PrintWebGLOutputRecursive(outFile, arr, arrName, N_dimensions, N_valuesPerItem, lo, high);
}

inline uint32_t FloatAsByteUNorm(float value)
{
	return uint32_t(255.f * value);
}

void SaveAsPPM(const std::vector<float>& arr, const std::string& fileName, const size_t dimensionSize[N_dimensions], size_t N_dimensions, size_t N_valuesPerItem)
{
	std::ofstream outfile(fileName);
	assert(N_dimensions >= 2);
	N_dimensions = 2; // ignore depth
	outfile << "P3" << std::endl << dimensionSize[0] << " " << dimensionSize[1] << std::endl << 255 << std::endl;
	const uint32_t pixCount = arr.size() / N_valuesPerItem;
	for (size_t i = 0; i < pixCount; ++i)
	{
		switch (N_valuesPerItem)
		{
		case 1: // monichromatic values
			outfile << FloatAsByteUNorm(arr[i]) << " "
				<< FloatAsByteUNorm(arr[i]) << " "
				<< FloatAsByteUNorm(arr[i]) << " ";
			break;
		case 2:
			outfile << FloatAsByteUNorm(arr[i * 2]) << " "
				<< FloatAsByteUNorm(arr[i * 2 + 1]) << " "
				<< 0 << " ";
			break;
		case 3:
			outfile << FloatAsByteUNorm(arr[i * 3]) << " "
				<< FloatAsByteUNorm(arr[i * 3 + 1]) << " "
				<< FloatAsByteUNorm(arr[i * 3 + 2]) << " ";
			break;
		default:
			assert(0);
		}
	}
	outfile.close();
}

inline float ComputeFinalScore(const std::vector<float>& arr, float distanceScore, size_t N_valuesPerItem, size_t ind1, size_t ind2)
{
	float valueSpaceScore = 0;
	for (size_t i = 0; i < N_valuesPerItem; ++i)
	{
		float val = (arr[ind1 * N_valuesPerItem + i] - arr[ind2 * N_valuesPerItem + i]);
		valueSpaceScore += val * val;
	}
	valueSpaceScore = FastPowSSEScalar(valueSpaceScore, (float)N_valuesPerItem / 2.0f);
	const float oneOverDistanceVarianceSq = 1.0f / (2.1f * 2.1f);
	return FastExp(-valueSpaceScore - distanceScore * oneOverDistanceVarianceSq);
}

inline __m128 ComputeFinalScoreSSE(const std::vector<float>& arr, float distanceScore[4], size_t N_valuesPerItem, size_t ind1, size_t ind2[4])
{
	SSERegister valueSpaceScore;
	valueSpaceScore.v = _mm_set_ps1(0.f);
	for (size_t i = 0; i < N_valuesPerItem; ++i)
	{
		float srcValue = arr[ind1 * N_valuesPerItem + i];
		float val0 = srcValue - arr[ind2[0] * N_valuesPerItem + i];
		float val1 = srcValue - arr[ind2[1] * N_valuesPerItem + i];
		float val2 = srcValue - arr[ind2[2] * N_valuesPerItem + i];
		float val3 = srcValue - arr[ind2[3] * N_valuesPerItem + i];
		__m128 val = _mm_set_ps(val0, val1, val2, val3);
		valueSpaceScore.v = Mad(val, val, valueSpaceScore.v);
	}
	valueSpaceScore.v = FastPowSSEVector(valueSpaceScore.v, (float)N_valuesPerItem / 2.0f);
	const float oneOverDistanceVarianceSq = 1.0f / (2.1f * 2.1f);
	return FastExpSSEVector(Negate(Mad(_mm_set_ps(distanceScore[0], distanceScore[1], distanceScore[2], distanceScore[3]), _mm_set_ps1(oneOverDistanceVarianceSq), valueSpaceScore.v)));
}

inline float ComputeDistanceScore(const int arr[], size_t Ndimensions)
{
	float distanceSq = 0;
	for (size_t i = 0; i < Ndimensions; ++i)
	{
		distanceSq += arr[i] * arr[i];
	}
	return distanceSq;
}

double sign(double v) { return (v >= 0.0) ? 1.0 : -1.0; }
double remap_tri( const double v )
{
	double r2 = 0.5 * v;
	double f1 = sqrt( r2 );
	double f2 = 1.0 - sqrt( r2 - 0.25 );
	return (v < 0.5) ? f1 : f2;
}

//note: splats to 24b
uint8_t* FloatDataToBytes(const std::vector<float>& arr, size_t dimensionSize, size_t N_valuesPerItem, bool do_remap_tri)
{
	uint8_t* bytes = new uint8_t[3 * arr.size() / N_valuesPerItem];

	for (size_t i = 0, n = arr.size() / N_valuesPerItem; i<n; ++i)
	{
		uint8_t lastChar = 0;
		for (size_t j = 0; j < N_valuesPerItem; j++)
		{
			double v_f = arr[i * N_valuesPerItem + j];

			if (do_remap_tri)
				v_f = remap_tri(v_f);

			int v_i = static_cast<int>(v_f * 256.0);
			v_i = (v_i < 0) ? 0 : v_i;
			v_i = (v_i > 255) ? 255 : v_i;

			lastChar = static_cast<char>(v_i);

			bytes[3 * i + j] = lastChar;
		}

		for (size_t j = N_valuesPerItem; j < 3; j++)
		{
			// splat
			if (N_valuesPerItem == 1)
			{
				bytes[3 * i + j] = lastChar;
			}
			else
			{
				bytes[3 * i + j] = 0;
			}
		}
	}

	return bytes;
}

//note: see http://gpuopen.com/vdr-follow-up-fine-art-of-film-grain/
void UnifyHistogram(std::vector<float>& arr, size_t N_valuesPerItem)
{
	for (size_t dim = 0; dim < N_valuesPerItem; dim++)
	{
		std::vector<std::pair<float, size_t> > entries(arr.size() / N_valuesPerItem);

		for (size_t i = 0, n = arr.size() / N_valuesPerItem; i < n; ++i)
		{
			entries[i] = std::make_pair(arr[i * N_valuesPerItem + dim], i);
		}

		std::sort(entries.begin(), entries.end(), [](const std::pair<float, size_t>& a, const std::pair<float, size_t>& b) -> bool
		{
			return a.first < b.first;
		});

		for (size_t i = 0, n = entries.size(); i<n; ++i)
		{
			float t = static_cast<float>(i) / static_cast<float>(n-1);
			size_t idx = entries[i].second;
			arr[idx * N_valuesPerItem + dim] = t;
		}
	}
}

inline size_t WrapDimension(size_t baseIndex, int offset, size_t dimSize)
{
	int posWrapped = (int)baseIndex + offset;
	if (posWrapped < 0)
	{
		posWrapped += dimSize;
	}
	if (posWrapped > (int)(dimSize - 1))
	{
		posWrapped -= dimSize;
	}

	return posWrapped;
}

struct KernelSample
{
	float	Weight;
	int32_t Distances[N_dimensions];
};

int main(int argc, char** argv)
{
	for (uint32_t dim = 0; dim <= N_dimensions; ++dim)
	{
		dimensionElementCount[dim] = ComputeElementCount(dim, dimensionSize);
	}
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dist(0, 1);
	// Note: we try to swap between 1 and 3 elements to try to jump over local minima
	const int maxSwapedElemCount = 3;
	std::uniform_int_distribution<> distInt(1, maxSwapedElemCount);
	std::uniform_int_distribution<> distSwap(0, (int)(totalElements - 1));

	std::vector<float> pattern[2] = { std::vector<float>(totalElements * N_valuesPerItem), std::vector<float>(totalElements * N_valuesPerItem) };
	uint32_t currentArray = 0;

	for (size_t i = 0; i < totalElements * N_valuesPerItem; ++i)
	{
		pattern[currentArray][i] = static_cast<float>(dist(gen));
	}
	UnifyHistogram(pattern[currentArray], N_valuesPerItem);
	//PrintCodeOutput("initialDist.txt", pattern[currentArray], "initialDist", true, dimensionSize, N_dimensions, N_valuesPerItem);

	if (chosenMethod == Method_SolidAngle)
	{
		std::vector<KernelSample> kernel;
		for (size_t i = 0; i < elementsToCheck; ++i)
		{
			KernelSample ks;
			float dist = 0.f;
			for (size_t d = 0; d < N_dimensions; ++d)
			{
				size_t dim = (i / IntPow(distanceToCheckBoth, d)) % distanceToCheckBoth;
				ks.Distances[d] = (int32_t)(dim - distanceToCheck);
				dist += float(ks.Distances[d] * ks.Distances[d]);
			}
			ks.Weight = ComputeDistanceScore(ks.Distances, N_dimensions);
			dist = powf(dist, 1.f / float(N_dimensions));
			//if (dist <= float(distanceToCheck) + 0.001) // keep kernel round to lessen the number of samples
			{
				kernel.push_back(ks);
			}
		}

		std::chrono::milliseconds time_start_ms = std::chrono::duration_cast<std::chrono::milliseconds >(
		std::chrono::system_clock::now().time_since_epoch());

		float bestScore = std::numeric_limits<float>::max();

		std::vector<bool>   touchedElemBits(totalElements, false);
		std::vector<size_t> touchedElemIndex;


		auto CoordToIndex = [](const int32_t srcCoord[N_dimensions]) -> int32_t
		{
			int32_t index = 0;
			for (size_t d = 0; d < N_dimensions; ++d)
			{
				index += dimensionElementCount[d] * srcCoord[d];
			}
			return index;
		};

		auto IndexToCoord = [](int32_t srcIndex, int32_t dstCoord[N_dimensions]) -> void
		{
			for (size_t d = 0; d < N_dimensions; ++d)
			{
				dstCoord[d] = FastModulo(FastDiv(srcIndex, dimensionElementCount[d]), dimensionSize[d]);
			}
		};

		auto ComputeScore = [&](int32_t srcCoord[N_dimensions], uint32_t currArray) -> float
		{
			float score = 0.f;
			int32_t srcElem = CoordToIndex(srcCoord);
			for (const KernelSample &ks : kernel)
			{
				size_t j = 0;
				for (size_t d = 0; d < N_dimensions; ++d)
				{
					j += WrapDimension(srcCoord[d], ks.Distances[d], dimensionSize[d]) * dimensionElementCount[d];
				}
				if (srcElem == j)
					continue;
				score += ComputeFinalScore(pattern[currArray], ks.Weight, N_valuesPerItem, srcElem, j);
			}
			assert(score >= 0.f);
			return score;
		};

		auto ComputeScoreSSE = [&](int32_t srcCoord[N_dimensions], uint32_t currArray) -> float
		{
			int32_t srcElem = CoordToIndex(srcCoord);
			SSERegister score;
			score.v = _mm_set_ps1(0.f);
			size_t neighOffsets[4] = { 0 };
			__declspec(align(16)) float  distWeights[4] = { 0 };
			uint32_t vIndex = 0;
			for (const KernelSample &ks : kernel)
			{			
				size_t j = 0;
				for (size_t d = 0; d < N_dimensions; ++d)
				{
					j += WrapDimension(srcCoord[d], ks.Distances[d], dimensionSize[d]) * dimensionElementCount[d];
				}
				if (srcElem == j)
					continue;
				neighOffsets[vIndex] = j;
				distWeights[vIndex] = ks.Weight;
				++vIndex;
				if (vIndex == 4)
				{
					score.v = _mm_add_ps(score.v, ComputeFinalScoreSSE(pattern[currArray], distWeights, N_valuesPerItem, srcElem, neighOffsets));
					vIndex = 0;
					distWeights[0] = 0.f; distWeights[1] = 0.f; distWeights[2] = 0.f; distWeights[3] = 0.f;
				}
			}
			if (vIndex != 0)
			{
				score.v = _mm_add_ps(score.v, ComputeFinalScoreSSE(pattern[currArray], distWeights, N_valuesPerItem, srcElem, neighOffsets));
			}			
			return score.s[0] + score.s[1] + score.s[2] + score.s[3];
		};

		auto MarkModifiedElems = [&](int32_t srcCoord[N_dimensions]) -> void
		{
			float score = 0.f;
			for (const KernelSample &ks : kernel)
			{
				size_t j = 0;
				for (size_t d = 0; d < N_dimensions; ++d)
				{
					j += WrapDimension(srcCoord[d], ks.Distances[d], dimensionSize[d]) * dimensionElementCount[d];
				}
				if (touchedElemBits[j] == false)
				{
					touchedElemIndex.push_back(j);
					touchedElemBits[j] = true;
				}
			}
		};

		if (useIncrementalUpdate)
		{
			pattern[currentArray ^ 1] = pattern[currentArray]; // both array start equal
			bestScore = 0.f;
			for (size_t i = 0; i < totalElements; ++i)
			{
				int32_t srcCoord[N_dimensions];
				IndexToCoord(i, srcCoord);
				bestScore += ComputeScoreSSE(srcCoord, 0);
			}
		}

		bool actuallyUseIncrementalUpdate = useIncrementalUpdate;
		// I removed this part, because, in my experiment it is not so easy to make a choice when size of each dimensions are not equals 
		// (this could be improved / automated  later) For example, big 3D texture with few slices. 
		// For now the choice to use incremental version is left to the user. 
		// As a rule of thumb : the bigger the texture, the higher odds are that the incremental version will be faster.
		// Below are the measurement I made that indicate that incremental version is faster (for square and cubic cases)
		// They are left as a reference.
		#if 0		
			switch(N_dimensions)
			{
				case 1:
					// inconsequential here, no-op
				break;
				case 2: 
					if (totalElements < IntPow(18, 2)) actuallyUseIncrementalUpdate = false;  // for a 18x18 texture or bigger, incremental version faster
				break;
				case 3:
				case 4:
					if (totalElements < IntPow(12, N_dimensions)) actuallyUseIncrementalUpdate = false; /// incremental version version faster for 12x12x12 (and bigger) or 12x12x12x12 (and bigger)
				break;
				default:
					assert(0); // case not handled
				break;
			}
		#endif

		for (size_t iter = 0; iter < numIterationsToFindDistribution; ++iter)
		{
			if (actuallyUseIncrementalUpdate)
			{
				uint32_t num_swaps = distInt(gen);
				size_t swapedElemIndex[maxSwapedElemCount * 2];
				for (size_t i = 0; i < num_swaps; ++i)
				{
					size_t from = distSwap(gen);
					size_t to = distSwap(gen);
					while (from == to)
						to = distSwap(gen);

					swapedElemIndex[2 * i] = from;
					swapedElemIndex[2 * i + 1] = to;

					// mark region where score must be recomputed

					int32_t toCoord[N_dimensions];
					IndexToCoord(to, toCoord);
					int32_t fromCoord[N_dimensions];
					IndexToCoord(from, fromCoord);

					MarkModifiedElems(toCoord);
					MarkModifiedElems(fromCoord);

					for (size_t vecDim = 0; vecDim < N_valuesPerItem; ++vecDim)
					{
						std::swap(pattern[0][from * N_valuesPerItem + vecDim], pattern[0][to * N_valuesPerItem + vecDim]);
					}
				}

				float scoreToRemove = 0.f;
				for (size_t elemIndex : touchedElemIndex)
				{
					int32_t elemCoord[N_dimensions];
					IndexToCoord(elemIndex, elemCoord);
					scoreToRemove += ComputeScoreSSE(elemCoord, 1); // remove score from previous distribution
				}

				float scoreToAdd = 0.f;
				for (size_t elemIndex : touchedElemIndex)
				{
					int32_t elemCoord[N_dimensions];
					IndexToCoord(elemIndex, elemCoord);
					scoreToAdd += ComputeScoreSSE(elemCoord, 0); // add score from current distribution
					touchedElemBits[elemIndex] = false;
				}

				float deltaScore = scoreToAdd - scoreToRemove;
				touchedElemIndex.clear();

				if (deltaScore < 0.f)
				{
					bestScore += deltaScore;
					// commit changes to other array
					for (uint32_t i = 0; i < num_swaps * 2; ++i)
					{
						const int modifiedIndex = swapedElemIndex[i];
						for (size_t vecDim = 0; vecDim < N_valuesPerItem; ++vecDim)
						{
							pattern[1][modifiedIndex * N_valuesPerItem + vecDim] = pattern[0][modifiedIndex * N_valuesPerItem + vecDim];
						}
					}
				}
				else
				{
					// rollback changes from other array
					for (uint32_t i = 0; i < num_swaps * 2; ++i)
					{
						const int modifiedIndex = swapedElemIndex[i];
						for (size_t vecDim = 0; vecDim < N_valuesPerItem; ++vecDim)
						{
							pattern[0][modifiedIndex * N_valuesPerItem + vecDim] = pattern[1][modifiedIndex * N_valuesPerItem + vecDim];
						}
					}
				}
			}
			else
			{
				// version with global score update
				// copy
				pattern[currentArray ^ 1] = pattern[currentArray];

				uint32_t num_swaps = distInt(gen);
				for (size_t i = 0; i < num_swaps; ++i)
				{
					size_t from = distSwap(gen);
					size_t to = distSwap(gen);
					while (from == to)
						to = distSwap(gen);

					for (size_t vecDim = 0; vecDim < N_valuesPerItem; ++vecDim)
					{
						std::swap(pattern[currentArray][from * N_valuesPerItem + vecDim], pattern[currentArray][to * N_valuesPerItem + vecDim]);
					}
				}
				float score = 0.0f;
				for (size_t i = 0; i < totalElements; ++i)
				{
					int32_t elemCoord[N_dimensions];
					IndexToCoord(i, elemCoord);
					score += ComputeScoreSSE(elemCoord, currentArray);
				}

				if (score < bestScore)
				{
					bestScore = score;
				}
				else
				{
					// swap back
					currentArray ^= 1;
				}
			}

			if (iter>0 && (iter % (numIterationsToFindDistribution / 100) == 0))
			{
				// tmp tmp
				float globalScoreDebug = 0.f;
				for (size_t i = 0; i < totalElements; ++i)
				{
					int32_t elemCoord[N_dimensions];
					IndexToCoord(i, elemCoord);
					globalScoreDebug += ComputeScoreSSE(elemCoord, currentArray);
				}
				std::cout << "Ref score : " << globalScoreDebug << std::endl;

				std::chrono::milliseconds time_ms = std::chrono::duration_cast<std::chrono::milliseconds >(
					std::chrono::system_clock::now().time_since_epoch());

				std::chrono::milliseconds elapsed = time_ms - time_start_ms;

				float pct = static_cast<float>(iter) / static_cast<float>(numIterationsToFindDistribution);
				float est_remain = static_cast<float>(elapsed.count()) / pct * (1 - pct);
				est_remain /= 1000.0f;

				std::cout << iter << "/" << numIterationsToFindDistribution << " best score: " << bestScore << " eta: " << static_cast<int>(est_remain) << "s";
				std::cout << "\t(" << (static_cast<int>(est_remain) / (60 * 60)) << " h " << (static_cast<int>(est_remain / 60) % 60) << " m " << (static_cast<int>(est_remain) % 60) << " s)";
				std::cout <<  std::endl;
			}
		}
	}
	else if (chosenMethod == Method_HighPass)
	{
		const float* convArr = nullptr;

		if (N_dimensions == 1)
			convArr = convWeights1D;
		else if (N_dimensions == 2)
			convArr = convWeights2D;
		else
			convArr = convWeights3D;

		for (size_t iter = 0; iter < 4; iter++)
		{
			// copy
			pattern[currentArray ^ 1] = pattern[currentArray];

			for (size_t i = 0; i < totalElements; ++i)
			{
				for (size_t vectorItem = 0; vectorItem < N_valuesPerItem; ++vectorItem)
				{
					float convSum = 0.0f;
					for (size_t elem = 0; elem < convSizeTotal; ++elem)
					{
						size_t j = 0;

						for (size_t d = 0; d < N_dimensions; ++d)
						{
							size_t sourceDim = (i / dimensionElementCount[d]) % dimensionSize[d];
							size_t offsetDim = (elem / IntPow(convSize, d)) % convSize;

							int offset = (int)offsetDim - convSize / 2;

							j += WrapDimension(sourceDim, offset, dimensionSize[d]) * dimensionElementCount[d];
						}

						convSum += pattern[currentArray ^ 1][j * N_valuesPerItem + vectorItem] * convArr[elem];
					}

					pattern[currentArray][i * N_valuesPerItem + vectorItem] = convSum;
				}
			}

			UnifyHistogram(pattern[currentArray], N_valuesPerItem);
		}
	}
	else
	{
		std::abort();
	}

	//PrintCodeOutput("finalDist.txt", pattern[currentArray], "finalDist", true, dimensionSize, N_dimensions, N_valuesPerItem);
	//PrintWebGLOutput("webgl.txt", pattern[currentArray], "finalDist", dimensionSize, N_dimensions, N_valuesPerItem, 0, totalElements);

	if (N_dimensions == 2)
	{
		char filename[512];
		memset(filename, 0, 512);
		sprintf(filename, "output_%dx%d_uni.bmp", (int)dimensionSize, (int)dimensionSize);
		uint8_t* bytedata = FloatDataToBytes(pattern[currentArray], dimensionSize, N_valuesPerItem, false);
		stbi_write_bmp(filename, dimensionSize, dimensionSize, 3, bytedata);
		std::cout << "wrote " << filename << std::endl;
		delete[] bytedata;

		memset(filename, 0, 512);
		sprintf(filename, "output_%dx%d_tri.bmp", (int)dimensionSize, (int)dimensionSize);
		bytedata = FloatDataToBytes(pattern[currentArray], dimensionSize, N_valuesPerItem, true);
		stbi_write_bmp(filename, dimensionSize, dimensionSize, 3, bytedata);
		std::cout << "wrote " << filename << std::endl;
		delete[] bytedata;
	}

	return 0;
}
