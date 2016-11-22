#include "blue_noise_generator.h"

#include <vector>
#include <atomic>
#include <thread>
#include <chrono>
#include <algorithm>
#include <random>

#include "config.h"
#include "blue_noise_generator_parameters.h"
#include "math_util.h"
#include "sse_util.h"


// private implementation for blue noise generator pImpl

struct KernelSample
{
	float	Weight;
	int32_t Distances[BlueNoiseGeneratorParameters::max_N_dimensions];
	int32_t DeltaIndex;
};

///////////
// PIMPL //
///////////

class BlueNoiseGeneratorImpl
{
public:
	BlueNoiseGeneratorImpl();
	BlueNoiseGenerator::EResult	GenerateBlueNoise(	const BlueNoiseGeneratorParameters	&generationParameters,
													std::vector<float>					&whiteNoiseResult,
													std::vector<float>					&blueNoiseResult,
													IBlueNoiseGenProgressMonitor		*progressMonitor = nullptr);
	static uint32_t GetMinTextureSize() { return uint32_t(_DistanceToCheckBoth); }
private:
	// utility / coordinate functions
	inline size_t	ComputeElementCount(size_t dimCount) const;
	inline size_t	WrapDimension(size_t baseIndex, int offset, size_t dimSize) const;
	inline int32_t	CoordToIndex(const int32_t srcCoord[]) const;
	inline void		IndexToCoord(int32_t srcIndex, int32_t dstCoord[]) const;

	// blue noise heuristic functions
	inline float	ComputeDistanceScore(const int arr[]) const;
	inline float	ComputeFinalLocalScore(const std::vector<float>& arr, float distanceScore, size_t ind1, size_t ind2) const;
	float			ComputeLocalScore(const int32_t srcCoord[], uint32_t currArray) const;
	float			ComputeLocalScoreScalar(const int32_t srcCoord[], uint32_t currArray) const;
	float			ComputeGlobalScore(uint32_t currArray) const;
	// sse specifics
	#ifdef USE_SSE
		inline __m128	ComputeFinalLocalScoreSSE(const std::vector<float>& arr, float distanceScore[4], size_t ind1, size_t ind2[4]) const;
		float			ComputeLocalScoreSSE(const int32_t srcCoord[], uint32_t currArray) const;
	#endif

	// main computation functions
	void			GenerateInitialWhiteNoise();
	void			PrecomputeKernel();
	void			ComputeBlueNoise(size_t numIter);
	void			MarkModifiedElems(int32_t srcCoord[], std::vector<bool> &touchedElemBits, std::vector<size_t> &touchedElemIndex);
	void			ComputeBlueNoiseIncremental(size_t numIter);
	void			ComputeBlueNoiseIncrementalMultiThreaded(size_t numIter);
	void			UnifyHistogram(std::vector<float>& arr);
	void			DoHighPass();
	void			NextIter(float deltaScore);

	// multi-threading handling
	uint64_t		ComputeMTRegionAcquisitionMask(size_t elemIndex) const;
	void			AcquireMTRegion(uint64_t claimedRegions);
	void			ReleaseMTRegion(uint64_t claimedRegions);
	uint32_t		ComputeAcquiredMTRegionSizeDivisorAsRShift() const;
	void			LockIterGuard();
	void			UnlockIterGuard();

private:
	static const size_t max_N_dimensions = BlueNoiseGeneratorParameters::max_N_dimensions;

	// internal working variables
	BlueNoiseGeneratorParameters		_GenParams;
	static const int					_MaxSwapedElemCount;
	std::vector<float>					_Pattern[2];
	uint32_t							_CurrentArray;
	std::vector<KernelSample>			_Kernel;
	size_t								_TotalElements;
	size_t								_DimensionElementCount[max_N_dimensions];
	float								_BestScore;
	volatile uint32_t					_IterTotal;
	bool								_ActuallyUseMultithreading;
	uint32_t							_AcquiredMTRegionSizeDivisorAsRShift;
	std::atomic<uint_least64_t>			_MTAcquiredRegions;
	std::atomic<int>					_IterGuard;
	IBlueNoiseGenProgressMonitor		*_ProgressMonitor;

	// Solid Angle method parameters
	static const size_t _DistanceToCheck;
	static const size_t _DistanceToCheckBoth;
	size_t				_ElementsToCheck;

	// Highpass filter parameters
	static const size_t _ConvSize;
	static const float  _ConvWeights1D[3];
	static const float  _ConvWeights2D[3 * 3];
	static const float  _ConvWeights3D[3 * 3 * 3];
	size_t			    _ConvSizeTotal;
};

/////////////////////
// STATICS MEMBERS //
/////////////////////

// Note: we try to swap between 1 and 3 elements to try to jump over local minimum
const int BlueNoiseGeneratorImpl::_MaxSwapedElemCount = 3u;

// Solid Angle method parameters
const size_t BlueNoiseGeneratorImpl::_DistanceToCheck	= 3; // original paper mentioned looking at whole array; however this is N^2 and super expensive, while exp(-4*4/2.2) ~= 0.000694216
const size_t BlueNoiseGeneratorImpl::_DistanceToCheckBoth = BlueNoiseGeneratorImpl::_DistanceToCheck * 2u + 1u; // in both directions

// Highpass filter parameters
const size_t BlueNoiseGeneratorImpl::_ConvSize = 3u;
const float BlueNoiseGeneratorImpl::_ConvWeights1D[3] = { -1, 2, -1 };
const float BlueNoiseGeneratorImpl::_ConvWeights2D[3 * 3] =
{
	-1, -2, -1,
	-2, 12, -2,
	-1, -2, -1
};
const float BlueNoiseGeneratorImpl::_ConvWeights3D[3 * 3 * 3] =
{
	-1, -2, -1,
	-2, -4, -2,
	-1, -2, -1,
	-2, -4, -2,
	-4, 56, -4,
	-2, -4, -2,
	-1, -2, -1,
	-2, -4, -2,
	-1, -2, -1,
};

////////////////////
// IMPLEMENTATION //
////////////////////

//===========================================================================================================================
BlueNoiseGeneratorImpl::BlueNoiseGeneratorImpl()
{
	_ConvSizeTotal             = 0u;
	_CurrentArray              = 0u;
	_TotalElements             = 0u;
	_ActuallyUseMultithreading = false;
	_ProgressMonitor           = nullptr;
	std::fill(std::begin(_DimensionElementCount), std::end(_DimensionElementCount), 0u);
	_BestScore = std::numeric_limits<float>::max();
}

//===========================================================================================================================
size_t BlueNoiseGeneratorImpl::ComputeElementCount(size_t dimCount) const
{
	size_t elemCount = 1;
	for (size_t currDim = 0; currDim < dimCount; ++currDim)
	{
		elemCount *= _GenParams.dimensionSize[currDim];
	}
	return elemCount;
}

//===========================================================================================================================
inline size_t BlueNoiseGeneratorImpl::WrapDimension(size_t baseIndex, int offset, size_t dimSize) const
{
	int posWrapped = (int)baseIndex + offset;
	if (posWrapped < 0)
	{
		posWrapped += int(dimSize);
	}
	if (posWrapped > (int)(dimSize - 1))
	{
		posWrapped -= int(dimSize);
	}
	return posWrapped;
}

//===========================================================================================================================
inline int32_t BlueNoiseGeneratorImpl::CoordToIndex(const int32_t srcCoord[]) const
{
	int32_t index = 0;
	for (size_t d = 0; d < _GenParams.N_dimensions; ++d)
	{
		index += int32_t(_DimensionElementCount[d] * srcCoord[d]);
	}
	return index;
}

//===========================================================================================================================
inline void BlueNoiseGeneratorImpl::IndexToCoord(int32_t srcIndex, int32_t dstCoord[]) const
{
	for (size_t d = 0; d < _GenParams.N_dimensions; ++d)
	{
		dstCoord[d] = FastModulo(FastDiv(srcIndex, int(_DimensionElementCount[d])), int(_GenParams.dimensionSize[d]));
	}
}

//===========================================================================================================================
inline float BlueNoiseGeneratorImpl::ComputeFinalLocalScore(const std::vector<float>& arr, float distanceScore, size_t ind1, size_t ind2) const
{
	const size_t N_valuesPerItem = _GenParams.N_valuesPerItem;
	float valueSpaceScore = 0;
	for (size_t i = 0; i < N_valuesPerItem; ++i)
	{
		float val = (arr[ind1 * N_valuesPerItem + i] - arr[ind2 * N_valuesPerItem + i]);
		valueSpaceScore += val * val;
	}
	valueSpaceScore = FastPowScalar(valueSpaceScore, (float)N_valuesPerItem / 2.f);
	const float oneOverDistanceVarianceSq = 1.0f / (2.1f * 2.1f);
	return FastExp(-valueSpaceScore - distanceScore * oneOverDistanceVarianceSq);
}

//===========================================================================================================================
#ifdef USE_SSE
	inline __m128 BlueNoiseGeneratorImpl::ComputeFinalLocalScoreSSE(const std::vector<float>& arr, float distanceScore[4], size_t ind1, size_t ind2[4]) const
	{
		const size_t N_valuesPerItem = _GenParams.N_valuesPerItem;
		SSERegister valueSpaceScore;
		valueSpaceScore.v = _mm_set_ps1(0.f);
		const float oneOverDistanceVarianceSq = 1.0f / (2.1f * 2.1f);
		__m128 mOneOverDistanceVarianceSq = _mm_set_ps1(oneOverDistanceVarianceSq);
		for (size_t i = 0; i < N_valuesPerItem; ++i)
		{
			float srcValue = arr[ind1 * N_valuesPerItem + i];
			float val0 = srcValue - arr[ind2[0] * N_valuesPerItem + i];
			float val1 = srcValue - arr[ind2[1] * N_valuesPerItem + i];
			float val3 = srcValue - arr[ind2[3] * N_valuesPerItem + i];
			float val2 = srcValue - arr[ind2[2] * N_valuesPerItem + i];
			__m128 val = _mm_set_ps(val0, val1, val2, val3);
			valueSpaceScore.v = Mad(val, val, valueSpaceScore.v);
		}
		if (_GenParams.N_valuesPerItem == 1)
		{
			valueSpaceScore.v = Abs(valueSpaceScore.v);
		}
		valueSpaceScore.v = FastPowSSEVector(valueSpaceScore.v, (float)N_valuesPerItem / 2.f);
		return FastExpSSEVector(Negate(Mad(_mm_set_ps(distanceScore[0], distanceScore[1], distanceScore[2], distanceScore[3]), mOneOverDistanceVarianceSq, valueSpaceScore.v)));
	}
#endif

//===========================================================================================================================
inline float BlueNoiseGeneratorImpl::ComputeDistanceScore(const int arr[]) const
{
	float distanceSq = 0;
	for (size_t i = 0; i < _GenParams.N_dimensions; ++i)
	{
		distanceSq += arr[i] * arr[i];
	}
	return distanceSq;
}

//===========================================================================================================================
float BlueNoiseGeneratorImpl::ComputeLocalScoreScalar(const int32_t srcCoord[], uint32_t currArray) const
{
	float score = 0.f;
	int32_t srcElem = CoordToIndex(srcCoord);
	for (const KernelSample &ks : _Kernel)
	{
		size_t j = 0;
		for (size_t d = 0; d < _GenParams.N_dimensions; ++d)
		{
			j += WrapDimension(srcCoord[d], ks.Distances[d], _GenParams.dimensionSize[d]) * _DimensionElementCount[d];
		}
		if (srcElem == j)
			continue;
		score += ComputeFinalLocalScore(_Pattern[currArray], ks.Weight, srcElem, j);
	}
	assert(score >= 0.f);
	return score;
};

//===========================================================================================================================
#ifdef USE_SSE
	float BlueNoiseGeneratorImpl::ComputeLocalScoreSSE(const int32_t srcCoord[], uint32_t currArray) const
	{
		int32_t srcElem = CoordToIndex(srcCoord);
		SSERegister score;
		score.v = _mm_set_ps1(0.f);
		size_t neighOffsets[4] = { 0 };
		__declspec(align(16)) float  distWeights[4] = { 0 }; // TODO : portable stuff here for alignment
		uint32_t vIndex = 0;
		bool needWrap = false;
		int32_t baseOffset = 0;
		for (size_t d = 0; d < _GenParams.N_dimensions; ++d)
		{
			if (srcCoord[d] - int32_t(_DistanceToCheck) < 0)
			{
				needWrap = true;
				break;
			}
			if (srcCoord[d] + int32_t(_DistanceToCheck) >= int32_t(_GenParams.dimensionSize[d]))
			{
				needWrap = true;
				break;
			}
			baseOffset += srcCoord[d] * _DimensionElementCount[d];
		}
		for (const KernelSample &ks : _Kernel)
		{
			size_t j = 0;
			if (needWrap)
			{
				for (size_t d = 0; d < _GenParams.N_dimensions; ++d)
				{
					j += WrapDimension(srcCoord[d], ks.Distances[d], _GenParams.dimensionSize[d]) * _DimensionElementCount[d];
				}
			}
			else
			{
				j = baseOffset + ks.DeltaIndex;
			}
			if (srcElem == j)
				continue;
			neighOffsets[vIndex] = j;
			distWeights[vIndex] = ks.Weight;
			++vIndex;
			if (vIndex == 4)
			{
				score.v = _mm_add_ps(score.v, ComputeFinalLocalScoreSSE(_Pattern[currArray], distWeights, srcElem, neighOffsets));
				vIndex = 0;
				distWeights[0] = 0.f; distWeights[1] = 0.f; distWeights[2] = 0.f; distWeights[3] = 0.f;
			}
		}
		if (vIndex != 0) // not a multiple of 4 ? handle remaining
		{
			score.v = _mm_add_ps(score.v, ComputeFinalLocalScoreSSE(_Pattern[currArray], distWeights, srcElem, neighOffsets));
		}
		return score.s[0] + score.s[1] + score.s[2] + score.s[3];
	}
#endif

//===========================================================================================================================
inline float BlueNoiseGeneratorImpl::ComputeLocalScore(const int32_t srcCoord[], uint32_t currArray) const
{
	#ifdef USE_SSE
		return ComputeLocalScoreSSE(srcCoord, currArray);
	#else
		return ComputeLocalScoreScalar(srcCoord, currArray);
	#endif
}

//===========================================================================================================================
float	BlueNoiseGeneratorImpl::ComputeGlobalScore(uint32_t currArray) const
{
	float score = 0.f;
	for (size_t i = 0; i < _TotalElements; ++i)
	{
		int32_t srcCoord[max_N_dimensions];
		IndexToCoord(int32_t(i), srcCoord);
#ifdef USE_SSE
		score += ComputeLocalScoreSSE(srcCoord, 0);
#else
		score += ComputeLocalScore(srcCoord, 0);
#endif
	}
	return score;
}

//===========================================================================================================================
//note: see http://gpuopen.com/vdr-follow-up-fine-art-of-film-grain/
void BlueNoiseGeneratorImpl::UnifyHistogram(std::vector<float>& arr)
{
	const size_t N_valuesPerItem = _GenParams.N_valuesPerItem;
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
		for (size_t i = 0, n = entries.size(); i < n; ++i)
		{
			float t = static_cast<float>(i) / static_cast<float>(n - 1);
			size_t idx = entries[i].second;
			arr[idx * N_valuesPerItem + dim] = t;
		}
	}
}

//===========================================================================================================================
void BlueNoiseGeneratorImpl::DoHighPass()
{
	const size_t N_valuesPerItem = _GenParams.N_valuesPerItem;
	_ConvSizeTotal = IntPow(_ConvSize, _GenParams.N_dimensions);
	const float* _ConvArr = nullptr;
	if (_GenParams.N_dimensions == 1)
		_ConvArr = _ConvWeights1D;
	else if (_GenParams.N_dimensions == 2)
		_ConvArr = _ConvWeights2D;
	else
		_ConvArr = _ConvWeights3D;
	for (size_t iter = 0; iter < 4; iter++)
	{
		// copy
		_Pattern[_CurrentArray ^ 1] = _Pattern[_CurrentArray];
		for (size_t i = 0; i < _TotalElements; ++i)
		{
			for (size_t vectorItem = 0; vectorItem < N_valuesPerItem; ++vectorItem)
			{
				float convSum = 0.0f;
				for (size_t elem = 0; elem < _ConvSizeTotal; ++elem)
				{
					size_t j = 0;
					for (size_t d = 0; d < _GenParams.N_dimensions; ++d)
					{
						size_t sourceDim = (i / _DimensionElementCount[d]) % _GenParams.dimensionSize[d];
						size_t offsetDim = (elem / IntPow(_ConvSize, d)) % _ConvSize;
						int offset = (int)offsetDim - _ConvSize / 2;
						j += WrapDimension(sourceDim, offset, _GenParams.dimensionSize[d]) * _DimensionElementCount[d];
					}
					convSum += _Pattern[_CurrentArray ^ 1][j * N_valuesPerItem + vectorItem] * _ConvArr[elem];
				}
				_Pattern[_CurrentArray][i * N_valuesPerItem + vectorItem] = convSum;
			}
		}
		UnifyHistogram(_Pattern[_CurrentArray]);
	}
}

//===========================================================================================================================
BlueNoiseGenerator::EResult BlueNoiseGeneratorImpl::GenerateBlueNoise(	const BlueNoiseGeneratorParameters &generationParams,
																		std::vector<float> &whiteNoiseResult,
																		std::vector<float> &blueNoiseResult,
																		IBlueNoiseGenProgressMonitor *progressMonitor)
{
	_GenParams = generationParams;
	// current limitation : each dimension must be bigger than the kernel size for the wrap to work properly
	for (uint32_t dim = 0; dim < _GenParams.N_dimensions; ++dim)
	{
		if (_GenParams.dimensionSize[dim] < _DistanceToCheckBoth) return BlueNoiseGenerator::Result_DimensionSmallerThanKernelSize;
	}
	_ProgressMonitor = progressMonitor;
	_TotalElements = ComputeElementCount(_GenParams.N_dimensions);
	// precompute elements per dimension for speed
	for (uint32_t dim = 0; dim <= _GenParams.N_dimensions; ++dim)
	{
		_DimensionElementCount[dim] = ComputeElementCount(dim);
	}

	_Pattern[0] = _Pattern[1] = std::vector<float>(_TotalElements * _GenParams.N_valuesPerItem);
	_CurrentArray = 0;
	_ElementsToCheck = IntPow(_DistanceToCheckBoth, _GenParams.N_dimensions);

	if (_ProgressMonitor) _ProgressMonitor->OnStartWhiteNoiseGeneration();

	GenerateInitialWhiteNoise();
	std::vector<float> whiteNoiseHolder = _Pattern[_CurrentArray];

	auto CommitResult = [&]() -> void
	{
		blueNoiseResult.swap(_Pattern[_CurrentArray]);
		whiteNoiseResult.swap(whiteNoiseHolder);
		_Pattern[0].clear();
		_Pattern[1].clear();
	};

	//////////////////////
	// high pass method //
	//////////////////////

	if (_GenParams.chosenMethod == BlueNoiseGeneratorParameters::Method_HighPass)
	{
		DoHighPass();
		CommitResult();
		return BlueNoiseGenerator::Result_OK;;
	}

	// TODO : possibly better to start from high pass distribution, add this as an enum ?
	 // DoHighPass();
	 // _Pattern[_CurrentArray ^ 1] = _Pattern[_CurrentArray];

	////////////////////////
	// solid angle method //
	////////////////////////

	assert(_GenParams.chosenMethod == BlueNoiseGeneratorParameters::Method_SolidAngle);
	PrecomputeKernel();
	if (_GenParams.useIncrementalUpdate)
	{
		// compute initial score
		_Pattern[_CurrentArray ^ 1] = _Pattern[_CurrentArray]; // both array start equal
		_BestScore = ComputeGlobalScore(0);
		_ActuallyUseMultithreading = _GenParams.useMultithreading && (_GenParams.N_dimensions == 2 || _GenParams.N_dimensions == 3);
		if (_ProgressMonitor) _ProgressMonitor->OnStartBlueNoiseGeneration();
		if (_ActuallyUseMultithreading)
		{
			ComputeBlueNoiseIncrementalMultiThreaded(_GenParams.numIterationsToFindDistribution);
		}
		else
		{
			ComputeBlueNoiseIncremental(_GenParams.numIterationsToFindDistribution);
		}
		CommitResult();
	}
	else
	{
		_BestScore = std::numeric_limits<float>::max();
		if (_ProgressMonitor) _ProgressMonitor->OnStartBlueNoiseGeneration();
		ComputeBlueNoise(_GenParams.numIterationsToFindDistribution);
		CommitResult();
	}
	return BlueNoiseGenerator::Result_OK;
}

//===========================================================================================================================
void BlueNoiseGeneratorImpl::GenerateInitialWhiteNoise()
{
	const size_t N_valuesPerItem = _GenParams.N_valuesPerItem;
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dist(0, 1);
	// white noise
	for (size_t i = 0; i < _TotalElements * N_valuesPerItem; ++i)
	{
		_Pattern[_CurrentArray][i] = static_cast<float>(dist(gen));
	}
	UnifyHistogram(_Pattern[_CurrentArray]);
}

//===========================================================================================================================
void BlueNoiseGeneratorImpl::PrecomputeKernel()
{
	_Kernel.clear();
	for (size_t i = 0; i < _ElementsToCheck; ++i)
	{
		KernelSample ks;
		float dist = 0.f;
		int32_t deltaIndex = 0;
		for (size_t d = 0; d < _GenParams.N_dimensions; ++d)
		{
			size_t dim = (i / IntPow(_DistanceToCheckBoth, d)) % _DistanceToCheckBoth;
			ks.Distances[d] = (int32_t)(dim - _DistanceToCheck);
			dist += float(ks.Distances[d] * ks.Distances[d]);
			deltaIndex += ks.Distances[d] * _DimensionElementCount[d];
		}
		ks.Weight = ComputeDistanceScore(ks.Distances);
		ks.DeltaIndex = deltaIndex;
		//float score = expf(-ks.Weight / (2.1f * 2.1f));
		//if (score > 0.008f)
		{
			_Kernel.push_back(ks);
		}
	}
}

//===========================================================================================================================
void BlueNoiseGeneratorImpl::MarkModifiedElems(int32_t srcCoord[], std::vector<bool> &touchedElemBits, std::vector<size_t> &touchedElemIndex)
{
	float score = 0.f;
	for (const KernelSample &ks : _Kernel)
	{
		size_t j = 0;
		for (size_t d = 0; d < _GenParams.N_dimensions; ++d)
		{
			j += WrapDimension(srcCoord[d], ks.Distances[d], _GenParams.dimensionSize[d]) * _DimensionElementCount[d];
		}
		if (touchedElemBits[j] == false)
		{
			touchedElemIndex.push_back(j);
			touchedElemBits[j] = true;
		}
	}
}

//===========================================================================================================================
void BlueNoiseGeneratorImpl::ComputeBlueNoise(size_t numIter)
{
	_IterTotal = 0;
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dist(0, 1);
	std::uniform_int_distribution<> distInt(1, _MaxSwapedElemCount);
	std::uniform_int_distribution<> distSwap(0, (int)(_TotalElements - 1));
	// version with global score update
	for (size_t iter = 0; iter < numIter; ++iter)
	{
		// copy
		_Pattern[_CurrentArray ^ 1] = _Pattern[_CurrentArray];
		uint32_t num_swaps = distInt(gen);
		for (size_t i = 0; i < num_swaps; ++i)
		{
			size_t from = distSwap(gen);
			size_t to = distSwap(gen);
			while (from == to)
				to = distSwap(gen);
			for (size_t vecDim = 0; vecDim < _GenParams.N_valuesPerItem; ++vecDim)
			{
				std::swap(_Pattern[_CurrentArray][from * _GenParams.N_valuesPerItem + vecDim], _Pattern[_CurrentArray][to * _GenParams.N_valuesPerItem + vecDim]);
			}
		}
		const float score = ComputeGlobalScore(_CurrentArray);
		if (score < _BestScore)
		{
			_BestScore = score;
		}
		else
		{
			// swap back
			_CurrentArray ^= 1;
		}
		++_IterTotal;
		if (_ProgressMonitor)
		{
			_ProgressMonitor->OnProgress(_IterTotal, _BestScore);
		}
	}
}

//===========================================================================================================================
void BlueNoiseGeneratorImpl::ComputeBlueNoiseIncrementalMultiThreaded(size_t numIter)
{
	 _MTAcquiredRegions.store(0ul);
	 _IterGuard.store(0);
	_IterTotal = 0;

	_AcquiredMTRegionSizeDivisorAsRShift = ComputeAcquiredMTRegionSizeDivisorAsRShift();
	size_t numThread = std::thread::hardware_concurrency();
	numThread = std::max(size_t(1u), numThread);

	std::vector<std::shared_ptr<std::thread> > threads;

	auto DoWork = [&](size_t numIter) -> void
	{
		this->ComputeBlueNoiseIncremental(numIter);
	};

	volatile bool finished = false;
	// monitoring thread
	auto MonitorThreadFunc = [&]() -> void
	{
		if (_ProgressMonitor)
		{
			while (!finished)
			{
				LockIterGuard();
				_ProgressMonitor->OnProgress(_IterTotal, _BestScore);
				UnlockIterGuard();
				std::this_thread::sleep_for(std::chrono::milliseconds(20)); // check once in a while
			}
		}
	};

	std::thread monitorThread(MonitorThreadFunc);

	for (size_t threadIndex = 0u; threadIndex < numThread; ++threadIndex)
	{
		threads.push_back(std::shared_ptr<std::thread>(new std::thread(DoWork, (numIter + numThread - 1) / numThread)));
	}
	for (auto& th : threads)
	{
		th->join();
	}
	finished = true;
	monitorThread.join();
	threads.clear();
}

//===========================================================================================================================
uint32_t BlueNoiseGeneratorImpl::ComputeAcquiredMTRegionSizeDivisorAsRShift() const
{
	size_t cubicSize = 0;
	uint32_t acquiredRegionSizeShift = 0u; // how many bit to shift a coord to create the region acquisition mask
	for (uint32_t dim = 0; dim < _GenParams.N_dimensions; ++dim)
	{
		cubicSize = std::max(cubicSize, _GenParams.dimensionSize[dim]);
	}
	cubicSize = NextPow2(uint32_t(cubicSize));
	if (_GenParams.N_dimensions == 2)
	{
		acquiredRegionSizeShift = (cubicSize <= 8) ? 0 : uint32_t(log2(float(cubicSize / 8)));
	}
	else
	{
		assert(_GenParams.N_dimensions == 3);
		acquiredRegionSizeShift = (cubicSize <= 4) ? 0 : uint32_t(log2(float(cubicSize / 4)));
	}
	return acquiredRegionSizeShift;
}

//===========================================================================================================================
uint64_t BlueNoiseGeneratorImpl::ComputeMTRegionAcquisitionMask(size_t elemIndex) const
{
	int32_t coord[max_N_dimensions];
	IndexToCoord(int32_t(elemIndex), coord);
	uint64_t claimedRegions(0u);
	for (const KernelSample &ks : _Kernel)
	{
		size_t bitIndex = 0;
		if (_GenParams.N_dimensions == 2)
		{
			// claim part of 8 x 8 grid
			bitIndex = (WrapDimension(coord[0], ks.Distances[0], _GenParams.dimensionSize[0]) >> _AcquiredMTRegionSizeDivisorAsRShift) +
				8 * (WrapDimension(coord[1], ks.Distances[1], _GenParams.dimensionSize[1]) >> _AcquiredMTRegionSizeDivisorAsRShift);
			assert(bitIndex < 64u);
			claimedRegions |= uint64_t(1u) << bitIndex;
		}
		else
		{
			assert(_GenParams.N_dimensions == 3);
			// claim part of 4 x 4 x 4 grid
			const size_t wrappedDim[] =
			{
				WrapDimension(coord[0], ks.Distances[0], _GenParams.dimensionSize[0]),
				WrapDimension(coord[1], ks.Distances[1], _GenParams.dimensionSize[1]),
				WrapDimension(coord[2], ks.Distances[2], _GenParams.dimensionSize[2])
			};
			bitIndex = (wrappedDim[0] >> _AcquiredMTRegionSizeDivisorAsRShift) +
						4 * (wrappedDim[1] >> _AcquiredMTRegionSizeDivisorAsRShift) +
						16 * (wrappedDim[2] >> _AcquiredMTRegionSizeDivisorAsRShift);

			assert(bitIndex < 64u);
			claimedRegions |= uint64_t(1u) << bitIndex;
		}
	}
	return claimedRegions;
}

//===========================================================================================================================
void BlueNoiseGeneratorImpl::AcquireMTRegion(uint64_t claimedRegions)
{
	while (true)
	{
		uint64_t previousBits = _MTAcquiredRegions.fetch_or(claimedRegions);
		if ((previousBits & claimedRegions) == uint64_t(0)) break; // all acquired at once ?
		uint64_t partialAcquiredBits = (~previousBits) & claimedRegions;
		// not all regions could be acquired simultaneously... release and spin lock
		_MTAcquiredRegions.fetch_and(~partialAcquiredBits);
	}
}

//===========================================================================================================================
void BlueNoiseGeneratorImpl::ReleaseMTRegion(uint64_t claimedRegions)
{
	uint64_t oldBits = _MTAcquiredRegions.fetch_and(~claimedRegions);
	assert((oldBits & claimedRegions) == claimedRegions);
}

//===========================================================================================================================
void BlueNoiseGeneratorImpl::LockIterGuard()
{
	while (true)
	{
		int oldValue = _IterGuard.fetch_add(1);
		if (oldValue == 0) break;
		_IterGuard.fetch_add(-1);
	}
}

//===========================================================================================================================
void BlueNoiseGeneratorImpl::UnlockIterGuard()
{
	_IterGuard.fetch_add(-1);
}

//===========================================================================================================================
void BlueNoiseGeneratorImpl::NextIter(float deltaScore)
{
	if (_ActuallyUseMultithreading)
	{
		LockIterGuard();
		_BestScore += deltaScore;
		++_IterTotal;
		UnlockIterGuard();
	}
	else
	{
		_BestScore += deltaScore;
		++_IterTotal;
	}
}

//===========================================================================================================================
void BlueNoiseGeneratorImpl::ComputeBlueNoiseIncremental(size_t numIter)
{
	_IterTotal = 0;

	std::vector<bool>   touchedElemBits(_TotalElements, false);
	std::vector<size_t> touchedElemIndex;

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dist(0, 1);
	std::uniform_int_distribution<> distInt(1, _MaxSwapedElemCount);
	std::uniform_int_distribution<> distSwap(0, (int)(_TotalElements - 1));

	for (size_t iter = 0; iter < numIter; ++iter)
	{
		uint32_t num_swaps = _ActuallyUseMultithreading ? 1u : distInt(gen);
		size_t swapedElemIndex[_MaxSwapedElemCount * 2];
		size_t from[_MaxSwapedElemCount];
		size_t to[_MaxSwapedElemCount];
		uint64_t claimedRegions;

		while (true)
		{
			claimedRegions = uint64_t(0u);
			// compute the swaps, possibly building a bitfield or the regions that need to be acquired (for multi-threaded version)
			for (size_t i = 0; i < num_swaps; ++i)
			{
				from[i] = distSwap(gen);
				to[i] = distSwap(gen);
				while (from[i] == to[i])
					to[i] = distSwap(gen);
				// update claimed region bitfield
				if (_ActuallyUseMultithreading)
				{
					claimedRegions |= ComputeMTRegionAcquisitionMask(from[i]);
					claimedRegions |= ComputeMTRegionAcquisitionMask(to[i]);
				}
			}
			if (!_ActuallyUseMultithreading) break;
			// because acquisition is very costly, diminish the chance of collision by doing a test
			// with no atomic
			if ((claimedRegions & _MTAcquiredRegions.load()) == uint64_t(0u)) break;
		}

		// acquire the region to be swapped
		if (_ActuallyUseMultithreading)
		{
			AcquireMTRegion(claimedRegions);
			// no they are acquired, be can do read-write op without colliding with another thread :)
		}

		// do the actual swaps
		for (size_t i = 0; i < num_swaps; ++i)
		{
			swapedElemIndex[2 * i] = from[i];
			swapedElemIndex[2 * i + 1] = to[i];
			// mark region where score must be recomputed
			int32_t toCoord[max_N_dimensions];
			IndexToCoord(int32_t(to[i]), toCoord);
			int32_t fromCoord[max_N_dimensions];
			IndexToCoord(int32_t(from[i]), fromCoord);
			MarkModifiedElems(toCoord, touchedElemBits, touchedElemIndex);
			MarkModifiedElems(fromCoord, touchedElemBits, touchedElemIndex);
			for (size_t vecDim = 0; vecDim < _GenParams.N_valuesPerItem; ++vecDim)
			{
				std::swap(_Pattern[0][from[i] * _GenParams.N_valuesPerItem + vecDim], _Pattern[0][to[i] * _GenParams.N_valuesPerItem + vecDim]);
			}
		}

		// compute score delta
		// we do it in local neighbourhood only, because  elements further than 3 / 4 cell from here do not contribute significantly

		float scoreToRemove = 0.f;
		for (size_t elemIndex : touchedElemIndex)
		{
			int32_t elemCoord[max_N_dimensions];
			IndexToCoord(int32_t(elemIndex), elemCoord);
			scoreToRemove += ComputeLocalScore(elemCoord, 1);
		}

		float scoreToAdd = 0.f;
		for (size_t elemIndex : touchedElemIndex)
		{
			int32_t elemCoord[max_N_dimensions];
			IndexToCoord(int32_t(elemIndex), elemCoord);
			scoreToAdd += ComputeLocalScore(elemCoord, 0); // add score from current distribution
			touchedElemBits[elemIndex] = false;
		}

		float deltaScore = scoreToAdd - scoreToRemove;
		touchedElemIndex.clear();

		auto CopyModifications = [&](uint32_t from, uint32_t to) -> void
		{
			for (size_t swapIndex = 0; swapIndex< num_swaps * 2; ++swapIndex)
			{
				size_t modifiedIndex = swapedElemIndex[swapIndex];
				for (size_t vecDim = 0; vecDim < _GenParams.N_valuesPerItem; ++vecDim)
				{
					_Pattern[to][modifiedIndex * _GenParams.N_valuesPerItem + vecDim] = _Pattern[from][modifiedIndex * _GenParams.N_valuesPerItem + vecDim];
				}
			}
		};

		if (deltaScore < 0.f)
		{
			NextIter(deltaScore);
			// commit changes to other array
			CopyModifications(0, 1);
		}
		else
		{
			NextIter(0.f);
			// rollback changes from other array
			CopyModifications(1, 0);
		}

		if (_ActuallyUseMultithreading)
		{
			ReleaseMTRegion(claimedRegions);
		}

		if (!_ActuallyUseMultithreading)
		{
			if (_ProgressMonitor)
			{
				_ProgressMonitor->OnProgress(_IterTotal, _BestScore);
			}
		}
	}
}

//===========================================================================================================================
BlueNoiseGenerator::BlueNoiseGenerator() : pImpl(new BlueNoiseGeneratorImpl)
{
}

//===========================================================================================================================
BlueNoiseGenerator::~BlueNoiseGenerator()
{
	delete pImpl;
}

//===========================================================================================================================
BlueNoiseGenerator::EResult BlueNoiseGenerator::GenerateBlueNoise(	const BlueNoiseGeneratorParameters &generationParameters,
																	std::vector<float> &whiteNoiseResult,
																	std::vector<float> &blueNoiseResult,
																	IBlueNoiseGenProgressMonitor *progressMonitor)
{
	return pImpl->GenerateBlueNoise(generationParameters, whiteNoiseResult, blueNoiseResult, progressMonitor);
}

//===========================================================================================================================
uint32_t BlueNoiseGenerator::GetMinTextureSize()
{
	return BlueNoiseGeneratorImpl::GetMinTextureSize();
}
