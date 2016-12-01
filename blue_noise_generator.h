#ifndef BLUE_NOISE_GENERATOR_H
#define BLUE_NOISE_GENERATOR_H

#include <memory>
#include <vector>

class BlueNoiseGeneratorParameters;


struct IBlueNoiseGenProgressMonitor
{
	virtual void OnStartWhiteNoiseGeneration() = 0;
	virtual void OnStartBlueNoiseGeneration() = 0;
	virtual void OnProgress(size_t iterCount, double bestScore, size_t swapCount, size_t swapAttempt) = 0;
	virtual void OnSliceGenerated(size_t sliceIndex, size_t totalSliceCount) = 0;
	virtual void OnSliceRefined(size_t sliceIndex, size_t totalSliceCount) = 0;
};


typedef float Real;

class BlueNoiseGenerator
{
public:
	enum EResult
	{
		Result_OK,
		Result_DimensionSmallerThanKernelSize,
	};
	BlueNoiseGenerator();
	~BlueNoiseGenerator();
	EResult GenerateBlueNoise(const BlueNoiseGeneratorParameters	&generationParameters,
								std::vector<float>					&whiteNoiseResult,
								std::vector<float>					&blueNoiseResult,
								IBlueNoiseGenProgressMonitor		*progressMonitor = nullptr);
	static uint32_t GetMinTextureSize();
	void GetCurrentBlueNoise(std::vector<float> &dest) const;
private:
	class BlueNoiseGeneratorImpl *pImpl;
};

#endif // BLUE_NOISE_GENERATOR_H
