#ifndef BLUE_NOISE_GENERATOR_H
#define BLUE_NOISE_GENERATOR_H

#include <memory>
#include <vector>

class BlueNoiseGeneratorParameters;


struct IBlueNoiseGenProgressMonitor
{
	virtual void OnStartWhiteNoiseGeneration() = 0;
	virtual void OnStartBlueNoiseGeneration() = 0;
	virtual void OnProgress(size_t iterCount, float bestScore) = 0;
};

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
private:
	class BlueNoiseGeneratorImpl *pImpl;
};

#endif // BLUE_NOISE_GENERATOR_H
