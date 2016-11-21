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
	BlueNoiseGenerator();
	~BlueNoiseGenerator();
	void GenerateBlueNoise(const BlueNoiseGeneratorParameters	&generationParameters,
						   std::vector<float>					&whiteNoiseResult,
						   std::vector<float>					&blueNoiseResult,
						   IBlueNoiseGenProgressMonitor			*progressMonitor = nullptr);
private:
	class BlueNoiseGeneratorImpl *pImpl;
};

#endif // BLUE_NOISE_GENERATOR_H
