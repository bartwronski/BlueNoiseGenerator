#define _CRT_SECURE_NO_WARNINGS
#include "debug_opt.h"

#include "blue_noise_generator.h"
#include "blue_noise_export_util.h"

#include <chrono>
#include <sstream>

class CBlueNoiseGenProgress : public IBlueNoiseGenProgressMonitor
{
public:
	CBlueNoiseGenProgress(size_t totalNumIter) : numIterationsToFindDistribution(totalNumIter), currPercent(0)
	{}
private:
	typedef std::chrono::milliseconds milliseconds;

	std::chrono::milliseconds	time_start_ms;
	size_t						numIterationsToFindDistribution;
	size_t                      currPercent;

private:
	virtual void OnProgress(size_t iterCount, float bestScore) override
	{
		size_t newPercent = size_t(double(iterCount) / (double(numIterationsToFindDistribution) / 100.0));
		if (iterCount > 0 && newPercent != currPercent)
		{
			currPercent = newPercent;
			milliseconds time_ms = GetCurrentTimeInMs();
			milliseconds elapsed = time_ms - time_start_ms;
			float pct = static_cast<float>(iterCount) / static_cast<float>(numIterationsToFindDistribution);
			float est_remain = static_cast<float>(elapsed.count()) / pct * (1 - pct);
			est_remain /= 1000.0f;
			// concatenate into a single stream to avoid multithreading problems
			std::ostringstream outString;
			outString << iterCount << "/" << numIterationsToFindDistribution << " best score: " << bestScore << " eta: " << static_cast<int>(est_remain) << " s. " << "elapsed = " << elapsed.count() / 1000.f << " s.";
			outString << " (" << (static_cast<int>(est_remain) / (60 * 60)) << " h " << (static_cast<int>(est_remain / 60) % 60) << " m " << (static_cast<int>(est_remain) % 60) << " s)";
			outString << std::endl;
			std::cout << outString.str();
		}
	}
	virtual void OnStartWhiteNoiseGeneration() override
	{
		// no-op
	}
	virtual void OnStartBlueNoiseGeneration() override
	{
		time_start_ms = GetCurrentTimeInMs();
	}
	static milliseconds GetCurrentTimeInMs()
	{
		return std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch());
	}
};


int main(int argc, char** argv)
{
	BlueNoiseGeneratorParameters params;
	params.chosenMethod = BlueNoiseGeneratorParameters::Method_SolidAngle;
	//params.chosenMethod = BlueNoiseGeneratorParameters::Method_HighPass;
	params.N_dimensions = 2u;
	params.dimensionSize[0] = 32;
	params.dimensionSize[1] = 32;
	params.dimensionSize[2] = 0;
	params.dimensionSize[3] = 0;
	params.N_valuesPerItem = 1u;
	params.useMultithreading = true;
	params.useIncrementalUpdate = true; // required for multithreading
	params.numIterationsToFindDistribution = 256 * 1024;

	std::vector<float> whiteNoise;
	std::vector<float> blueNoise;

	CBlueNoiseGenProgress progress(params.numIterationsToFindDistribution);
	BlueNoiseGenerator generator;
	generator.GenerateBlueNoise(params, whiteNoise, blueNoise, &progress);


	BlueNoiseExportUtil::SaveAsPPM(whiteNoise, "white_noise.ppm", params.dimensionSize, params.N_dimensions, params.N_valuesPerItem);
	//BlueNoiseExportUtil::SaveAsBMP(whiteNoise, "white_noise.bmp", params.dimensionSize, params.N_valuesPerItem);
	BlueNoiseExportUtil::SaveAsPPM(blueNoise, "blue_noise.ppm", params.dimensionSize, params.N_dimensions, params.N_valuesPerItem);

	// debugging code for output 3D texture
	// front
/*	for (size_t z = 0u; z < params.dimensionSize[2]; ++z)
	{
		char buf[512];
		sprintf(buf, "blue_noise_z%u.ppm", z);
		BlueNoiseExportUtil::SaveAsPPM(blueNoise, buf, params.dimensionSize, params.N_dimensions, params.N_valuesPerItem, uint32_t(z));
	}
	// side
	std::vector<float> sidePattern(blueNoise.size());
	for (size_t z = 0u; z < params.dimensionSize[2]; ++z)
	{
		for (size_t y = 0u; y < params.dimensionSize[1]; ++y)
		{
			for (size_t x = 0u; x < params.dimensionSize[0]; ++x)
			{
				sidePattern[z * (params.dimensionSize[0] * params.dimensionSize[1]) + y * params.dimensionSize[0] + x] = blueNoise[x * (params.dimensionSize[0] * params.dimensionSize[1]) + y * params.dimensionSize[0] + z];
			}
		}
	}

	for (size_t z = 0u; z < params.dimensionSize[2]; ++z)
	{
		char buf[512];
		sprintf(buf, "blue_noise_side_z%u.ppm", z);
		BlueNoiseExportUtil::SaveAsPPM(sidePattern, buf, params.dimensionSize, params.N_dimensions, params.N_valuesPerItem, uint32_t(z));
	}*/

}
