#define _CRT_SECURE_NO_WARNINGS
#include "debug_opt.h"

#include "blue_noise_generator.h"
#include "blue_noise_export_util.h"

#include <chrono>
#include <sstream>
#include <iomanip>

BlueNoiseGeneratorParameters params;
BlueNoiseGenerator generator;


std::vector<float> InvertTex3DXZ(const std::vector<float> &tex, const BlueNoiseGeneratorParameters &params);
void OutputDebugTex3D(const std::vector<float> &blueNoise, const BlueNoiseGeneratorParameters &params, const std::string &filenameFront, const std::string &filenameSide);

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
	static std::string SecondsToHumanReadable(float seconds);
	static milliseconds GetCurrentTimeInMs();
	//
	virtual void OnProgress(size_t iterCount, double bestScore, size_t swapCount, size_t swapAttempt) override;
	virtual void OnStartWhiteNoiseGeneration() override;
	virtual void OnStartBlueNoiseGeneration() override;
	virtual void OnSliceGenerated(size_t sliceIndex, size_t totalSliceCount) override;
	virtual void OnSliceRefined(size_t sliceIndex, size_t totalSliceCount) override;
};

int main(int argc, char** argv)
{

	// 2D Texture example
	params.chosenMethod = BlueNoiseGeneratorParameters::Method_SolidAngle;
	params.N_dimensions = 2u;
	params.dimensionSize[0] = 128;
	params.dimensionSize[1] = 128;
	params.dimensionSize[2] = 0;
	params.dimensionSize[3] = 0;
	params.N_valuesPerItem = 1u;
	params.useMultithreading = true;
	params.useIncrementalUpdate = true; // required for multithreading
	params.numIterationsToFindDistribution = 4 * 256 * 1024;


	// 3D Texture example : use no wrap for Z and progressive approach with BlueNoiseGeneratorParameters::Method_IndependantSlices (otherwise very slow, and do not converge well)
	/*params.chosenMethod = BlueNoiseGeneratorParameters::Method_IndependantSlices;
	params.N_dimensions = 3u;
	params.dimensionSize[0] = 128;
	params.dimensionSize[1] = 128;
	params.dimensionSize[2] = 8;
	params.dimensionSize[3] = 0;
	params.N_valuesPerItem = 1u;
	params.useMultithreading = true;
	params.useIncrementalUpdate = true; // required for multithreading
	params.numIterationsToFindDistribution = 256 * 1024;*/

	std::vector<float> whiteNoise;
	std::vector<float> blueNoise;
	CBlueNoiseGenProgress progress(params.numIterationsToFindDistribution);
	BlueNoiseGenerator::EResult result = generator.GenerateBlueNoise(params, whiteNoise, blueNoise, &progress);
	switch (result)
	{
	case BlueNoiseGenerator::Result_DimensionSmallerThanKernelSize:
		std::cout << "Texture dimension must be greater or equal than " << generator.GetMinTextureSize() << " texels." << std::endl;
		getchar();
		return -1;
		break;
	}
	BlueNoiseExportUtil::SaveAsPPM(whiteNoise, "mt_white_noise.ppm", params.dimensionSize, params.N_dimensions, params.N_valuesPerItem);
	BlueNoiseExportUtil::SaveAsPPM(blueNoise, "mt_blue_noise.ppm", params.dimensionSize, params.N_dimensions, params.N_valuesPerItem);
	BlueNoiseExportUtil::PrintCodeOutput("noise.txt", blueNoise, "BlueNoise", false, params.dimensionSize, params.N_dimensions, params.N_valuesPerItem);
}

std::vector<float> InvertTex3DXZ(const std::vector<float> &tex, const BlueNoiseGeneratorParameters &params)
{
	size_t		swappedDimensionSize[BlueNoiseGeneratorParameters::max_N_dimensions];
	swappedDimensionSize[0] = params.dimensionSize[2];
	swappedDimensionSize[1] = params.dimensionSize[1];
	swappedDimensionSize[2] = params.dimensionSize[0];

	std::vector<float> sidePattern(tex.size());
	for (size_t z = 0u; z < swappedDimensionSize[2]; ++z)
	{
		for (size_t y = 0u; y < swappedDimensionSize[1]; ++y)
		{
			for (size_t x = 0u; x < swappedDimensionSize[0]; ++x)
			{
				for (size_t itemValueIndex = 0; itemValueIndex < params.N_valuesPerItem; ++itemValueIndex)
				{
					const float srcValue = tex[params.N_valuesPerItem * (x * (params.dimensionSize[0] * params.dimensionSize[1]) + y * params.dimensionSize[0] + z) + itemValueIndex];
					sidePattern[params.N_valuesPerItem * (z * (swappedDimensionSize[0] * swappedDimensionSize[1]) + y * swappedDimensionSize[0] + x) + itemValueIndex] = srcValue;
				}
			}
		}
	}
	return sidePattern;
}

void OutputDebugTex3D(const std::vector<float> &blueNoise, const BlueNoiseGeneratorParameters &params, const std::string &filenameFront, const std::string &filenameSide)
{
	// debugging code for output 3D texture
	// front
	if (params.N_dimensions == 3)
	{
		for (size_t z = 0u; z < params.dimensionSize[2]; ++z)
		{
			char buf[512];
			sprintf(buf, "%s_z%u.ppm", filenameFront.c_str(), z);
			BlueNoiseExportUtil::SaveAsPPM(blueNoise, buf, params.dimensionSize, params.N_dimensions, params.N_valuesPerItem, uint32_t(z));
		}
		// side
		size_t		swappedDimensionSize[BlueNoiseGeneratorParameters::max_N_dimensions];
		swappedDimensionSize[0] = params.dimensionSize[2];
		swappedDimensionSize[1] = params.dimensionSize[1];
		swappedDimensionSize[2] = params.dimensionSize[0];

		std::vector<float> sidePattern = InvertTex3DXZ(blueNoise, params);

		for (size_t z = 0u; z < params.dimensionSize[2]; ++z)
		{
			char buf[512];
			sprintf(buf, "%s_z%u.ppm", filenameSide.c_str(), z);
			BlueNoiseExportUtil::SaveAsPPM(sidePattern, buf, swappedDimensionSize, params.N_dimensions, params.N_valuesPerItem, uint32_t(z));
		}
	}
}

std::string CBlueNoiseGenProgress::SecondsToHumanReadable(float seconds)
{
	std::ostringstream outString;
	int im = static_cast<int>(seconds);
	outString << " (" << (im / (60 * 60)) << " h " << ((im / 60) % 60) << " m " << (im % 60) << " s)";
	return outString.str();
}
void CBlueNoiseGenProgress::OnProgress(size_t iterCount, double bestScore, size_t swapCount, size_t swapAttempt)
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
		outString << SecondsToHumanReadable(est_remain);
		//outString << "swap ratio : " << (double)swapCount / swapAttempt;
		outString << std::endl;
		std::cout << outString.str();
		//
		/*{
		std::vector<float> blueNoise;
		generator.GetCurrentBlueNoise(blueNoise);
		static int count = 0;
		std::ostringstream name;
		name << "mt_blue_noise_" << std::setfill('0') << std::setw(5) << count << ".ppm";
		BlueNoiseExportUtil::SaveAsPPM(blueNoise, name.str(), params.dimensionSize, params.N_dimensions, params.N_valuesPerItem);
		++count;
		//OutputDebugTex3D(blueNoise, params, "mt_blue_noise_front", "mt_blue_noise_side");
		}*/
	}
}
void CBlueNoiseGenProgress::OnStartWhiteNoiseGeneration()
{
	std::cout << "Blue noise generation started ..." << std::endl;
	// no-op
}
void CBlueNoiseGenProgress::OnStartBlueNoiseGeneration()
{
	time_start_ms = GetCurrentTimeInMs();
}
CBlueNoiseGenProgress::milliseconds CBlueNoiseGenProgress::GetCurrentTimeInMs()
{
	return std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch());
}
void CBlueNoiseGenProgress::OnSliceGenerated(size_t sliceIndex, size_t totalSliceCount)
{
	std::cout << sliceIndex << "/" << totalSliceCount << std::endl;
}
void CBlueNoiseGenProgress::OnSliceRefined(size_t sliceIndex, size_t totalSliceCount)
{
	// no-op
	if (sliceIndex != 0)
	{
		milliseconds time_ms = GetCurrentTimeInMs();
		milliseconds elapsed = time_ms - time_start_ms;
		float pct = static_cast<float>(sliceIndex) / static_cast<float>(totalSliceCount);
		float est_remain = static_cast<float>(elapsed.count()) / pct * (1 - pct);
		est_remain /= 1000.0f;
		// concatenate into a single stream to avoid multithreading problems
		std::ostringstream outString;
		outString << "(" << sliceIndex << "/" << totalSliceCount << ") ";
		outString << " eta: " << static_cast<int>(est_remain) << " s. " << "elapsed = " << elapsed.count() / 1000.f << " s.";
		outString << SecondsToHumanReadable(est_remain);
		outString << std::endl;
		std::cout << outString.str();
	}
	/*{
		std::ostringstream name;
		name << "mt_blue_noise_slice_" << std::setfill('0') << std::setw(5) << sliceIndex << ".ppm";
		std::vector<float> blueNoise;
		generator.GetCurrentBlueNoise(blueNoise);
		BlueNoiseExportUtil::SaveAsPPM(blueNoise, name.str(), params.dimensionSize, params.N_dimensions, params.N_valuesPerItem, sliceIndex);
	}
	{
		std::ostringstream name;
		name << "mt_blue_noise_slice_side_" << std::setfill('0') << std::setw(5) << sliceIndex << ".ppm";
		std::vector<float> blueNoise;
		generator.GetCurrentBlueNoise(blueNoise);

		size_t		swappedDimensionSize[BlueNoiseGeneratorParameters::max_N_dimensions];
		swappedDimensionSize[0] = params.dimensionSize[2];
		swappedDimensionSize[1] = params.dimensionSize[1];
		swappedDimensionSize[2] = params.dimensionSize[0];

		BlueNoiseExportUtil::SaveAsPPM(InvertTex3DXZ(blueNoise, params), name.str(), swappedDimensionSize, params.N_dimensions, params.N_valuesPerItem, 0);
	}*/

	/*std::ostringstream name;
	name << "mt_blue_noise_slice_" << std::setfill('0') << std::setw(5) << sliceIndex << ".ppm";
	std::vector<float> blueNoise;
	generator.GetCurrentBlueNoise(blueNoise);
	BlueNoiseExportUtil::SaveAsPPM(blueNoise, name.str(), params.dimensionSize, params.N_dimensions, params.N_valuesPerItem);*/
}