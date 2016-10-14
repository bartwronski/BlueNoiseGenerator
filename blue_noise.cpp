#include <random>
#include <string>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <chrono>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

size_t IntPow(size_t base, size_t exp)
{
	size_t ret = 1;
	while (exp--)
	{
		ret *= base;
	}
	return ret;
}

enum EMethod
{
	Method_SolidAngle,
	Method_HighPass,
};

const size_t dimensionSize = 32;
const size_t N_dimensions = 2;
const size_t N_valuesPerItem = 1;
const size_t totalElements = IntPow(dimensionSize, N_dimensions);

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

void PrintCodeOutput(const std::string& fileName, const std::vector<float>& arr, const std::string& arrName, bool mathematica, size_t dimensionSize, size_t N_dimensions, size_t N_valuesPerItem)
{
	std::ofstream outFile;
	outFile.open(fileName, std::ios::out | std::ios::trunc);

	const size_t arrSize = arr.size();

	if (!mathematica)
	{
		outFile << "static const float " << arrName;
		for (size_t d = 0; d < N_dimensions; ++d)
		{
			outFile << "[" << dimensionSize << "]";
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
			size_t dim = (i / IntPow(dimensionSize, d)) % dimensionSize;

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
			size_t dim = (i / IntPow(dimensionSize, d)) % dimensionSize;

			if (dim == dimensionSize - 1)
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

void PrintWebGLOutputRecursive(std::ofstream& outFile, const std::vector<float>& arr, const std::string& arrName, size_t dimensionSize, size_t N_dimensions, size_t N_valuesPerItem, size_t lo, size_t high)
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
		PrintWebGLOutputRecursive(outFile, arr, arrName, dimensionSize, N_dimensions, N_valuesPerItem, lo, mid);
		outFile << "} else {" << std::endl;
		PrintWebGLOutputRecursive(outFile, arr, arrName, dimensionSize, N_dimensions, N_valuesPerItem, mid, high);
		outFile << std::endl << "}";
	}
};

void PrintWebGLOutput(const std::string& fileName, const std::vector<float>& arr, const std::string& arrName, size_t dimensionSize, size_t N_dimensions, size_t N_valuesPerItem, size_t lo, size_t high)
{
	std::ofstream outFile;
	outFile.open(fileName, std::ios::out | std::ios::trunc);

	PrintWebGLOutputRecursive(outFile, arr, arrName, dimensionSize, N_dimensions, N_valuesPerItem, lo, high);
}

inline float ComputeFinalScore(const std::vector<float>& arr, float distanceScore, size_t N_valuesPerItem, size_t ind1, size_t ind2)
{
	float valueSpaceScore = 0;
	for (size_t i = 0; i < N_valuesPerItem; ++i)
	{
		float val = (arr[ind1 * N_valuesPerItem + i] - arr[ind2 * N_valuesPerItem + i]);
		valueSpaceScore += val * val;
	}

	valueSpaceScore = powf(valueSpaceScore, (float)N_valuesPerItem / 2.0f);
	const float oneOverDistanceVarianceSq = 1.0f / (2.1f * 2.1f);

	return expf(-valueSpaceScore - distanceScore * oneOverDistanceVarianceSq);
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

float sign(float v) { return (v >= 0) ? 1.0f : -1.0f; }
//note: remaps uniform input to triangular value
//      see https://www.shadertoy.com/view/4t2SDh
float remap_tri(float v)
{
	// Convert uniform distribution into triangle-shaped distribution.
	float orig = v*2.0f - 1.0f;
	v = orig / sqrtf(fabs(orig));
	v = (v < -1.0f) ? -1.0f : v; // Nerf the NaN generated by 0*rsqrt(0). Thanks @FioraAeterna!
	v = v - sign(orig) + 0.5f;

	// result is range [-0.5,1.5] which is useful for actual dithering.
	// convert to [0,1] for output
	return (v + 0.5f) * 0.5f;
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
			float v_f = arr[i * N_valuesPerItem + j];

			if (do_remap_tri)
				v_f = remap_tri(v_f);

			int v_i = static_cast<int>(v_f * 255.0f + 0.5f);
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

		std::sort(entries.begin(), entries.end(), [](const auto& a, const auto & b) -> bool
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

int main(int argc, char** argv)
{
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dist(0, 1);
	// Note: we try to swap between 1 and 3 elements to try to jump over local minima
	std::uniform_int_distribution<> distInt(1, 3);
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
		std::vector<float> distanceWeights(elementsToCheck);

		for (size_t i = 0; i < elementsToCheck; ++i)
		{
			int distances[N_dimensions];
			for (size_t d = 0; d < N_dimensions; ++d)
			{
				size_t dim = (i / IntPow(distanceToCheckBoth, d)) % distanceToCheckBoth;
				distances[d] = (int)dim - distanceToCheck;
			}

			distanceWeights[i] = ComputeDistanceScore(distances, N_dimensions);
		}

		std::chrono::milliseconds time_start_ms = std::chrono::duration_cast<std::chrono::milliseconds >(
			std::chrono::system_clock::now().time_since_epoch());

		float bestScore = std::numeric_limits<float>::max();
		for (size_t iter = 0; iter < numIterationsToFindDistribution; ++iter)
		{
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
				for (size_t elem = 0; elem < elementsToCheck; ++elem)
				{
					size_t j = 0;

					for (size_t d = 0; d < N_dimensions; ++d)
					{
						size_t sourceDim = (i / IntPow(dimensionSize, d)) % dimensionSize;
						size_t offsetDim = (elem / IntPow(distanceToCheckBoth, d)) % distanceToCheckBoth;

						int offset = (int)offsetDim - distanceToCheck;

						j += WrapDimension(sourceDim, offset, dimensionSize) * IntPow(dimensionSize, d);
					}

					if (i == j)
						continue;

					score += ComputeFinalScore(pattern[currentArray], distanceWeights[elem], N_valuesPerItem, i, j);
				}
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


			if (iter>0 && (iter % (numIterationsToFindDistribution / 100) == 0))
			{
				std::chrono::milliseconds time_ms = std::chrono::duration_cast<std::chrono::milliseconds >(
					std::chrono::system_clock::now().time_since_epoch());

				std::chrono::milliseconds elapsed = time_ms - time_start_ms;

				float pct = static_cast<float>(iter) / static_cast<float>(numIterationsToFindDistribution);
				float est_remain = static_cast<float>(elapsed.count()) / pct * (1 - pct);
				est_remain /= 1000.0f;

				std::cout << iter << "/" << numIterationsToFindDistribution << " best score: " << bestScore << " eta: " << static_cast<int>(est_remain) << "s" << std::endl;
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
							size_t sourceDim = (i / IntPow(dimensionSize, d)) % dimensionSize;
							size_t offsetDim = (elem / IntPow(convSize, d)) % convSize;

							int offset = (int)offsetDim - convSize / 2;

							j += WrapDimension(sourceDim, offset, dimensionSize) * IntPow(dimensionSize, d);
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
