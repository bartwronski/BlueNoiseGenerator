#include <random>
#include <string>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <iomanip>
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

const size_t dimensionSize = 32;
const size_t N_dimensions = 2;
const size_t N_valuesPerItem = 1;
const size_t totalElements = IntPow(dimensionSize, N_dimensions);
const size_t distanceToCheck = 3; // original paper mentioned looking at whole array; however this is N^2 and super expensive, while exp(-4*4/2.2) ~= 0.000694216
const size_t distanceToCheckBoth = distanceToCheck * 2 + 1; // in both directions

const size_t elementsToCheck = IntPow(distanceToCheckBoth, N_dimensions);

const size_t numIterationsToFindDistribution = 256 * 1024;


void PrintCodeOutput(const std::vector<float>& arr, const std::string& arrName, bool mathematica, size_t dimensionSize, size_t N_dimensions, size_t N_valuesPerItem)
{
	const size_t arrSize = arr.size();

	if (mathematica)
	{
		std::cout << arrName << " = " << std::endl;
	}
	else
	{
		std::cout << "static const float " << arrName;
		for (size_t d = 0; d < N_dimensions; ++d)
		{
			std::cout << "[" << dimensionSize << "]";
		}
		if (N_valuesPerItem > 1)
		{
			std::cout << "[" << N_valuesPerItem << "]";
		}
		std::cout << " = " << std::endl;
	}


	for (size_t i = 0; i < arrSize / N_valuesPerItem; ++i)
	{
		for (size_t d = 0; d < N_dimensions; ++d)
		{
			size_t dim = (i / IntPow(dimensionSize, d)) % dimensionSize;

			if (dim == 0)
			{
				std::cout << "{";
			}
			else
			{
				break;
			}
		}

		if (N_valuesPerItem == 1)
		{
			std::cout << arr[i];
		}
		else
		{
			std::cout << "{";
			for (size_t v = 0; v < N_valuesPerItem; ++v)
			{
				std::cout << arr[i * N_valuesPerItem + v];
				if (v < N_valuesPerItem - 1)
				{
					std::cout << ", ";
				}
			}
			std::cout << "}";
		}

		for (size_t d = 0; d < N_dimensions; ++d)
		{
			size_t dim = (i / IntPow(dimensionSize, d)) % dimensionSize;

			if (dim == dimensionSize - 1)
			{
				std::cout << "}";
			}
			else
			{
				std::cout << ",";
				if (d > 0)
				{
					std::cout << std::endl;
				}
				break;
			}
		}
	}
	if (!mathematica)
	{
		std::cout << ";";
	}

	std::cout << std::endl << std::endl;
}

std::string dimNames[4] = { "x", "y", "z", "w" };

void PrintWebGLOutputRecursive(const std::vector<float>& arr, const std::string& arrName, size_t dimensionSize, size_t N_dimensions, size_t N_valuesPerItem, size_t lo, size_t high)
{
	if (high - lo == 1)
	{
		std::cout << "return ";

		if (N_valuesPerItem > 1)
		{
			std::cout << "vec" << N_valuesPerItem << "(";
			for (size_t v = 0; v < N_valuesPerItem; ++v)
			{
				std::cout << arr[lo * N_valuesPerItem + v];
				if (v < N_valuesPerItem - 1)
				{
					std::cout << ", ";
				}
			}

			std::cout << ");";
		}
		else
		{
			std::cout << arr[lo] << ";";
		}
	}
	else
	{
		size_t mid = (lo + high) / 2;

		std::cout << "if(" << arrName << " < " << mid << ") " << std::endl << "{" << std::endl;
		PrintWebGLOutputRecursive(arr, arrName, dimensionSize, N_dimensions, N_valuesPerItem, lo, mid);
		std::cout << "} else {" << std::endl;
		PrintWebGLOutputRecursive(arr, arrName, dimensionSize, N_dimensions, N_valuesPerItem, mid, high);
		std::cout << std::endl << "}";
	}
};

void PrintWebGLOutput(const std::vector<float>& arr, const std::string& arrName, size_t dimensionSize, size_t N_dimensions, size_t N_valuesPerItem)
{
	for (size_t i = 0; i < arr.size() / N_valuesPerItem; ++i)
	{
		std::cout << "if (";

		for (size_t d = 0; d < N_dimensions; ++d)
		{
			size_t dim = (i / IntPow(dimensionSize, d)) % dimensionSize;
			std::cout << arrName << "." << dimNames[d] << " == " << std::setw(3) << dim;
			if (d < N_dimensions - 1)
			{
				std::cout << " && ";
			}
		}
		std::cout << ") return ";

		if (N_valuesPerItem > 1)
		{
			std::cout << "vec" << N_valuesPerItem << "(";
			for (size_t v = 0; v < N_valuesPerItem; ++v)
			{
				std::cout << std::setw(6) << arr[i * N_valuesPerItem + v];
				if (v < N_valuesPerItem - 1)
				{
					std::cout << ", ";
				}
			}

			std::cout << ");" << std::endl;
		}
		else
		{
			std::cout << arr[i] << ";" << std::endl;
		}
	}
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
	uint8_t* bytes = new uint8_t[3 * arr.size()];

	for (size_t i = 0, n = arr.size(); i<n; ++i)
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
			float t = static_cast<float>(i) / static_cast<float>(n);
			size_t idx = entries[i].second;
			arr[idx * N_valuesPerItem + dim] = t;
		}
	}
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

	//PrintCodeOutput(pattern[currentArray], "initialDist", true, dimensionSize, N_dimensions, N_valuesPerItem);

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

					int finalPositionCompare = (int)(sourceDim + offset);
					if (finalPositionCompare < 0)
					{
						finalPositionCompare += dimensionSize;
					}
					if (finalPositionCompare > dimensionSize - 1)
					{
						finalPositionCompare -= dimensionSize;
					}

					j += finalPositionCompare * IntPow(dimensionSize, d);
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

	//PrintCodeOutput(pattern[currentArray], "finalDist", true, dimensionSize, N_dimensions, N_valuesPerItem);
	//PrintWebGLOutputRecursive(pattern[currentArray], "finalDist", dimensionSize, N_dimensions, N_valuesPerItem, 0, totalElements);

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
