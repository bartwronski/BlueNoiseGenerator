#include <random>
#include <string>
#include <iostream>
#include <cmath>
#include <iomanip>

size_t IntPow(size_t base, size_t exp)
{
	size_t ret = 1;
	while (exp--)
	{
		ret *= base;
	}
	return ret;
}


const size_t dimensionSize = 16;
const size_t N_dimensions = 2;
const size_t N_valuesPerItem = 1;
const size_t totalElements = IntPow(dimensionSize, N_dimensions);
const size_t distanceToCheck = 3; // original paper mentioned looking at whole array; however this is N^2 and super expensive, while exp(-4*4/2.2) ~= 0.000694216
const size_t distanceToCheckBoth = distanceToCheck * 2 + 1; // in both directions

const size_t elementsToCheck = IntPow(distanceToCheckBoth, N_dimensions);

const size_t numAttemptToFindNoBias = 8;
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

			if (dim == dimensionSize-1)
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

std::string dimNames[4] = {"x", "y", "z", "w"};

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

inline float ComputeSimilarityScore(const std::vector<float>& arr, size_t N_valuesPerItem, size_t ind1, size_t ind2)
{
	float score = 0;
	for (size_t i = 0; i < N_valuesPerItem; ++i)
	{
		float val = (arr[ind1 * N_valuesPerItem + i] - arr[ind2 * N_valuesPerItem + i]);
		score += val * val;
	}
	// Orig paper formula 1.0f/sqrt(distance ^ vector_size)
	// 1.0f/sqrt(sqrt(distanceSq) ^ vector_size)
	// 1.0f/(distanceSq^(vector_size/4)
	return powf(score, -(float)N_valuesPerItem / 4.0f);
}

inline float ComputeDistanceScore(const int arr[], size_t Ndimensions)
{
	float distanceSq = 0;

	for (size_t i = 0; i < Ndimensions; ++i)
	{
		distanceSq += arr[i] * arr[i];
	}
	// Note 1.1 - variance parameter recommended in original paper
	return expf(-distanceSq / (2.0f * 1.1f));
}

int main(int argc, char** argv)
{
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dist(0, 1);
	// Note: we try to swap between 1 and 3 elements to try to jump over local minima
	std::uniform_int_distribution<> distInt(1, 3);
	std::uniform_int_distribution<> distSwap(0, (int)(totalElements - 1));

	std::vector<float> pattern[2] = {std::vector<float>(totalElements * N_valuesPerItem), std::vector<float>(totalElements * N_valuesPerItem)};
	uint32_t currentArray = 0;

	float currentSum = 1000000.0f;
	for (size_t it = 0; it < numAttemptToFindNoBias; ++it)
	{
		// copy
		pattern[currentArray ^ 1] = pattern[currentArray];

		float sum = 0.0f;
		for (size_t i = 0; i < totalElements * N_valuesPerItem; ++i)
		{
			float val = static_cast<float>(dist(gen));
			pattern[currentArray][i] = val;
			sum += val;
		}
		const float desiredSum = totalElements * N_valuesPerItem * 0.5f;
		if (fabs(currentSum - desiredSum) < fabs(sum - desiredSum))
		{
			// swap back
			currentArray ^= 1;
		}
		else
		{
			// std::cout << sum << std::endl;
			currentSum = sum;
		}
	}
	PrintCodeOutput(pattern[currentArray], "initialDist", true, dimensionSize, N_dimensions, N_valuesPerItem);

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

	float bestScore = 100000.0f;
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
					size_t sourceDim = (i	 / IntPow(dimensionSize, d))		% dimensionSize;
					size_t offsetDim = (elem / IntPow(distanceToCheckBoth, d))	% distanceToCheckBoth;

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

				float similarityScore = ComputeSimilarityScore(pattern[currentArray], N_valuesPerItem, i, j);
				score += similarityScore * distanceWeights[elem];
			}
		}
		
		if (score < bestScore)
		{
			// std::cout << score << " " << bestScore << " " << iter << std::endl;
			bestScore = score;
		}
		else
		{
			// swap back
			currentArray ^= 1;
		}
	}

	PrintCodeOutput(pattern[currentArray], "finalDist", true, dimensionSize, N_dimensions, N_valuesPerItem);
	PrintWebGLOutput(pattern[currentArray], "finalDist", dimensionSize, N_dimensions, N_valuesPerItem);

	return 0;
}