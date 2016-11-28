#ifndef BLUE_NOISE_GENERATOR_PARAMETERS_H
#define BLUE_NOISE_GENERATOR_PARAMETERS_H

#include "config.h"
#include <stdint.h>

class BlueNoiseGeneratorParameters
{
public:
	static const size_t max_N_dimensions = 4;
	enum EMethod
	{
		Method_WhiteNoise,
		Method_SolidAngle,
		Method_IndependantSlices,
		Method_HighPass
	};
	//
	EMethod			chosenMethod;
	size_t			N_dimensions;
	size_t			dimensionSize[max_N_dimensions];
	size_t			N_valuesPerItem;
	bool			useMultithreading;
	bool			useIncrementalUpdate; // required for multithreading
	size_t			numIterationsToFindDistribution;
	int32_t			refineSpecificSlice;
public:
	// ctor with default values
	BlueNoiseGeneratorParameters();
};


#endif // BLUE_NOISE_GENERATOR_PARAMETERS_H
