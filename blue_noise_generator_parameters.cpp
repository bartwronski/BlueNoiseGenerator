#include "blue_noise_generator_parameters.h"

BlueNoiseGeneratorParameters::BlueNoiseGeneratorParameters()
{
	chosenMethod = Method_SolidAngle;
	N_dimensions						= 2u;
	dimensionSize[0] = 64;
	dimensionSize[1] = 64;
	dimensionSize[2] = 0;
	dimensionSize[3] = 0;
	N_valuesPerItem						= 1u;
	useMultithreading					= true;
	useIncrementalUpdate				= true;
	numIterationsToFindDistribution		= 256 * 1024;
}

