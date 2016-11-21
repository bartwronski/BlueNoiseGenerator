#ifndef BLUE_NOISE_EXPORT_UTIL_H
#define BLUE_NOISE_EXPORT_UTIL_H

#define _CRT_SECURE_NO_DEPRECATE // avoid security warning on fopen and the like on windows

#include <string>
#include <vector>
#include <stdint.h>
#include <iostream>
#include <fstream>

#include "config.h"
#include "blue_noise_generator_parameters.h"

class BlueNoiseExportUtil
{
public:
	enum { max_N_dimensions = BlueNoiseGeneratorParameters::max_N_dimensions };

	static void PrintCodeOutput(const std::string&			fileName,
								const std::vector<float>&	arr,
								const std::string&			arrName,
								bool						mathematica,
								const size_t				dimensionSize[max_N_dimensions],
								size_t						N_dimensions,
								size_t						N_valuesPerItem);

	static void PrintWebGLOutput(const std::string&			fileName,
								 const std::vector<float>&	arr,
								 const std::string&			arrName,
								 size_t						N_dimensions,
								 size_t						N_valuesPerItem,
								 size_t						lo,
								 size_t						high);

	static void SaveAsPPM(const std::vector<float>&			arr,
						  const std::string&				fileName,
						  const size_t						dimensionSize[max_N_dimensions],
						  size_t							N_dimensions,
						  size_t							N_valuesPerItem,
						  uint32_t							slice = 0u);

	static void SaveAsBMP(const std::vector<float>&			arr,
						  const std::string&				fileName,
						  const size_t						dimensionSize[2],
						  size_t							N_valuesPerItem);
};









#endif
