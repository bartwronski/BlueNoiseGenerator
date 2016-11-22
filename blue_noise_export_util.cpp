#include "blue_noise_export_util.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#include <iomanip>


// some utility functions
namespace
{
	//===========================================================================================================================
	inline uint32_t FloatAsByteUNorm(float value)
	{
		return uint32_t(255.f * value);
	}

	//===========================================================================================================================
	double remap_tri(const double v)
	{
		double r2 = 0.5 * v;
		double f1 = sqrt(r2);
		double f2 = 1.0 - sqrt(r2 - 0.25);
		return (v < 0.5) ? f1 : f2;
	}

	//===========================================================================================================================
	size_t ComputeElementCount(const size_t	dimensionSize[BlueNoiseExportUtil::max_N_dimensions], size_t dimCount)
	{
		size_t elemCount = 1;
		for (size_t currDim = 0; currDim < dimCount; ++currDim)
		{
			elemCount *= dimensionSize[currDim];
		}
		return elemCount;
	}

	//===========================================================================================================================
	//note: splats to 24b
	uint8_t* FloatDataToBytes(const std::vector<float>& arr, size_t N_valuesPerItem, bool do_remap_tri)
	{
		uint8_t* bytes = new uint8_t[3 * arr.size() / N_valuesPerItem];
		for (size_t i = 0, n = arr.size() / N_valuesPerItem; i < n; ++i)
		{
			uint8_t lastChar = 0;
			for (size_t j = 0; j < N_valuesPerItem; j++)
			{
				double v_f = arr[i * N_valuesPerItem + j];
				if (do_remap_tri)
					v_f = remap_tri(v_f);
				int v_i = static_cast<int>(v_f * 256.0);
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

	//===========================================================================================================================
	void PrintWebGLOutputRecursive(std::ofstream&				outFile,
									const std::vector<float>&	arr,
									const std::string&			arrName,
									size_t						N_dimensions,
									size_t						N_valuesPerItem,
									size_t						lo,
									size_t						high)
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
			PrintWebGLOutputRecursive(outFile, arr, arrName, N_dimensions, N_valuesPerItem, lo, mid);
			outFile << "} else {" << std::endl;
			PrintWebGLOutputRecursive(outFile, arr, arrName, N_dimensions, N_valuesPerItem, mid, high);
			outFile << std::endl << "}";
		}
	}

}

//===========================================================================================================================
void BlueNoiseExportUtil::PrintCodeOutput (const std::string&			fileName,
											const std::vector<float>&	arr,
											const std::string&			arrName,
											bool						mathematica,
											const size_t				dimensionSize[max_N_dimensions],
											size_t						N_dimensions,
											size_t						N_valuesPerItem)
{
	std::ofstream outFile;
	outFile.open(fileName, std::ios::out | std::ios::trunc);
	const size_t arrSize = arr.size();
	if (!mathematica)
	{
		outFile << "static const float " << arrName;
		for (size_t d = 0; d < N_dimensions; ++d)
		{
			outFile << "[" << dimensionSize[d] << "]";
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
			size_t dim = (i / ComputeElementCount(dimensionSize, d)) % dimensionSize[d];
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
			size_t dim = (i / ComputeElementCount(dimensionSize, d)) % dimensionSize[d];

			if (dim == dimensionSize[d] - 1)
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

//===========================================================================================================================
void BlueNoiseExportUtil::PrintWebGLOutput(const std::string&			fileName,
											const std::vector<float>&	arr,
											const std::string&			arrName,
											size_t						N_dimensions,
											size_t						N_valuesPerItem,
											size_t						lo,
											size_t						high)
{
	std::ofstream outFile;
	outFile.open(fileName, std::ios::out | std::ios::trunc);
	PrintWebGLOutputRecursive(outFile, arr, arrName, N_dimensions, N_valuesPerItem, lo, high);
}

//===========================================================================================================================
void BlueNoiseExportUtil::SaveAsPPM(const std::vector<float>&			arr,
									 const std::string&					fileName,
									 const size_t						dimensionSize[max_N_dimensions],
									 size_t								N_dimensions,
									 size_t								N_valuesPerItem,
									 uint32_t							slice)
{
	std::ofstream outfile(fileName);
	assert(N_dimensions >= 2);
	size_t srcOffset = 0;
	if (N_dimensions == 3)
	{
		srcOffset += slice * (dimensionSize[0] * dimensionSize[1]);
	}
	N_dimensions = 2; // ignore depth
	outfile << "P3" << std::endl << dimensionSize[0] << " " << dimensionSize[1] << std::endl << 255 << std::endl;
	const uint32_t pixCount = uint32_t(dimensionSize[0] * dimensionSize[1]);
	for (size_t i = srcOffset; i < pixCount + srcOffset; ++i)
	{
		switch (N_valuesPerItem)
		{
		case 1: // monichromatic values
			outfile << FloatAsByteUNorm(arr[i]) << " "
				<< FloatAsByteUNorm(arr[i]) << " "
				<< FloatAsByteUNorm(arr[i]) << " ";
			break;
		case 2:
			outfile << FloatAsByteUNorm(arr[i * 2]) << " "
				<< FloatAsByteUNorm(arr[i * 2 + 1]) << " "
				<< 0 << " ";
			break;
		case 3:
			outfile << FloatAsByteUNorm(arr[i * 3]) << " "
				<< FloatAsByteUNorm(arr[i * 3 + 1]) << " "
				<< FloatAsByteUNorm(arr[i * 3 + 2]) << " ";
			break;
		default:
			assert(0);
		}
	}
	outfile.close();
}

//===========================================================================================================================
void BlueNoiseExportUtil::SaveAsBMP(const std::vector<float>&			arr,
									 const std::string&					fileName,
									 const size_t						dimensionSize[2],
									 size_t								N_valuesPerItem)
{
	uint8_t* bytedata = FloatDataToBytes(arr, N_valuesPerItem, false);
	stbi_write_bmp(fileName.c_str(), int(dimensionSize[0]), int(dimensionSize[1]), 3, bytedata);
	delete[] bytedata;
}
