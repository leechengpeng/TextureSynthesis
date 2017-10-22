#pragma once
#include <ctime>
#include <vector>
#include <opencv2\core.hpp>

namespace TS
{
	class TextureSynthesis
	{
	public:
		TextureSynthesis();
		~TextureSynthesis();

		void synthesize(const cv::Mat& vSampleTexture, cv::Mat& vTexture);

	private:
		void __genNoiseTexture(cv::Mat& vTexture);
		void __calSampleNeighborhoods(const cv::Mat& vSample, std::vector<std::vector<const uchar*>>& voSampleNeighborhoods);
		void __calNeighborhoods(const cv::Mat& vImage, const CvPoint& vPoint, int NeighborhoodSize, std::vector<const uchar*>& voNeighborhoods);

		double __calEuclideanDistance(const uchar *vVectorsA, const uchar *vVectorsB, unsigned int vDim);
		const uchar* __getMostSimilarPixel(const std::vector<std::vector<const uchar*>>& vSampleNeighborhoods, const std::vector<const uchar*>& vNeighborhoods, unsigned int vDim);

		int m_NeighborhoodSize;
		int m_ResolutionLevel;
	};

	TextureSynthesis::TextureSynthesis() : m_NeighborhoodSize(9), m_ResolutionLevel(1)
	{

	}

	TextureSynthesis::~TextureSynthesis()
	{
		m_NeighborhoodSize = 9;
		m_ResolutionLevel  = 1;
	}

	// *********************************************************
	// Function:
	void TextureSynthesis::synthesize(const cv::Mat& vSampleTexture, cv::Mat& vTexture)
	{
		_ASSERT(!vTexture.empty());

		__genNoiseTexture(vTexture);

		std::vector<std::vector<const uchar*>> SampleNeighborhoods;
		__calSampleNeighborhoods(vSampleTexture, SampleNeighborhoods);

		for (int RowIndex = 0, Channels = vTexture.channels(); RowIndex < vTexture.rows; ++RowIndex)
		{
			auto pRows = vTexture.ptr<uchar>(RowIndex);
			for (int ColIndex = 0; ColIndex < vTexture.cols; ++ColIndex)
			{
				std::vector<const uchar*> Neighborhoods;
				__calNeighborhoods(vTexture, CvPoint(ColIndex, RowIndex), m_NeighborhoodSize, Neighborhoods);

				auto MostSimilarPixel = __getMostSimilarPixel(SampleNeighborhoods, Neighborhoods, vSampleTexture.channels());

				_ASSERT(MostSimilarPixel);
				memcpy(pRows + ColIndex * vSampleTexture.channels(), MostSimilarPixel, vSampleTexture.channels());
			}
		}
	}

	// *********************************************************
	// Function:
	const uchar* TextureSynthesis::__getMostSimilarPixel(const std::vector<std::vector<const uchar*>>& vSampleNeighborhoods, const std::vector<const uchar*>& vNeighborhoods, unsigned int vDim)
	{
		double MinLength = 0.0;
		const uchar* MostSimilarPixel = nullptr;

		auto Iter = vSampleNeighborhoods.begin();
		for (size_t i = 0; i < Iter->size(); ++i)
		{
			MinLength += __calEuclideanDistance((*Iter)[i], vNeighborhoods[i], vDim);
		}

		while (Iter != vSampleNeighborhoods.end())
		{
			double Length = 0.0;
			for (size_t i = 0; i < Iter->size(); ++i)
			{
				Length += __calEuclideanDistance((*Iter)[i], vNeighborhoods[i], vDim);
			}

			if (Length < MinLength)
			{
				MinLength = Length;
				MostSimilarPixel = Iter->back();
			}

			++Iter;
		}

		return MostSimilarPixel;
	}

	// *********************************************************
	// Function:
	double TextureSynthesis::__calEuclideanDistance(const uchar *vVectorsA, const uchar *vVectorsB, unsigned int vDim)
	{
		_ASSERT(vVectorsA && vVectorsB);

		double Sum = 0.0;
		for (unsigned int Index = 0; Index < vDim; ++Index)
		{
			Sum += (vVectorsA[Index] - vVectorsB[Index]) * (vVectorsA[Index] - vVectorsB[Index]);
		}

		return Sum;
	}

	// *********************************************************
	// Function:
	void TextureSynthesis::__genNoiseTexture(cv::Mat& vTexture)
	{
		srand((unsigned)time(NULL));

		for (int RowIndex = 0, Channels = vTexture.channels(); RowIndex < vTexture.rows; ++RowIndex)
		{
			auto pRows = vTexture.ptr<uchar>(RowIndex);
			for (int ColIndex = 0; ColIndex < vTexture.cols; ++ColIndex)
			{
				for (int ChannelIndex = 0; ChannelIndex < Channels; ++ChannelIndex)
				{
					pRows[ColIndex * Channels + ChannelIndex] = rand() % 255;
				}
			}
		}
	}

	// *********************************************************
	// Function:
	void TextureSynthesis::__calSampleNeighborhoods(const cv::Mat& vSample, std::vector<std::vector<const uchar*>>& voSampleNeighborhoods)
	{
		for (int RowIndex = 0; RowIndex < vSample.rows; ++RowIndex)
		{
			for (int ColIndex = 0; ColIndex < vSample.cols; ++ColIndex)
			{
				std::vector<const uchar*> Neighborhoods;
				__calNeighborhoods(vSample, CvPoint(ColIndex, RowIndex), m_NeighborhoodSize, Neighborhoods);

				voSampleNeighborhoods.push_back(Neighborhoods);
			}
		}
	}

	// *********************************************************
	// Function:
	void TextureSynthesis::__calNeighborhoods(const cv::Mat& vImage, const CvPoint& vPoint, int NeighborhoodSize, std::vector<const uchar*>& voNeighborhoods)
	{
		_ASSERT(NeighborhoodSize > 2 && NeighborhoodSize % 2 != 0);
	
		int NeighborhoodCols     = NeighborhoodSize;
		int NeighborhoodRows     = (NeighborhoodSize + 1) / 2;
		int HalfNeighborhoodSize = (NeighborhoodSize - 1) / 2;

		for (int RowIndex = vPoint.y - HalfNeighborhoodSize; RowIndex < vPoint.y; ++RowIndex)
		{
			for (int ColIndex = vPoint.x - HalfNeighborhoodSize; ColIndex <= vPoint.x + HalfNeighborhoodSize; ++ColIndex)
			{
				CvPoint Point;
				RowIndex < 0 ? Point.y = vImage.rows + RowIndex : Point.y = RowIndex;
				ColIndex < 0 ? Point.x = vImage.cols + ColIndex : ColIndex >= vImage.cols ? Point.x = ColIndex - vImage.cols : Point.x = ColIndex;

				voNeighborhoods.push_back(vImage.ptr<uchar>(Point.y) + Point.x * vImage.channels());
			}
		}

		auto pRows = vImage.ptr<uchar>(vPoint.y);
		for (int ColIndex = vPoint.x - HalfNeighborhoodSize; ColIndex <= vPoint.x; ++ColIndex)
		{
			int x = (ColIndex < 0 ? vImage.cols + ColIndex : ColIndex);
			voNeighborhoods.push_back(pRows + x * vImage.channels());
		}
	}
}