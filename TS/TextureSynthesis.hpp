#pragma once
#include <ctime>
#include <vector>
#include <opencv2\core.hpp>
#include "TSVQ.hpp"

namespace TS
{
	class TextureSynthesis
	{
	public:
		TextureSynthesis();
		~TextureSynthesis();

		void synthesize(const cv::Mat& vSampleTexture, cv::Mat& vTexture, bool vAccelerate = false);

	private:
		void __genNoiseTexture(cv::Mat& vTexture);
		void __calSampleNeighborhoods(const cv::Mat& vSample, std::vector<uchar*>& voSampleNeighborhoods);
		uchar* __calNeighborhoods(const cv::Mat& vImage, const CvPoint& vPoint);

		double __calEuclideanDistance(const uchar *vVectorsA, const uchar *vVectorsB, unsigned int vDim);
		const uchar* __getMostSimilarPixel(const std::vector<uchar*>& vSampleNeighborhoods, uchar* vNeighborhoods, unsigned int vDim);

		int m_Level;
		int m_NeighborhoodSize;
		int m_TemplateSize;
		int m_HalfTemplateSize;
	};

	TextureSynthesis::TextureSynthesis() : m_TemplateSize(9), m_Level(1), m_NeighborhoodSize(0), m_HalfTemplateSize(0)
	{

	}

	TextureSynthesis::~TextureSynthesis()
	{
		m_Level = 1;
		m_TemplateSize = 9;
	}

	// *********************************************************
	// Function:
	void TextureSynthesis::synthesize(const cv::Mat& vSampleTexture, cv::Mat& vTexture, bool vAccelerate)
	{
		_ASSERT(!vTexture.empty() && m_TemplateSize > 2 && m_TemplateSize % 2 != 0);

		m_HalfTemplateSize = (m_TemplateSize - 1) / 2;
		m_NeighborhoodSize = (m_HalfTemplateSize * m_TemplateSize + m_HalfTemplateSize + 1) * vSampleTexture.channels();

		__genNoiseTexture(vTexture);

		std::vector<uchar*> SampleNeighborhoods;
		__calSampleNeighborhoods(vSampleTexture, SampleNeighborhoods);

		LLL::TSVQ<uchar> Accelerator;
		if (vAccelerate)
		{
			Accelerator.build(SampleNeighborhoods, m_NeighborhoodSize, 100);	//FIXME: magic number
		}

		for (int RowIndex = 0, Channels = vTexture.channels(); RowIndex < vTexture.rows; ++RowIndex)
		{
			auto pRows = vTexture.ptr<uchar>(RowIndex);
			for (int ColIndex = 0; ColIndex < vTexture.cols; ++ColIndex)
			{
				auto pNeighborhood = __calNeighborhoods(vTexture, CvPoint(ColIndex, RowIndex));

				const uchar* MostSimilarPixel = nullptr;
				if (vAccelerate)
				{
					MostSimilarPixel = Accelerator.quantizeVector(pNeighborhood);
				}
				else
				{
					MostSimilarPixel = __getMostSimilarPixel(SampleNeighborhoods, pNeighborhood, vSampleTexture.channels());
				}
				MostSimilarPixel += m_NeighborhoodSize - vSampleTexture.channels();

				_ASSERT(MostSimilarPixel);
				memcpy(pRows + ColIndex * vSampleTexture.channels(), MostSimilarPixel, vSampleTexture.channels());
			}
		}
	}

	// *********************************************************
	// Function:
	const uchar* TextureSynthesis::__getMostSimilarPixel(const std::vector<uchar*>& vSampleNeighborhoods, uchar* vNeighborhoods, unsigned int vDim)
	{
		_ASSERT(!vSampleNeighborhoods.empty());

		auto Iter = vSampleNeighborhoods.begin();

		uchar* MostSimilarPixel = *Iter;
		double MinLength = __calEuclideanDistance(*Iter, vNeighborhoods, m_NeighborhoodSize);

		while (Iter != vSampleNeighborhoods.end())
		{
			double Length = __calEuclideanDistance(*Iter, vNeighborhoods, m_NeighborhoodSize);
			if (Length < MinLength)
			{
				MinLength = Length;
				MostSimilarPixel = *Iter;
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
			double Difference = vVectorsA[Index] - vVectorsB[Index];
			Sum += Difference * Difference;
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
	void TextureSynthesis::__calSampleNeighborhoods(const cv::Mat& vSample, std::vector<uchar*>& voSampleNeighborhoods)
	{
		for (int RowIndex = 0; RowIndex < vSample.rows; ++RowIndex)
		{
			for (int ColIndex = 0; ColIndex < vSample.cols; ++ColIndex)
			{
				voSampleNeighborhoods.push_back(__calNeighborhoods(vSample, CvPoint(ColIndex, RowIndex)));
			}
		}
	}

	// *********************************************************
	// Function:
	uchar* TextureSynthesis::__calNeighborhoods(const cv::Mat& vImage, const CvPoint& vPoint)
	{
		int Index = 0;
		uchar* pNeighborhood = new uchar[m_NeighborhoodSize]();

		for (int RowIndex = vPoint.y - m_HalfTemplateSize; RowIndex < vPoint.y; ++RowIndex)
		{
			for (int ColIndex = vPoint.x - m_HalfTemplateSize; ColIndex <= vPoint.x + m_HalfTemplateSize; ++ColIndex)
			{
				CvPoint Point;
				RowIndex < 0 ? Point.y = vImage.rows + RowIndex : Point.y = RowIndex;
				ColIndex < 0 ? Point.x = vImage.cols + ColIndex : ColIndex >= vImage.cols ? Point.x = ColIndex - vImage.cols : Point.x = ColIndex;

				memcpy(pNeighborhood + Index, vImage.ptr<uchar>(Point.y) + Point.x * vImage.channels(), vImage.channels());
				Index += vImage.channels();

				_ASSERT(Index <= m_NeighborhoodSize);
			}
		}

		auto pRows = vImage.ptr<uchar>(vPoint.y);
		for (int ColIndex = vPoint.x - m_HalfTemplateSize; ColIndex <= vPoint.x; ++ColIndex)
		{
			int x = (ColIndex < 0 ? vImage.cols + ColIndex : ColIndex);

			memcpy(pNeighborhood + Index, pRows + x * vImage.channels(), vImage.channels());
			Index += vImage.channels();

			_ASSERT(Index <= m_NeighborhoodSize);
		}

		return pNeighborhood;
	}
}