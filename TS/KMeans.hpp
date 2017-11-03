#pragma once
#include <vector>

namespace LT
{
	enum KMeansFlags
	{
		KMEANS_RANDOM_CENTERS,
		KMEANS_UNIFORM,
		// TODO: Use kmeans++ center initialization by Arthur and Vassilvitskii[Arthur2007].
		KMEANS_PP_CENTERS,
	};

	template <typename T>
	class KMeans
	{
	public:
		KMeans() : m_MaxIters(10), m_Means() {}

		void cluster(const std::vector<T*>& vDataSet, unsigned vDimension, unsigned vK, std::vector<int>& voLabels, KMeansFlags vFlag = KMEANS_RANDOM_CENTERS);

		const std::vector<T*> getKMeans() const { return m_Means; }

	private:
		void __init(const std::vector<T*>& vDataSet, unsigned vDimension, unsigned vK, KMeansFlags vFlag = KMEANS_RANDOM_CENTERS);
		void __cluster(const std::vector<T*>& vDataSet, unsigned vDimension, unsigned vK, std::vector<int>& voLabels);

		double __distance(const T *vVectorA, const T *vVectorB, unsigned vDim) const;

		unsigned m_MaxIters;
		std::vector<T*> m_Means;
	};

	template <typename T>
	void KMeans<T>::cluster(const std::vector<T*>& vDataSet, unsigned vDimension, unsigned vK, std::vector<int>& voLabels, KMeansFlags vFlag)
	{
		_ASSERT(!vDataSet.empty() && vDimension > 0 && vK > 0 && vDataSet.size() == voLabels.size() && vK < vDataSet.size());

		__init(vDataSet, vDimension, vK, vFlag);

		for (size_t i = 0; i < m_MaxIters; ++i)
		{
			// cluster
			for (size_t k = 0; k < vDataSet.size(); ++k)
			{
				unsigned Index = 0;

				auto MinIndex = Index;
				auto MinDistance = __distance(vDataSet[k], m_Means[Index], vDimension);

				for (++Index; Index < vK; ++Index)
				{
					auto Distance = __distance(vDataSet[k], m_Means[Index], vDimension);

					if (Distance < MinDistance)
					{
						MinIndex = Index;
						MinDistance = Distance;
					}
				}
				voLabels[k] = MinIndex;
			}

			// adjust center
			std::vector<std::vector<double>> Center(vK, std::vector<double>(vDimension));
			std::vector<int> Counts(vK, 0);
			for (size_t k = 0; k < vDataSet.size(); ++k)
			{
				for (size_t Index = 0; Index < vDimension; ++Index)
				{
					Center[voLabels[k]][Index] += vDataSet[k][Index];
				}
				++Counts[voLabels[k]];
			}

			for (size_t k = 0; k < vK; ++k)
			{
				for (size_t Index = 0; Index < vDimension; ++Index)
				{
					_ASSERT(Counts[k] > 0);
					m_Means[k][Index] = static_cast<T>(Center[k][Index] / Counts[k]);
				}
			}
		}

		int i = 0;
	}

	template <typename T>
	void KMeans<T>::__cluster(const std::vector<T*>& vDataSet, unsigned vDimension, unsigned vK, std::vector<int>& voLabels)
	{

	}

	//******************************************************************************
	//FUNCTION:
	template <typename T>
	void KMeans<T>::__init(const std::vector<T*>& vDataSet, unsigned vDimension, unsigned vK, KMeansFlags vFlag)
	{
		if (!m_Means.empty())
		{
			for (auto Data : m_Means) { delete[] Data; }
			m_Means.clear();
		}
		m_Means.resize(vK);
		for (size_t i = 0; i < vK; ++i)
		{
			m_Means[i] = new T[vDimension]();
		}

		if (vFlag == KMEANS_RANDOM_CENTERS)
		{
			for (size_t i = 0; i < vK; ++i)
			{
				for (size_t k = 0; k < vDimension; ++k)
				{
					m_Means[i][k] = vDataSet[i][k];
				}
			}
		}
		else if (vFlag == KMEANS_UNIFORM)
		{
			for (size_t k = 0; k < vK; ++k)
			{
				int Select = k * vDataSet.size() / vK;

				for (size_t i = 0; i < vDimension; ++i)
				{
					m_Means[k][i] = vDataSet[Select][i];
				}

				int i = 0;
			}
		}
		else if (vFlag == KMEANS_PP_CENTERS)
		{
			// TODO
		}
	}

	//******************************************************************************
	//FUNCTION:
	template <typename T>
	double KMeans<T>::__distance(const T *vVectorA, const T *vVectorB, unsigned vDim) const
	{
		_ASSERT(vVectorA && vVectorB);

		double Sum = 0.0;
		for (unsigned Index = 0; Index < vDim; ++Index)
		{
			Sum += (vVectorA[Index] - vVectorB[Index]) * (vVectorA[Index] - vVectorB[Index]);
		}

		return sqrt(Sum);
	}
}