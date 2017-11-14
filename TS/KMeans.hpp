#pragma once
#include <ctime>
#include <vector>
#include "Math.hpp"

namespace LT
{
	// TODO: Improving centers select algorithm
	enum KMeansFlags
	{
		KMEANS_RANDOM_CENTERS,
		KMEANS_UNIFORM_CENTERS,
		KMEANS_PP_CENTERS,
	};

	template <typename T>
	class KMeans
	{
	public:
		KMeans();
		~KMeans();

		void cluster(const std::vector<T*>& vDataSet, unsigned vDimension, unsigned vK, std::vector<int>& voLabels, KMeansFlags vFlag = KMEANS_UNIFORM_CENTERS);

		unsigned getK() const;
		const T* getMean(unsigned vIndex) const;

		void setMaxIters(unsigned vMaxIters);

	private:
		void __init(const std::vector<T*>& vDataSet, unsigned vDimension, unsigned vK, KMeansFlags vFlag);
		void __initCenters(const std::vector<T*>& vDataSet, unsigned vDimension, unsigned vK, KMeansFlags vFlag);

		void __cluster(const std::vector<T*>& vDataSet, unsigned vDimension, unsigned vK, std::vector<int>& voLabels);
		void __adjustCenters(const std::vector<T*>& vDataSet, unsigned vDimension, unsigned vK, std::vector<int>& voLabels);

		unsigned m_K;
		unsigned m_MaxIters;
		std::vector<T*> m_Means;
	};

	template <typename T>
	KMeans<T>::KMeans() : m_K(), m_MaxIters(10), m_Means()
	{
	}

	template <typename T>
	KMeans<T>::~KMeans()
	{
		for (size_t i = 0; i < m_Means.size(); ++i)
		{
			if (m_Means[i]) delete[] m_Means[i];
		}
	}

	//******************************************************************************
	//FUNCTION:
	template <typename T>
	void KMeans<T>::cluster(const std::vector<T*>& vDataSet, unsigned vDimension, unsigned vK, std::vector<int>& voLabels, KMeansFlags vFlag)
	{
		_ASSERT(!vDataSet.empty() && vDimension > 0 && vK > 0 && vDataSet.size() == voLabels.size() && vK < vDataSet.size());

		__init(vDataSet, vDimension, vK, vFlag);

		for (size_t Iter = 0; Iter < m_MaxIters; ++Iter)
		{
			__cluster(vDataSet, vDimension, vK, voLabels);
			__adjustCenters(vDataSet, vDimension, vK, voLabels);
		}
	}

	//******************************************************************************
	//FUNCTION:
	template <typename T>
	inline unsigned KMeans<T>::getK() const
	{
		return m_K;
	}

	//******************************************************************************
	//FUNCTION:
	template <typename T>
	inline const T* KMeans<T>::getMean(unsigned vIndex) const
	{
		_ASSERT(vIndex > 0 && vIndex < m_K);
		return m_Means[vIndex];
	}

	//******************************************************************************
	//FUNCTION:
	template <typename T>
	inline void KMeans<T>::setMaxIters(unsigned vMaxIters)
	{
		_ASSERT(vMaxIters > 0);
		m_MaxIters = vMaxIters;
	}

	//******************************************************************************
	//FUNCTION:
	template <typename T>
	void KMeans<T>::__init(const std::vector<T*>& vDataSet, unsigned vDimension, unsigned vK, KMeansFlags vFlag)
	{
		KMeans<T>::~KMeans();
		m_Means.resize(vK);

		for (size_t i = 0; i < vK; ++i)
		{
			m_Means[i] = new T[vDimension]();
		}

		__initCenters(vDataSet, vDimension, vK, vFlag);
	}

	//******************************************************************************
	//FUNCTION:
	template <typename T>
	void KMeans<T>::__initCenters(const std::vector<T*>& vDataSet, unsigned vDimension, unsigned vK, KMeansFlags vFlag)
	{
		if (vFlag == KMEANS_UNIFORM_CENTERS)
		{
			for (size_t k = 0; k < vK; ++k)
			{
				int Select = k * vDataSet.size() / vK;

				for (size_t i = 0; i < vDimension; ++i)
				{
					m_Means[k][i] = vDataSet[Select][i];
				}
			}
		}
		else if (vFlag == KMEANS_RANDOM_CENTERS)
		{
			// TODO
		}
		else if (vFlag == KMEANS_PP_CENTERS)
		{
			// TODO
		}
	}

	//******************************************************************************
	//FUNCTION:
	template <typename T>
	void KMeans<T>::__cluster(const std::vector<T*>& vDataSet, unsigned vDimension, unsigned vK, std::vector<int>& voLabels)
	{
		for (size_t i = 0; i < vDataSet.size(); ++i)
		{
			size_t MinIndex = 0;
			double MinDistance = Math::distance(vDataSet[i], m_Means[MinIndex], vDimension);

			for (size_t k = 1; k < vK; ++k)
			{
				double Distance = Math::distance(vDataSet[i], m_Means[k], vDimension);
				if (Distance < MinDistance)
				{
					MinDistance = Distance;
					MinIndex = k;
				}
			}

			voLabels[i] = MinIndex;
		}
	}

	//******************************************************************************
	//FUNCTION:
	template <typename T>
	void KMeans<T>::__adjustCenters(const std::vector<T*>& vDataSet, unsigned vDimension, unsigned vK, std::vector<int>& voLabels)
	{
		std::vector<std::vector<double>> Centers(vK, std::vector<double>(vDimension));
		std::vector<int> Counts(vK, 0);
		for (size_t k = 0; k < vDataSet.size(); ++k)
		{
			for (size_t Index = 0; Index < vDimension; ++Index)
			{
				Centers[voLabels[k]][Index] += vDataSet[k][Index];
			}
			++Counts[voLabels[k]];
		}

		for (size_t k = 0; k < vK; ++k)
		{
			for (size_t Index = 0; Index < vDimension; ++Index)
			{
				_ASSERT(Counts[k] > 0);
				m_Means[k][Index] = static_cast<T>(Centers[k][Index] / Counts[k]);
			}
		}
	}
}