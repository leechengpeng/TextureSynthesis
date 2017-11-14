#pragma once
#include <cassert>

namespace LT
{
	namespace Math
	{
		constexpr double ERROR = 0.0000001;

		//******************************************************************************
		//FUNCTION:
		template <typename T>
		double distance(const T* vVector1, const T* vVertor2, const unsigned vDimension)
		{
			assert(vVector1 && vVertor2);

			double Distance = 0.0;
			for (size_t i = 0; i < vDimension; ++i)
			{
				Distance += (vVector1[i] - vVertor2[i]) * (vVector1[i] - vVertor2[i]);
			}

			return sqrt(Distance);
		}
	}
}