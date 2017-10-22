#include <iostream>
#include <opencv2\core.hpp>
#include <opencv2\highgui\highgui.hpp>
#include "TextureSynthesis.hpp"

int main()
{
	cv::Mat Sample = cv::imread("Sample.jpg");
	if (Sample.empty())
	{
		std::cerr << "Loading Sample failed..." << std::endl;
		return -1;
	}

	cv::Mat SynthesisTexture(64, 64, Sample.type());
	TS::TextureSynthesis TS;
	TS.synthesize(Sample, SynthesisTexture);

	cv::imshow("SynthesisTexture", SynthesisTexture);
	cvWaitKey();

	return 0;
}