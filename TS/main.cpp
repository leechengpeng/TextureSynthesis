#include <string>
#include <iostream>
#include <opencv2\core.hpp>
#include <opencv2\imgproc.hpp>
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

	unsigned TextureSize = 128;
	cv::Mat SynthesisTexture(TextureSize, TextureSize, Sample.type());

	clock_t Start = clock();

	TS::TextureSynthesis TS;
	TS.synthesize(Sample, SynthesisTexture);

	clock_t End = clock();
	std::cout << "Time-consuming is: " << (End - Start) / 1000.0 << "s" << std::endl;

	// cv::imshow("SynthesisTexture", SynthesisTexture);
	// TODO: add time to file name
	cv::imwrite("Result/SynthesisTexture.jpg", SynthesisTexture);
	cvWaitKey();

	return 0;
}