/*
 * Tests.cpp
 *
 *  Created on: 29-12-2011
 *      Author: Michal
 */

#define BOOST_TEST_MODULE WaveletTransformCompressionTest
#include <boost/test/included/unit_test.hpp>

#include "HaarWaveletTransform.h"

const double DOUBLE_CLOSE = 0.0001;

BOOST_AUTO_TEST_SUITE( HaarTransformTests )

BOOST_AUTO_TEST_CASE( ShiftMultiplicationTest ) {
	unsigned int i = 1;
	i <<= 1;
	BOOST_CHECK_EQUAL(i, 2u);
	i <<= 1;
	BOOST_CHECK_EQUAL(i, 4u);
	i <<= 1;
	BOOST_CHECK_EQUAL(i, 8u);
	i <<= 1;
	BOOST_CHECK_EQUAL(i, 16u);
}

BOOST_AUTO_TEST_CASE( SimpleHaarTransformTest ) {
	double testData[16];
	double inputData = 0.0;
	for(int i = 0; i < 16; ++i) {
		inputData += 1.0;
		testData[i] = inputData;
	}

	HaarWaveletTransform haarTransform(testData, 4, 4);
	haarTransform.transform();

	BOOST_CHECK_CLOSE(haarTransform.getItem(0, 0), 17.0 / 2.0, DOUBLE_CLOSE);
	BOOST_CHECK_CLOSE(haarTransform.getItem(0, 1), -1.0, DOUBLE_CLOSE);
	BOOST_CHECK_CLOSE(haarTransform.getItem(0, 2), -1.0 / 2.0, DOUBLE_CLOSE);
	BOOST_CHECK_CLOSE(haarTransform.getItem(0, 3), -1.0 / 2.0, DOUBLE_CLOSE);
	BOOST_CHECK_CLOSE(haarTransform.getItem(1, 0), -4.0, DOUBLE_CLOSE);
	BOOST_CHECK_CLOSE(haarTransform.getItem(1, 1), 0.0, DOUBLE_CLOSE);
	BOOST_CHECK_CLOSE(haarTransform.getItem(1, 2), 0.0, DOUBLE_CLOSE);
	BOOST_CHECK_CLOSE(haarTransform.getItem(1, 3), 0.0, DOUBLE_CLOSE);
	BOOST_CHECK_CLOSE(haarTransform.getItem(2, 0), -2.0, DOUBLE_CLOSE);
	BOOST_CHECK_CLOSE(haarTransform.getItem(2, 1), 0.0, DOUBLE_CLOSE);
	BOOST_CHECK_CLOSE(haarTransform.getItem(2, 2), 0.0, DOUBLE_CLOSE);
	BOOST_CHECK_CLOSE(haarTransform.getItem(2, 3), 0.0, DOUBLE_CLOSE);
	BOOST_CHECK_CLOSE(haarTransform.getItem(3, 0), -2.0, DOUBLE_CLOSE);
	BOOST_CHECK_CLOSE(haarTransform.getItem(3, 1), 0.0, DOUBLE_CLOSE);
	BOOST_CHECK_CLOSE(haarTransform.getItem(3, 2), 0.0, DOUBLE_CLOSE);
	BOOST_CHECK_CLOSE(haarTransform.getItem(3, 3), 0.0, DOUBLE_CLOSE);
}

BOOST_AUTO_TEST_CASE( SimpleHaarInverseTransformTest ) {
	unsigned int testWidth = 64;
	unsigned int testHeight = 64;
	double* testData = new double[testWidth * testHeight];
	double inputData = 0.0;
	for(unsigned int i = 0; i < testWidth * testHeight; ++i) {
		inputData += 1.0;
		testData[i] = inputData;
	}

	HaarWaveletTransform haarTransform(testData, testWidth, testHeight);
	haarTransform.transform();
	haarTransform.inverseTransform();

	for(unsigned int row = 0; row < testHeight; ++row) {
		for(unsigned int column = 0; column < testWidth; ++column) {
			BOOST_CHECK_CLOSE(testData[(row * testWidth) + column], haarTransform.getItem(row, column), DOUBLE_CLOSE);
		}
	}
	delete[] testData;
}

BOOST_AUTO_TEST_CASE( SimpleHaarInverseTransformDifferentWidthHeightTest ) {
	unsigned int testWidth = 64;
	unsigned int testHeight = 128;
	double* testData = new double[testWidth * testHeight];
	double inputData = 0.0;
	for(unsigned int i = 0; i < testWidth * testHeight; ++i) {
		inputData += 1.0;
		testData[i] = inputData;
	}

	HaarWaveletTransform haarTransform(testData, testWidth, testHeight);
	haarTransform.transform();
	haarTransform.inverseTransform();

	for(unsigned int row = 0; row < testHeight; ++row) {
		for(unsigned int column = 0; column < testWidth; ++column) {
			BOOST_CHECK_CLOSE(testData[(row * testWidth) + column], haarTransform.getItem(row, column), DOUBLE_CLOSE);
		}
	}
}

BOOST_AUTO_TEST_SUITE_END()
