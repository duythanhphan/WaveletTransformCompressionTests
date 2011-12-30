/*
 * Tests.cpp
 *
 *  Created on: 29-12-2011
 *      Author: Michal
 */

#define BOOST_TEST_MODULE WaveletTransformCompressionTest
#include <boost/test/included/unit_test.hpp>

#include "HaarWaveletTransform.h"
#include "UnsignedInteger.h"

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

BOOST_AUTO_TEST_SUITE( UnsignedIntegerTests )

BOOST_AUTO_TEST_CASE( isPowerOfTwoTest ) {
	unsigned int powersOfTwo[32];
	powersOfTwo[0] = 1;
	unsigned int nextPow = 1;
	for(unsigned int i = 0; i < 31; ++i) {
		nextPow <<= 1;
		powersOfTwo[i + 1] = nextPow;
	}

	for(unsigned int i = 0; i < 32; ++i) {
		BOOST_CHECK(UnsignedInteger::isPowerOfTwo(powersOfTwo[i]) == true);
	}

	unsigned int notPowersOfTwo[10] = {0, 3, 130, 6, 97, 2000, 10000, 700, 100000, 25};
	for(unsigned int i = 0; i < 10; ++i) {
		BOOST_CHECK(UnsignedInteger::isPowerOfTwo(notPowersOfTwo[i]) == false);
	}
}

BOOST_AUTO_TEST_CASE( getClosestPowerOfTwoTest ) {
	BOOST_CHECK_EQUAL(UnsignedInteger::getClosestPowerOfTwo(0), 1u);
	BOOST_CHECK_EQUAL(UnsignedInteger::getClosestPowerOfTwo(3), 4u);
	BOOST_CHECK_EQUAL(UnsignedInteger::getClosestPowerOfTwo(5), 8u);
	BOOST_CHECK_EQUAL(UnsignedInteger::getClosestPowerOfTwo(7), 8u);
	BOOST_CHECK_EQUAL(UnsignedInteger::getClosestPowerOfTwo(120), 128u);
	BOOST_CHECK_EQUAL(UnsignedInteger::getClosestPowerOfTwo(16300), 16384u);
	BOOST_CHECK_EQUAL(UnsignedInteger::getClosestPowerOfTwo(16385), 32768u);

	BOOST_CHECK_EQUAL(UnsignedInteger::getClosestPowerOfTwo(2), 2u);
	BOOST_CHECK_EQUAL(UnsignedInteger::getClosestPowerOfTwo(1024), 1024u);
	BOOST_CHECK_EQUAL(UnsignedInteger::getClosestPowerOfTwo(256), 256u);
	BOOST_CHECK_EQUAL(UnsignedInteger::getClosestPowerOfTwo(4), 4u);
}

BOOST_AUTO_TEST_SUITE_END()
