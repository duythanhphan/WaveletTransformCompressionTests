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
#include "Quantizer.h"
#include "Heap.h"
#include "HuffmanCoding.h"
#include "RLE.h"
#include "Encoder.h"
#include "HuffmanDecoder.h"
#include "RLEDecoder.h"

const double DOUBLE_CLOSE = 0.0001;

unsigned int getMask(unsigned int size) {
	unsigned int mask = 0;
	for(unsigned int i = 0; i < size; ++i) {
		UnsignedInteger::setBitFromLeft(&mask, i);
	}
	return mask;
}

bool isPrefix(unsigned int prefix, unsigned int code, unsigned int prefixSize) {
	unsigned int mask = getMask(prefixSize);
	unsigned int test = code & mask;
	if(prefix == test) {
		return true;
	}

	return false;
}

BOOST_AUTO_TEST_SUITE( SimpleTests )

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

BOOST_AUTO_TEST_CASE( ShiftTest ) {
	unsigned int toShift = 0;
	toShift |= 1;
	toShift |= 1 << 1;
	toShift |= 1 << 4;
	toShift |= 1 << 6;
	BOOST_CHECK_EQUAL(toShift, 83u);
	//toShift
	//6  5  4  3  2  1  0
	//1  0  1  0  0  1  1
	//31 30 29 28 27 26 25

	unsigned int size = sizeof(unsigned int) * 8 - 1;
	unsigned int shift = size - 6;
	toShift <<= shift;

	unsigned int test = 0;
	test |= 1 << 31;
	test |= 1 << 29;
	test |= 1 << 26;
	test |= 1 << 25;

	BOOST_CHECK_EQUAL(test, toShift);

	unsigned int a = 3;
	a <<= 2;
	BOOST_CHECK_EQUAL(a, 12u);
}

BOOST_AUTO_TEST_CASE( CreateOffMaskTest ) {
	std::map<unsigned int, unsigned int> m_offMasks;
	unsigned int m_iMaxBitPosition = (sizeof(unsigned int) * 8) - 1;

	unsigned int mask = 0;
	unsigned int offIndex = m_iMaxBitPosition;
	for(unsigned int i = 0; i < m_iMaxBitPosition; ++i) {
		mask |= 1 << i;
		m_offMasks.insert(std::pair<unsigned int, unsigned int>(offIndex, mask));
		--offIndex;
	}

	unsigned int test = 0;
	test = m_offMasks[31];
	BOOST_CHECK_EQUAL(test, 1u);
	test = m_offMasks[30];
	BOOST_CHECK_EQUAL(test, 3u);
}

BOOST_AUTO_TEST_CASE( MaskDecodeTest ) {
	unsigned int testData[2];
	UnsignedInteger::setBitFromRight(&testData[0], 0);
	UnsignedInteger::setBitFromRight(&testData[0], 2);
	UnsignedInteger::setBitFromLeft(&testData[1], 0);
	UnsignedInteger::setBitFromLeft(&testData[1], 1);
	UnsignedInteger::setBitFromLeft(&testData[1], 3);

	unsigned int code = 0;
	unsigned int codeSize = 7;
	UnsignedInteger::setBitFromLeft(&code, 0);
	UnsignedInteger::setBitFromLeft(&code, 2);
	UnsignedInteger::setBitFromLeft(&code, 3);
	UnsignedInteger::setBitFromLeft(&code, 4);
	UnsignedInteger::setBitFromLeft(&code, 6);

	unsigned int pos = 29;
	unsigned int decoded = 0;

	unsigned int lmaskIndex = 32 - pos;
	unsigned int rmaskIndex = codeSize - lmaskIndex;
	BOOST_CHECK_EQUAL(lmaskIndex, 3u);
	BOOST_CHECK_EQUAL(rmaskIndex, 4u);

	unsigned int lMask = 0;
	UnsignedInteger::setBitFromRight(&lMask, 0);
	UnsignedInteger::setBitFromRight(&lMask, 1);
	UnsignedInteger::setBitFromRight(&lMask, 2);
	unsigned int rMask = 0;
	UnsignedInteger::setBitFromLeft(&rMask, 0);
	UnsignedInteger::setBitFromLeft(&rMask, 1);
	UnsignedInteger::setBitFromLeft(&rMask, 2);
	UnsignedInteger::setBitFromLeft(&rMask, 3);

	decoded  = testData[0] & lMask;
	BOOST_CHECK_EQUAL(decoded, 5u);
	decoded <<= pos;
	unsigned int testDecoded1 = 0;
	UnsignedInteger::setBitFromLeft(&testDecoded1, 0);
	UnsignedInteger::setBitFromLeft(&testDecoded1, 2);
	BOOST_CHECK_EQUAL(decoded, testDecoded1);

	decoded |= (testData[1] & rMask) >> lmaskIndex;
	BOOST_CHECK_EQUAL(decoded, code);
}

BOOST_AUTO_TEST_CASE( CodeLessOperator ) {
	std::map<HuffmanCoding<double>::Code, double > testMap;
	HuffmanCoding<double>::Code code;
	code.code = 8;
	code.size = 4;
	testMap.insert(std::pair<HuffmanCoding<double>::Code, double >(code, 4.0));
	code.code = 2;
	code.size = 2;
	testMap.insert(std::pair<HuffmanCoding<double>::Code, double >(code, 2.0));
	code.code = 1;
	code.size = 1;
	testMap.insert(std::pair<HuffmanCoding<double>::Code, double >(code, 1.0));
	code.code = 4;
	code.size = 3;
	testMap.insert(std::pair<HuffmanCoding<double>::Code, double >(code, 3.0));

	std::map<HuffmanCoding<double>::Code, double >::iterator it = testMap.begin();
	BOOST_REQUIRE(it != testMap.end());
	BOOST_CHECK_EQUAL(it->first.code, 1u);
	BOOST_CHECK_EQUAL(it->first.size, 1u);
	BOOST_CHECK_EQUAL(it->second, 1.0);

	++it;
	BOOST_REQUIRE(it != testMap.end());
	BOOST_CHECK_EQUAL(it->first.code, 2u);
	BOOST_CHECK_EQUAL(it->first.size, 2u);
	BOOST_CHECK_EQUAL(it->second, 2.0);

	++it;
	BOOST_REQUIRE(it != testMap.end());
	BOOST_CHECK_EQUAL(it->first.code, 4u);
	BOOST_CHECK_EQUAL(it->first.size, 3u);
	BOOST_CHECK_EQUAL(it->second, 3.0);

	++it;
	BOOST_REQUIRE(it != testMap.end());
	BOOST_CHECK_EQUAL(it->first.code, 8u);
	BOOST_CHECK_EQUAL(it->first.size, 4u);
	BOOST_CHECK_EQUAL(it->second, 4.0);
}

BOOST_AUTO_TEST_CASE( ShiftFilledZeroTest ) {
	unsigned int test = 0;
	unsigned int mask = 0;
	for(unsigned int i = 0; i < UnsignedInteger::NUMBER_OF_BITS; ++i) {
		UnsignedInteger::setBitFromRight(&test, i);
		UnsignedInteger::setBitFromRight(&mask, i);
	}
	test <<= 30;

	unsigned int result = 0;
	UnsignedInteger::setBitFromLeft(&result, 0);
	UnsignedInteger::setBitFromLeft(&result, 1);

	BOOST_CHECK_EQUAL(test, result);
}

BOOST_AUTO_TEST_CASE( isPrefixTest ) {
	unsigned int prefix = 0;
	UnsignedInteger::setBitFromLeft(&prefix, 1);
	unsigned int code = 0;
	UnsignedInteger::setBitFromLeft(&code, 1);
	UnsignedInteger::setBitFromLeft(&code, 3);

	BOOST_CHECK(isPrefix(prefix, code, 2u));

	unsigned int otherCode = 0;
	UnsignedInteger::setBitFromLeft(&otherCode, 2);
	UnsignedInteger::setBitFromLeft(&otherCode, 3);
	BOOST_CHECK(isPrefix(otherCode, code, 4u) == false);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( HaarTransformTests )

//After change division in Haar from 2 to sqrt(2) test is no logner valid.
//BOOST_AUTO_TEST_CASE( SimpleHaarTransformTest ) {
//	double testData[16];
//	double inputData = 0.0;
//	for(int i = 0; i < 16; ++i) {
//		inputData += 1.0;
//		testData[i] = inputData;
//	}
//
//	HaarWaveletTransform haarTransform(testData, 4, 4);
//	haarTransform.transform();
//
//	BOOST_CHECK_CLOSE(haarTransform.getItem(0, 0), 17.0 / 2.0, DOUBLE_CLOSE);
//	BOOST_CHECK_CLOSE(haarTransform.getItem(0, 1), -1.0, DOUBLE_CLOSE);
//	BOOST_CHECK_CLOSE(haarTransform.getItem(0, 2), -1.0 / 2.0, DOUBLE_CLOSE);
//	BOOST_CHECK_CLOSE(haarTransform.getItem(0, 3), -1.0 / 2.0, DOUBLE_CLOSE);
//	BOOST_CHECK_CLOSE(haarTransform.getItem(1, 0), -4.0, DOUBLE_CLOSE);
//	BOOST_CHECK_CLOSE(haarTransform.getItem(1, 1), 0.0, DOUBLE_CLOSE);
//	BOOST_CHECK_CLOSE(haarTransform.getItem(1, 2), 0.0, DOUBLE_CLOSE);
//	BOOST_CHECK_CLOSE(haarTransform.getItem(1, 3), 0.0, DOUBLE_CLOSE);
//	BOOST_CHECK_CLOSE(haarTransform.getItem(2, 0), -2.0, DOUBLE_CLOSE);
//	BOOST_CHECK_CLOSE(haarTransform.getItem(2, 1), 0.0, DOUBLE_CLOSE);
//	BOOST_CHECK_CLOSE(haarTransform.getItem(2, 2), 0.0, DOUBLE_CLOSE);
//	BOOST_CHECK_CLOSE(haarTransform.getItem(2, 3), 0.0, DOUBLE_CLOSE);
//	BOOST_CHECK_CLOSE(haarTransform.getItem(3, 0), -2.0, DOUBLE_CLOSE);
//	BOOST_CHECK_CLOSE(haarTransform.getItem(3, 1), 0.0, DOUBLE_CLOSE);
//	BOOST_CHECK_CLOSE(haarTransform.getItem(3, 2), 0.0, DOUBLE_CLOSE);
//	BOOST_CHECK_CLOSE(haarTransform.getItem(3, 3), 0.0, DOUBLE_CLOSE);
//}

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
	unsigned int testWidth = 128;
	unsigned int testHeight = 128;
	double* testData = new double[testWidth * testHeight];
	double inputData = 0.0;
	for(unsigned int i = 0; i < testWidth * testHeight; ++i) {
		inputData += 1.0;
		testData[i] = inputData;
	}

	HaarWaveletTransform haarTransform(testData, testWidth);
	haarTransform.transform();
	haarTransform.inverseTransform();

	for(unsigned int row = 0; row < testHeight; ++row) {
		for(unsigned int column = 0; column < testWidth; ++column) {
			BOOST_CHECK_CLOSE(testData[(row * testWidth) + column], haarTransform.getItem(row, column), DOUBLE_CLOSE);
		}
	}

	delete[] testData;
}

//After chagne in WaveletTransform test not applicable
//BOOST_AUTO_TEST_CASE( HaarTransformSizeDifferentFromPowerOfTwo ) {
//	unsigned int testWidth = 65;
//	unsigned int testHeight = 129;
//	double* testData = new double[testWidth * testHeight];
//	double inputData = 0.0;
//	for(unsigned int i = 0; i < testWidth * testHeight; ++i) {
//		inputData += 1.0;
//		testData[i] = inputData;
//	}
//
//	HaarWaveletTransform haarTransform(testData, testWidth, testHeight);
//	haarTransform.transform();
//	haarTransform.inverseTransform();
//
//	for(unsigned int row = 0; row < testHeight; ++row) {
//		for(unsigned int column = 0; column < testWidth; ++column) {
//			BOOST_CHECK_CLOSE(testData[(row * testWidth) + column], haarTransform.getItem(row, column), DOUBLE_CLOSE);
//		}
//	}
//
//	delete[] testData;
//}

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

BOOST_AUTO_TEST_CASE( setBitFromRightTest ) {
	unsigned int testInt = 0;

	UnsignedInteger::setBitFromRight(&testInt, 0);
	BOOST_CHECK_EQUAL(testInt, 1u);

	UnsignedInteger::setBitFromRight(&testInt, 1);
	BOOST_CHECK_EQUAL(testInt, 3u);

	UnsignedInteger::setBitFromRight(&testInt, 2);
	BOOST_CHECK_EQUAL(testInt, 7u);

	UnsignedInteger::setBitFromRight(&testInt, 3);
	BOOST_CHECK_EQUAL(testInt, 15u);

	testInt = 0;
	UnsignedInteger::setBitFromRight(&testInt, 16);
	BOOST_CHECK_EQUAL(testInt, 65536u);
}

BOOST_AUTO_TEST_CASE( setBitFromLeftTest ) {
	unsigned int testInt = 0;

	UnsignedInteger::setBitFromLeft(&testInt, 31);
	BOOST_CHECK_EQUAL(testInt, 1u);

	UnsignedInteger::setBitFromLeft(&testInt, 30);
	BOOST_CHECK_EQUAL(testInt, 3u);

	UnsignedInteger::setBitFromLeft(&testInt, 29);
	BOOST_CHECK_EQUAL(testInt, 7u);

	UnsignedInteger::setBitFromLeft(&testInt, 28);
	BOOST_CHECK_EQUAL(testInt, 15u);
}

BOOST_AUTO_TEST_CASE( reverseIntTest ) {
	unsigned int toReverse = 0;

	UnsignedInteger::setBitFromLeft(&toReverse, 0);
	UnsignedInteger::setBitFromLeft(&toReverse, 1);
	UnsignedInteger::setBitFromLeft(&toReverse, 3);
	UnsignedInteger::setBitFromLeft(&toReverse, 4);
	UnsignedInteger::setBitFromLeft(&toReverse, 5);
	UnsignedInteger::setBitFromLeft(&toReverse, 8);

	unsigned int reversed = UnsignedInteger::reverse(toReverse);
	BOOST_CHECK_EQUAL(reversed, 315u);

	BOOST_CHECK_EQUAL(UnsignedInteger::reverse(reversed), toReverse);

	toReverse = 0;
	UnsignedInteger::setBitFromRight(&toReverse, 0);
	UnsignedInteger::setBitFromRight(&toReverse, 1);
	reversed = UnsignedInteger::reverse(toReverse);
	unsigned int check = 0;
	UnsignedInteger::setBitFromLeft(&check, 0);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( QuantizerTests )

BOOST_AUTO_TEST_CASE( QuantizerRoundTest ) {
	BOOST_CHECK_EQUAL(0.0, Quantizer::round(-0.3));
	BOOST_CHECK_EQUAL(0.0, Quantizer::round(0.3));
	BOOST_CHECK_EQUAL(-1.0, Quantizer::round(-0.6));
	BOOST_CHECK_EQUAL(1.0, Quantizer::round(0.6));
	BOOST_CHECK_EQUAL(0.0, Quantizer::round(0.00000001));
	BOOST_CHECK_EQUAL(0.0, Quantizer::round(-0.00000001));
	BOOST_CHECK_EQUAL(0.0, Quantizer::round(0.5));
	BOOST_CHECK_EQUAL(-1.0, Quantizer::round(-0.5));
	BOOST_CHECK_EQUAL(10.0, Quantizer::round(10.123456789123456789));
	BOOST_CHECK_EQUAL(-10.0f, Quantizer::round(-10.123456789123456789));
}

BOOST_AUTO_TEST_CASE( QuantizerGetApproximationTest ) {
	BOOST_CHECK_EQUAL(0.0, Quantizer::getApproximation(0.3, 1.0));
	BOOST_CHECK_EQUAL(0.0, Quantizer::getApproximation(-0.3, 1.0));
	BOOST_CHECK_EQUAL(-1.0, Quantizer::getApproximation(-0.5, 1.0));
	BOOST_CHECK_EQUAL(0.0, Quantizer::getApproximation(0.5, 1.0));

	BOOST_CHECK_EQUAL(0.0, Quantizer::getApproximation(0.1 / 2.0, 0.1));
	BOOST_CHECK_EQUAL(-0.1, Quantizer::getApproximation(-0.1 / 2.0, 0.1));
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( HeapTests )

BOOST_AUTO_TEST_CASE( HeapDequeue ) {
	int heapData[10] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
	int result[10] = {10, 9, 8, 7, 6, 5, 4, 3, 2, 1};
	Heap<int> testHeap(heapData, 10);
	testHeap.buildHeap();

	int counter = 0;
	while(testHeap.getHeapSize() > 0) {
		BOOST_REQUIRE(counter < 10);
		BOOST_CHECK_EQUAL(testHeap.dequeue(), result[counter]);
		++counter;
	}
}

BOOST_AUTO_TEST_CASE( HeapEnqueue ) {
	int heapData[20];
	Heap<int> testHeap(heapData, 0, 20);
	int result[20] = {20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1};
	int inputData[20] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20};

	for(int i = 0; i < 20; ++i) {
		testHeap.enqueue(inputData[i]);
	}

	int counter = 0;
	while(testHeap.getHeapSize() > 0) {
		BOOST_REQUIRE(counter < 20);
		BOOST_CHECK_EQUAL(testHeap.dequeue(), result[counter]);
		++counter;
	}
}

BOOST_AUTO_TEST_CASE( QueueItemTest ) {
	HuffmanCoding<char>::Leaf pLeafs[7];
	pLeafs[0].value = 'a';
	pLeafs[0].count = 7;
	pLeafs[1].value = 'b';
	pLeafs[1].count = 10;
	pLeafs[2].value = 'c';
	pLeafs[2].count = 2;
	pLeafs[3].value = 'd';
	pLeafs[3].count = 30;
	pLeafs[4].value = 'e';
	pLeafs[4].count = 11;
	pLeafs[5].value = 'f';
	pLeafs[5].count = 6;
	pLeafs[6].value = 'g';
	pLeafs[6].count = 25;

	HuffmanCoding<char>::QueueItem pQueueItems[7];
	for(int i = 0; i < 7; ++i) {
		pQueueItems[i].node = &pLeafs[i];
	}

	Heap<HuffmanCoding<char>::QueueItem > heap(pQueueItems, 7);
	heap.buildHeap();

	HuffmanCoding<char>::QueueItem queueItem;
	int result[7] = {2, 6, 7, 10, 11, 25, 30};
	int counter = 0;
	while(heap.getHeapSize() > 0) {
		BOOST_REQUIRE(counter < 7);
		queueItem = heap.dequeue();
		BOOST_CHECK_EQUAL(queueItem.node->count, result[counter]);
		++counter;
	}
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( HuffmanCodingTests )

BOOST_AUTO_TEST_CASE( QueueItemTest ) {
	HuffmanCoding<double>::QueueItem item1, item2;
	HuffmanCoding<double>::Node node, node2;
	node.count = 120;
	node2.count = 10;
	HuffmanCoding<double>::Node* nullNode = 0;

	BOOST_CHECK_EQUAL(item1.node, nullNode);

	item1.node = &node;

	item2 = item1;
	BOOST_CHECK_EQUAL(item2.node, &node);

	HuffmanCoding<double>::QueueItem item3(item2);
	BOOST_CHECK_EQUAL(item3.node, &node);

	HuffmanCoding<double>::QueueItem item4;
	item4.node = &node2;
	BOOST_CHECK_EQUAL(item4 > item1, true);
	BOOST_CHECK_EQUAL(item1 > item4, false);
	BOOST_CHECK_EQUAL(item4 < item1, false);
	BOOST_CHECK_EQUAL(item1 < item4, true);
}

BOOST_AUTO_TEST_CASE( HuffmanCodingSimpleTest ) {
	HuffmanCoding<char>::Leaf* pLeafs = new HuffmanCoding<char>::Leaf[7];
	pLeafs[0].value = 'a';
	pLeafs[0].count = 7;
	pLeafs[1].value = 'b';
	pLeafs[1].count = 10;
	pLeafs[2].value = 'c';
	pLeafs[2].count = 2;
	pLeafs[3].value = 'd';
	pLeafs[3].count = 30;
	pLeafs[4].value = 'e';
	pLeafs[4].count = 11;
	pLeafs[5].value = 'f';
	pLeafs[5].count = 6;
	pLeafs[6].value = 'g';
	pLeafs[6].count = 25;
	HuffmanCoding<char> huffmanCoding(pLeafs, 7);

	huffmanCoding.createCodeTable();

	std::map<char, HuffmanCoding<char>::Code > codeTable;
	bool createResult = huffmanCoding.getTable(codeTable);
	BOOST_CHECK(createResult);

	std::map<char, HuffmanCoding<char>::Code >::iterator it;
	it = codeTable.find('f');
	BOOST_CHECK(it != codeTable.end());

	unsigned int testCode = 0;
	UnsignedInteger::setBitFromLeft(&testCode, 2);
	UnsignedInteger::setBitFromLeft(&testCode, 3);
	BOOST_CHECK_EQUAL(testCode, it->second.code);
	BOOST_CHECK_EQUAL(4u, it->second.size);

	it = codeTable.find('a');
	BOOST_CHECK(it != codeTable.end());
	BOOST_CHECK_EQUAL(0u, it->second.code);

	it = codeTable.find('c');
	BOOST_CHECK(it != codeTable.end());
	testCode = 0;
	UnsignedInteger::setBitFromLeft(&testCode, 2);
	BOOST_CHECK_EQUAL(testCode, it->second.code);

	it = codeTable.find('b');
	BOOST_CHECK(it != codeTable.end());
	testCode = 0;
	UnsignedInteger::setBitFromLeft(&testCode, 1);
	BOOST_CHECK_EQUAL(testCode, it->second.code);

	it = codeTable.find('e');
	BOOST_CHECK(it != codeTable.end());
	testCode = 0;
	UnsignedInteger::setBitFromLeft(&testCode, 1);
	UnsignedInteger::setBitFromLeft(&testCode, 2);
	BOOST_CHECK_EQUAL(testCode, it->second.code);

	it = codeTable.find('g');
	BOOST_CHECK(it != codeTable.end());
	testCode = 0;
	UnsignedInteger::setBitFromLeft(&testCode, 0);
	BOOST_CHECK_EQUAL(testCode, it->second.code);

	it = codeTable.find('d');
	BOOST_CHECK(it != codeTable.end());
	testCode = 0;
	UnsignedInteger::setBitFromLeft(&testCode, 0);
	UnsignedInteger::setBitFromLeft(&testCode, 1);
	BOOST_CHECK_EQUAL(testCode, it->second.code);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( RLETests )

BOOST_AUTO_TEST_CASE( RLEConstructionTest ) {
	int data[10] = {1, 1, 1, 1, 5, 5, 5, 5, 2, 2};
	RLE<int> rle(data, 10);

	BOOST_CHECK_EQUAL(rle.getEncodedDataSize(), 0u);
	BOOST_CHECK(rle.getData() != 0);

	rle.encode();
	BOOST_REQUIRE_EQUAL(rle.getEncodedDataSize(), 3u);

	RLE<int>::Run* pData = rle.getData();
	BOOST_CHECK_EQUAL(pData[0].value, 1);
	BOOST_CHECK_EQUAL(pData[0].run, 4u);

	BOOST_CHECK_EQUAL(pData[1].value, 5);
	BOOST_CHECK_EQUAL(pData[1].run, 4u);

	BOOST_CHECK_EQUAL(pData[2].value, 2);
	BOOST_CHECK_EQUAL(pData[2].run, 2u);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( EncoderTests )

struct CodesFixture {
	unsigned int code1;
	unsigned int size1;
	unsigned int code2;
	unsigned int size2;
	unsigned int code3;
	unsigned int size3;

	CodesFixture() {
		code1 = 0;
		UnsignedInteger::setBitFromLeft(&code1, 1);
		UnsignedInteger::setBitFromLeft(&code1, 2);
		size1 = 4;

		code2 = 0;
		UnsignedInteger::setBitFromLeft(&code2, 0);
		UnsignedInteger::setBitFromLeft(&code2, 1);
		UnsignedInteger::setBitFromLeft(&code2, 5);
		UnsignedInteger::setBitFromLeft(&code2, 6);
		UnsignedInteger::setBitFromLeft(&code2, 7);
		size2 = 8;

		code3 = 0;
		UnsignedInteger::setBitFromLeft(&code3, 2);
		size3 = 3;
	}
	~CodesFixture() { }
};

BOOST_FIXTURE_TEST_CASE(EncodeBitsTest, CodesFixture) {
	Encoder encoder(5);

	encoder.encode(code1, size1);
	encoder.encode(code1, size1);
	encoder.encode(code1, size1);
	encoder.encode(code2, size2);
	encoder.encode(code1, size1);
	encoder.encode(code3, size3);
	encoder.encode(code2, size2);

	BOOST_CHECK_EQUAL(encoder.encodedSize(), 2u);
	unsigned int testValue = 0;
	UnsignedInteger::setBitFromLeft(&testValue, 0);
	UnsignedInteger::setBitFromLeft(&testValue, 1);
	UnsignedInteger::setBitFromLeft(&testValue, 2);
	BOOST_CHECK_EQUAL(encoder.getData()[1], testValue);

	encoder.encode(code3, size3);
	UnsignedInteger::setBitFromLeft(&testValue, 5);
	BOOST_CHECK_EQUAL(encoder.getData()[1], testValue);

	encoder.encode(code3, size3);
	UnsignedInteger::setBitFromLeft(&testValue, 8);
	BOOST_CHECK_EQUAL(encoder.getData()[1], testValue);
}

BOOST_FIXTURE_TEST_CASE(EncoderReallocateTest, CodesFixture) {
	Encoder encoder(1);
	BOOST_CHECK_EQUAL(encoder.getSize(), 1u);

	encoder.encode(code1, size1);
	encoder.encode(code1, size1);
	encoder.encode(code1, size1);
	encoder.encode(code2, size2);
	encoder.encode(code1, size1);
	encoder.encode(code3, size3);
	encoder.encode(code2, size2);

	BOOST_CHECK_EQUAL(encoder.getSize(), 2u);

	BOOST_CHECK_EQUAL(encoder.encodedSize(), 2u);
	unsigned int testValue = 0;
	UnsignedInteger::setBitFromLeft(&testValue, 0);
	UnsignedInteger::setBitFromLeft(&testValue, 1);
	UnsignedInteger::setBitFromLeft(&testValue, 2);
	BOOST_CHECK_EQUAL(encoder.getData()[1], testValue);

	encoder.encode(code2, size2);
	encoder.encode(code2, size2);
	encoder.encode(code2, size2);
	encoder.encode(code3, size3);
	encoder.encode(code2, size2);

	testValue = 0;
	UnsignedInteger::setBitFromLeft(&testValue, 3);
	UnsignedInteger::setBitFromLeft(&testValue, 4);
	UnsignedInteger::setBitFromLeft(&testValue, 5);
	BOOST_REQUIRE_EQUAL(encoder.encodedSize(), 3u);
	BOOST_CHECK_EQUAL(encoder.getData()[2], testValue);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( HuffmanDecoderTests )

BOOST_AUTO_TEST_CASE( HuffmanDecoderTest ) {
	HuffmanCoding<double>::Leaf* pLeafs = new HuffmanCoding<double>::Leaf[5];
	pLeafs[0].value = 1.0;
	pLeafs[0].count = 5;
	pLeafs[1].value = 2.0;
	pLeafs[1].count = 2;
	pLeafs[2].value = 3.0;
	pLeafs[2].count = 3;
	pLeafs[3].value = 4.0;
	pLeafs[3].count = 3;
	pLeafs[4].value = 5.0;
	pLeafs[4].count = 1;

	HuffmanCoding<double>::Leaf* nullLeaf = 0;
	BOOST_REQUIRE_EQUAL(pLeafs[0].left, nullLeaf);
	BOOST_REQUIRE_EQUAL(pLeafs[0].right, nullLeaf);
	BOOST_REQUIRE_EQUAL(pLeafs[0].parent, nullLeaf);
	BOOST_REQUIRE_EQUAL(pLeafs[1].left, nullLeaf);
	BOOST_REQUIRE_EQUAL(pLeafs[1].right, nullLeaf);
	BOOST_REQUIRE_EQUAL(pLeafs[1].parent, nullLeaf);
	BOOST_REQUIRE_EQUAL(pLeafs[2].left, nullLeaf);
	BOOST_REQUIRE_EQUAL(pLeafs[2].right, nullLeaf);
	BOOST_REQUIRE_EQUAL(pLeafs[2].parent, nullLeaf);
	BOOST_REQUIRE_EQUAL(pLeafs[3].left, nullLeaf);
	BOOST_REQUIRE_EQUAL(pLeafs[3].right, nullLeaf);
	BOOST_REQUIRE_EQUAL(pLeafs[3].parent, nullLeaf);
	BOOST_REQUIRE_EQUAL(pLeafs[4].left, nullLeaf);
	BOOST_REQUIRE_EQUAL(pLeafs[4].right, nullLeaf);
	BOOST_REQUIRE_EQUAL(pLeafs[4].parent, nullLeaf);

	double data[14] = {1.0, 2.0, 5.0, 2.0, 1.0, 3.0, 4.0, 1.0, 3.0, 4.0, 1.0, 4.0, 3.0, 1.0};
	HuffmanCoding<double> huffmanCoding(pLeafs, 5);
	huffmanCoding.createCodeTable();
	std::map<double, HuffmanCoding<double>::Code > codeTable;
	huffmanCoding.getTable(codeTable);

	unsigned int testCode = 0;
	BOOST_CHECK_EQUAL(testCode, codeTable[4.0].code);
	BOOST_CHECK_EQUAL(2u, codeTable[4.0].size);

	UnsignedInteger::setBitFromLeft(&testCode, 0);
	UnsignedInteger::setBitFromLeft(&testCode, 1);
	BOOST_CHECK_EQUAL(testCode, codeTable[1.0].code);
	BOOST_CHECK_EQUAL(2u, codeTable[1.0].size);

	testCode = 0;
	UnsignedInteger::setBitFromLeft(&testCode, 0);
	BOOST_CHECK_EQUAL(testCode, codeTable[3.0].code);
	BOOST_CHECK_EQUAL(2u, codeTable[3.0].size);

	testCode = 0;
	UnsignedInteger::setBitFromLeft(&testCode, 1);
	UnsignedInteger::setBitFromLeft(&testCode, 2);
	BOOST_CHECK_EQUAL(testCode, codeTable[2.0].code);
	BOOST_CHECK_EQUAL(3u, codeTable[2.0].size);

	testCode = 0;
	UnsignedInteger::setBitFromLeft(&testCode, 1);
	BOOST_CHECK_EQUAL(testCode, codeTable[5.0].code);
	BOOST_CHECK_EQUAL(3u, codeTable[5.0].size);

	std::map<HuffmanCoding<double>::Code, double > decodeTable;
	std::map<double, HuffmanCoding<double>::Code>::iterator it;
	for(it = codeTable.begin(); it != codeTable.end(); ++it) {
		decodeTable.insert(std::pair<HuffmanCoding<double>::Code, double >(it->second, it->first));
	}

	Encoder encoder(5);
	for(int i = 0; i < 14; ++i) {
		it = codeTable.find(data[i]);
		BOOST_REQUIRE(it != codeTable.end());
		encoder.encode(it->second.code, it->second.size);
	}

	unsigned int expectedEncodedBits = 0;
	for(unsigned int i = 0; i < 14; ++i) {
		expectedEncodedBits += codeTable[data[i]].size;
	}
	unsigned int expectedEncodedBytes = (unsigned int)ceil((double)expectedEncodedBits / (double)UnsignedInteger::NUMBER_OF_BITS);
	BOOST_CHECK_EQUAL(expectedEncodedBytes, encoder.encodedSize());

	HuffmanDecoder<double> huffmanDecoder(encoder.getData(), encoder.encodedSize(), &decodeTable, 14);
	huffmanDecoder.decode();

	double* pDecodedData = huffmanDecoder.getDecodedData();
	BOOST_REQUIRE_EQUAL(huffmanDecoder.getDecodedDataSize(), 14u);
	for(unsigned int i = 0; i < huffmanDecoder.getDecodedDataSize(); ++i) {
		BOOST_CHECK_EQUAL(pDecodedData[i], data[i]);
	}
}

BOOST_AUTO_TEST_CASE( HuffmanDecoderSimpleLongCode ) {
	std::map<HuffmanCoding<char>::Code, char > decodeTable;

	unsigned int codeA = 0;
	UnsignedInteger::setBitFromLeft(&codeA, 28);
	HuffmanCoding<char>::Code code;
	code.code = codeA;
	code.size = 29;
	decodeTable.insert(std::pair<HuffmanCoding<char>::Code, char >(code, 'a'));

	unsigned int codeB = 0;
	UnsignedInteger::setBitFromLeft(&codeB, 0);
	UnsignedInteger::setBitFromLeft(&codeB, 2);
	UnsignedInteger::setBitFromLeft(&codeB, 3);
	UnsignedInteger::setBitFromLeft(&codeB, 5);
	UnsignedInteger::setBitFromLeft(&codeB, 6);
	code.code = codeB;
	code.size = 7;
	decodeTable.insert(std::pair<HuffmanCoding<char>::Code, char >(code, 'b'));

	Encoder encoder(10);
	encoder.encode(codeA, 29);
	encoder.encode(codeB, 7);
	encoder.encode(codeB, 7);

	HuffmanDecoder<char> huffmanDecoder(encoder.getData(), encoder.encodedSize(), &decodeTable, 3);
	huffmanDecoder.decode();

	BOOST_REQUIRE_EQUAL(huffmanDecoder.getDecodedDataSize(), 3u);
	BOOST_CHECK_EQUAL(huffmanDecoder.getDecodedData()[0], 'a');
	BOOST_CHECK_EQUAL(huffmanDecoder.getDecodedData()[1], 'b');
	BOOST_CHECK_EQUAL(huffmanDecoder.getDecodedData()[2], 'b');
}

BOOST_AUTO_TEST_CASE( HuffmanDecoderCode32 ) {
	std::map<HuffmanCoding<char>::Code, char > decodeTable;

	unsigned int codeA = 0;
	UnsignedInteger::setBitFromLeft(&codeA, 3);
	HuffmanCoding<char>::Code code;
	code.code = codeA;
	code.size = 4;
	decodeTable.insert(std::pair<HuffmanCoding<char>::Code, char >(code, 'a'));

	unsigned int codeB = 0;
	UnsignedInteger::setBitFromLeft(&codeB, 0);
	UnsignedInteger::setBitFromLeft(&codeB, 2);
	code.code = codeB;
	code.size = 3;
	decodeTable.insert(std::pair<HuffmanCoding<char>::Code, char >(code, 'b'));

	Encoder encoder(10);
	for(int i = 0; i < 8; ++i) {
		encoder.encode(codeA, 4);
	}
	encoder.encode(codeB, 3);
	encoder.encode(codeB, 3);
	encoder.encode(codeB, 3);

	HuffmanDecoder<char> huffmanDecoder(encoder.getData(), encoder.encodedSize(), &decodeTable, 11);
	huffmanDecoder.decode();

	BOOST_REQUIRE_EQUAL(huffmanDecoder.getDecodedDataSize(), 11u);
	for(int i = 0; i < 8; ++i) {
		BOOST_CHECK_EQUAL(huffmanDecoder.getDecodedData()[i], 'a');
	}
	for(int i = 8; i < 11; ++i) {
		BOOST_CHECK_EQUAL(huffmanDecoder.getDecodedData()[i], 'b');
	}
}

BOOST_AUTO_TEST_CASE( HuffmanDecoderLongCode ) {
	HuffmanCoding<char>::Leaf* pLeafs = new HuffmanCoding<char>::Leaf[8];
	pLeafs[0].value = 'a';
	pLeafs[0].count = 8;
	pLeafs[1].value = 'b';
	pLeafs[1].count = 5;
	pLeafs[2].value = 'c';
	pLeafs[2].count = 10;
	pLeafs[3].value = 'd';
	pLeafs[3].count = 12;
	pLeafs[4].value = 'e';
	pLeafs[4].count = 4;
	pLeafs[5].value = 'f';
	pLeafs[5].count = 3;
	pLeafs[6].value = 'g';
	pLeafs[6].count = 7;
	pLeafs[7].value = 'h';
	pLeafs[7].count = 15;

	HuffmanCoding<char> huffmanCoding(pLeafs, 8);
	huffmanCoding.createCodeTable();
	std::map<char, HuffmanCoding<char>::Code > codeTable;
	huffmanCoding.getTable(codeTable);

	std::map<char, HuffmanCoding<char>::Code >::iterator it, itPrefix;
	std::map<HuffmanCoding<char>::Code, char > decodeTable;
	for(it = codeTable.begin(); it != codeTable.end(); ++it) {
		decodeTable.insert(std::pair<HuffmanCoding<char>::Code, char >(it->second, it->first));
	}

	for(it = codeTable.begin(); it != codeTable.end(); ++it) {
//		printf("%c: %d\n", it->first, it->second.size);

		for(itPrefix = codeTable.begin(); itPrefix != codeTable.end(); ++itPrefix) {
			if(it == itPrefix) {
				continue;
			}

			BOOST_CHECK_MESSAGE(isPrefix(itPrefix->second.code,  it->second.code, itPrefix->second.size) == false,
					"code for: " << itPrefix->first << " is prefix of: " << it->first);
		}
	}

//	printf("code for d: %d\n", codeTable['d'].code);

	char data[64] = {
			'a', 'a', 'b', 'b', 'b', 'c', 'b', 'b',
			'a', 'a', 'a', 'c', 'c', 'c', 'c', 'a',
			'c', 'c', 'a', 'a', 'c', 'h', 'h', 'h',
			'h', 'c', 'h', 'c', 'h', 'h', 'h', 'h',
			'f', 'h', 'f', 'h', 'f', 'h', 'h', 'h',
			'h', 'd', 'd', 'e', 'd', 'g', 'g', 'g',
			'e', 'd', 'd', 'g', 'd', 'd', 'd', 'g',
			'd', 'd', 'g', 'd', 'e', 'g', 'e', 'd'
	};
	Encoder encoder(10);

//	printf("\n\n");
	for(int i = 0; i < 64; ++i) {
		it = codeTable.find(data[i]);
		BOOST_REQUIRE(it != codeTable.end());
		encoder.encode(it->second.code, it->second.size);
	}

	HuffmanDecoder<char> huffmanDecoder(encoder.getData(), encoder.encodedSize(), &decodeTable, 64);
	huffmanDecoder.decode();

//	for(unsigned int i = 0; i < huffmanDecoder.getDecodedDataSize(); ++i) {
//		printf("%c ", huffmanDecoder.getDecodedData()[i]);
//		if((i + 1) % 8 == 0) {
//			printf("\n");
//		}
//	}
//	printf("\n");

	BOOST_REQUIRE_EQUAL(huffmanDecoder.getDecodedDataSize(), 64u);
	char* decodedData = huffmanDecoder.getDecodedData();
	for(unsigned int i = 0; i < huffmanDecoder.getDecodedDataSize(); ++i) {
		BOOST_CHECK_EQUAL(decodedData[i], data[i]);
	}
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( RLEDecoderTests )

BOOST_AUTO_TEST_CASE( RLEDecoderDecodeTest) {
	RLE<double>::Run runs[5];
	runs[0].run = 5;
	runs[0].value = 1.0;
	runs[1].run = 1;
	runs[1].value = 2.0;
	runs[2].run = 3;
	runs[2].value = 3.0;
	runs[3].run = 6;
	runs[3].value = 4.0;
	runs[4].run = 2;
	runs[4].value = 5.0;

	double decodeMemory[17];
	RLEDecoder<double> rleDecoder(runs, 5);
	bool decodeResult = rleDecoder.decode(decodeMemory, 17);
	BOOST_CHECK(decodeResult);

	double result[17] = {1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 3.0, 3.0, 3.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 5.0, 5.0};
	for(int i = 0; i < 17; ++i) {
		BOOST_CHECK_EQUAL(decodeMemory[i], result[i]);
	}

	RLEDecoder<double> rleDecoder2(runs, 5);
	double decode1[6];
	double decode2[9];
	rleDecoder2.decode(decode1, 6);
	rleDecoder2.decode(decode2, 9);

	double decode1Result[6] = {1.0, 1.0, 1.0, 1.0, 1.0, 2.0};
	double decode2Result[9] = {3.0, 3.0, 3.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0};

	for(int i = 0; i < 6; ++i) {
		BOOST_CHECK_EQUAL(decode1[i], decode1Result[i]);
	}
	for(int i = 0; i < 9; ++i) {
		BOOST_CHECK_EQUAL(decode2[i], decode2Result[i]);
	}
}

BOOST_AUTO_TEST_CASE( RLEDecoderEncodeDecodeTest ) {
	double data[20] = {
			1.0, 1.0, 2.0, 3.0, 3.0,
			5.0, 5.0, 5.0, 5.0, 5.0,
			4.0, 4.0, 0.0, 0.0, 0.0,
			0.0, 0.0, 0.0, 0.0, 0.0
	};

	RLE<double> rleCoder(data, 20);
	rleCoder.encode();

	double result[20] = {0.0};
	RLEDecoder<double> rleDecoder(rleCoder.getData(), rleCoder.getEncodedDataSize());
	bool decodeResult = rleDecoder.decode(result, 20);
	BOOST_CHECK(decodeResult);

	for(int i = 0; i < 20; ++i) {
		BOOST_CHECK_EQUAL(data[i], result[i]);
	}
}

BOOST_AUTO_TEST_SUITE_END()
