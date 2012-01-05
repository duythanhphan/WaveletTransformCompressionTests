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

const double DOUBLE_CLOSE = 0.0001;

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

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( HaarTransformTests )

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

	delete[] testData;
}

BOOST_AUTO_TEST_CASE( HaarTransformSizeDifferentFromPowerOfTwo ) {
	unsigned int testWidth = 65;
	unsigned int testHeight = 129;
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

	delete[] pLeafs;
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
