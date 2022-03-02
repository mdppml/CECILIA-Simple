#include <cstdlib>
#include <iostream>
#include <deque>
#include <chrono>
#include <iomanip>
#include <assert.h>
#include "../../core/cnn.h"
using namespace std;


void print1DArray(string const &str1, double* x, uint32_t size) {
    cout << "======================= " << str1 << " =======================" << endl;
    for(uint32_t i = 0; i < size; i++) {
        cout << x[i] << "\t";
    }
    cout << endl;
    cout << "==============================================================" << endl;
}

int main(int argc, char* argv[]) {
    uint8_t role = atoi(argv[1]);
    uint16_t cport = atoi(argv[2]);
    string caddress(argv[3]);
    uint16_t hport = atoi(argv[4]);
    string haddress(argv[5]);

    uint16_t window_size = atoi(argv[6]);
    string matrix_string = argv[7];

    uint32_t matrix_col_size = atoi(argv[8]);
    uint32_t matrix_row_size = atoi(argv[9]);

    //ensure ports are not 0 for helper or client
    cout << "Setting ports for helper/client...";
    if (cport != 0) {
        assert(cport < 1 << (sizeof(uint16_t) * 8));
    }
    if (hport != 0) {
        assert(hport < 1 << (sizeof(uint16_t) * 8));
    }

    Party *proxy;
    cout << "Creating Proxy...";
    if (role==0)
        proxy = new Party(P1,hport, haddress, cport, caddress);
    else
        proxy = new Party(P2,hport, haddress, cport, caddress);

    //init
    uint16_t matrix_size = matrix_row_size * matrix_col_size;
    double *doubleMatrix = new double [matrix_size];

    //fill matrix (as vector) with params from input
    stringstream mss(matrix_string);
    uint32_t i = 0;
    while(mss.good() && i<matrix_size){
        string elements;
        getline(mss, elements, ',');
        double element = std::atof(elements.c_str());
        doubleMatrix[i] = element;
        i++;
    }
    //convert matrix elements
    uint64_t *matrix = convert2uint64(doubleMatrix, matrix_size, 15);
    uint64_t *resortedMatrix = new uint64_t [matrix_size];
    //                                                     here same param twice but could be different sizes
    RST(matrix, matrix_col_size, matrix_row_size, window_size, window_size, resortedMatrix);

    //helper variable to create shares of v1 and v2
    uint64_t *mTmp = new uint64_t [matrix_size]; //pick random value in ring --> shareOfMatrix is calculated from that
    uint64_t *shareOfMatrix;
    
    for(uint32_t i = 0; i < matrix_size; i++){
        mTmp[i] = proxy->generateCommonRandom();
    }
    cout << "Start creating share for proxy " << to_string(role) << endl;
    //create shares
    if(role == P1) {
        for(uint32_t i = 0; i < matrix_size; i++){
            shareOfMatrix[i] = resortedMatrix[i] - mTmp[i];
        }
    }
    else {
        shareOfMatrix = mTmp;
    }

    cout << "___MMAXPOOL___: " << endl;
    cout << "Ground Truth: " << endl;
    print1DMatrixByWindows("matrix for MAX", doubleMatrix, matrix_row_size, matrix_col_size, 1, 1);
    print1DMatrixByWindows("matrix for MMAX", doubleMatrix, matrix_row_size, matrix_col_size, window_size, window_size);
    delete []matrix;

    proxy->SendBytes(CNN_MAX, matrix_size);
    uint64_t maxElement = MAX(proxy, shareOfMatrix, matrix_size);
    //print1DMatrixByWindows("resorted matrix", resortedMatrix, matrix_row_size, matrix_col_size, window_size, window_size);

    uint16_t window_length = window_size*window_size;
    uint16_t numberOfWins = static_cast<uint16_t >(floor(matrix_size / window_length));
    cout << "Final result MAX_VAL: " << REC(proxy, maxElement) << " which is " << convert2double(maxElement) << endl;


    uint64_t *maxElements;
    uint32_t mmaxParams;
    //contains matrix_size bitwise at 16 MSBs and window size at the 16 LSBs
    mmaxParams = ((uint64_t)matrix_row_size << 48) + ((uint64_t)matrix_col_size << 32) + (uint64_t) window_size;
    proxy->SendBytes(CNN_MMAX, mmaxParams);
    maxElements = MAX(proxy, shareOfMatrix, matrix_row_size, matrix_col_size, window_size);

    print1DMatrixByWindows("Final result MMAX", convert2double(REC(proxy, maxElements, numberOfWins), numberOfWins, 0), 1, numberOfWins, 1, 1);


    // CNN inference pipeline, call functions sequantially for inference (Matrix MUL, MAXPOOL, RELU, vs...)


    proxy->SendBytes(CORE_END);
    proxy->PrintBytes();
    return 0;
}

