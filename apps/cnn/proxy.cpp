#include <cstdlib>
#include <iostream>
#include <deque>
#include <chrono>
#include <iomanip>
#include <assert.h>
#include "../../core/cnn.h"
#include <sstream>
using namespace std;
uint64_t*** CL(Party* proxy, uint64_t** input, uint64_t i_size, uint64_t** kernel, uint32_t k_size, uint32_t k_number, uint8_t stride);
uint64_t*** CL(Party* proxy, uint64_t*** input, uint64_t i_dim, uint32_t i_number, uint64_t* kernel, uint32_t k_dim, uint8_t stride);

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

    uint32_t matrix_size = atoi(argv[8]);
    uint32_t stride = atoi(argv[9]);

    //ensure ports are not 0 for helper or client
    cout << "Setting ports for helper/client..." << endl;
    if (cport != 0) {
        assert(cport < 1 << (sizeof(uint16_t) * 8));
    }
    if (hport != 0) {
        assert(hport < 1 << (sizeof(uint16_t) * 8));
    }

    Party *proxy;
    cout << "Creating Proxy..." << endl;
    if (role==0)
        proxy = new Party(P1,hport, haddress, cport, caddress);
    else
        proxy = new Party(P2,hport, haddress, cport, caddress);

    //init TODO reading data at this point...
    double** doubleMatrix = new double *[matrix_size];

    //fill matrix (as vector) with params from input
    stringstream mss(matrix_string);
    for (uint32_t j = 0; j < matrix_size; j++){
        doubleMatrix[j] = new double [matrix_size];
        uint32_t i = 0;
        while(mss.good() && i<matrix_size){
            string elements;
            getline(mss, elements, ',');
            doubleMatrix[j][i] = std::atof(elements.c_str());
        }
        i++;
    }
    //convert matrix elements
    uint64_t** matrix = convert2uint64(doubleMatrix, matrix_size, matrix_size, 15);
    cout << "created matrix" << endl;
    // init input parameter
    uint64_t *data = new uint64_t [128];
    uint64_t *weights = new uint64_t [128];
    for (uint8_t i = 0; i < 128; i++){
        weights[i] = 1; //TODO where do we get the weights of the layers from? adjustable via input file or fixed?
    }
    cout << "init weights done" << endl;
    uint64_t **kernel = new uint64_t*[5];
    for (uint8_t i = 0; i < 5; i++) {
        kernel[i] = new uint64_t[25];
        for (uint8_t j = 0; j < 25; j++) {
            kernel[i][j] = proxy->generateRandom();
        }
    }
    cout << "init kernel done" << endl;
    //helper variable to create shares of v1 and v2
    uint64_t **mTmp = new uint64_t *[matrix_size]; //pick random value in ring --> shareOfMatrix is calculated from that
    uint64_t **shareOfMatrix = new uint64_t *[matrix_size];

    for(uint32_t i = 0; i < matrix_size; i++){
        mTmp[i] = new uint64_t[matrix_size];
        for(uint32_t j = 0; j < matrix_size; j++) {
            mTmp[i][j] = proxy->generateCommonRandom();
        }
    }
    cout << "Start creating share for proxy " << to_string(role) << endl;
    //create shares
    if(role == P1) {
        for(uint32_t k = 0; k < matrix_size; k++){
            shareOfMatrix[k] = new uint64_t[matrix_size];
            for(uint32_t j = 0; j < matrix_size; j++) {
                shareOfMatrix[k][j] = matrix[k][j] - mTmp[k][j];
            }
        }
    }
    else {
        shareOfMatrix = mTmp;
    }
    cout << "finished creating shares... Call Convolutional layer" << endl;

    // CNN inference pipeline, call functions sequantially for inference (Matrix MUL, MAXPOOL, RELU, vs...)

    // convolutional layer
    uint64_t mmaxParams[4];
    uint64_t kernelNum = 5;
    mmaxParams[0] = matrix_size;
    mmaxParams[1] = 5; // kernel size
    mmaxParams[2] = kernelNum; // kernel number
    mmaxParams[3] = stride; // stride

    // PERFORMING CONVOLUTION
    proxy->SendBytes(CNN_CL1);

    unsigned char *ptr_out = proxy->getBuffer1();
    addVal2CharArray(mmaxParams, &ptr_out, 4);
    Send(proxy->getSocketHelper(), proxy->getBuffer1(), 4 * 8);

    uint64_t*** conv1 = CL(proxy, matrix, matrix_size, kernel, 5, kernelNum, stride);
    //          size after conv                   after maxpool divide by 2
    matrix_size = ((matrix_size - 5)/stride + 1) / 2;

    // convolutional layer 2
    uint64_t *kernel2 = new uint64_t[25];
    for (uint8_t j = 0; j < 25; j++) {
        kernel2[j] = proxy->generateRandom();
    }

    mmaxParams[0] = matrix_size;
    mmaxParams[1] = kernelNum; // matrix number is same as previous kernel number
    mmaxParams[2] = 5; // kernel size
    mmaxParams[3] = stride; // stride

    // PERFORMING CONVOLUTION
    proxy->SendBytes(CNN_CL2);

    unsigned char *ptr_out2 = proxy->getBuffer1();
    addVal2CharArray(mmaxParams, &ptr_out2, 4);
    Send(proxy->getSocketHelper(), proxy->getBuffer1(), 4 * 8);

    uint64_t*** conv2 = CL(proxy, conv1, matrix_size, kernelNum, kernel2, 5, stride);
    matrix_size = ((matrix_size - 5)/stride + 1) / 2;
    delete conv1;

    // fully connected layer:
    uint32_t nodes_out = 100;
    uint64_t **weights2 = new uint64_t*[nodes_out];

    for (uint32_t n = 0; n<nodes_out; n++){
        weights2[n] = new uint64_t [matrix_size*matrix_size];
        for (uint32_t j = 0; j < matrix_size*matrix_size; j++) {
            weights2[n][j] = proxy->generateRandom();
        }
    }

    mmaxParams[0] = matrix_size;
    mmaxParams[1] = kernelNum; // input number is same as previous kernel number
    mmaxParams[2] = nodes_out;

    // PERFORMING CONVOLUTION
    proxy->SendBytes(CNN_FCL);

    unsigned char *ptr_out3 = proxy->getBuffer1();
    addVal2CharArray(mmaxParams, &ptr_out3, 3);
    Send(proxy->getSocketHelper(), proxy->getBuffer1(), 3 * 8);
    uint64_t * output = FCL(proxy, conv2, matrix_size, kernelNum, weights2, nodes_out);

    for(uint64_t n = 0; n<nodes_out; n++){
        cout << "Node " << n << ": " << output[n] << endl;
    }


    proxy->SendBytes(CORE_END);
    proxy->PrintBytes();
    return 0;
}

