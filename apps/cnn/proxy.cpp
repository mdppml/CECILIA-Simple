#include <cstdlib>
#include <iostream>
#include <fstream>
#include <deque>
#include <chrono>
#include <iomanip>
#include <assert.h>
#include "../../core/cnn.h"
#include <sstream>
//#include "onnx_pb.h"+
//#include "../../submodules/onnxruntime"

using namespace std;
//uint64_t*** CL(Party* proxy, uint64_t** input, uint64_t i_size, uint64_t** kernel, uint32_t k_size, uint32_t k_number, uint8_t stride);
//uint64_t*** CL(Party* proxy, uint64_t*** input, uint64_t i_dim, uint32_t i_number, uint64_t* kernel, uint32_t k_dim, uint8_t stride);

void print1DArray(string const &str1, double* x, uint32_t size) {
    cout << "======================= " << str1 << " =======================" << endl;
    for(uint32_t i = 0; i < size; i++) {
        cout << x[i] << "\t";
    }
    cout << endl;
    cout << "==============================================================" << endl;
}

int main(int argc, char* argv[]) {
    if (argc < 6){
        cout << "Calling proxy without specifying role (1), port (2), address (3), helpers port (4) and helpers adress (5) is not possible." << endl;
        return -1;
    }
    uint8_t role = atoi(argv[1]);
    uint16_t cport = atoi(argv[2]);
    string caddress(argv[3]);
    uint16_t hport = atoi(argv[4]);
    string haddress(argv[5]);

    string file_path = argv[6];

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

    // init input parameter
    uint64_t matrix_size = 784;
    uint64_t** data = new uint64_t* [28]; //TODO this are the images
    double** model_weights;
    double** weight_dims;

    //read weights of each layer from model:
    ifstream weights_file;
    string line;
    uint64_t pos = 0;
    uint64_t number_of_layer = 0;
    string layer_dim; // dimension of the weights matrix of one layer
    string delim = ",";
    weights_file.open(file_path);
    if (weights_file.is_open()){
        while (getline(weights_file,line)){
            // line starts with L= and contains number of layers afterwards.
            if (line.find("L=") == 0 && line.length() >= 2){
                number_of_layer = stoi(line.substr(2));
            }
            // line starts with * and contains dimensions of this layer weights as k,l
            if (line[0] == '*'){
                // parse dimensions of layer:
                while ((pos = line.find(delim)) != string::npos) {
                    layer_dim = line.substr(0, pos);
                    line.erase(0, pos + delim.length()); // remove the part which was parsed
                }
            }
        }
        weights_file.close();
    }

    //convert matrix elements

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
                shareOfMatrix[k][j] = data[k][j] - mTmp[k][j];
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
    mmaxParams[3] = 1; // stride

    // PERFORMING CONVOLUTION
    proxy->SendBytes(CNN_CL);

    unsigned char *ptr_out = proxy->getBuffer1();
    addVal2CharArray(mmaxParams, &ptr_out, 4);
    Send(proxy->getSocketHelper(), proxy->getBuffer1(), 4 * 8);

    uint64_t*** conv1 = CL(proxy, data, matrix_size, kernel, 5, kernelNum, 1);
    //          size after conv                   after maxpool divide by 2
    matrix_size = ((matrix_size - 5) + 1) / 2;

    // convolutional layer 2
    uint64_t *kernel2 = new uint64_t[25];
    for (uint8_t j = 0; j < 25; j++) {
        kernel2[j] = proxy->generateRandom();
    }

    mmaxParams[0] = matrix_size;
    mmaxParams[1] = kernelNum; // matrix number is same as previous kernel number
    mmaxParams[2] = 5; // kernel size
    mmaxParams[3] = 1; // stride

    // PERFORMING CONVOLUTION
    proxy->SendBytes(CNN_CL2);
    cout << "CL2 is " << CNN_CL2 << endl;

    unsigned char *ptr_out2 = proxy->getBuffer1();
    addVal2CharArray(mmaxParams, &ptr_out2, 4);
    Send(proxy->getSocketHelper(), proxy->getBuffer1(), 4 * 8);

    uint64_t*** conv2 = CL(proxy, conv1, matrix_size, kernelNum, kernel2, 5, 1);
    matrix_size = ((matrix_size - 5)/1 + 1) / 2;
    cout << "finished CL2" << endl;
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
    cout << "params for FCL sent" << endl;
    uint64_t * output = FCL(proxy, conv2, matrix_size, kernelNum, weights2, nodes_out);
    cout << "FCL done" << endl;

    for(uint64_t n = 0; n<nodes_out; n++){
        cout << "Node " << n << ": " << output[n] << endl;
    }


    proxy->SendBytes(CORE_END);
    proxy->PrintBytes();
    return 0;
}

