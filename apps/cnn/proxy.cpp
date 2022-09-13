#include <cstdlib>
#include <iostream>
#include <deque>
#include <chrono>
#include <assert.h>
#include "../../core/cnn.h"
#include "model_parser.h"
#include "mnist/mnist_reader.hpp"

using namespace std;
const string MODEL_DIR = "../../apps/cnn/model_files/";
const string CHAMELEON_MODEL_FILES = "Chameleon_CNN/", LENET_NN_MODEL_FILES = "LeNet_trained/", MINIONN_MODEL_FILE = "MiniONN.txt";
const string MNIST_DATA = "../../apps/cnn/mnist/";

uint64_t layer_number;
uint32_t i_number, i_channel, i_width, i_height;
uint32_t k_number, k_dim;
uint32_t stride, padding;
uint32_t divisor;
uint32_t maxpool_window_dim;
uint32_t curr_layer;

void initParams(uint32_t mode) {
    curr_layer = 0;
    i_number = 1, i_channel = 1, i_width = 28, i_height = 28;
    k_dim = 5;
    stride = 1;
    padding = 0;
    maxpool_window_dim = 2;
    switch (mode) {
        case 0:{
            layer_number = 3;
            stride = 2;
            padding = 2;
            maxpool_window_dim = 0;
            // i_width and i_height are adjusted after padding is performed; all other parameters are not modified.
            break;
        }
        case 1:{
            layer_number = 4;
            break;
        }
        case 2:{
            cout << "TODO single cell network" << endl; //TODO
            break;
        }
        default: {
            layer_number = 4;
            break;
        }
    }
    divisor = stride;
    if(maxpool_window_dim > 0){
        divisor *= maxpool_window_dim;
    }
}

void resetParams() {
    i_channel = 1;
    curr_layer = 0;

    if(padding > 0){
        i_height = 28 + 2*padding;
        i_width =  28 + 2*padding;
    }
}

int main(int argc, char* argv[]) {
    if (argc < 6){
        cout << "Calling proxy without specifying role (1), port (2), address (3), helpers port (4) and helpers adress (5) is not possible." << endl;
        cout << "Specify the Network Mode with the 6th parameter: 0 = CHAMELEON, 1 = LeNet5" << endl;
        return 1;
    }
    uint8_t role = atoi(argv[1]);
    uint16_t cport = atoi(argv[2]);
    string caddress(argv[3]);
    uint16_t hport = atoi(argv[4]);
    string haddress(argv[5]);
    uint32_t nn_mode = atoi(argv[6]);

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

    // parse model parameters
    initParams(nn_mode);
    double**** model_weights = new double ***[layer_number];
    uint32_t* bias_dimensions = new uint32_t [layer_number];
    switch (nn_mode) {
        case 0:{
            cout << "CHAMELEON" << endl;
            model_weights = getChameleonParameters(MODEL_DIR + CHAMELEON_MODEL_FILES, 5);
            bias_dimensions[0] = 5;
            bias_dimensions[1] = 100;
            bias_dimensions[2] = 10;
            break;
        }
        case 1:{
            cout << "LeNet by SecureNN" << endl;
            model_weights = getLeNetParameters(MODEL_DIR + LENET_NN_MODEL_FILES, true);
            bias_dimensions[0] = 20;
            bias_dimensions[1] = 50;
            //bias_dimensions[2] = 800;
            bias_dimensions[2] = 800;
            bias_dimensions[3] = 500;
            break;
        }
        case 2:{
            //TODO Single cell network
        }
        default: {
            cout << "No neural network mode is matching the supported one." << endl;
            return -1;
        }
    }
    //read in MNIST data
    cout << "Reading input data... ";
    mnist::MNIST_dataset<vector, vector<uint8_t>, uint8_t> dataset =
            mnist::read_dataset<vector, vector, uint8_t, uint8_t>(MNIST_DATA);
    vector<vector<unsigned char>> training = dataset.test_images;
    vector<unsigned char> label = dataset.test_labels;

    uint64_t*** data = new uint64_t** [training.size()];
    for (uint32_t i = 0; i < i_number; ++i) {
        data[i] = new uint64_t *[i_height];
        for (uint32_t r = 0; r < i_height; ++r) {
            data[i][r] = new uint64_t [i_width];
            for (uint32_t c = 0; c < i_width; ++c) {
                double pixelValue = training.at(i).at(r*i_width + c);
                data[i][r][c] = proxy->createShare(pixelValue);                 // store directly as secret shares
            }
        }
        if(padding > 0){
            data[i] = PAD(data[i], i_height, i_width, 0, padding);
        }
    }
    cout << "Creating secret shares: " << endl;
    // secret shares of model parameter
    uint64_t **bias = new uint64_t *[layer_number];
    for (int layer = 0; layer < layer_number; ++layer) {
        cout << " - Bias for layer " << layer << " of length " << bias_dimensions[layer] << endl;
        bias[layer] = proxy->createShare(model_weights[layer + layer_number][0][0], bias_dimensions[layer]);
        //print1DArray("BIAS: ", convert2double(REC(proxy, bias[layer], bias_dimensions[layer]), bias_dimensions[layer]), bias_dimensions[layer]);
        delete[] model_weights[layer + layer_number][0][0];
        delete[] model_weights[layer + layer_number][0];
        delete[] model_weights[layer + layer_number];
    }
    k_number = bias_dimensions[0];
    uint64_t ***kernel = new uint64_t**[k_number]; // weights of all kernel for first CL
    for (uint32_t i = 0; i < k_number; i++) {
        //cout << "Kernel " << i << endl;
        kernel[i] = proxy->createShare(model_weights[0][i], i_channel, k_dim*k_dim);
        //print2DArray("KERNEL: ", convert2double(REC(proxy, kernel[i], i_channel, k_dim*k_dim), i_channel, k_dim*k_dim), i_channel, k_dim*k_dim);
        delete[] model_weights[0][i][0]; // if more than one channel: delete for each channel
        delete[] model_weights[0][i];
    }
    delete[] model_weights[0];

    double* prediction = new double [i_number];
    uint32_t correct = 0, incorrect = 0;
    // CNN INFERENCE PIPELINE
    for (uint32_t image = 0; image < i_number; ++image) {
        cout << "START INFERENCE PIPELINE" << endl;
        k_number = bias_dimensions[0];
        resetParams();
        auto *** input = new uint64_t **[i_channel];
        input[0] = data[image]; // currently only 1 channel for input supported

        // PERFORMING CONVOLUTION
        //print2DArray("input", convert2double(REC(proxy, input[0], i_height, i_width), i_height, i_width), i_height, i_width);
        //print2DArray("KERNEL: ", convert2double(REC(proxy, kernel[0], i_channel, k_dim*k_dim), i_channel, k_dim*k_dim), i_channel, k_dim*k_dim);
        //print1DArray("BIAS: ", convert2double(REC(proxy, bias[curr_layer], bias_dimensions[curr_layer]), bias_dimensions[curr_layer]), bias_dimensions[curr_layer]);
        uint64_t*** conv = CL(proxy, input, i_channel, i_height, i_width, kernel, k_dim, k_number, stride, maxpool_window_dim, bias[curr_layer]);
        i_channel = k_number;
        //          size after conv
        i_height = ((i_height - k_dim) + 1) / divisor;
        i_width = ((i_width - k_dim) + 1) / divisor;
        curr_layer++;

        //init what the architectures share for FCL:
        uint32_t nodes_out;
        uint32_t nodes_in;

        uint64_t* prev_layer_res;
        uint64_t **weights;

        switch (nn_mode) {                              // from here on network architectures differ
            case 0:{ // Chameleon
                nodes_in = i_height*i_width*i_channel;
                nodes_out = bias_dimensions[curr_layer];
                // FULLY CONNECTED LAYER
                weights = proxy->createShare(model_weights[curr_layer][0], nodes_out, nodes_in);
                if(image == i_number-1) {
                    delete[] model_weights[curr_layer][0];
                    delete[] model_weights[curr_layer];
                }
                uint64_t *flattened = FLT(conv, i_height, i_width, i_channel);
                prev_layer_res = FCL(proxy, flattened, nodes_in, weights, nodes_out, bias[curr_layer]);

                nodes_in = nodes_out;
                break;
            }
            case 1: { // SecureNN
                // PERFORMING CONVOLUTION
                k_number = bias_dimensions[curr_layer];
                i_channel =  bias_dimensions[curr_layer-1];
                cout << "k: " << k_number << ", c: " << i_channel << endl;
                cout << "w: " << i_width << ", h: " << i_height << endl;
                kernel = new uint64_t**[k_number];
                for (uint32_t i = 0; i < k_number; i++) {
                    kernel[i] = proxy->createShare(model_weights[curr_layer][i], i_channel, k_dim*k_dim);
                    if(image == i_number-1) {
                        delete[] model_weights[curr_layer][i];
                    }
                }
                if(image == i_number-1) {
                    delete[] model_weights[curr_layer];
                }
                cout << "call CL2 (LeNet_NN)" << endl;
                conv = CL(proxy, conv, i_channel, i_height, i_width, kernel, k_dim, k_number, stride, maxpool_window_dim, bias[curr_layer]);
                i_channel = k_number;
                i_height = ((i_height - k_dim) + 1) / divisor;
                i_width = ((i_width - k_dim) + 1) / divisor;
                cout << "finished CL2 (LeNet_NN)" << endl;
                curr_layer++;

                cout << "curr layer: " << curr_layer << ", w: " << i_width << ", h: " << i_height << ", i channel: " << i_channel << endl;
                // fully connected layer:
                nodes_out = bias_dimensions[curr_layer+1];
                nodes_in = i_height * i_width * i_channel; // equals bias_dimensions[curr_layer]
                cout << "create shares..." << nodes_in << " x " << nodes_out << endl;
                weights = proxy->createShare(model_weights[curr_layer][0], nodes_out, nodes_in);
                cout << "done: delete. " << endl;
                if(image == i_number-1) {
                    delete[] model_weights[curr_layer][0];
                    delete[] model_weights[curr_layer];
                }
                uint64_t * flattened = FLT(conv, i_height, i_width, i_channel);
                prev_layer_res = FCL(proxy, flattened, nodes_in, weights, nodes_out, bias[curr_layer]);
                cout << "finished FCL1 (LeNet_NN)" << endl;

                nodes_in = nodes_out;
                break;
            }
            case 2:{
                //TODO single cell
                nodes_in = i_height*i_width*i_channel;
            }
        }
        // last FULLY CONNECTED LAYER
        nodes_out = 10;
        curr_layer++;
        weights = proxy->createShare(model_weights[curr_layer][0], nodes_out, nodes_in);
        if(image == i_number-1) {
            delete[] model_weights[curr_layer][0];
            delete[] model_weights[curr_layer];
        }

        uint64_t * output = FCL(proxy, prev_layer_res, nodes_in, weights, nodes_out, bias[curr_layer]);
        switch (nn_mode) {                              // from here on network architectures differ
            case 0:{ // Chameleon
                prediction[image] = ARGMAX(proxy, output, nodes_out);
                if (prediction[image] == label[image]){
                    correct++;
                }
                else{
                    incorrect++;
                }
                break;
            }
            case 1: { // SecureNN
                //TODO ASM(output[i]) = RELU(output[i])/RELU(output, nodes_out))
                break;
            }
            case 2:{

            }
        }
    }
    cout << "accuracy: " << correct/i_number << endl;
    proxy->PrintBytes();
    return 0;
}


