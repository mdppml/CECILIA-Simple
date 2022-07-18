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
const string MINIONN_MODEL_FILE = "MiniONN.txt", CHAMELEON_MODEL_FILES = "Chameleon_CNN/";
const string MNIST_DATA = "../../apps/cnn/mnist/";

uint64_t layer_number;
uint32_t i_number, i_channel, i_width, i_height;
uint32_t k_number, k_dim;
uint32_t stride, padding;
uint32_t divisor;
bool doMaxpool;
uint32_t curr_layer;

void resetParams(uint32_t mode, vector<int> params);
void initParams(uint32_t mode);

int main(int argc, char* argv[]) {
    if (argc < 6){
        cout << "Calling proxy without specifying role (1), port (2), address (3), helpers port (4) and helpers adress (5) is not possible." << endl;
        cout << "Specify the Network Mode with the 6th parameter: 0 = CHAMELEON, 1 = MINIONN" << endl;
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
    vector<int> dimensions;
    initParams(nn_mode);
    double**** model_weights = new double ***[layer_number];
    uint32_t* bias_dimensions = new uint32_t [layer_number];
    switch (nn_mode) {
        case 0:{
            cout << "CHAMELEON" << endl;
            model_weights = getChameleonParameters(MODEL_DIR + CHAMELEON_MODEL_FILES, 5, i_channel);
            bias_dimensions[0] = 5;
            bias_dimensions[1] = 100;
            bias_dimensions[2] = 10;

            break;
        }
        case 1:{
            cout << "MINIONN" << endl;
            model_weights[0] = getONNXParameters(MODEL_DIR + MINIONN_MODEL_FILE, dimensions, layer_number);
            uint64_t * weight_dims = new uint64_t [dimensions.at(0)];
            for (uint64_t i = 0; i < dimensions.at(0); i++){
                int w = dimensions[i + 1];
                weight_dims[i] = convert2uint64(w);
            }
            //TODO weights are not clear
            //TODO bias are unknown
            print1DArray("weight dimensions", convert2double(weight_dims, dimensions.at(0)), dimensions.at(0));
            //print2DArray("model weights - layer 2: ", model_weights[0], dimensions.at(1), dimensions.at(2));
            print2DArray("model weights - layer 3: ", model_weights[0][1], dimensions.at(3), 1, 0);
            resetParams(nn_mode, dimensions);
            break;
        }
        default:{
            cout << "No or unknown network has been specified. Use random weights to perform CNN inference for 3 convolutional layer and one fully connected layer" << endl;
            // fill weights for convolutional layer
            for(uint32_t layer = 0; layer<3; layer++){
                model_weights[layer] = new double **[k_number];
                for(uint64_t kernel = 0; kernel< k_number; kernel++){
                    model_weights[layer][kernel] = random_2D_data(proxy, k_dim-layer, k_dim-layer, -5.0, 5.0);
                }
                k_number -= (layer + 1);
                // layer 0 has 5 kernel: 5x5    layer 1 has 4 kernel: 4x4 and   layer 2 has 2 kernel: 3x3
            }
            uint32_t in_fcl = (((((i_height - k_dim + 1) / divisor) - k_dim) / divisor) - k_dim - 1) / divisor;
            uint32_t out_fcl = 10;
            // fill weights for fully connected layer
            model_weights[3] = new double **[1];
            model_weights[3][0] = random_2D_data(proxy, in_fcl, out_fcl, -5.0, 5.0);
            break;
        }
    }
    //read in MNIST data
    mnist::MNIST_dataset<vector, vector<uint8_t>, uint8_t> dataset =
            mnist::read_dataset<vector, vector, uint8_t, uint8_t>(MNIST_DATA);

    vector<vector<unsigned char>> training = dataset.test_images;
    uint64_t*** data = new uint64_t** [training.size()];
    vector<unsigned char> label = dataset.test_labels;
    for (uint32_t i = 0; i < training.size(); ++i) {
        data[i] = new uint64_t *[i_height];
        for (uint32_t r = 0; r < i_height; ++r) {
            data[i][r] = new uint64_t [i_width];
            for (uint32_t c = 0; c < i_width; ++c) {
                double pixelValue = training.at(i).at(r*i_width + c);
                data[i][r][c] = proxy->createShare(pixelValue);                 // store directly as secret shares
            }
        }
    }
    // secret shares of model parameter
    uint64_t **bias = new uint64_t *[layer_number];
    for (int layer = 0; layer < layer_number; ++layer) {
        bias[layer] = proxy->createShare(model_weights[layer + 3][0][0], bias_dimensions[layer]);
        //print1DArray("BIAS: ", convert2double(REC(proxy, bias[layer], bias_dimensions[layer]), bias_dimensions[layer]), bias_dimensions[layer]);
    }
    uint64_t ***kernel = new uint64_t**[k_number]; // weights of all kernel for first CL
    for (uint32_t i = 0; i < k_number; i++) {
        kernel[i] = proxy->createShare(model_weights[0][i], i_channel, k_dim*k_dim);
    }

    // CNN INFERENCE PIPELINE
    for (uint32_t image = 0; image < i_number; ++image) {
        resetParams(nn_mode, dimensions);
        uint64_t *** input = new uint64_t **[i_channel];
        if(padding > 0){
            for(uint32_t i = 0; i<i_channel; i++){
                input[i] = PAD(data[image], i_height, i_width, 0, padding);
            }
            i_height += 2*padding;
            i_width += 2*padding;
        }

        // PERFORMING CONVOLUTION
        uint64_t*** conv = CL(proxy, input, i_channel, i_height, i_width, kernel, k_dim, k_number, stride, doMaxpool, bias[curr_layer]);
        i_channel = k_number;
        //          size after conv
        i_height = ((i_height - k_dim) + 1) / divisor;
        i_width = ((i_width - k_dim) + 1) / divisor;
        curr_layer++;

        //init what the architectures share for FCL:
        uint32_t nodes_out;
        uint32_t nodes_in = i_height*i_width*i_channel;

        uint64_t* prev_layer_res;
        uint64_t **weights;

        switch (nn_mode) {                              // from here on network architectures differ
            case 0:{
                k_number = 1;
                nodes_out = 100;
                // FULLY CONNECTED LAYER
                weights = proxy->createShare(model_weights[curr_layer][0], nodes_out, nodes_in);
                uint64_t *flattened = FLT(conv, i_height, i_width, i_channel);
                prev_layer_res = FCL(proxy, flattened, nodes_in, weights, nodes_out, bias[curr_layer]);

                curr_layer++;
                nodes_in = nodes_out;
                i_channel = 1;
                nodes_out = 10;
                break;
            }
            case 1:{
                // PERFORMING CONVOLUTION
                kernel = new uint64_t**[k_number];
                for (uint32_t i = 0; i < k_number; i++) {
                    kernel[i] = proxy->createShare(model_weights[curr_layer][i], 1, k_dim*k_dim);
                }
                //TODO might be problematic if conv changes size...
                conv = CL(proxy, conv, i_channel, i_height, i_width, kernel, k_dim, k_number, stride, doMaxpool, bias[curr_layer]);
                i_channel = k_number;
                i_height = ((i_height - k_dim) + 1) / (stride*2);
                i_width = ((i_width - k_dim) + 1) / (stride*2);
                cout << "finished CL2 (MiniONN)" << endl;

                curr_layer++;
                // fully connected layer:
                nodes_out = 100;
                nodes_in = i_height * i_width * i_channel;
                weights = proxy->createShare(model_weights[curr_layer][0], nodes_out, nodes_in);
                prev_layer_res = FLT(conv, i_height, i_width, i_channel);
                break;
            }
            default:{
                // random network: 2 more convolutions, then fully connected layer:
                for (uint32_t c = 0; c<2; c++){
                    k_number -= (c+1);
                    k_dim -= 1;
                    // PERFORMING CONVOLUTION
                    kernel = new uint64_t**[k_number];
                    for (uint32_t i = 0; i < k_number; i++) {
                        kernel[i] = proxy->createShare(model_weights[curr_layer][i], i_channel, k_dim*k_dim);
                    }
                    //TODO might be problematic if conv changes size...
                    conv = CL(proxy, conv, i_channel, i_height, i_width, kernel, k_dim, k_number, stride, doMaxpool, bias[curr_layer]);
                    i_channel = k_number;
                    i_height = ((i_height - k_dim) + 1) / (stride*2);
                    i_width = ((i_width - k_dim) + 1) / (stride*2);
                    curr_layer++;
                    cout << "finished CL" << curr_layer << " (random mode)"<< endl;
                }
                // fully connected layer:
                nodes_out = 10;
                nodes_in = i_height * i_width * i_channel;
                weights = proxy->createShare(model_weights[curr_layer][0], nodes_out, nodes_in);
                prev_layer_res = FLT(conv, i_height, i_width, i_channel);
                break;
            }
        }
        // FULLY CONNECTED LAYER
        uint64_t * output = FCL(proxy, prev_layer_res, nodes_in, weights, nodes_out, bias[curr_layer]);
        uint64_t maxNode = MAX(proxy, output, nodes_out);
        //TODO if MINIONN: do ASM(output[i]) = RELU(output[i])/RELU(output, nodes_out))
        cout << "Recognized number " << convert2double(REC(proxy, maxNode)) << ": " << endl;
        print1DArray("output", convert2double(REC(proxy, output, nodes_out), nodes_out), nodes_out);
        cout << "LABEL = " << to_string(label[image]) << endl;
    }

    proxy->PrintBytes();
    return 0;
}

void initParams(uint32_t mode) {
    layer_number = 4;
    curr_layer = 0;
    i_number = 10, i_channel = 1, i_width = 28, i_height = 28;
    k_number = 5, k_dim = 5;
    stride = 1;
    padding = 0;
    doMaxpool = true;
    switch (mode) {
        case 0:{
            layer_number = 3;
            stride = 2;
            padding = 2;
            doMaxpool = false;
            // i_width and i_height are adjusted after padding is performed; all other parameters are not modified.
            break;
        }
        default:{
            k_number = 16;
            break;
        }
    }
    divisor = stride;
    if(doMaxpool){
        divisor *= 2;
    }
}

void resetParams(uint32_t mode, vector<int> params) {
    i_channel = 1, i_width = 28, i_height = 28;
    k_number = 5;
    curr_layer = 0;
    switch (mode) {
        case 0: {
            // i_width and i_height are adjusted after padding; all other parameters are not modified.
            break;
        }
        case 1: {
            uint32_t input_param_start = params.at(0) + 1;
            i_width = params.at(input_param_start + 2);
            i_height = params.at(input_param_start + 2);
            break;
        }
        default: {
            k_number = 16;
            break;
        }
    }
}

