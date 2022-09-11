//
// Created by debora on 08.07.22.
//
#include <iostream>
#include "../../core/core.h"
#include "../../core/cnn.h"
#include "model_parser.h"

const string MODEL_DIR = "../../apps/cnn/model_files/";
const string MINIONN_MODEL_FILE = "MiniONN.txt";

uint32_t i_number, i_channel, i_width, i_height;
uint32_t k_number, k_dim;
uint32_t stride, padding;
uint32_t maxpool_window_dim;
uint32_t divisor;

void initParams(uint32_t mode) {
    i_number = 1, i_channel = 1, i_width = 28, i_height = 28;
    k_dim = 5;
    stride = 1;
    padding = 0;
    maxpool_window_dim = 2;
    switch (mode) {
        case 0:{
            stride = 2;
            padding = 2;
            maxpool_window_dim = 0;
            // i_width and i_height are adjusted after padding is performed; all other parameters are not modified.
            break;
        }
        case 2:{
            k_number = 20;
        }
        case 3:{
            //TODO single cell
            break;
        }
        // for all other cases nothing to change
    }
    divisor = stride;
    if(maxpool_window_dim > 0){
        divisor *= 2;
    }
}

void resetParams(uint32_t mode) {
    i_channel = 1;
    switch (mode) {
        case 0:{
            k_number = 5;
            break;
        }
        case 1:{
            k_number = 20;
            break;
        }
        case 2:{
            //TODO single cell
            break;
        }
    }
    if(padding > 0){
        i_height = 28 + 2*padding;
        i_width =  28 + 2*padding;
    }
}

int main(int argc, char* argv[]) {
    if (argc != 4) {
        clog << "Error: The program requires exactly three arguments; the IP address of the helper, the port it is listening on and the network mode." << endl;
        return 1;
    }
    string address(argv[1]);
    uint16_t port = strtol(argv[2], nullptr, 10);
    uint32_t nn_mode = atoi(argv[3]);

    auto *helper = new Party(HELPER,port,address);

    // parse model parameters
    initParams(nn_mode);

    for (int image = 0; image < i_number; ++image) {
        resetParams(nn_mode);

        // CNN INFERENCE PIPELINE
        // PERFORMING CONVOLUTION
        CL(helper, nullptr, i_channel, i_height, i_width, nullptr, k_dim, k_number, stride, maxpool_window_dim, nullptr);
        i_channel = k_number;
        //          size after conv         after maxpool divide by 2
        i_height = ((i_height - k_dim) + 1) / divisor;
        i_width = ((i_width - k_dim) + 1) / divisor;
        //init what the architectures share for FCL:
        uint32_t nodes_out;
        uint32_t nodes_in = i_height * i_width * i_channel;

        // from here on network architectures differ
        switch (nn_mode) {
            case 0: {
                nodes_out = 100;
                // FULLY CONNECTED LAYER
                FCL(helper, nullptr, nodes_in, nullptr, nodes_out, nullptr);
                nodes_in = nodes_out;
                break;
            }
            case 1: {
                // PERFORMING CONVOLUTION
                k_number = 50;
                cout << "k: " << k_number << ", c: " << i_channel << endl;
                cout << "w: " << i_width << ", h: " << i_height << endl;
                CL(helper, nullptr, i_channel, i_height, i_width, nullptr, k_dim, k_number, stride, maxpool_window_dim,
                   nullptr);
                i_channel = k_number;
                i_height = ((i_height - k_dim) + 1) / divisor;
                i_width = ((i_width - k_dim) + 1) / divisor;

                cout << "k: " << k_number << ", c: " << i_channel << endl;
                cout << "w: " << i_width << ", h: " << i_height << endl;
                // fully connected layer:
                nodes_out = 500;
                nodes_in = i_height * i_width * i_channel;
                cout << "call FCL.." << endl;
                FCL(helper, nullptr, nodes_in, nullptr, nodes_out, nullptr);
                cout << "finished FCL1 (LeNet_NN)" << endl;

                nodes_in = nodes_out;
                break;
            }
            case 2:{
                //TODO single cell
            }
            default: {
                cout << "no neural network mode is matching the supported one." << endl;
                return -1;
            }
        }
        nodes_out = 10;
        // FULLY CONNECTED LAYER
        FCL(helper, nullptr, nodes_in, nullptr, nodes_out, nullptr);
        MAX(helper, nullptr, nodes_out);
    }
    helper->PrintBytes();
    return 0;
}
