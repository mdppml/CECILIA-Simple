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
bool doMaxpool;
uint32_t divisor;

void resetParams(uint32_t mode, vector<int> params);
void initParams(uint32_t mode);

int main(int argc, char* argv[]) {
    if (argc != 4) {
        clog << "Error: The program requires exactly two arguments; the IP address of the helper and the port it is listening on." << endl;
        return 1;
    }
    string address(argv[1]);
    uint16_t port = strtol(argv[2], nullptr, 10);
    uint32_t nn_mode = atoi(argv[3]);

    auto *helper = new Party(HELPER,port,address);

    // parse model parameters
    vector<int> dimensions;
    initParams(nn_mode);

    for (int image = 0; image < i_number; ++image) {
        resetParams(nn_mode, dimensions);

        // CNN INFERENCE PIPELINE
        // PERFORMING CONVOLUTION
        CL(helper, nullptr, i_channel, i_height, i_width, nullptr, k_dim, k_number, stride, doMaxpool, nullptr);
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
                k_number = 1;
                nodes_out = 100;
                // FULLY CONNECTED LAYER
                FCL(helper, nullptr, nodes_in, nullptr, nodes_out, nullptr);

                nodes_in = nodes_out;
                i_channel = 1;
                nodes_out = 10;
                break;
            }
            case 1: {
                // PERFORMING CONVOLUTION
                CL(helper, nullptr, i_channel, i_height, i_width, nullptr, k_dim, k_number, stride, doMaxpool, nullptr);
                i_channel = k_number;
                i_height = ((i_height - k_dim) + 1) / (stride * 2);
                i_width = ((i_width - k_dim) + 1) / (stride * 2);

                // fully connected layer:
                nodes_out = 100;
                nodes_in = i_height * i_width * i_channel;
                break;
            }
            default: {
                // random network: 2 more convolutions, then fully connected layer:
                for (uint32_t c = 0; c < 2; c++) {
                    k_number -= (c + 1);
                    k_dim -= 1;
                    // PERFORMING CONVOLUTION
                    CL(helper, nullptr, i_channel, i_height, i_width, nullptr, k_dim, k_number, stride, doMaxpool,
                       nullptr);
                    i_channel = k_number;
                    i_height = ((i_height - k_dim) + 1) / (stride * 2);
                    i_width = ((i_width - k_dim) + 1) / (stride * 2);
                    cout << "finished CL" << c + 2 << " (random mode)" << endl;
                }
                // fully connected layer:
                nodes_out = 10;
                nodes_in = i_height * i_width * i_channel;
                break;
            }
        }
        // FULLY CONNECTED LAYER
        FCL(helper, nullptr, nodes_in, nullptr, nodes_out, nullptr);
        MAX(helper, nullptr, nodes_out);
    }
    helper->PrintBytes();
    return 0;
}

void initParams(uint32_t mode) {
    i_number = 10, i_channel = 1, i_width = 28, i_height = 28;
    k_number = 5, k_dim = 5;
    stride = 1;
    padding = 0;
    doMaxpool = true;
    switch (mode) {
        case 0:{
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
    k_number = 5, k_dim = 5;
    stride = 1;

    switch (mode) {
        case 0:{
            i_height += 2 * padding;
            i_width += 2 * padding;
            stride = 2;
            break;
        }
        case 1:{
            uint32_t input_param_start = params.at(0)+1;
            i_width = params.at(input_param_start + 2);
            i_height = params.at(input_param_start + 2);
            break;
        }
        default:{
            k_number = 16;
            break;
        }
    }
}
