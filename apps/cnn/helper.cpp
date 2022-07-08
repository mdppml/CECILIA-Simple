//
// Created by debora on 08.07.22.
//
#include <iostream>
#include "../../core/core.h"
#include "../../core/cnn.h"
#include "model_parser.h"

const string MODEL_DIR = "../../apps/cnn/model_files/";
const string MINIONN_MODEL_FILE = "MiniONN.txt";

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
    vector<int> d;
    // init input parameter --> all supported networks expect the same input shape.
    uint32_t i_number = 10;
    uint32_t i_channel = 1;
    uint32_t i_width = 28;
    uint32_t i_height = 28;
    uint32_t k_number = 5;
    uint32_t k_dim = 5;
    uint32_t stride = 1;
    uint32_t padding = 0;
    bool doMaxpool = true;

    uint32_t divisor = stride;
    if(doMaxpool){
        divisor *= 2;
    }
    switch (nn_mode) {
        case 0:{
            stride = 2;
            divisor = stride;
            doMaxpool = false;
            padding = 2;
            // this is a well defined CNN architecture trained for the MNIST dataset
            break;
        }
        case 1:{
            cout << "MINIONN" << endl;
            getONNXParameters(MODEL_DIR + MINIONN_MODEL_FILE, d, 4);
            uint32_t input_param_start = d.at(0)+1;
            i_number = d.at(input_param_start);
            i_channel = d.at(input_param_start + 1);
            i_width = d.at(input_param_start + 2);
            i_height = d.at(input_param_start + 2);
            k_number = 5; //TODO
            k_dim = 5;
            break;
        }
        default:{
            cout << "No or unknown network has been specified. Use random weights to perform CNN inference for 3 convolutional layer and one fully connected layer" << endl;
            k_number = 16;
            uint32_t in_fcl = (((((i_height - k_dim + 1) / divisor) - k_dim) / divisor) - k_dim - 1) / divisor;
            uint32_t out_fcl = 10;
            break;
        }
    }
    //create shares
    if(padding > 0){
        i_height += 2*padding;
        i_width += 2*padding;
    }

    // CNN INFERENCE PIPELINE

    // PERFORMING CONVOLUTION
    CL(helper, nullptr, i_channel, i_height, i_width, nullptr, k_dim, k_number, stride, doMaxpool);
    i_channel = k_number;
    //          size after conv         after maxpool divide by 2
    i_height = ((i_height - k_dim) + 1) / divisor;
    i_width = ((i_width - k_dim) + 1) / divisor;
    //init what the architectures share for FCL:
    uint32_t nodes_out;
    uint32_t nodes_in = i_height*i_width*i_channel;

    // from here on network architectures differ
    switch (nn_mode) {
        case 0:{
            k_number = 1;
            nodes_out = 100;
            // FULLY CONNECTED LAYER
            FCL(helper, nullptr, nodes_in, nullptr, nodes_out);

            nodes_in = nodes_out;
            i_channel = 1;
            nodes_out = 10;
            break;
        }
        case 1:{
            // PERFORMING CONVOLUTION
            CL(helper, nullptr, i_channel, i_height, i_width, nullptr, k_dim, k_number, stride, doMaxpool);
            i_channel = k_number;
            i_height = ((i_height - k_dim) + 1) / (stride*2);
            i_width = ((i_width - k_dim) + 1) / (stride*2);

            // fully connected layer:
            nodes_out = 100;
            nodes_in = i_height * i_width * i_channel;
            break;
        }
        default:{
            // random network: 2 more convolutions, then fully connected layer:
            for (uint32_t c = 0; c<2; c++){
                k_number -= (c+1);
                k_dim -= 1;
                // PERFORMING CONVOLUTION
                CL(helper, nullptr, i_channel, i_height, i_width, nullptr, k_dim, k_number, stride, doMaxpool);
                i_channel = k_number;
                i_height = ((i_height - k_dim) + 1) / (stride*2);
                i_width = ((i_width - k_dim) + 1) / (stride*2);
                cout << "finished CL" << c+2 << " (random mode)"<< endl;
            }
            // fully connected layer:
            nodes_out = 10;
            nodes_in = i_height * i_width * i_channel;
            break;
        }
    }
    // FULLY CONNECTED LAYER
    FCL(helper, nullptr, nodes_in, nullptr, nodes_out);

    helper->PrintBytes();
    return 0;
}