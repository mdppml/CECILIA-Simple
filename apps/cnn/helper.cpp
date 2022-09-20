//
// Created by debora on 08.07.22.
//
#include <iostream>
#include "../../core/core.h"
#include "../../core/cnn.h"
#include "model_parser.h"
#include "../../core/auc.h"

const string MODEL_DIR = "../../apps/cnn/model_files/";
const string MINIONN_MODEL_FILE = "MiniONN.txt";

uint32_t i_number, i_channel, i_width, i_height;
uint32_t k_number, k_dim;
uint32_t stride, padding;
uint32_t maxpool_window_dim;
uint32_t divisor;
uint32_t curr_layer;
bool trained_for_MNIST;

uint32_t nodes_out;
uint32_t nodes_in;

void initParams(uint32_t mode) {
    //MNIST data
    trained_for_MNIST = true;
    i_number = 10, i_channel = 1, i_width = 28, i_height = 28;
    k_dim = 5, k_number = 20;
    stride = 1;
    padding = 0;
    maxpool_window_dim = 2;
    switch (mode) {
        case 0:{ //CHAMELEON
            k_number = 5;
            stride = 2;
            padding = 2;
            maxpool_window_dim = 0;
            // i_width and i_height are adjusted after padding is performed; all other parameters are not modified.
            break;
        }
        //nothing to change for LeNet
        case 2:{ //CellCNN
            cout << "init params for single cell network" << endl; //TODO
            // Multi-cell data
            trained_for_MNIST = false;
            k_number = 8;
            i_channel = 64; // number of samples
            i_width = 35;   // number of markers
            i_height = 2000;// number of cells
            k_dim = 35; // 1x35
            maxpool_window_dim = 2000; // 2000x1
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
    if(trained_for_MNIST){
        i_channel = 1;
        i_height = 28;
        i_width =  28;
    }
    else{
        //i_channel = 64; // number of samples
        i_width = 35;   // number of markers
        i_height = 2000;// number of cells
    }
    if(padding > 0){
        i_height += 2*padding;
        i_width += 2*padding;
    }
    if(padding > 0) {
        i_height += 2 * padding;
        i_width += 2 * padding;
    }
}

void updateParamsForFCL(uint32_t n_out, bool firstFCL){
    if(firstFCL){
        nodes_in = i_height*i_width*i_channel;
    }
    else{
        nodes_in = nodes_out;
        curr_layer++; //if firstFCL: curr_layer has been updated by updateParamsAfterCL
    }
    nodes_out = n_out;
}

void updateParamsAfterCL(){
    if(trained_for_MNIST) {
        i_height = (i_height - k_dim + 1) / divisor;
        i_width = (i_width - k_dim + 1) / divisor;
    }
    else{
        i_height = i_channel;
        i_width = k_number; //TODO ??
    }
    i_channel = k_number;
    curr_layer++;
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
        resetParams();
        // CNN INFERENCE PIPELINE
        switch (nn_mode) {
            case 0: {
                CL(helper, nullptr, i_channel, i_height, i_width, nullptr, k_dim, k_number, stride, maxpool_window_dim,
                          nullptr, true);
                updateParamsAfterCL();
                updateParamsForFCL(100, true);
                // FULLY CONNECTED LAYER
                FCL(helper, nullptr, nodes_in, nullptr, nodes_out, nullptr);
                updateParamsForFCL(10, false);
                break;
            }
            case 1: {
                cout << i_channel << " " << i_height << " " << i_width << " (c x h x w) " << k_dim << " x " << k_number << " s= " << stride << " " << maxpool_window_dim << endl;
                CL(helper, nullptr, i_channel, i_height, i_width, nullptr, k_dim, k_number, stride, maxpool_window_dim,
                   nullptr, false);
                updateParamsAfterCL();
                k_number = 50;
                cout << i_channel << " " << i_height << " " << i_width << " (c x h x w) " << k_dim << " x " << k_number << " s= " << stride << " " << maxpool_window_dim << endl;
                CL(helper, nullptr, i_channel, i_height, i_width, nullptr, k_dim, k_number, stride, maxpool_window_dim,
                   nullptr, true);
                updateParamsAfterCL();

                // fully connected layer:
                updateParamsForFCL(500, true);
                FCL(helper, nullptr, nodes_in, nullptr, nodes_out, nullptr);

                updateParamsForFCL(10, false);
                break;
            }
            case 2:{
                for (int k = 0; k < k_number; ++k) {
                    MATVECMUL(helper, nullptr, nullptr, 1, i_height, i_width);
                }
                //TODO max with asymmetric size
                updateParamsAfterCL();
                updateParamsForFCL(2, true);
            }
            default: {
                cout << "no neural network mode is matching the supported one." << endl;
                return -1;
            }
        }
        // FULLY CONNECTED LAYER
        FCL(helper, nullptr, nodes_in, nullptr, nodes_out, nullptr);
        switch (nn_mode) {                              // from here on network architectures differ
            case 0:{ // Chameleon
                ARGMAX(helper, nullptr, nodes_out);
                break;
            }
            case 1: { // SecureNN
                RELU(helper, nullptr, nodes_out);
                //DIVISION(helper, nullptr, nullptr); // vectorized DIV
                for (int i = 0; i < nodes_out; ++i) {
                    DIVISION(helper, 0, 0);
                }
                //TODO if DIV by sum... this is factor --> only perform argmax...
                convert2double(REC(helper, ARGMAX(helper, nullptr, nodes_out)));
                //print1DArray("Input to ARGMAX:", convert2double(REC(helper, output, nodes_out), nodes_out), nodes_out);
                break;
            }
            case 2:{

            }
        }
    }
    helper->PrintBytes();
    return 0;
}
