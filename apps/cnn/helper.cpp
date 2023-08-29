//
// Created by debora on 08.07.22.
//
#include <iostream>
#include "../../core/core.h"
#include "../../core/cnn.h"
#include "model_parser.h"
#include "../../core/auc.h"

const string MODEL_DIR = "../../apps/cnn/model_files/"; // FIXME this is outdated
const string MINIONN_MODEL_FILE = "MiniONN.txt";

uint32_t i_number, i_channel, i_width, i_height;
uint32_t k_number, k_dim;
uint32_t stride, padding;
uint32_t max_win_width, max_win_height;
uint32_t h_divisor, w_divisor;
uint32_t curr_layer;
bool trained_for_MNIST;

uint32_t nodes_out;
uint32_t nodes_in;

/**
 * Initializes all parameters defining the structure of the used networks (Chameleon (mode 0), (LeNet mode 1), CellCnn (mode 2) but not yet usable, and potentially others.)
 * Those parameters include number of images for which inference shall be computed, dimensions of input images, dimensions of kernel and others.
 * @param mode the mode specifying the network for which the parameters are initialized.
 */
void initParams(uint32_t mode) {
    //MNIST data
    trained_for_MNIST = true;
    i_number = 10000, i_channel = 1, i_width = 28, i_height = 28;
    k_dim = 5, k_number = 20;
    stride = 1;
    padding = 0;
    max_win_width = 2, max_win_height = 2;
    switch (mode) {
        case 0:{ //CHAMELEON
            k_number = 5;
            stride = 2;
            padding = 2;
            max_win_width = 0;
            max_win_height = 0;
            // i_width and i_height are adjusted after padding is performed; all other parameters are not modified.
            break;
        }
        //nothing to change for LeNet
        case 2:{ //CellCNN
            cout << "init params for single cell network" << endl; //TODO
            // Multi-cell data
            trained_for_MNIST = false;
            k_number = 8;
            i_channel = 1; // number of samples
            i_width = 35;   // number of markers
            i_height = 2000;// number of cells
            k_dim = 35; // 1x35
            max_win_width = 1;
            max_win_height = 2000; // 2000x1
            break;
        }
    }
    w_divisor = stride;
    if(max_win_width > 0){
        w_divisor *= max_win_width;
    }
    h_divisor = stride;
    if(max_win_height > 0){
        h_divisor *= max_win_height;
    }
}
/**
 * Reset those network parameters that have been changed during inference computation, so that the state is equal to that after initParams() have been called.
 */
void resetParams() {
    i_channel = 1;
    if(trained_for_MNIST){
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
}
/**
 * Update parameters as preparation for a subsequent Fully Connected Layer
 * @param bias_dimensions vector defining the number of bias values needed per layer.
 * @param firstFCL is true if the subsequent Fully Connected Layer (FCL) is the first FCL within the according network.
 */
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
/**
 * Update parameters after Convolutional Layer.
 */
void updateParamsAfterCL(){
    if(trained_for_MNIST) {
        i_height = (i_height - k_dim + 1) / h_divisor;
        i_width = (i_width - k_dim + 1) / w_divisor;
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
        cout << "INFERENCE " << image << endl;
        // CNN INFERENCE PIPELINE
        switch (nn_mode) {
            case 0: {
                cout << "call CL for CHAMELEON" << endl;
                k_number = 5;
                ConvolutionalLayer(helper, nullptr, i_channel, i_height, i_width, nullptr, k_dim, k_number, stride, max_win_height, max_win_width,
                          nullptr, true);
                updateParamsAfterCL();
                updateParamsForFCL(100, true);
                // FULLY CONNECTED LAYER
                FullyConnectedLayer(helper, nullptr, nodes_in, nullptr, nodes_out, nullptr);
                ReLU(helper, nullptr, nodes_out);
                updateParamsForFCL(10, false);
                break;
            }
            case 1: {
                k_number = 20;
                cout << "call CL for LeNet" << endl;
                ConvolutionalLayer(helper, nullptr, i_channel, i_height, i_width, nullptr, k_dim, k_number, stride, max_win_height, max_win_width,
                   nullptr, false);
                updateParamsAfterCL();

                k_number = 50;
                ConvolutionalLayer(helper, nullptr, i_channel, i_height, i_width, nullptr, k_dim, k_number, stride, max_win_height, max_win_width,
                   nullptr, true);
                updateParamsAfterCL();

                // fully connected layer:
                updateParamsForFCL(500, true);
                FullyConnectedLayer(helper, nullptr, nodes_in, nullptr, nodes_out, nullptr);
                ReLU(helper, nullptr, nodes_out);

                updateParamsForFCL(10, false);
                break;
            }
            case 2:{
                k_number = 8;
                for (int k = 0; k < k_number; ++k) {
                    MatrixVectorMultiply(helper, nullptr, nullptr, 0, i_height*i_width, 0);
                }
                updateParamsAfterCL();
                updateParamsForFCL(2, true);
                break;
            }
            default: {
                cout << "no neural network mode is matching the supported one." << endl;
                return -1;
            }
        }
        // FULLY CONNECTED LAYER
        FullyConnectedLayer(helper, nullptr, nodes_in, nullptr, nodes_out, nullptr);
        ArgMax(helper, nullptr, nodes_out);
    }
    helper->PrintBytes();
    return 0;
}
