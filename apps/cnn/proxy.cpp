#include <cstdlib>
#include <iostream>
#include <deque>
#include <chrono>
#include <assert.h>
#include "../../core/cnn.h"
#include "../../core/auc.h"
#include "model_parser.h"
#include "mnist_data/mnist_reader.hpp"

using namespace std;
const string MODEL_DIR = "../../apps/cnn/model_files/";
const string CHAMELEON_MODEL_FILES = "Chameleon_CNN/", LENET_NN_MODEL_FILES = "LeNet_trained/LeNetNN_b128_n15/",
CELL_CNN_MODEL_FILES = "Cell_CNN/", MINIONN_MODEL_FILE = "MiniONN.txt";
const string MNIST_DATA = "../../apps/cnn/mnist_data/";

uint64_t layer_number;
uint32_t i_number, i_channel, i_width, i_height;
uint32_t k_number, k_dim;
uint32_t stride, padding;
uint32_t h_divisor, w_divisor;
uint32_t max_win_width, max_win_height;
uint32_t curr_layer;
bool trained_for_MNIST;

uint32_t nodes_out;
uint32_t nodes_in;

void initParams(uint32_t mode) {
    curr_layer = 0;
    //MNIST data
    trained_for_MNIST = true;
    i_number = 10, i_channel = 1, i_width = 28, i_height = 28;
    k_dim = 5;
    stride = 1;
    padding = 0;
    max_win_width = 2, max_win_height = 2;
    switch (mode) {
        case 0:{ //CHAMELEON
            layer_number = 3;
            stride = 2;
            padding = 2;
            max_win_width = 0, max_win_height = 0;
            // i_width and i_height are adjusted after padding is performed; all other parameters are not modified.
            break;
        }
        case 1:{ //LeNet
            layer_number = 4;
            break;
        }
        case 2:{ //CellCNN
            cout << "init params for single cell network" << endl; //TODO
            layer_number = 2;
            // Multi-cell data
            trained_for_MNIST = false;
            //i_channel = 64; // number of samples
            i_width = 35;   // number of markers
            i_height = 2000;// number of cells
            k_dim = 35; // 1x35
            max_win_height = 2000;
            max_win_width = 1;
            break;
        }
    }
    h_divisor = stride;
    if(max_win_height > 0){
        h_divisor *= max_win_height;
    }
    w_divisor = stride;
    if(max_win_width > 0){
        w_divisor *= max_win_width;
    }
}

void resetParams() {
    curr_layer = 0;
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
}

void updateParamsForFCL(uint32_t* bias_dimensions, bool firstFCL){
    if(firstFCL){
        nodes_in = i_height*i_width*i_channel;
    }
    else{
        nodes_in = nodes_out;
        curr_layer++; //if firstFCL: curr_layer has been updated by updateParamsAfterCL
    }
    nodes_out = bias_dimensions[curr_layer];
}

void updateParamsAfterCL(){
    if(trained_for_MNIST) {
        i_height = (i_height - k_dim + 1) / h_divisor;
        i_width = (i_width - k_dim + 1) / w_divisor;
    }
    else{
        i_height = i_channel;
        i_width = i_height; //both dimensions are eliminated to 1 after maxpool
    }
    i_channel = k_number;
    curr_layer++;
}

int main(int argc, char* argv[]) {
    if (argc < 6){
        cout << "Calling proxy without specifying role (1), port (2), address (3), helpers port (4) and helpers adress (5) is not possible." << endl;
        cout << "Specify the Network Mode with the 6th parameter: 0 = CHAMELEON, 1 = LeNet5, 2 = Cell_CNN" << endl;
        return 1;
    }
    uint8_t role = atoi(argv[1]);
    uint16_t cport = atoi(argv[2]);
    string caddress(argv[3]);
    uint16_t hport = atoi(argv[4]);
    string haddress(argv[5]);
    uint32_t nn_mode = atoi(argv[6]);

    // parse model parameters
    initParams(nn_mode);
    double**** model_weights;
    auto* bias_dimensions = new uint32_t [layer_number];
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
            bias_dimensions[2] = 500;
            bias_dimensions[3] = 10;
            break;
        }
        case 2:{
            cout << "CellCNN" << endl;
            model_weights = getCellCnnParameters(MODEL_DIR + CELL_CNN_MODEL_FILES, 8);
            bias_dimensions[0] = 8;
            bias_dimensions[1] = 2;
            break;
        }
        default: {
            cout << "No neural network mode is matching the supported one." << endl;
            return -1;
        }
    }
    cout << "Reading input data... ";
    vector<vector<unsigned char>> test_set;
    vector<unsigned char> test_label;
    if(trained_for_MNIST){
        //read in MNIST data
        mnist::MNIST_dataset<vector, vector<uint8_t>, uint8_t> dataset =
                mnist::read_dataset<vector, vector, uint8_t, uint8_t>(MNIST_DATA);
        test_set = dataset.test_images;
        test_label = dataset.test_labels;
    }
    else{
        //TODO read multi cell data
        /*double* random_label = random_1D_data(proxy, i_number, 10.0);
        for (int i = 0; i < i_number; ++i) {
            test_label.push_back(random_label[i]);
        }
        double** random_values = random_2D_data(proxy, i_number, i_height*i_width, -10.0, 10.0);
        for (int i = 0; i < i_number; ++i) {
            vector<unsigned char> values;
            for (int j = 0; j < i_height*i_width; ++j) {
                values.push_back(random_values[i][j]);
            }
            test_set.push_back(values);
        }*/
    }

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

    auto*** data = new uint64_t** [test_set.size()];
    for (uint32_t i = 0; i < i_number; ++i) {
        data[i] = new uint64_t *[i_height];
        for (uint32_t r = 0; r < i_height; ++r) {
            data[i][r] = new uint64_t [i_width];
            for (uint32_t c = 0; c < i_width; ++c) {
                double pixelValue = test_set.at(i).at(r * i_width + c);
                data[i][r][c] = proxy->createShare(pixelValue);                 // store directly as secret shares
            }
        }
        if(padding > 0){
            uint64_t padding_value = proxy->createShare(0);
            data[i] = PAD(data[i], i_height, i_width, padding_value, padding);
        }
    }
    cout << "Creating secret shares: " << endl;
    // secret shares: bias
    auto **bias = new uint64_t *[layer_number];
    for (int layer = 0; layer < layer_number; ++layer) {
        //cout << " - Bias for layer " << layer << " of length " << bias_dimensions[layer] << endl;
        bias[layer] = proxy->createShare(model_weights[layer + layer_number][0][0], bias_dimensions[layer]);
        delete[] model_weights[layer + layer_number][0][0];
        delete[] model_weights[layer + layer_number][0];
        delete[] model_weights[layer + layer_number];
    }

    // secret shares: kernel
    k_number = bias_dimensions[0];
    auto ***kernel = new uint64_t**[i_channel]; // weights of all kernel for first CL
    for (int i = 0; i < i_channel; ++i) {
        kernel[i] = new uint64_t *[k_number];
    }
    uint32_t k_size = nn_mode == 2 ? k_dim : k_dim*k_dim;
    for (uint32_t k = 0; k < k_number; k++) {
        uint64_t ** kernel_weights = proxy->createShare(model_weights[0][k], i_channel, k_size);
        for (int i = 0; i < i_channel; ++i) {
            kernel[i][k] = kernel_weights[i];
        }
        delete[] model_weights[0][k][0];
        delete[] model_weights[0][k];       // if more than one channel: delete for each channel
    }
    delete[] model_weights[0];

    auto* prediction = new double [i_number];
    double correct = 0, incorrect = 0;
    // CNN INFERENCE PIPELINE
    for (uint32_t image = 0; image < i_number; ++image) {
        cout << "INFERENCE PIPELINE " << image << endl;
        k_number = bias_dimensions[0];
        resetParams();
        auto *** input = new uint64_t **[i_channel];
        input[0] = data[image]; // currently only 1 channel for input supported

        uint64_t*** conv;
        uint64_t* prev_layer_res;
        uint64_t **weights;
        switch (nn_mode) {                              // from here on network architectures differ
            case 0:{ // Chameleon
                conv = CL(proxy, input, i_channel, i_height, i_width, kernel, k_dim, k_number, stride, max_win_height, max_win_width, bias[curr_layer], true);
                updateParamsAfterCL();
                updateParamsForFCL(bias_dimensions, true);
                // FULLY CONNECTED LAYER
                weights = proxy->createShare(model_weights[curr_layer][0], nodes_out, nodes_in);
                if(image == i_number-1) {
                    delete[] model_weights[curr_layer][0];
                    delete[] model_weights[curr_layer];
                }
                prev_layer_res = FCL(proxy, conv[0][0], nodes_in, weights, nodes_out, bias[curr_layer]);
                updateParamsForFCL(bias_dimensions, false);
                break;
            }
            case 1: { // SecureNN
                conv = CL(proxy, input, i_channel, i_height, i_width, kernel, k_dim, k_number, stride, max_win_height, max_win_width, bias[curr_layer], false);
                updateParamsAfterCL();
                // PERFORMING CONVOLUTION
                k_number = bias_dimensions[curr_layer];
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
                conv = CL(proxy, conv, i_channel, i_height, i_width, kernel, k_dim, k_number, stride, max_win_height, max_win_width, bias[curr_layer], true);
                updateParamsAfterCL();

                // fully connected layer:
                updateParamsForFCL(bias_dimensions, true);
                weights = proxy->createShare(model_weights[curr_layer][0], nodes_out, nodes_in);
                if(image == i_number-1) {
                    delete[] model_weights[curr_layer][0];
                    delete[] model_weights[curr_layer];
                }
                //uint64_t * flattened = FLT(conv, i_height, i_width, i_channel);
                prev_layer_res = FCL(proxy, conv[0][0], nodes_in, weights, nodes_out, bias[curr_layer]);
                updateParamsForFCL(bias_dimensions, false);
                break;
            }
            case 2:{
                // CONVOLUTIONAL LAYER (adapted to non-symmetric filter size
                conv = new uint64_t **[k_number];
                for (int k = 0; k < k_number; ++k) {
                    conv[k] = MATVECMUL(proxy, input, kernel[k], 1, i_height, i_width);
                    conv[k][0][0] = MAX(proxy, conv[k][0], i_height);
                }
                updateParamsAfterCL();

                prev_layer_res = FLT(conv, i_height, i_width, i_channel);
                updateParamsForFCL(bias_dimensions, true);
            }
        }
        delete[] conv;
        // last FULLY CONNECTED LAYER
        weights = proxy->createShare(model_weights[curr_layer][0], nodes_out, nodes_in);
        if(image == i_number-1) {
            delete[] model_weights[curr_layer][0];
            delete[] model_weights[curr_layer];
        }

        uint64_t * output = FCL(proxy, prev_layer_res, nodes_in, weights, nodes_out, bias[curr_layer]);
        switch (nn_mode) {                              // from here on network architectures differ
            case 0:{ // Chameleon (treated same as SecureNN)
            }
            case 1: { // SecureNN
                prediction[image] = convert2double(REC(proxy, ARGMAX(proxy, output, nodes_out)));
                //print1DArray("Input to ARGMAX:", convert2double(REC(proxy, output, nodes_out), nodes_out), nodes_out);
                cout << ": predicted " << prediction[image] << ", correct is " << int(test_label[image]) << endl;
                break;
            }
            case 2:{
                prediction[image] = convert2double(REC(proxy, ARGMAX(proxy, output, nodes_out)));
                cout << ": predicted " << prediction[image] << ", correct is " << int(test_label[image]) << endl;
                break;
            }
        }
        delete[] output;
        if (prediction[image] == int(test_label[image])){
            correct++;
        }
        else{
            incorrect++;
        }
    }
    delete[] data[0];
    delete[] data;
    cout << "accuracy: " << printf("%.3f", correct/i_number) << " (" << correct << "/" << i_number << ")" << endl;
    print1DArray("Prediction: ", prediction, i_number);
    string s_correct = "";
    for (int i = 0; i < i_number; ++i)
        s_correct += to_string(test_label[i]) + "\t";

    cout << "Correct label: " << endl << s_correct << endl;
    proxy->PrintBytes();
    return 0;
}


