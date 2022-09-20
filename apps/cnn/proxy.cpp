#include <cstdlib>
#include <iostream>
#include <deque>
#include <chrono>
#include <assert.h>
#include "../../core/cnn.h"
#include "../../core/auc.h"
#include "model_parser.h"
#include "mnist/mnist_reader.hpp"

using namespace std;
const string MODEL_DIR = "../../apps/cnn/model_files/";
const string CHAMELEON_MODEL_FILES = "Chameleon_CNN/", LENET_NN_MODEL_FILES = "LeNet_trained/LeNetNN_b4_n10/", MINIONN_MODEL_FILE = "MiniONN.txt";
const string MNIST_DATA = "../../apps/cnn/mnist/";

uint64_t layer_number;
uint32_t i_number, i_channel, i_width, i_height;
uint32_t k_number, k_dim;
uint32_t stride, padding;
uint32_t divisor;
uint32_t maxpool_window_dim;
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
    maxpool_window_dim = 2;
    switch (mode) {
        case 0:{ //CHAMELEON
            layer_number = 3;
            stride = 2;
            padding = 2;
            maxpool_window_dim = 0;
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
    double**** model_weights;
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
            bias_dimensions[2] = 500;
            bias_dimensions[3] = 10;
            break;
        }
        case 2:{
            //TODO Single cell network
            cout << "CellCNN" << endl;
            uint32_t k = 8;
            model_weights = new double ***[layer_number*2];
            model_weights[0] = new double **[k]; // weights for cl
            for(uint32_t c = 0; c < k; c++){
                model_weights[0][c] = new double*[i_width];
                model_weights[0][c] = random_2D_data(proxy, i_channel, k_dim, -254.0, 255.0);
            }
            model_weights[1] = new double **[1]; // weights for fcl
            model_weights[1][0] = random_2D_data(proxy, 2, k, -254.0, 255.0);

            model_weights[2] = new double **[1]; // bias for cl
            model_weights[2][0] = random_2D_data(proxy, 1, 8, -5.0, 5.0);

            model_weights[3] = new double **[1]; // bias for fcl
            model_weights[3][0] = random_2D_data(proxy, 1, 2, -5.0, 5.0);

            bias_dimensions[0] = 8;
            bias_dimensions[1] = 2;
            break;
        }
        default: {
            cout << "No neural network mode is matching the supported one." << endl;
            return -1;
        }
    }
    vector<vector<unsigned char>> test_set;
    vector<unsigned char> test_label;
    if(trained_for_MNIST){
        //read in MNIST data
        cout << "Reading input data... ";
        mnist::MNIST_dataset<vector, vector<uint8_t>, uint8_t> dataset =
                mnist::read_dataset<vector, vector, uint8_t, uint8_t>(MNIST_DATA);
        test_set = dataset.test_images;
        test_label = dataset.test_labels;
    }
    else{
        cout << "generate random labels...";
        //TODO read multi cell data
        double* random_label = random_1D_data(proxy, i_channel, 10.0);
        for (int i = 0; i < i_channel; ++i) {
            test_label.push_back(random_label[i]);
        }
        cout << "generate random input data...";
        double** random_values = random_2D_data(proxy, i_channel, i_height*i_width, -10.0, 10.0);
        for (int i = 0; i < i_channel; ++i) {
            vector<unsigned char> values;
            for (int j = 0; j < i_height*i_width; ++j) {
                values.push_back(random_values[i][j]);
            }
            test_set.push_back(values);
        }
    }

    uint64_t*** data = new uint64_t** [test_set.size()];
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
    uint32_t k_size = nn_mode == 2 ? k_dim : k_dim*k_dim;
    for (uint32_t i = 0; i < k_number; i++) {
        cout << "Kernel " << i << endl;
        kernel[i] = proxy->createShare(model_weights[0][i], i_channel, k_size);
        //print2DArray("KERNEL: ", convert2double(REC(proxy, kernel[i], i_channel, k_dim*k_dim), i_channel, k_dim*k_dim), i_channel, k_dim*k_dim);
        delete[] model_weights[0][i][0]; // if more than one channel: delete for each channel
        delete[] model_weights[0][i];
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

        // PERFORMING CONVOLUTION
        //print2DArray("input", convert2double(REC(proxy, input[0], i_height, i_width), i_height, i_width), i_height, i_width);
        //print2DArray("KERNEL: ", convert2double(REC(proxy, kernel[0], i_channel, k_dim*k_dim), i_channel, k_dim*k_dim), i_channel, k_dim*k_dim);
        //print1DArray("BIAS: ", convert2double(REC(proxy, bias[curr_layer], bias_dimensions[curr_layer]), bias_dimensions[curr_layer]), bias_dimensions[curr_layer]);
        uint64_t*** conv;
        uint64_t* prev_layer_res;
        uint64_t **weights;
        switch (nn_mode) {                              // from here on network architectures differ
            case 0:{ // Chameleon
                conv = CL(proxy, input, i_channel, i_height, i_width, kernel, k_dim, k_number, stride, maxpool_window_dim, bias[curr_layer], true);
                updateParamsAfterCL();
                updateParamsForFCL(bias_dimensions, true);
                // FULLY CONNECTED LAYER
                weights = proxy->createShare(model_weights[curr_layer][0], nodes_out, nodes_in);
                if(image == i_number-1) {
                    delete[] model_weights[curr_layer][0];
                    delete[] model_weights[curr_layer];
                }
                //uint64_t *flattened = FLT(conv, i_height, i_width, i_channel);
                prev_layer_res = FCL(proxy, conv[0][0], nodes_in, weights, nodes_out, bias[curr_layer]);
                updateParamsForFCL(bias_dimensions, false);
                break;
            }
            case 1: { // SecureNN
                conv = CL(proxy, input, i_channel, i_height, i_width, kernel, k_dim, k_number, stride, maxpool_window_dim, bias[curr_layer], false);
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
                conv = CL(proxy, conv, i_channel, i_height, i_width, kernel, k_dim, k_number, stride, maxpool_window_dim, bias[curr_layer], true);
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
                }
                //TODO MAX with assymetric window size...
                cout << "CL for CellCNN done, maxpool not yet." << endl;
                updateParamsAfterCL();

                // FULLY CONNECTED LAYER
                /*updateParamsForFCL(bias_dimensions, true);
                weights = proxy->createShare(model_weights[curr_layer][0], nodes_out, nodes_in);
                if(image == i_number-1) {
                    delete[] model_weights[curr_layer][0];
                    delete[] model_weights[curr_layer];
                }
                //uint64_t *flattened = FLT(conv, i_height, i_width, i_channel);
                prev_layer_res = FCL(proxy, conv[0][0], nodes_in, weights, nodes_out, bias[curr_layer]);*/
                updateParamsForFCL(bias_dimensions, true);
            }
        }
        // last FULLY CONNECTED LAYER
        cout << "prepare kernel for second FCL" << endl;
        weights = proxy->createShare(model_weights[curr_layer][0], nodes_out, nodes_in);
        if(image == i_number-1) {
            delete[] model_weights[curr_layer][0];
            delete[] model_weights[curr_layer];
        }

        cout << "second FCL" << endl;
        uint64_t * output = FCL(proxy, prev_layer_res, nodes_in, weights, nodes_out, bias[curr_layer]);
        switch (nn_mode) {                              // from here on network architectures differ
            case 0:{ // Chameleon
                prediction[image] = convert2double(REC(proxy, ARGMAX(proxy, output, nodes_out)));
                //print1DArray("Input to ARGMAX:", convert2double(REC(proxy, output, nodes_out), nodes_out), nodes_out);
                cout << ": predicted " << prediction[image] << ", correct is " << int(test_label[image]) << endl;
                break;
            }
            case 1: { // SecureNN
                uint64_t * out_relu = RELU(proxy, output, nodes_out);
                uint64_t sum_out_relu = out_relu[0];
                for (int i = 0; i < nodes_out; ++i) {
                    sum_out_relu += out_relu[i];
                }
                uint64_t v_norm_val[nodes_out];
                // memcpy(v_norm_val, &sum_out_relu, nodes_out*sizeof(sum_out_relu));
                //output = DIVISION(proxy, out_relu, sum_out_relu); // vectorized DIV
                for (int i = 0; i < nodes_out; ++i) {
                    output[i] = DIVISION(proxy, out_relu[i], sum_out_relu);
                }
                //TODO if DIV by sum... this is factor --> only perform argmax...
                prediction[image] = convert2double(REC(proxy, ARGMAX(proxy, output, nodes_out)));
                //print1DArray("Input to ARGMAX:", convert2double(REC(proxy, output, nodes_out), nodes_out), nodes_out);
                cout << ": predicted " << prediction[image] << ", correct is " << int(test_label[image]) << endl;
                print1DArray("output nodes:", convert2double(REC(proxy, output, nodes_out), nodes_out), nodes_out);
                break;
            }
            case 2:{

            }
        }
        if (prediction[image] == int(test_label[image])){
            correct++;
        }
        else{
            incorrect++;
        }
    }
    cout << "accuracy: " << printf("%.3f", correct/i_number) << " (" << correct << "/" << i_number << ")" << endl;
    proxy->PrintBytes();
    return 0;
}


