#include <cstdlib>
#include <iostream>
#include <deque>
#include <chrono>
#include <assert.h>
#include "../../core/cnn.h"
#include "model_parser.h"

using namespace std;
const string MODEL_DIR = "../../apps/cnn/model_files/";
const string MINIONN_MODEL_FILE = "MiniONN.txt", CHAMELEON_MODEL_FILES = "Chameleon_CNN/";

int main(int argc, char* argv[]) {
    if (argc < 6){
        cout << "Calling proxy without specifying role (1), port (2), address (3), helpers port (4) and helpers adress (5) is not possible." << endl;
        cout << "Specify the Network Mode with the 6th parameter: 0 = CHAMELEON, 1 = MINIONN" << endl;
        return 1;
    }
    cout << "start proxy" << endl;
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
    vector<int> d;
    double**** model_weights;

    // init input parameter --> all supported networks expect the same input shape.
    uint32_t i_number = 10;
    uint32_t i_channel = 1;
    uint32_t i_width = 28;
    uint32_t i_height = 28;
    uint32_t k_number = 5;
    uint32_t k_dim = 5;
    switch (nn_mode) {
        case 0:{
            cout << "CHAMELEON" << endl;
            // this is a well defined CNN architecture trained for the MNIST dataset
            model_weights = getChameleonParameters(MODEL_DIR + CHAMELEON_MODEL_FILES, k_number, 1);
            break;
        }
        case 1:{
            cout << "MINIONN" << endl;
            model_weights[0] = getONNXParameters(MODEL_DIR + MINIONN_MODEL_FILE, d, 4); //TODO I dont know which weights belong to which kernel in one layer
            uint64_t * weight_dims = new uint64_t [d.at(0)];
            for (uint64_t i = 0; i<d.at(0); i++){
                int w = d[i+1];
                weight_dims[i] = convert2uint64(w);
            }
            print1DArray("weight dimensions", convert2double(weight_dims, d.at(0)), d.at(0));
            //print2DArray("model weights - layer 2: ", model_weights[0], d.at(1), d.at(2));
            print2DArray("model weights - layer 3: ", model_weights[0][1], d.at(3), 1, 0);

            uint32_t input_param_start = d.at(0)+1;
            i_number = d.at(input_param_start);
            i_channel = d.at(input_param_start + 1);
            i_width = d.at(input_param_start + 2);
            i_height = d.at(input_param_start + 2);
            k_number = 5; //TODO
            k_dim = 5;
            break;
        }
        default:
            cout << "No or unknown network has been specified. Use random data." << endl;
            //TODO generate random weights
            break;
    }

    cout << "parsed parameters: # of input = " << i_number << ", i_channel = " << i_channel << ", widht = " << i_width << ", height = " << i_height << endl;
    uint64_t matrix_size = i_height*i_width;
    uint64_t*** data = new uint64_t** [i_channel]; //TODO this are the images

    //convert matrix elements
    cout << "init weights done" << endl;
    uint64_t ***kernel = new uint64_t**[k_number];
    print2DArray("kernel 0", model_weights[0][0], i_channel, k_dim*k_dim);
    for (uint32_t i = 0; i < k_number; i++) {
        kernel[i] = proxy->createShare(model_weights[0][i], i_channel, k_dim*k_dim);
    }
    cout << "init kernel done" << endl;
    //helper variable to create shares of v1 and v2
    uint64_t ***mTmp = new uint64_t **[i_channel]; //pick random value in ring --> shareOfMatrix is calculated from that
    uint64_t ***shareOfMatrix = new uint64_t **[i_channel];

    for(uint32_t i = 0; i < i_channel; i++){
        data[i] = convert2uint64(random_2D_data(proxy, i_height, i_width, 0.0, 255.0), i_height, i_width); //TODO remove when read from file
        mTmp[i] = convert2uint64(random_2D_data(proxy, i_height, i_width, 0.0, 255.0), i_height, i_width);
    }
    cout << "Start creating share for proxy " << to_string(role) << endl;
    //create shares
    if(role == P1) {
        for(uint32_t k = 0; k < i_channel; k++){
            shareOfMatrix[k] = new uint64_t*[i_height];
            for(uint32_t j = 0; j < i_height; j++) {
                shareOfMatrix[k][j] = new uint64_t [i_width];
                for(uint32_t l = 0; l < i_width; l++) {
                    shareOfMatrix[k][j][l] = data[k][j][l] - mTmp[k][j][l];
                }
            }
        }
    }
    else {
        shareOfMatrix = mTmp;
    }
    cout << "finished creating shares... Call Convolutional layer" << endl;

    // CNN inference pipeline
    // convolutional layer
    uint32_t mmaxParams[6];
    mmaxParams[0] = i_channel;
    mmaxParams[1] = i_height;
    mmaxParams[2] = i_width;
    mmaxParams[3] = k_dim; // kernel size
    mmaxParams[4] = k_number; // kernel number is also output channel
    mmaxParams[5] = 1; // stride
    proxy->SendBytes(CNN_CL, mmaxParams, 6);

    // PERFORMING CONVOLUTION
    uint64_t*** conv = CL(proxy, data, i_channel, i_height, i_width, kernel, k_dim, k_number, 1);
    cout << "finished CL1" << endl;
    //          size after conv                   after maxpool divide by 2
    i_channel = k_number;
    i_height = ((i_height - k_dim) + 1) / 2;
    i_width = ((i_width - k_dim) + 1) / 2;
    cout << "new params: i_channel = " << i_channel << ", h = " << i_height << ", w = " << i_width << endl;
    //what architectures share for FCL:
    uint32_t nodes_out = 100;
    uint32_t nodes_in = i_height*i_width*i_channel;
    uint64_t **weights2;
    if(nn_mode == 0){
        // from here on network architectures differ
        k_number = 1;
        // fully connected layer:
        weights2 = proxy->createShare(model_weights[1][0], nodes_out, nodes_in);
    }
    else{
        // PERFORMING CONVOLUTION
        //TODO might be problematic if conv changes size...
        conv = CL(proxy, conv, i_channel, i_height, i_width, kernel, k_dim, k_number, 1);
        i_channel = k_number;
        i_height = ((i_height - k_dim) + 1) / 2;
        i_width = ((i_width - k_dim) + 1) / 2;
        cout << "finished CL2" << endl;

        // fully connected layer:
        weights2 = new uint64_t*[nodes_in];

        for (uint32_t n = 0; n<nodes_in; n++){
            weights2[n] = new uint64_t [nodes_out];
            for (uint32_t j = 0; j < matrix_size*matrix_size; j++) {
                weights2[n][j] = proxy->generateRandom();
            }
        }
    }

    mmaxParams[0] = i_height;
    mmaxParams[1] = k_number; // input number is same as previous kernel number
    mmaxParams[2] = nodes_out;

    // PERFORMING CONVOLUTION
    proxy->SendBytes(CNN_FCL, mmaxParams, 3);
    cout << "params for FCL sent" << endl;
    uint64_t * output = FCL(proxy, conv, i_height, k_number, weights2, nodes_out);
    cout << "FCL done" << endl;

    for(uint64_t n = 0; n<nodes_out; n++){
        cout << "Node " << n << ": " << output[n] << endl;
    }

    proxy->SendBytes(CORE_END);
    proxy->PrintBytes();
    return 0;
}

