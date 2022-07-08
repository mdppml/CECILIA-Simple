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
    uint32_t stride = 1;
    uint32_t padding = 0;
    bool doMaxpool = true;

    uint32_t divisor = stride;
    if(doMaxpool){
        divisor *= 2;
    }
    switch (nn_mode) {
        case 0:{
            cout << "CHAMELEON" << endl;
            stride = 2;
            divisor = stride;
            doMaxpool = false;
            padding = 2;
            // this is a well defined CNN architecture trained for the MNIST dataset
            model_weights = getChameleonParameters(MODEL_DIR + CHAMELEON_MODEL_FILES, k_number, i_channel);
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
        default:{
            cout << "No or unknown network has been specified. Use random weights to perform CNN inference for 3 convolutional layer and one fully connected layer" << endl;
            k_number = 16;
            model_weights = new double ***[4];
            // fill weights for convolutional layer
            for(uint32_t layer = 0; layer<3; layer++){
                model_weights[layer] = new double **[k_number];
                for(uint64_t kernel = 0; kernel< k_number; kernel++){
                    model_weights[layer][kernel] = random_2D_data(proxy, k_dim-layer, k_dim-layer, -5.0, 5.0);
                }
                k_number -= (layer + 1);
                // layer 0 has 5 kernel: 5x5
                // layer 1 has 4 kernel: 4x4 and
                // layer 2 has 2 kernel: 3x3
            }
            uint32_t in_fcl = (((((i_height - k_dim + 1) / divisor) - k_dim) / divisor) - k_dim - 1) / divisor;
            uint32_t out_fcl = 10;
            // fill weights for fully connected layer
            model_weights[3] = new double **[1];
            model_weights[3][0] = random_2D_data(proxy, in_fcl, out_fcl, -5.0, 5.0);
            break;
        }
    }
    //create shares
    uint64_t*** data = new uint64_t** [i_channel]; //TODO read in real images
    for(uint64_t i = 0; i<i_channel; i++){
        data[i] = proxy->createShare(random_2D_data(proxy, i_height, i_width, 0.0, 255.0), i_height, i_width);
    }
    if(padding > 0){
        for(uint32_t i = 0; i<i_channel; i++){
            data[i] = PAD(data[i], i_height, i_width, 0, padding);
        }
        i_height += 2*padding;
        i_width += 2*padding;
    }

    uint64_t ***kernel = new uint64_t**[k_number];
    for (uint32_t i = 0; i < k_number; i++) {
        kernel[i] = proxy->createShare(model_weights[0][i], i_channel, k_dim*k_dim);
    }

    // CNN INFERENCE PIPELINE
    /*uint32_t params[7];
    params[0] = i_channel;
    params[1] = i_height;
    params[2] = i_width;
    params[3] = k_dim; // kernel size
    params[4] = k_number; // kernel number is also output channel
    params[5] = stride; // stride
    params[6] = doMaxpool;
    proxy->SendBytes(CNN_CL, params, 7);*/

    // PERFORMING CONVOLUTION
    uint64_t*** conv = CL(proxy, data, i_channel, i_height, i_width, kernel, k_dim, k_number, stride, doMaxpool);
    i_channel = k_number;
    //          size after conv         after maxpool divide by 2
    i_height = ((i_height - k_dim) + 1) / divisor;
    i_width = ((i_width - k_dim) + 1) / divisor;
    //init what the architectures share for FCL:
    uint32_t nodes_out;
    uint32_t nodes_in = i_height*i_width*i_channel;

    uint64_t* prev_layer_res;
    uint64_t **weights = nullptr;
    // from here on network architectures differ
    switch (nn_mode) {
        case 0:{
            k_number = 1;
            nodes_out = 100;
            // FULLY CONNECTED LAYER
            weights = proxy->createShare(model_weights[1][0], nodes_out, nodes_in);
            uint64_t *flattened = FLT(conv, i_height, i_width, i_channel);
            /*params[0] = nodes_in;
            params[1] = nodes_out;
            proxy->SendBytes(CNN_FCL, params, 2);*/
            prev_layer_res = FCL(proxy, flattened, nodes_in, weights, nodes_out);

            nodes_in = nodes_out;
            i_channel = 1;
            nodes_out = 10;
            break;
        }
        case 1:{
            // PERFORMING CONVOLUTION
            kernel = new uint64_t**[k_number];
            for (uint32_t i = 0; i < k_number; i++) {
                kernel[i] = proxy->createShare(model_weights[1][i], i_channel, k_dim*k_dim);
            }
            //TODO might be problematic if conv changes size...
            conv = CL(proxy, conv, i_channel, i_height, i_width, kernel, k_dim, k_number, stride, doMaxpool);
            i_channel = k_number;
            i_height = ((i_height - k_dim) + 1) / (stride*2);
            i_width = ((i_width - k_dim) + 1) / (stride*2);
            cout << "finished CL2 (MiniONN)" << endl;

            // fully connected layer:
            nodes_out = 100;
            nodes_in = i_height * i_width * i_channel;
            weights = proxy->createShare(model_weights[2][0], nodes_out, nodes_in);
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
                    kernel[i] = proxy->createShare(model_weights[c+1][i], i_channel, k_dim*k_dim);
                }
                //TODO might be problematic if conv changes size...
                conv = CL(proxy, conv, i_channel, i_height, i_width, kernel, k_dim, k_number, stride, doMaxpool);
                i_channel = k_number;
                i_height = ((i_height - k_dim) + 1) / (stride*2);
                i_width = ((i_width - k_dim) + 1) / (stride*2);
                cout << "finished CL" << c+2 << " (random mode)"<< endl;
            }
            // fully connected layer:
            nodes_out = 10;
            nodes_in = i_height * i_width * i_channel;
            weights = proxy->createShare(model_weights[3][0], nodes_out, nodes_in);
            prev_layer_res = FLT(conv, i_height, i_width, i_channel);
            break;
        }
    }
    // FULLY CONNECTED LAYER
    /*params[0] = nodes_in;
    params[1] = nodes_out;
    proxy->SendBytes(CNN_FCL, params, 2);*/

    uint64_t * output = FCL(proxy, prev_layer_res, nodes_in, weights, nodes_out);

    for(uint64_t n = 0; n<nodes_out; n++){
        cout << "Node " << n << ": " << output[n] << endl;
    }

 //   proxy->SendBytes(CORE_END);
    proxy->PrintBytes();
    return 0;
}

