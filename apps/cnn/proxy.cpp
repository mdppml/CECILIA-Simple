#include <cstdlib>
#include "stdio.h"
#include <iostream>
#include <deque>
#include <chrono>
#include <assert.h>
#include "../../core/cnn.h"
#include "model_parser.h"
#include "mnist_data/mnist_reader.hpp"

#include <fstream>

using namespace std;
const string MODEL_DIR = "../apps/cnn/model_files/";
const string CHAMELEON_MODEL_FILES = MODEL_DIR + "Chameleon_CNN/", LENET_NN_MODEL_FILES = MODEL_DIR + "LeNet_trained/LeNetNN_b4_n5/",
        CELL_CNN_MODEL_FILES = MODEL_DIR + "Cell_CNN/", MINIONN_MODEL_FILE = "MiniONN.txt";
const string LENET_CORRECT_PATH = LENET_NN_MODEL_FILES + "correctness/",
        CHAMELEON_CORRECT_PATH = CHAMELEON_MODEL_FILES + "correctness/";
const string MNIST_DATA = "../apps/cnn/mnist_data/";
const double PIXEL_MAX = 255; // can be used for normalization if desired.
const bool eval_correctness = true;

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

/**
 * Initializes all parameters defining the structure of the used networks (Chameleon (mode 0), (LeNet mode 1), CellCnn (mode 2) but not yet usable, and potentially others.)
 * Those parameters include number of images for which inference shall be computed, dimensions of input images, dimensions of kernel and others.
 * @param mode the mode specifying the network for which the parameters are initialized.
 */
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
        case 0: { //CHAMELEON
            layer_number = 3;
            stride = 2;
            padding = 2;
            max_win_width = 0, max_win_height = 0;
            break;
        }
        case 1: { //LeNet
            layer_number = 4;
            break;
        }
        case 2: { //CellCNN
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
        default: {
            cout << "No neural network mode is matching " << mode << ". Parameters can not be initialized." << endl;
        }
    }
    h_divisor = stride;
    if (max_win_height > 0) {
        h_divisor *= max_win_height;
    }
    w_divisor = stride;
    if (max_win_width > 0) {
        w_divisor *= max_win_width;
    }
}
/**
 * Reset those network parameters that have been changed during inference computation, so that the state is equal to that after initParams() have been called.
 */
void resetParams() {
    curr_layer = 0;
    i_channel = 1;
    if (trained_for_MNIST) {
        i_channel = 1;
        i_height = 28;
        i_width = 28;
    } else {
        //i_channel = 64; // number of samples
        i_width = 35;   // number of markers
        i_height = 2000;// number of cells
    }
    if (padding > 0) {
        i_height += 2 * padding;
        i_width += 2 * padding;
    }
}
/**
 * Update parameters as preparation for a subsequent Fully Connected Layer
 * @param bias_dimensions vector defining the number of bias values needed per layer.
 * @param firstFCL is true if the subsequent Fully Connected Layer (FCL) is the first FCL within the according network.
 */
void updateParamsForFCL(const uint32_t *bias_dimensions, bool firstFCL) {
    if (firstFCL) {
        nodes_in = i_height * i_width * i_channel;
    } else {
        nodes_in = nodes_out;
        curr_layer++; //if firstFCL: curr_layer has been updated by updateParamsAfterCL
    }
    nodes_out = bias_dimensions[curr_layer];
}
/**
 * Update parameters after Convolutional Layer.
 */
void updateParamsAfterCL() {
    if (trained_for_MNIST) {
        i_height = (i_height - k_dim + 1) / h_divisor;
        i_width = (i_width - k_dim + 1) / w_divisor;
    } else {
        i_height = i_channel;
        i_width = i_height; //both dimensions are eliminated to 1 after maxpool
    }
    i_channel = k_number;
    curr_layer++;
}
/**
 * Deletion routine for secret shares of images.
 * @param images 4D matrix where all secret shared images are stored.
 * @param image_to_delete position of the image which should be deleted to free up storage.
 */
void delete_image_share(uint64_t**** images, uint32_t image_to_delete){
    for (int j = 0; j < i_channel; ++j) {
        for (uint32_t r = 0; r < i_height; ++r) {
            delete[] images[image_to_delete][j][r];
        }
        delete[] images[image_to_delete][j];
    }
    delete[] images[image_to_delete];
}
/**
 * Deletion routine for the network model weights.
 * @param model_weights matrix containing kernel and bias weights for each layer
 * @param i_layer specifies the position of the layer for which kernel weights or bias weights shall be deleted.
 * @param i_channel number of channels for the specified weights; only needed if kernel weights are deleted and therefore i_layer < than the number of layers in the network
 * @param number_of_kernels number of kernels for the specified weights; only needed if kernel weights are deleted and therefore i_layer < than the number of layers in the network
 */
void delete_model_weights(double ****model_weights, uint32_t i_layer, uint32_t i_channel, uint32_t number_of_kernels){
    if (i_layer < layer_number){
        // kernel weights:
        for (int j = 0; j < i_channel; ++j) {
            for (int k = 0; k < number_of_kernels; ++k) {
                delete[] model_weights[i_layer][j][k];
            }
            delete[] model_weights[i_layer][j];
        }
    }
    else{
        // bias weights are 1d vectors stored in the same matrix
        delete[] model_weights[i_layer][0][0];
        delete[] model_weights[i_layer][0];
    }
    delete[] model_weights[i_layer];
}

int main(int argc, char *argv[]) {
    if (argc < 6) {
        cout
                << "Calling proxy without specifying role (1), port (2), address (3), helpers port (4) and helpers adress (5) is not possible."
                << endl;
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
    double ****model_weights;
    auto *bias_dimensions = new uint32_t[layer_number];
    switch (nn_mode) {
        case 0: {
            cout << "CHAMELEON" << endl;
            model_weights = getChameleonParameters(CHAMELEON_MODEL_FILES, 5);
            bias_dimensions[0] = 5;
            bias_dimensions[1] = 100;
            bias_dimensions[2] = 10;
            break;
        }
        case 1: {
            cout << "LeNet by SecureNN" << endl;
            model_weights = getLeNetParameters(LENET_NN_MODEL_FILES, true);
            bias_dimensions[0] = 20;
            bias_dimensions[1] = 50;
            bias_dimensions[2] = 500;
            bias_dimensions[3] = 10;
            break;
        }
        case 2: {
            cout << "CellCNN" << endl;
            model_weights = getCellCnnParameters(CELL_CNN_MODEL_FILES, 8);
            bias_dimensions[0] = 8;
            bias_dimensions[1] = 2;
            break;
        }
        default: {
            cout << "No neural network mode is matching the supported one." << endl;
            return -1;
        }
    }
    cout << "Reading input data... " << endl;
    vector<vector<unsigned char>> test_set;
    vector<unsigned char> test_label;
    if (trained_for_MNIST) {
        //read in MNIST data
        mnist::MNIST_dataset<vector, vector<uint8_t>, uint8_t> dataset =
                mnist::read_dataset<vector, vector, uint8_t, uint8_t>(MNIST_DATA);
        test_set = dataset.test_images;
        test_label = dataset.test_labels;
        if (eval_correctness and role == 0) {
            for (int i = 0; i < i_number; ++i) {
                ofstream image_file;
                string path;
                switch (nn_mode) {
                    case 0:
                        path = CHAMELEON_CORRECT_PATH;
                        break;
                    case 1:
                        path = LENET_CORRECT_PATH;
                }
                path += "eval_secure_" + to_string(i) + ".txt";
                image_file.open(path, std::ios::out);
                if (!image_file) {
                    // create file first
                    cout << "create file: " << path << endl;
                    image_file.open(path, std::ios::app);
                }

                for (int r = 0; r < i_height; ++r) {
                    for (int c = 0; c < i_width; ++c) {
                        image_file << static_cast<float>(test_set.at(i).at(r * i_width + c)/PIXEL_MAX);
                        if (c < (i_width - 1)) {
                            image_file << ",";
                        }
                    }
                    image_file << endl;
                }
                image_file << endl;
                image_file.close();

                cout << "close file: " << path << endl;
            }
        }
    } else {
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
    if (role == 0)
        proxy = new Party(P1, hport, haddress, cport, caddress);
    else
        proxy = new Party(P2, hport, haddress, cport, caddress);

    cout << "Creating secret shares: " << endl;
    // secret shares: data
    auto ****input = new uint64_t ***[i_number];
    for (int i = 0; i < i_number; ++i) {
        input[i] = new uint64_t **[i_channel];
        for (int j = 0; j < i_channel; ++j) {
            input[i][j] = new uint64_t *[i_height];
            for (uint32_t r = 0; r < i_height; ++r) {
                input[i][j][r] = new uint64_t[i_width];
                for (uint32_t c = 0; c < i_width; ++c) {
                    double pixelValue = test_set.at(i).at(r * i_width + c)/PIXEL_MAX;
                    input[i][j][r][c] = proxy->createShare(pixelValue); // store directly as secret shares
                }
            }
            if (padding > 0) {
                uint64_t padding_value = proxy->createShare(0.0);
                input[i][j] = PAD(input[i][j], i_height, i_width, padding_value, padding);
            }
        }

    }

    // secret shares: bias
    auto **bias = new uint64_t *[layer_number];
    for (int layer = 0; layer < layer_number; ++layer) {
        bias[layer] = proxy->createShare(model_weights[layer + layer_number][0][0], bias_dimensions[layer]);
        delete_model_weights(model_weights, layer + layer_number, i_channel, bias_dimensions[layer]);
    }

    // secret shares: kernel
    k_number = bias_dimensions[0];
    uint32_t k_size = nn_mode == 2 ? k_dim : k_dim * k_dim;
    auto ***kernel = new uint64_t **[i_channel]; // weights of all kernel for first CL

    auto *prediction = new double[i_number];
    double correct = 0, incorrect = 0;
    // CNN INFERENCE PIPELINE
    for (uint32_t image = 0; image < i_number; ++image) {
        cout << "INFERENCE PIPELINE " << image << endl;
        k_number = bias_dimensions[0];
        resetParams();
        uint64_t ***conv;
        uint64_t *prev_layer_res;
        uint64_t **weights;
        //must be set for each inference iteration because kernel is set for each layer.
        for (int i = 0; i < i_channel; ++i) {
            kernel[i] = proxy->createShare(model_weights[curr_layer][i], k_number, k_size);
        }
        if (image == i_number - 1) {
            delete_model_weights(model_weights, curr_layer, i_channel, bias_dimensions[curr_layer]);
        }

        switch (nn_mode) {                              // from here on network architectures differ
            case 0: { // Chameleon
                print2DArray("kernel values", convert2double(REC(proxy, kernel[0], k_number, k_dim * k_dim), k_number, k_dim * k_dim), k_number, k_dim * k_dim);
                print2DArray("input", convert2double(REC(proxy, input[image][0], i_height, i_width), i_height, i_width), i_height, i_width);
                conv = CL(proxy, input[image], i_channel, i_height, i_width, kernel, k_dim, k_number, stride, max_win_height,
                          max_win_width, bias[curr_layer], true);
                delete_image_share(input, image);
                updateParamsAfterCL();

                updateParamsForFCL(bias_dimensions, true);
                // FULLY CONNECTED LAYER
                weights = proxy->createShare(model_weights[curr_layer][0], nodes_out, nodes_in);
                if (image == i_number - 1) {
                    delete_model_weights(model_weights, curr_layer, i_channel, bias_dimensions[curr_layer]);
                }
                prev_layer_res = FCL(proxy, conv[0][0], nodes_in, weights, nodes_out, bias[curr_layer]);
                prev_layer_res = RELU(proxy, prev_layer_res, nodes_out);
                updateParamsForFCL(bias_dimensions, false);
                break;
            }
            case 1: { // SecureNN
                conv = CL(proxy, input[image], i_channel, i_height, i_width, kernel, k_dim, k_number, stride, max_win_height,
                          max_win_width, bias[curr_layer], false);
                delete_image_share(input, image);
                updateParamsAfterCL();
                // PERFORMING CONVOLUTION
                k_number = bias_dimensions[curr_layer];
                kernel = new uint64_t **[i_channel];
                for (uint32_t i = 0; i < i_channel; i++) {
                    kernel[i] = proxy->createShare(model_weights[curr_layer][i], k_number, k_dim * k_dim);
                }
                if (image == i_number - 1) {
                    delete_model_weights(model_weights, curr_layer, i_channel, bias_dimensions[curr_layer]);
                }

                conv = CL(proxy, conv, i_channel, i_height, i_width, kernel, k_dim, k_number, stride, max_win_height,
                          max_win_width, bias[curr_layer], true);
                updateParamsAfterCL();

                // fully connected layer:
                updateParamsForFCL(bias_dimensions, true);
                weights = proxy->createShare(model_weights[curr_layer][0], nodes_out, nodes_in);
                if (image == i_number - 1) {
                    delete_model_weights(model_weights, curr_layer, 0, nodes_out);
                }
                prev_layer_res = FCL(proxy, conv[0][0], nodes_in, weights, nodes_out, bias[curr_layer]);
                prev_layer_res = RELU(proxy, prev_layer_res, nodes_out);
                updateParamsForFCL(bias_dimensions, false);
                break;
            }
            case 2: {
                // CONVOLUTIONAL LAYER (adapted to non-symmetric filter size
                conv = new uint64_t **[k_number];
                for (int k = 0; k < k_number; ++k) {
                    conv[k] = MATVECMUL(proxy, input[image], kernel[k], 1, i_height, i_width);
                    conv[k][0][0] = MAX(proxy, conv[k][0], i_height);
                }
                delete_image_share(input, image);
                updateParamsAfterCL();

                prev_layer_res = FLT(conv, i_height, i_width, i_channel);
                updateParamsForFCL(bias_dimensions, true);
            }
        }
        // last FULLY CONNECTED LAYER
        weights = proxy->createShare(model_weights[curr_layer][0], nodes_out, nodes_in);
        if (image == i_number - 1) {
            delete_model_weights(model_weights, curr_layer, 0, nodes_out);
        }

        uint64_t *output = FCL(proxy, prev_layer_res, nodes_in, weights, nodes_out, bias[curr_layer]);
        prediction[image] = convert2double(REC(proxy, ARGMAX(proxy, output, nodes_out)));
        print1DArray("Input to ARGMAX:", convert2double(REC(proxy, output, nodes_out), nodes_out), nodes_out);
        cout << ": predicted " << prediction[image] << ", correct is " << int(test_label[image]) << endl;

        if (prediction[image] == int(test_label[image])) {
            correct++;
        } else {
            incorrect++;
        }
        double *inference_res = convert2double(REC(proxy, output, nodes_out), nodes_out);
        if (trained_for_MNIST and eval_correctness and proxy->getPRole() == P1) {
            ofstream image_file;
            string path;
            switch (nn_mode) {
                case 0:
                    path = CHAMELEON_CORRECT_PATH;
                    break;
                case 1:
                    path = LENET_CORRECT_PATH;
            }
            path += "eval_secure_" + to_string(image) + ".txt";
            image_file.open(path, std::ios::app);
            if (!image_file) {
                cout << "Error opening file for appending prediction at " << path << endl;

                return 1;
            }
            for (int v = 0; v < nodes_out; ++v) {
                image_file << inference_res[v] << "\t";
            }
            image_file << endl;
            image_file << prediction[image] << endl;
            image_file.close();
        }
        // deleting routine:
        delete[] inference_res;
        delete[] conv;
        delete[] prev_layer_res;
        delete[] output;

        for (int o = 0; o < nodes_out; ++o) {
            delete[] weights[o];
        }
        delete[] weights;
    }
    delete[] model_weights;
    proxy->PrintBytes();
    proxy->piK();

    print1DArray("Prediction: ", prediction, i_number);
    double* ground_truth = new double [i_number];
    for (int i = 0; i < i_number; ++i)
        ground_truth[i] = test_label[i];
    print1DArray("Ground Truth: ", ground_truth, i_number);

    cout << "Accuracy: " << printf("%.2f", correct / i_number) << " (" << printf("%.1f", correct) << "/" << i_number << ")" << endl;
    return 0;
}