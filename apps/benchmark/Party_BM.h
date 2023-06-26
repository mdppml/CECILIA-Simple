//
// Created by noah on 08/07/22.
//

#ifndef CECILIA_PARTY_BM_H
#define CECILIA_PARTY_BM_H

#include <atomic>
#include "../../core/Party.h"
#include "../../utils/connection.h"
#include "../../core/cnn.h"
#include "../../core/rkn.h"
#include "../../core/auc.h"
#include <thread>
#include <iostream>
#include <stdexcept>
#include <chrono>
#include <string>
#include <map>

using namespace std::chrono;

constexpr int MIN_VAL = -100;
constexpr int MAX_VAL = 100;

class Party_BM : public Party {
public:
     /**
      *
      * @param role either P1 or P2
      * @param helper_port
      * @param helper_ip
      * @param p1_port
      * @param p1_ip
      * @param vector_length the vector length for operations involving vectors
      * @param matrix_x the column count for matrices in matrix operations
      * @param matrix_y the row count for matrices in matrix operations
      * @param window_size
      * @param kernel_size
      * @param kernel_count
      * @param repeats
      */
    explicit Party_BM(
            role role,
            uint16_t helper_port,
            const string &helper_ip,
            uint16_t p1_port,
            const string &p1_ip,
            int vector_length,
            int matrix_x,
            int matrix_y,
            int window_size,
            int kernel_size,
            int kernel_count,
            int repeats
            ) : Party(role, helper_port, helper_ip, p1_port, p1_ip) {
        this->_matrix_x = matrix_x;
        this->_matrix_y = matrix_y;
        this->_vector_length = vector_length;
        _flattened_matrix_size = _matrix_x * _matrix_y;
        this->_flattened_3d_matrix_size = _matrix_x * _matrix_y * _vector_length;
        this->_repeats = repeats;
        this->_window_size = window_size;
        this->_kernel_count = kernel_count;
        this->_kernel_size = kernel_size;
        if (role == HELPER) {
            throw std::invalid_argument("Tried to call the P1 and P2 constructor with a helper role");
        }
        //generate random data:
        _vector = init_vector(_vector_length);
        _number = _vector[0];
        _vector_of_matrices = init_3d_matrix(_matrix_x, _matrix_y, _vector_length);
        _matrix_MATMATMUL = init_3d_matrix(_matrix_y, _matrix_x, _vector_length);
        _vector_of_gram_matrices = new uint64_t ** [_vector_length];
        for (int i = 0; i < _vector_length; i++) {
            _vector_of_gram_matrices[i] = randomOrthogonalMatrix(this, matrix_x);
        }
        _flattened_matrix = init_vector(_flattened_matrix_size);
        _matrix_MATVECMUL = init_matrix(_matrix_x, _vector_length);
        _double_matrix = random_2D_data(this, _matrix_y, _matrix_x, 0, 255);
        _matrix_kernel = init_3d_matrix(_kernel_size*_kernel_size, _kernel_count, _vector_length);
        _bias_vector = init_vector(_kernel_count);
        _flattened_3d_matrix = init_vector(_flattened_3d_matrix_size);
        _matrix_fcl_weights = init_matrix(_flattened_3d_matrix_size, _kernel_count);
    }



    /**
      *
      * @param helper_port
      * @param helper_ip
      * @param p1_port
      * @param p1_ip
      * @param vector_length the vector length for operations involving vectors
      * @param matrix_x the column count for matrices in matrix operations
      * @param matrix_y the row count for matrices in matrix operations
      * @param window_size
      * @param kernel_size
      * @param kernel_count
      * @param repeats
      */
    explicit Party_BM(
            uint16_t helper_port,
            const string &helper_ip,
            int vector_length,
            int matrix_x,
            int matrix_y,
            int window_size,
            int kernel_size,
            int kernel_count,
            int repeats
    ) : Party(HELPER, helper_port, helper_ip, 0, "") {
        _matrix_x = matrix_x;
        _matrix_y = matrix_y;
        _vector_length = vector_length;
        _window_size = window_size;
        _kernel_size = kernel_size;
        _kernel_count = kernel_count;
        _repeats = repeats;
        _flattened_matrix_size = _matrix_x * _matrix_y;
        _flattened_3d_matrix_size = _matrix_x * _matrix_y * _vector_length;
        _number = 0;
        _vector = nullptr;
        _vector_of_matrices = nullptr;
        _vector_of_gram_matrices = nullptr;
        _matrix_MATVECMUL = nullptr;
        _flattened_matrix = nullptr;
        _double_matrix = nullptr;
        _matrix_MATMATMUL = nullptr;
        _matrix_kernel = nullptr;
        _bias_vector = nullptr;
        _flattened_3d_matrix = nullptr;
        _matrix_fcl_weights = nullptr;

    }

    ~Party_BM() {
        delete[] _vector;
        delete[] _vector_of_matrices;
        delete[] _vector_of_gram_matrices;
        delete[] _matrix_MATMATMUL;
        delete[] _matrix_kernel;
        delete[] _bias_vector;
        delete[] _flattened_3d_matrix;
    }

    tuple<double, double> benchmark(const string& function) {
        if (! function_mappings.count(function)) {
            if (getPRole() == P1) {
                cerr << "\nUnknown function: " << function << endl;
            }
            return make_tuple(-1, -1);
        }
        void(Party_BM::*function_pointer)() = function_mappings[function];
        clock_t cpu_start = clock();
        time_point time_start = steady_clock::now();
        for (int i = 0; i < _repeats; i++) {
            (this->*function_pointer)();
        }
        clock_t cpu_end = clock();
        time_point time_end = steady_clock::now();
        long duration = duration_cast<milliseconds>(time_end - time_start).count();
        double cpu_time = ((cpu_end - cpu_start) / 1000.0) / _repeats;
        double real_time = duration / static_cast<double>(_repeats);
        return make_tuple(cpu_time, real_time);
    }


    tuple<int, string*> get_all_function_names() {
        int count = 0;
        for (const auto &pair : function_mappings) {
            count +=1;
        }
        auto function_names = new string[count];
        int i = 0;
        for (const auto &pair : function_mappings) {
            function_names[i] = pair.first;
            i++;
        }
        return make_tuple(count, function_names);
    }

    void partial_summation() { //TODO ask about d
        PartialSum(this, _vector, _vector_length, 3);
    }

    void multiplication() {
        Multiply(this, _vector, _vector, _vector_length);
    }

    void division_cnn() {
        Divide(this, _number, _number);
    }

    /*void division_auc() {
        MDIVISION(this, _vector, _vector, _vector_length);
    }*/

    void mmux() {
        Multiplex(this, _vector, _vector, _vector, _vector_length);
    }

    void most_significant_bit() {
        MostSignificantBit(this, _vector, _vector_length);
    }

    /*void most_significant_bit_auc() {
        AUCMSB(this, _vector, _vector_length);
    }*/

    void modular_conversion() {
        ModularConversion(this, _vector, _vector_length);
    }

    void compare() {
        Compare(this, _vector, _vector, _vector_length);
    }

    void exponential() {
        Exp(this, _vector, _vector_length);
    }

    void dot_product() {
        DotProduct(this, _flattened_matrix, _flattened_matrix, _flattened_matrix_size, _matrix_x);
    }

    void matrix_matrix_multiplication() {
        if (this->getPRole() == HELPER) {
            MatrixMatrixMultiply(this, nullptr, nullptr, _vector_length * _matrix_y * _matrix_x * _matrix_y, 0, 0);
        } else {
            MatrixMatrixMultiply(this, _vector_of_matrices, _matrix_MATMATMUL, _vector_length, _matrix_y, _matrix_x,
                                 _matrix_y);
        }

    }

    void matrix_vector_multiplication() {
        if (this->getPRole() == HELPER) {
            MatrixVectorMultiply(this, nullptr, nullptr, 0, _matrix_y * _vector_length * _matrix_x, 0);
        } else {
            MatrixVectorMultiply(this, _vector_of_matrices, _matrix_MATVECMUL, _vector_length, _matrix_y, _matrix_x);
        }
    }

    void reconstruct() {
        if (getPRole() != HELPER) {
            Reconstruct(this, _vector, _vector_length);
        }
    }

    void modular_inverse() {
        ModularInverse(this, _number);
    }

    void round() {
        MRound(this, _vector, _vector_length);
    }

    void max() {
        Max(this, _flattened_matrix, _flattened_matrix_size);
    }

    void max_per_window() {
        Max(this, _flattened_matrix, _matrix_y, _matrix_x, _window_size, _window_size);
    }

    void argmax() {
        ArgMax(this, _flattened_matrix, _flattened_matrix_size);
    }

    void relu() {
        ReLU(this, _flattened_matrix, _flattened_matrix_size);
    }

    void derivative_relu() {
        DerivativeReLU(this, _flattened_matrix, _flattened_matrix_size);
    }

    void conv_layer() {
        ConvolutionalLayer(this, _vector_of_matrices, _vector_length, _matrix_y, _matrix_x, _matrix_kernel,
                           _kernel_size, _kernel_count, 2, _window_size, _window_size, _bias_vector, false);
    }

    void fully_connected_layer() {
        FullyConnectedLayer(this, _flattened_3d_matrix, _flattened_3d_matrix_size, _matrix_fcl_weights, _kernel_count,
                            _bias_vector);
    }
    void inverse_square_root() {
        InverseSqrt(this, _vector_of_gram_matrices, _vector_length, _matrix_x);
    }

    void createShare_BM() {
        if (getPRole() != HELPER) {
            createShare(_double_matrix, _matrix_y, _matrix_x);
        }
    }

    void gram_matrix_to_kernel_matrix() {
        GaussianKernel(this, _vector_of_gram_matrices, 1, _vector_length, _matrix_x);
    }

private:
    uint64_t* init_vector(size_t size) {
        double* original_vector = random_1D_data(this, size, 255, false);
        uint64_t* vector = createShare(
                original_vector,
                size);
        delete[] original_vector;
        return vector;
    }

    uint64_t** init_matrix(size_t x, size_t y) {
        double** vector_of_vectors = random_2D_data(this, y, x, 0, 255);
        uint64_t** matrix = createShare(
                vector_of_vectors,
                y,
                x
        );
        delete[] vector_of_vectors;
        return matrix;
    }

    uint64_t*** init_3d_matrix(size_t x, size_t y, size_t z) {
        auto*** matrix_3d = new uint64_t ** [z];
        double** matrix = random_2D_data(this, y, x, 0, 255);
        for (size_t i = 0; i < z; i++) {
            matrix_3d[i] = createShare(
                    matrix,
                    y,
                    x
            );
        }
        delete[] matrix;
        return(matrix_3d);
    }

    static std::map<std::string, void(Party_BM::*)()> create_map() {
        std::map<std::string, void(Party_BM::*)()> map;
        map["PartialSum"] = &Party_BM::partial_summation;
        map["Multiply"] = &Party_BM::multiplication;
        map["Divide"] = &Party_BM::division_cnn;
        //map["MDIVISION"] = &Party_BM::division_auc;
        map["Multiplex"] = &Party_BM::mmux;
        map["MostSignificantBit"] = &Party_BM::most_significant_bit;
        //map["AUCMSB"] = &Party_BM::most_significant_bit_auc;
        map["ModularConversion"] = &Party_BM::modular_conversion;
        map["Compare"] = &Party_BM::compare;
        map["Exp"] = &Party_BM::exponential;
        map["DotProduct"] = &Party_BM::dot_product;
        map["MatrixMatrixMultiply"] = &Party_BM::matrix_matrix_multiplication;
        map["MatrixVectorMultiply"] = &Party_BM::matrix_vector_multiplication;
        map["Reconstruct"] = &Party_BM::reconstruct;
        map["ModularInverse"] = &Party_BM::modular_inverse;
        //map["MRound"] = &Party_BM::round;
        map["ReLU"] = &Party_BM::relu;
        map["DerivativeReLU"] = &Party_BM::derivative_relu;
        map["ConvolutionalLayer"] = &Party_BM::conv_layer;
        map["FullyConnectedLayer"] = &Party_BM::fully_connected_layer;
        map["InverseSqrt"] = &Party_BM::inverse_square_root;
        map["createShare"] = &Party_BM::createShare_BM;
        map["GaussianKernel"] = &Party_BM::gram_matrix_to_kernel_matrix;
        map["ArgMax"] = &Party_BM::argmax;
        map["Max"] = &Party_BM::max;
        map["MAXPOOL"] = &Party_BM::max_per_window;
        return map;
    }
    std::map<std::string, void(Party_BM::*)()> function_mappings = create_map();
    int _matrix_x, _matrix_y, _window_size, _kernel_size, _kernel_count, _vector_length, _repeats;
    int _flattened_matrix_size, _flattened_3d_matrix_size;
    uint64_t _number;
    uint64_t *_vector, *_flattened_matrix, *_bias_vector, *_flattened_3d_matrix;
    uint64_t **_matrix_MATVECMUL, **_matrix_fcl_weights;
    uint64_t ***_vector_of_matrices, ***_matrix_MATMATMUL, ***_vector_of_gram_matrices, ***_matrix_kernel;
    double** _double_matrix;
};




#endif //CECILIA_PARTY_BM_H
