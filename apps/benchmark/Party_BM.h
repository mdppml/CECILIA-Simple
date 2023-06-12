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
      * @param p0_port
      * @param p0_ip
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
            uint16_t p0_port,
            const string &p0_ip,
            int vector_length,
            int matrix_x,
            int matrix_y,
            int window_size,
            int kernel_size,
            int kernel_count,
            int repeats
            ) : Party(role, helper_port, helper_ip, p0_port, p0_ip) {
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
            throw std::invalid_argument("Tried to call the P0 and P1 constructor with a helper role");
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
        PSUM(this, _vector, _vector_length, 3);
    }

    void multiplication() {
        MUL(this, _vector, _vector, _vector_length);
    }

    void division_cnn() {
        DIV(this, _number, _number);
    }

    /*void division_auc() {
        MDIVISION(this, _vector, _vector, _vector_length);
    }*/

    void mmux() {
        MUX(this, _vector, _vector, _vector, _vector_length);
    }

    void most_significant_bit() {
        MSB(this, _vector, _vector_length);
    }

    /*void most_significant_bit_auc() {
        AUCMSB(this, _vector, _vector_length);
    }*/

    void modular_conversion() {
        MOC(this, _vector, _vector_length);
    }

    void compare() {
        CMP(this, _vector, _vector, _vector_length);
    }

    void exponential() {
        EXP(this, _vector, _vector_length);
    }

    void dot_product() {
        DP(this, _flattened_matrix, _flattened_matrix, _flattened_matrix_size, _matrix_x);
    }

    void matrix_matrix_multiplication() {
        if (this->getPRole() == HELPER) {
            MATMATMUL(this, nullptr, nullptr, _vector_length * _matrix_y * _matrix_x * _matrix_y, 0, 0);
        } else {
            MATMATMUL(this, _vector_of_matrices, _matrix_MATMATMUL, _vector_length, _matrix_y, _matrix_x, _matrix_y);
        }

    }

    void matrix_vector_multiplication() {
        if (this->getPRole() == HELPER) {
            MATVECMUL(this, nullptr, nullptr, 0, _matrix_y * _vector_length * _matrix_x, 0);
        } else {
            MATVECMUL(this, _vector_of_matrices, _matrix_MATVECMUL, _vector_length, _matrix_y, _matrix_x);
        }
    }

    void reconstruct() {
        if (getPRole() != HELPER) {
            REC(this, _vector, _vector_length);
        }
    }

    void modular_inverse() {
        MDI(this, _number);
    }

    void round() {
        MRound(this, _vector, _vector_length);
    }

    void max() {
        MAX(this, _flattened_matrix,_flattened_matrix_size);
    }

    void max_per_window() {
        MAX(this, _flattened_matrix, _matrix_y, _matrix_x, _window_size, _window_size);
    }

    void argmax() {
        ARGMAX(this, _flattened_matrix, _flattened_matrix_size);
    }

    void relu() {
        RELU(this, _flattened_matrix, _flattened_matrix_size);
    }

    void derivative_relu() {
        DRELU(this, _flattened_matrix, _flattened_matrix_size);
    }

    void conv_layer() {
        CL(this, _vector_of_matrices, _vector_length, _matrix_y, _matrix_x, _matrix_kernel, _kernel_size, _kernel_count, 2, _window_size, _window_size, _bias_vector, false);
    }

    void fully_connected_layer() {
        FCL(this, _flattened_3d_matrix, _flattened_3d_matrix_size, _matrix_fcl_weights, _kernel_count, _bias_vector);
    }
    void inverse_square_root() {
        INVSQRT(this, _vector_of_gram_matrices, _vector_length, _matrix_x);
    }

    void createShare_BM() {
        if (getPRole() != HELPER) {
            createShare(_double_matrix, _matrix_y, _matrix_x);
        }
    }

    void gram_matrix_to_kernel_matrix() {
        GM2KM(this, _vector_of_gram_matrices, 1, _vector_length, _matrix_x);
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
        map["PSUM"] = &Party_BM::partial_summation;
        map["MUL"] = &Party_BM::multiplication;
        map["DIV"] = &Party_BM::division_cnn;
        //map["MDIVISION"] = &Party_BM::division_auc;
        map["MUX"] = &Party_BM::mmux;
        map["MSB"] = &Party_BM::most_significant_bit;
        //map["AUCMSB"] = &Party_BM::most_significant_bit_auc;
        map["MOC"] = &Party_BM::modular_conversion;
        map["CMP"] = &Party_BM::compare;
        map["EXP"] = &Party_BM::exponential;
        map["DP"] = &Party_BM::dot_product;
        map["MATMATMUL"] = &Party_BM::matrix_matrix_multiplication;
        map["MATVECMUL"] = &Party_BM::matrix_vector_multiplication;
        map["REC"] = &Party_BM::reconstruct;
        map["MDI"] = &Party_BM::modular_inverse;
        //map["MRound"] = &Party_BM::round;
        map["RELU"] = &Party_BM::relu;
        map["DRELU"] = &Party_BM::derivative_relu;
        map["CL"] = &Party_BM::conv_layer;
        map["FCL"] = &Party_BM::fully_connected_layer;
        map["INVSQRT"] = &Party_BM::inverse_square_root;
        map["createShare"] = &Party_BM::createShare_BM;
        map["GM2KM"] = &Party_BM::gram_matrix_to_kernel_matrix;
        map["ARGMAX"] = &Party_BM::argmax;
        map["MAX"] = &Party_BM::max;
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
