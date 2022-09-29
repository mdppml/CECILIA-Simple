//
// Created by noah on 08/07/22.
//

#ifndef CECILIA_PARTY_BM_H
#define CECILIA_PARTY_BM_H

#include <atomic>
#include "../../core/Party.h"
#include "../../utils/connection.h"
#include "../../core/auc.h"
#include "../../core/cnn.h"
#include "../../core/rkn.h"
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
    /** This constructor is used for P1 and P2.
     *
     * @param role either P1 or P2
     * @param count the vector length for all vectorised operations (i.e. how often to apply the operation)
     * @param vector_length the vector length for all operations that work on vectors
     * @param matrix_x the column count for matrices in matrix operations
     * @param matrix_y the row count for matrices in matrix operations
     * @param helper_port
     * @param helper_ip
     * @param p1_port
     * @param p1_ip
     */
    explicit Party_BM(
            role role,
            uint16_t helper_port,
            const string &helper_ip,
            uint16_t p1_port,
            const string &p1_ip,
            int count,
            int vector_length,
            int matrix_x,
            int matrix_y,
            int window_size,
            int cycle_count
            ) : Party(role, helper_port, helper_ip, p1_port, p1_ip) {
        this->_matrix_x = matrix_x;
        this->_matrix_y = matrix_y;
        _flattened_matrix_size = _matrix_x * _matrix_y;
        this->_count = count;
        this->_cycle_count = cycle_count;
        this->_vector_length = vector_length;
        this->_window_size = window_size;
        if (role == HELPER) {
            throw std::invalid_argument("Tried to call the P1 and P2 constructor with a helper role");
        }
        //generate random data:
        double* vector = random_1D_data(this, _count, 255, false);
        _vector = createShare(
                vector,
                _vector_length);
        _number = _vector[0];
        double** vector_of_vectors = random_2D_data(this, _count, _vector_length, 0, 255);
        _vector_of_vectors = createShare(
                vector_of_vectors,
                _count,
                _vector_length
            );
        _vector_of_matrices = new uint64_t ** [_count];
        for (int i = 0; i < _count; i++) {
            _vector_of_matrices[i] = createShare(
                    random_2D_data(this, _matrix_y, _matrix_x, 0, 255),
                    _matrix_y,
                    _matrix_x
                );
        }
        _vector_of_gram_matrices = new uint64_t ** [_count];
        for (int i = 0; i < _count; i++) {
            _vector_of_gram_matrices[i] = randomOrthogonalMatrix(this, matrix_x);
        }
        _flattened_vector_of_vectors_size = _vector_length * _count;
        delete[] vector;
        vector = random_1D_data(this, _flattened_vector_of_vectors_size, 255, false);
        _flattened_vector_of_vectors = createShare(
                vector,
                _flattened_vector_of_vectors_size
            );
        delete[] vector;
        vector = random_1D_data(this, _flattened_matrix_size, 255, false);
        _flattened_matrix = createShare(
                vector,
                _flattened_matrix_size
            );
        double** matrix = random_2D_data(this, _count, matrix_y, 0, 255);
        _matrix_MATVECMUL = createShare(
                matrix,
                _count,
                matrix_x
            );
        _double_matrix = random_2D_data(this, _matrix_y, _matrix_x, 0, 255);
        delete[] vector;
        delete[] vector_of_vectors;
        delete[] matrix;
    }

    /** This constructor is exclusively used for the helper.
     *
     * @param count
     * @param vector_length
     * @param matrix_x
     * @param matrix_y
     * @param helper_port
     * @param helper_ip
     */
    explicit Party_BM(
            uint16_t helper_port,
            const string &helper_ip,
            int count,
            int vector_length,
            int matrix_x,
            int matrix_y,
            int window_size,
            int cycle_count
    ) : Party(HELPER, helper_port, helper_ip, 0, "") {
        this->_matrix_x = matrix_x;
        this->_matrix_y = matrix_y;
        this->_count = count;
        this->_vector_length = vector_length;
        this->_window_size = window_size;
        this->_cycle_count = cycle_count;
        _flattened_vector_of_vectors_size = _vector_length * _count;
        _flattened_matrix_size = _matrix_x * _matrix_y;
        _number = 0;
        _vector = nullptr;
        _vector_of_vectors = nullptr;
        _vector_of_matrices = nullptr;
        _vector_of_gram_matrices = nullptr;
        _flattened_vector_of_vectors = nullptr;
        _matrix_MATVECMUL = nullptr;
        _flattened_matrix = nullptr;
        _double_matrix = nullptr;
    }

    ~Party_BM() {
        delete[] _vector;
        delete[] _vector_of_vectors;
        delete[] _vector_of_matrices;
        delete[] _flattened_vector_of_vectors;
        delete[] _vector_of_gram_matrices;
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
        for (int i = 0; i < _cycle_count; i++) {
            (this->*function_pointer)();
        }
        clock_t cpu_end = clock();
        time_point time_end = steady_clock::now();
        long duration = duration_cast<milliseconds>(time_end - time_start).count();
        double cpu_time = ((cpu_end - cpu_start) / 1000.0) / _cycle_count;
        double real_time = duration / static_cast<double>(_cycle_count);
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
        PSUM(this, _vector, _count, 3);
    }

    void multiplication() {
        MUL(this, _vector, _vector, _count);
    }

    void partial_multiplication() {
        PMUL(this, _vector, _vector, _count);
    }

    void division_cnn() {
        DIV(this, _number, _number);
    }

    void division_auc() {
        MDIVISION(this, _vector, _vector, _count);
    }

    void mmux() {
        MUX(this, _vector, _vector, _vector, _count);
    }

    void most_significant_bit() {
        MSB(this, _vector, _count);
    }

    void most_significant_bit2() {
        MSBv2(this, _vector, _count);
    }

    void most_significant_bit_auc() {
        AUCMSB(this, _vector, _count);
    }

    void modular_conversion() {
        MOC(this, _vector, _count);
    }

    void compare() {
        CMP(this, _vector, _vector, _count);
    }

    void exponential() {
        EXP(this, _vector, _count);
    }

    void dot_product() {
        DP(this, _flattened_vector_of_vectors, _flattened_vector_of_vectors, _flattened_vector_of_vectors_size, _vector_length);
    }

    void matrix_matrix_multiplication() {
        MATMATMUL(this, _vector_of_matrices, _vector_of_matrices, _count, _matrix_y, _matrix_x, _matrix_x);
    }

    void matrix_vector_multiplication() {
        MATVECMUL(this, _vector_of_matrices, _matrix_MATVECMUL, _count, _matrix_y, _matrix_x);
    }

    void reconstruct() {
        if (getPRole() != HELPER) {
            REC(this, _vector, _count);
        }
    }

    void modular_inverse() {
        MDI(this, _number);
    }

    void round() {
        MRound(this, _vector, _count);
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
        RELU(this, _vector, _count);
    }

    void derivative_relu() {
        DRELU(this, _vector, _count);
    }

    void inverse_square_root() {
        INVSQRT(this, _vector_of_gram_matrices, _count, _matrix_x);
    }

    void createShare_BM() {
        if (getPRole() != HELPER) {
            createShare(_double_matrix, _matrix_y, _matrix_x);
        }
    }

    void gram_matrix_to_kernel_matrix() {
        GM2KM(this, _vector_of_gram_matrices, 1, _count, _matrix_x);
    }

private:
    static std::map<std::string, void(Party_BM::*)()> create_map() {
        std::map<std::string, void(Party_BM::*)()> map;
        map["PSUM"] = &Party_BM::partial_summation;
        map["MUL"] = &Party_BM::multiplication;
        map["PMUL"] = &Party_BM::partial_multiplication;
        map["DIV"] = &Party_BM::division_cnn;
        map["MDIVISION"] = &Party_BM::division_auc;
        map["MUX"] = &Party_BM::mmux;
        map["MSB"] = &Party_BM::most_significant_bit;
        map["MSBv2"] = &Party_BM::most_significant_bit2;
        map["AUCMSB"] = &Party_BM::most_significant_bit_auc;
        map["MOC"] = &Party_BM::modular_conversion;
        map["CMP"] = &Party_BM::compare;
        map["EXP"] = &Party_BM::exponential;
        map["DP"] = &Party_BM::dot_product;
        map["MATMATMUL"] = &Party_BM::matrix_matrix_multiplication;
        map["MATVECMUL"] = &Party_BM::matrix_vector_multiplication;
        map["REC"] = &Party_BM::reconstruct;
        //map["MDI"] = &Party_BM::modular_inverse;
        map["MRound"] = &Party_BM::round;
        map["RELU"] = &Party_BM::relu;
        map["DRELU"] = &Party_BM::derivative_relu;
        map["INVSQRT"] = &Party_BM::inverse_square_root;
        map["createShare"] = &Party_BM::createShare_BM;
        map["GM2KM"] = &Party_BM::gram_matrix_to_kernel_matrix;
        map["ARGMAX"] = &Party_BM::argmax;
        map["MAX"] = &Party_BM::max;
        map["MAXPOOL"] = &Party_BM::max_per_window;
        return map;
    }
    std::map<std::string, void(Party_BM::*)()> function_mappings = create_map();
    int _matrix_x, _matrix_y, _window_size, _count, _vector_length, _flattened_vector_of_vectors_size, _cycle_count;
    int _flattened_matrix_size;
    uint64_t _number;
    uint64_t *_vector, *_flattened_vector_of_vectors, *_flattened_matrix;
    uint64_t **_vector_of_vectors, **_matrix_MATVECMUL;
    uint64_t ***_vector_of_matrices, ***_vector_of_gram_matrices;
    double** _double_matrix;
};




#endif //CECILIA_PARTY_BM_H