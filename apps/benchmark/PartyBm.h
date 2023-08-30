//
// Created by noah on 08/07/22.
//

#ifndef CECILIA_PARTYBM_H
#define CECILIA_PARTYBM_H

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

constexpr int kMinVal = -100;
constexpr int kMaxVal = 100;

class PartyBm : public Party {
public:
     /**
      *
      * @param role either proxy1 or proxy2
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
    explicit PartyBm(
             Role role,
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
        this->matrix_x_ = matrix_x;
        this->matrix_y_ = matrix_y;
        this->vector_length_ = vector_length;
         flattened_matrix_size_ = matrix_x_ * matrix_y_;
        this->flattened_3d_matrix_size_ = matrix_x_ * matrix_y_ * vector_length_;
        this->repeats_ = repeats;
        this->window_size_ = window_size;
        this->kernel_count_ = kernel_count;
        this->kernel_size_ = kernel_size;
        if (role == helper) {
            throw std::invalid_argument("Tried to call the proxy1 and proxy2 constructor with a helper Role");
        }
        //generate random data:
        vector_ = InitialiseVector(vector_length_);
         number_ = vector_[0];
         vector_of_matrices_ = Initialise3dMatrix(matrix_x_, matrix_y_, vector_length_);
         matrix_for_matrix_multiplication_ = Initialise3dMatrix(matrix_y_, matrix_x_, vector_length_);
         vector_of_gram_matrices_ = new uint64_t ** [vector_length_];
        for (int i = 0; i < vector_length_; i++) {
            vector_of_gram_matrices_[i] = RandomOrthogonalMatrix(this, matrix_x);
        }
         flattened_matrix_ = InitialiseVector(flattened_matrix_size_);
         matrix_for_vector_multiplication_ = InitialiseMatrix(matrix_x_, vector_length_);
         double_matrix_ = Random2dData(this, matrix_y_, matrix_x_, 0, 255);
         matrix_kernel_ = Initialise3dMatrix(kernel_size_ * kernel_size_, kernel_count_, vector_length_);
         bias_vector_ = InitialiseVector(kernel_count_);
         flattened_3d_matrix_ = InitialiseVector(flattened_3d_matrix_size_);
         matrix_fcl_weights_ = InitialiseMatrix(flattened_3d_matrix_size_, kernel_count_);
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
    explicit PartyBm(
            uint16_t helper_port,
            const string &helper_ip,
            int vector_length,
            int matrix_x,
            int matrix_y,
            int window_size,
            int kernel_size,
            int kernel_count,
            int repeats
    ) : Party(helper, helper_port, helper_ip, 0, "") {
        matrix_x_ = matrix_x;
        matrix_y_ = matrix_y;
        vector_length_ = vector_length;
        window_size_ = window_size;
        kernel_size_ = kernel_size;
        kernel_count_ = kernel_count;
        repeats_ = repeats;
        flattened_matrix_size_ = matrix_x_ * matrix_y_;
        flattened_3d_matrix_size_ = matrix_x_ * matrix_y_ * vector_length_;
        number_ = 0;
        vector_ = nullptr;
        vector_of_matrices_ = nullptr;
        vector_of_gram_matrices_ = nullptr;
        matrix_for_vector_multiplication_ = nullptr;
        flattened_matrix_ = nullptr;
        double_matrix_ = nullptr;
        matrix_for_matrix_multiplication_ = nullptr;
        matrix_kernel_ = nullptr;
        bias_vector_ = nullptr;
        flattened_3d_matrix_ = nullptr;
        matrix_fcl_weights_ = nullptr;

    }

    ~PartyBm() {
        delete[] vector_;
        delete[] vector_of_matrices_;
        delete[] vector_of_gram_matrices_;
        delete[] matrix_for_matrix_multiplication_;
        delete[] matrix_kernel_;
        delete[] bias_vector_;
        delete[] flattened_3d_matrix_;
    }

    tuple<double, double> Benchmark(const string& function) {
        if (! function_mappings.count(function)) {
            if (GetPRole() == proxy1) {
                cerr << "\nUnknown function: " << function << endl;
            }
            return make_tuple(-1, -1);
        }
        void(PartyBm::*function_pointer)() = function_mappings[function];
        clock_t cpu_start = clock();
        time_point time_start = steady_clock::now();
        for (int i = 0; i < repeats_; i++) {
            (this->*function_pointer)();
        }
        clock_t cpu_end = clock();
        time_point time_end = steady_clock::now();
        long duration = duration_cast<milliseconds>(time_end - time_start).count();
        double cpu_time = ((cpu_end - cpu_start) / 1000.0) / repeats_;
        double real_time = duration / static_cast<double>(repeats_);
        return make_tuple(cpu_time, real_time);
    }


    tuple<int, string*> GetAllFunctionNames() {
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

    void BenchmarkPartialSum() { //TODO ask about d
        PartialSum(this, vector_, vector_length_, 3);
    }

    void BenchmarkMultiply() {
        Multiply(this, vector_, vector_, vector_length_);
    }

    void BenchmarkAucDivide() {
        AucDivide(this, number_, number_);
    }

    void BenchmarkMultiplex() {
        Multiplex(this, vector_, vector_, vector_, vector_length_);
    }

    void BenchmarkMostSignificantBit() {
        MostSignificantBit(this, vector_, vector_length_);
    }

    void BenchmarkModularConversion() {
        ModularConversion(this, vector_, vector_length_);
    }

    void BenchmarkCompare() {
        Compare(this, vector_, vector_, vector_length_);
    }

    void BenchmarkExp() {
        Exp(this, vector_, vector_length_);
    }

    void BenchmarkDotProduct() {
        DotProduct(this, flattened_matrix_, flattened_matrix_, flattened_matrix_size_, matrix_x_);
    }

    void BenchmarkMatrixMatrixMultiply() {
        if (this->GetPRole() == helper) {
            MatrixMatrixMultiply(this, nullptr, nullptr, vector_length_ * matrix_y_ * matrix_x_ * matrix_y_, 0, 0);
        } else {
            MatrixMatrixMultiply(this, vector_of_matrices_, matrix_for_matrix_multiplication_, vector_length_, matrix_y_, matrix_x_,
                                 matrix_y_);
        }

    }

    void BenchmarkMatrixVectorMultiply() {
        if (this->GetPRole() == helper) {
            MatrixVectorMultiply(this, nullptr, nullptr, 0, matrix_y_ * vector_length_ * matrix_x_, 0);
        } else {
            MatrixVectorMultiply(this, vector_of_matrices_, matrix_for_vector_multiplication_, vector_length_, matrix_y_, matrix_x_);
        }
    }

    void BenchmarkReconstruct() {
        if (GetPRole() != helper) {
            Reconstruct(this, vector_, vector_length_);
        }
    }

    void BenchmarkModularInverse() {
        ModularInverse(this, number_);
    }

    void BenchmarkRound() {
        Round(this, vector_, vector_length_);
    }

    void BenchmarkMax() {
        Max(this, flattened_matrix_, flattened_matrix_size_);
    }

    void BenchmarkMaxPerWindow() {
        Max(this, flattened_matrix_, matrix_y_, matrix_x_, window_size_, window_size_);
    }

    void BenchmarkArgMax() {
        ArgMax(this, flattened_matrix_, flattened_matrix_size_);
    }

    void BenchmarkReLU() {
        Relu(this, flattened_matrix_, flattened_matrix_size_);
    }

    void BenchmarkDerivateReLU() {
        DerivativeRelu(this, flattened_matrix_, flattened_matrix_size_);
    }

    void BenchmarkConvolutionalLayer() {
        ConvolutionalLayer(this, vector_of_matrices_, vector_length_, matrix_y_, matrix_x_, matrix_kernel_,
                           kernel_size_, kernel_count_, 2, window_size_, window_size_, bias_vector_, false);
    }

    void BenchmarkFullyConnectedLayer() {
        FullyConnectedLayer(this, flattened_3d_matrix_, flattened_3d_matrix_size_, matrix_fcl_weights_, kernel_count_,
                            bias_vector_);
    }
    void BenchmarkInverseSqrt() {
        InverseSqrt(this, vector_of_gram_matrices_, vector_length_, matrix_x_);
    }

    void BenchmarkGaussianKernel() {
        GaussianKernel(this, vector_of_gram_matrices_, 1, vector_length_, matrix_x_);
    }

private:
    uint64_t* InitialiseVector(size_t size) {
        double* original_vector = Random1dData(this, size, 255, false);
        uint64_t* vector = CreateShare(
                original_vector,
                size);
        delete[] original_vector;
        return vector;
    }

    uint64_t** InitialiseMatrix(size_t x, size_t y) {
        double** vector_of_vectors = Random2dData(this, y, x, 0, 255);
        uint64_t** matrix = CreateShare(
                vector_of_vectors,
                y,
                x
        );
        delete[] vector_of_vectors;
        return matrix;
    }

    uint64_t*** Initialise3dMatrix(size_t x, size_t y, size_t z) {
        auto*** matrix_3d = new uint64_t ** [z];
        double** matrix = Random2dData(this, y, x, 0, 255);
        for (size_t i = 0; i < z; i++) {
            matrix_3d[i] = CreateShare(
                    matrix,
                    y,
                    x
            );
        }
        delete[] matrix;
        return(matrix_3d);
    }

    static std::map<std::string, void(PartyBm::*)()> CreateMap() {
        std::map<std::string, void(PartyBm::*)()> map;
        map["PartialSum"] = &PartyBm::BenchmarkPartialSum;
        map["Multiply"] = &PartyBm::BenchmarkMultiply;
        map["aucDivide"] = &PartyBm::BenchmarkAucDivide;
        map["Multiplex"] = &PartyBm::BenchmarkMultiplex;
        map["MostSignificantBit"] = &PartyBm::BenchmarkMostSignificantBit;
        map["ModularConversion"] = &PartyBm::BenchmarkModularConversion;
        map["Compare"] = &PartyBm::BenchmarkCompare;
        map["Exp"] = &PartyBm::BenchmarkExp;
        map["DotProduct"] = &PartyBm::BenchmarkDotProduct;
        map["MatrixMatrixMultiply"] = &PartyBm::BenchmarkMatrixMatrixMultiply;
        map["MatrixVectorMultiply"] = &PartyBm::BenchmarkMatrixVectorMultiply;
        map["Reconstruct"] = &PartyBm::BenchmarkReconstruct;
        map["ModularInverse"] = &PartyBm::BenchmarkModularInverse;
        map["Relu"] = &PartyBm::BenchmarkReLU;
        map["DerivativeRelu"] = &PartyBm::BenchmarkDerivateReLU;
        map["ConvolutionalLayer"] = &PartyBm::BenchmarkConvolutionalLayer;
        map["FullyConnectedLayer"] = &PartyBm::BenchmarkFullyConnectedLayer;
        map["InverseSqrt"] = &PartyBm::BenchmarkInverseSqrt;
        map["GaussianKernel"] = &PartyBm::BenchmarkGaussianKernel;
        map["ArgMax"] = &PartyBm::BenchmarkArgMax;
        map["Max"] = &PartyBm::BenchmarkMax;
        map["MAXPOOL"] = &PartyBm::BenchmarkMaxPerWindow;
        return map;
    }
    std::map<std::string, void(PartyBm::*)()> function_mappings = CreateMap();
    int matrix_x_, matrix_y_, window_size_, kernel_size_, kernel_count_, vector_length_, repeats_;
    int flattened_matrix_size_, flattened_3d_matrix_size_;
    uint64_t number_;
    uint64_t *vector_, *flattened_matrix_, *bias_vector_, *flattened_3d_matrix_;
    uint64_t **matrix_for_vector_multiplication_, **matrix_fcl_weights_;
    uint64_t ***vector_of_matrices_, ***matrix_for_matrix_multiplication_, ***vector_of_gram_matrices_, ***matrix_kernel_;
    double** double_matrix_;
};




#endif //CECILIA_PARTYBM_H
