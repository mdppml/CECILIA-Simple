//
// Created by Ali Burak on 29.09.23.
//



#include <cstdlib>
#include <iostream>
#include <fstream>
#include <chrono>
#include "../../core/core.h"

using namespace std;

void printHead(double *a, int size) {
    for(int i = 0; i < size; i++) {
        cout << a[i] << endl;
    }
}

void MNWE_test(Party *proxy) {
    ofstream txt;


    cout << "hello world" << endl;

    auto n_matrices = 2;
    auto n1 = 5;
    auto m1 = 3; //100--- creation works, loop fails
    auto n2 = m1;
    auto m2 = 4;


    auto headsize = 5;


    double **mat1 = new double *[n1];
    double **mat2 = new double *[n1];


    double **mat3 = new double *[m1];
    double **mat4 = new double *[m1];

    double *vec1 = new double[m1];
    double *vec2 = new double[m1];

    // matrices
    for (int row = 0; row < n1; ++row) {
        mat1[row] = new double[m1];
        for (int col = 0; col < m1; ++col) {
            mat1[row][col] = 1;
        }
    }

    for (int row = 0; row < n1; ++row) {
        mat2[row] = new double[m1];
        for (int col = 0; col < m1; ++col) {
            mat2[row][col] = 2.5;
        }
    }

    for (int row = 0; row < m1; ++row) {
        mat3[row] = new double[m2];
        for (int col = 0; col < m2; ++col) {
            mat3[row][col] = 3.6;
        }
    }

    for (int row = 0; row < m1; ++row) {
        mat4[row] = new double[m2];
        for (int col = 0; col < m2; ++col) {
            mat4[row][col] = 3;
        }
    }


    // vectors
    for (int row = 0; row < m1; ++row) {
        vec1[row] = 2.5;
    }
    for (int row = 0; row < m1; ++row) {
        vec2[row] = 3.5;
    }


    uint64_t a = proxy->CreateShare(2.5);
    cout << "a: " << a << endl;
    cout << "Rec a: " << ConvertToDouble(Reconstruct(proxy, a)) << endl;
    uint64_t c = proxy->CreateShare(2.5);
    cout << "c: " << c << endl;
    cout << "Rec c: " << ConvertToDouble(Reconstruct(proxy, c)) << endl;
    uint64_t b = ConvertToUint64(3.5);
    cout << "b in our ring: " << b << endl;
    cout << "b in real number domain: " << ConvertToDouble(b) << endl;
    cout << "a * b = " << ConvertToDouble(Reconstruct(proxy, LocalMultiply(a, b))) << endl;
    proxy->SendBytes(coreMultiply);
    cout << "a * c = " << ConvertToDouble(Reconstruct(proxy, Multiply(proxy, a, c))) << endl;


    cout << "Created data" << endl;

    auto mat1_share = proxy->CreateShare(mat1, n1, m1);
    auto mat2_share = proxy->CreateShare(mat2, n1, m1);
    auto mat3_share = proxy->CreateShare(mat3, n2, m2);
    auto mat4_share = proxy->CreateShare(mat4, n2, m2);

    auto vec1_share = proxy->CreateShare(vec1, m1);
    auto vec2_share = proxy->CreateShare(vec2, m1);

    auto rec_mat1 = Reconstruct(proxy, mat1_share, n1, m1);
    auto rec_mat2 = Reconstruct(proxy, mat2_share, n2, m2);

    auto uint_rec_vec1 = Reconstruct(proxy, vec1_share, m1);
    auto rec_vec1 = ConvertToDouble(uint_rec_vec1, m1);
    auto rec_vec2 = ConvertToDouble(Reconstruct(proxy, vec2_share, m2), m2);

    for(int i = 0; i < m1; i++) {
        cout << vec1_share[i] << " ";
    }
    cout << endl;

    for(int i = 0; i < m1; i++) {
        cout << uint_rec_vec1[i] << " ";
    }
    cout << endl;

    for(int i = 0; i < m1; i++) {
        cout << rec_vec1[i] << " ";
    }
    cout << endl;

    uint64_t ***mats1 = new uint64_t **[2];
    mats1[0] = mat1_share;
    mats1[1] = mat2_share;

    uint64_t ***mats2 = new uint64_t **[2];
    mats2[0] = mat3_share;
    mats2[1] = mat4_share;

    uint64_t **vecs = new uint64_t *[2];
    vecs[0] = vec1_share;
    vecs[1] = vec2_share;

// **************************************
//    cout << "MatrixMatrixMultiply test" << endl;
//    if(m1 != n2) {
//        cout << "m1 cannot be different from n2 for MatrixMatrixMultiply" << endl;
//        return;
//    }
//    uint32_t *params = new uint32_t[3];
//    uint32_t size = n1 * n2 * m2; // size of the vector ---- what is this?
//    params[0] = n1;
//    params[1] = n2;
//    params[2] = m2;
//    proxy->SendBytes(coreMatrixMatrixMultiply, params, 3);
//    uint64_t **mul_share = MatrixMatrixMultiply(proxy, mat1_share, mat2_share, n1, n2, m2);
//
//    cout << "Multiplication successfull" << endl;
//
//    auto recov_gram_tmp = Reconstruct(proxy, mul_share, n1, m2);
//    double **recov_gram = ConvertToDouble(recov_gram_tmp, n1, m2);
//
//    cout << "Recovered mul (head upper left): " << endl;
//    for (int row=0; row < ((n1 < headsize) ? n1 : headsize); ++row){
//        for (int col=0; col < ((m2 < headsize) ? m2 : headsize); ++col){
//            cout << recov_gram[row][col]<< " ";
//        }
//        cout << endl;
//    }
//
//    cout << "Ground truth" << endl;
//    auto gt_matmatmul = ConvertToDouble(LocalMatrixMatrixMultiply(rec_mat1, rec_mat2, n1, m1, m2), n1, m2);
//    for (int row=0; row < ((n1 < headsize) ? n1 : headsize); ++row){
//        for (int col=0; col < ((m2 < headsize) ? m2 : headsize); ++col){
//            cout << gt_matmatmul[row][col]<< " ";
//        }
//        cout << endl;
//    }

// **************************************
    cout << "VectorizedMatrixMatrixMultiply test" << endl;
    if(m1 != n2) {
        cout << "m1 cannot be different from n2 for MatrixMatrixMultiply" << endl;
        return;
    }
    uint32_t *params = new uint32_t[3];
    uint32_t size = n1 * n2 * m2; // size of the vector ---- what is this?
    params[0] = 2;
    params[1] = n1;
    params[2] = m1;
    params[3] = m2;
    proxy->SendBytes(coreVectorisedMatrixMatrixMultiply, params, 4);
    uint64_t ***mul_share = MatrixMatrixMultiply(proxy, mats1, mats2, 2, n1, m1, m2);

    cout << "Multiplication successfull" << endl;

    auto recov_gram_tmp = Reconstruct(proxy, mul_share, 2, n1, m2);
    double ***recov_gram = new double**[2];
    recov_gram[0] = ConvertToDouble(recov_gram_tmp[0], n1, m2);
    recov_gram[1] = ConvertToDouble(recov_gram_tmp[1], n1, m2);

    cout << "Recovered mul (head upper left): " << endl;
    for(int n = 0; n < 2; n++) {
        for (int row=0; row < ((n1 < headsize) ? n1 : headsize); ++row){
            for (int col=0; col < ((m2 < headsize) ? m2 : headsize); ++col){
                cout << recov_gram[n][row][col]<< " ";
            }
            cout << endl;
        }
        cout << "=================" << endl;
    }

//    cout << "Ground truth" << endl;
//    auto gt_matmatmul = ConvertToDouble(LocalMatrixMatrixMultiply(rec_mat1, rec_mat2, n1, m1, m2), n1, m2);
//    for (int row=0; row < ((n1 < headsize) ? n1 : headsize); ++row){
//        for (int col=0; col < ((m2 < headsize) ? m2 : headsize); ++col){
//            cout << gt_matmatmul[row][col]<< " ";
//        }
//        cout << endl;
//    }

// **************************************
//    cout << "Single MatrixVectorMultiply test" << endl;
//    uint32_t *params = new uint32_t[2];
//    uint32_t size = n1 * n2 * m2; // size of the vector ---- what is this?
//    params[0] = n1;
//    params[1] = m1;
//    proxy->SendBytes(coreMatrixVectorMultiply, params, 2);
//    uint64_t *mul_share = MatrixVectorMultiply(proxy, mat1_share, vec1_share, n1, m1);
//
//    cout << "Multiplication successfull" << endl;
//
//    auto recov_vec_tmp = Reconstruct(proxy, mul_share,m1);
//    double *recov_vec = ConvertToDouble(recov_vec_tmp, m1);
//
//    cout << "Recovered mul (head): " << endl;
//    for (int row=0; row < ((m1 < headsize) ? m1 : headsize); ++row){
//        cout << recov_vec[row] << endl;
//    }

// **************************************
//    cout << "Vectorized MatrixVectorMultiply test" << endl;
//    uint32_t *params = new uint32_t[3];
//    uint32_t size = n1 * n2 * m2; // size of the vector ---- what is this?
//    params[0] = 2;
//    params[1] = n1;
//    params[2] = m1;
//    proxy->SendBytes(coreVectorisedMatrixVectorMultiply, params, 3);
//    uint64_t **mul_share = MatrixVectorMultiply(proxy, mats, vecs, n_matrices, n1, m1);
//
//    cout << "Multiplication successfull" << endl;
//
//    auto recov_vec_tmp = Reconstruct(proxy, mul_share, n_matrices, n1);
//    double **recov_vecs = ConvertToDouble(recov_vec_tmp, n_matrices, n1);
//
//    cout << "Recovered mul (head upper left): " << endl;
//    for(int mat = 0; mat < n_matrices; mat++) {
//        for (int row=0; row < ((n1 < headsize) ? n1 : headsize); ++row){
//            cout << recov_vecs[mat][row] << endl;
//        }
//        cout << endl;
//    }




//    cout << "Ground truth" << endl;
//    auto gt_matmatmul = ConvertToDouble(LocalMatrixVectorMultiply(rec_mat1, rec_mat2, n1, m1, m2), n1, m2);
//    for (int row=0; row < ((n1 < headsize) ? n1 : headsize); ++row){
//        for (int col=0; col < ((m2 < headsize) ? m2 : headsize); ++col){
//            cout << gt_matmatmul[row][col]<< " ";
//        }
//        cout << endl;
//    }

//    delete [] mat1;
//    delete [] mat2;
//    delete [] mat1_share;
//    delete [] mat2_share;
//    delete [] recov_gram_tmp;
//    delete [] recov_gram;
}



int main(int argc, char* argv[]) {
    uint8_t role = atoi(argv[1]);
    uint16_t cport = atoi(argv[2]);
    string caddress(argv[3]);
    uint16_t hport = atoi(argv[4]);
    string haddress(argv[5]);

    cout << "Starting something" << endl;

    Party *proxy;
    if (role == 0)
        proxy = new Party(proxy1, hport, haddress, cport, caddress);
    else
        proxy = new Party(proxy2, hport, haddress, cport, caddress);

    cout << "Party creation is done" << endl;

    auto start = chrono::high_resolution_clock::now();
    MNWE_test(proxy);
    auto end = chrono::high_resolution_clock::now();
    double time_taken =
            chrono::duration_cast<chrono::nanoseconds>(end - start).count();
    time_taken *= 1e-9;
    cout<<time_taken<<endl;


    proxy->SendBytes(coreEnd);
    //proxy->PrintBytes();

    return 0;
}
