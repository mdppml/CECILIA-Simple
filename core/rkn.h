//
// Created by aburak on 24.03.22.
//

#ifndef CECILIA_RKN_H
#define CECILIA_RKN_H

#include "core.h"
#include "../utils/flib.h"
//#include <Eigen/Eigen>

//using namespace Eigen;

uint64_t** randomOrthogonalMatrix(Party* proxy, uint32_t size) {
    /*
     * Generate a size-by-size random unit-length orthogonal matrix in P1 and matrix of zeros in P2
     * The solution is based on this: https://scicomp.stackexchange.com/a/34974
     *
     * Input(s)
     * size: the dimension of the matrix, which is size-by-size
     *
     * Output(s)
     * Returns a size-by-size matrix
     */
/*    if( proxy->getPRole() == P1 ||  proxy->getPRole() == P2) {
        // generate a vector of random values
        double *rand_vals = new double[size * size];
        for (uint32_t i = 0; i < size * size; i++) {
            // I guess there is no need for both proxies to know the orthogonal matrix
    //        rand_vals[i] = convert2double((generateCommonRandom() & MAXRAND));
            rand_vals[i] = convert2double(proxy->generateCommonRandom() & ORTHMASK);
        }

        // create an orthogonal matrix by using the generated random array in Eigen - the result is a Matrix
        Map<Matrix<double, Dynamic, Dynamic, RowMajor>> X(rand_vals, size, size);
        Matrix<double, Dynamic, Dynamic, RowMajor> XtX = X.transpose() * X;
        SelfAdjointEigenSolver<Matrix<double, Dynamic, Dynamic, RowMajor>> es(XtX);
        Matrix<double, Dynamic, Dynamic, RowMajor> S = es.operatorInverseSqrt();
        Matrix<double, Dynamic, Dynamic, RowMajor> orth = X * S;

        // convert a Matrix of double into two dimensional uint64_t array
        double *tmp = new double[size * size];
        Map<Matrix<double, Dynamic, Dynamic, RowMajor>>(tmp, size, size) = orth;
        uint64_t **M = new uint64_t *[size];
        double** debugging_M = new double*[size];
        for (int i = 0; i < size; i++) {
            M[i] = new uint64_t[size];
            debugging_M[i] = new double[size];
            for (int j = 0; j < size; j++) {
                M[i][j] = convert2uint64(tmp[i * size + j]);
                debugging_M[i][j] = tmp[i * size + j];
            }
        }
        return M;
    }
    return NULL;*/
}

void EIG(Party *proxy, uint32_t size, double epsilon = 0.01) {
    /*
     * Perform eigenvalue decomposition of a single Gram matrix
     */
/*    int p_role = proxy->getPRole();
    int socket_p1 = proxy->getSocketP1();
    int socket_p2 = proxy->getSocketP2();

    if (p_role == HELPER) {
        // receive the shares of the masked Gram matrix
        thread thr1 = thread(Receive, socket_p1, proxy->getBuffer1(), size * size * 8);
        thread thr2 = thread(Receive, socket_p2, proxy->getBuffer2(), size * size * 8);
        thr1.join();
        thr2.join();

        unsigned char *ptr = proxy->getBuffer1();
        unsigned char *ptr2 = proxy->getBuffer2();

        uint64_t **G1;
        uint64_t **G2;
        convert22DArray(&ptr, G1, size, size);
        convert22DArray(&ptr2, G2, size, size);

        // perform eigenvalue decomposition
        // 1. convert to double array
        double *masked_G = new double[size * size];
        for (uint32_t i = 0; i < size * size; i++) {
            masked_G[i] = convert2double(G1[i / size][i % size] + G2[i / size][i % size]);
        }

        // 2. convert to Matrix
        Map<Matrix<double, Dynamic, Dynamic, RowMajor>> masked_matrix_G(masked_G, size,
                                                                        size); // make sure that this conversion preserves the same orientation - row-major vs column-major!!!

        // 3. perform eigenvalue decomposition
        EigenSolver<Matrix<double, Dynamic, Dynamic, RowMajor>> eig_solver;
        eig_solver.compute(masked_matrix_G);
        ptr = proxy->getBuffer1();
        ptr2 = proxy->getBuffer2();

        // eigenvector - size-by-size
        if(DEBUG_FLAG >= 2)
            cout << "Step 3.1: eigenvectors" << endl;
        Matrix<double, Dynamic, Dynamic, RowMajor> eig_vectors = eig_solver.eigenvectors().real();
        double *matrix2double = new double[size * size];
        Map<Matrix<double, Dynamic, Dynamic, RowMajor>>(matrix2double, size,
                size) = eig_vectors; // convert from Matrix to double*
        uint64_t **M1 = new uint64_t *[size];
        uint64_t **M2 = new uint64_t *[size];
        uint64_t tmp_share;
        for (int i = 0; i < size; i++) {
            M1[i] = new uint64_t[size];
            M2[i] = new uint64_t[size];
            for (int j = 0; j < size; j++) {
                tmp_share = proxy->generateRandom();
                M1[i][j] = convert2uint64(matrix2double[i * size + j]) -
                           tmp_share; // make sure that this conversion preserves the same orientation - row-major vs column-major!!!
                M2[i][j] = tmp_share;
                addVal2CharArray(M1[i][j], &ptr);
                addVal2CharArray(M2[i][j], &ptr2);
            }
        }

        // eigenvalues - size-by-1
        double *double_eig_values = new double[size];
        Matrix<double, Dynamic, 1> vec_eig_values = eig_solver.eigenvalues().real();
        Map<Matrix<double, Dynamic, 1>>(double_eig_values, size,
                1) = vec_eig_values; // Since it is a vector, do we need to specify that it is row-major?

//        if(DEBUG_FLAG >= 3)
        print1DArray("Lambda + s", double_eig_values, size);

//        double* delta = new double[size]; // a vector of random values to mask eigenvalues
        uint64_t *masked_eig_values = new uint64_t[size];
        double alpha = MIN_DELTA + ((double) rand() / RAND_MAX) * (MAX_DELTA - MIN_DELTA); // a scalar random value to mask eigenvalues
        double* delta = new double[size];
        for (uint32_t i = 0; i < size; i++) {
            delta[i] = MIN_DELTA + ((double) rand() / RAND_MAX) * (MAX_DELTA - MIN_DELTA);
            masked_eig_values[i] = convert2uint64((double_eig_values[i] + epsilon) * delta[i] + alpha);
            addVal2CharArray(masked_eig_values[i], &ptr);
            addVal2CharArray(convert2uint64(delta[i]), &ptr2);
        }

        addVal2CharArray(convert2uint64(alpha), &ptr2); // add scalar alpha to buffer2

        // 4. send
        //  - the share of eigenvectors (size * size * 8 bits) and masked eigenvalues (size * 8 bits) to P1
        //  - the share of eigenvectors (size * size * 8 bits), the mask delta (size * 8 bits) and scalar alpha (8 bits) to P2
        thr1 = thread(Send, socket_p1, proxy->getBuffer1(), size * size * 8 + size * 8);
        thr2 = thread(Send, socket_p2, proxy->getBuffer2(), size * size * 8 + size * 8 + 8);
        thr1.join();
        thr2.join();

        for (uint32_t i = 0; i < size; i++) {
            delete[] M1[i];
            delete[] M2[i];
        }
        delete[] M1;
        delete[] M2;
        delete[] double_eig_values;
        delete[] masked_eig_values;
        delete[] matrix2double;
        delete[] masked_G;
        delete[] delta;
    }
//    else if(p_role == P1 || p_role == P2) {
//        proxy->SendBytes(RKN_EIG, size);
//    }*/
}

void EIG(Party* proxy, uint32_t n_gms, uint32_t size, double epsilon = 0.01) {
    /*
     * Perform eigenvalue decomposition of a single Gram matrix
     */
//    cout << "EIG:Size: " << size << endl;
//    cout << "EIG:n_gms: " << n_gms << endl;
    int p_role = proxy->getPRole();
    int socket_p1 = proxy->getSocketP1();
    int socket_p2 = proxy->getSocketP2();

    if (p_role == HELPER) {
        // receive the shares of the masked Gram matrix
        thread thr1 = thread(Receive, socket_p1, proxy->getBuffer1(), n_gms * size * size * 8);
        thread thr2 = thread(Receive, socket_p2, proxy->getBuffer2(), n_gms * size * size * 8);
        thr1.join();
        thr2.join();

        unsigned char *ptr = proxy->getBuffer1();
        unsigned char *ptr2 = proxy->getBuffer2();

        uint64_t ***G1;
        uint64_t ***G2;
        convert23DArray(&ptr, G1, n_gms, size, size);
        convert23DArray(&ptr2, G2, n_gms, size, size);

        // perform eigenvalue decomposition
        // 1. convert to double array
        double **masked_G = new double*[n_gms];
        for(uint32_t g = 0; g < n_gms; g++) {
            masked_G[g] = new double[size * size];
            for (uint32_t i = 0; i < size * size; i++) {
                masked_G[g][i] = convert2double(G1[g][i / size][i % size] + G2[g][i / size][i % size]);
            }
        }

        // initialize pointers to the beginning of the buffers to send eigenvalue- and eigenvector-related things
        ptr = proxy->getBuffer1();
        ptr2 = proxy->getBuffer2();

        // pre-allocation
        double *matrix2double = new double[size * size];
        uint64_t **M1 = new uint64_t *[size];
        uint64_t **M2 = new uint64_t *[size];
        for (int i = 0; i < size; i++) {
            M1[i] = new uint64_t[size];
            M2[i] = new uint64_t[size];
        }
        double *double_eig_values = new double[size];
        uint64_t *masked_eig_values = new uint64_t[size];
        double* delta = new double[size];
        uint64_t** dbg_masked_eig_values = new uint64_t*[n_gms];

        // for each gram matrix
        for(uint32_t g = 0; g < n_gms; g++) {
            // 2. convert to Matrix
            Map<Matrix<double, Dynamic, Dynamic, RowMajor>> masked_matrix_G(masked_G[g], size,
                                                                            size); // make sure that this conversion preserves the same orientation - row-major vs column-major!!!

            // 3. perform eigenvalue decomposition
            EigenSolver<Matrix<double, Dynamic, Dynamic, RowMajor>> eig_solver;
            eig_solver.compute(masked_matrix_G);

            // eigenvector - size-by-size
            Matrix<double, Dynamic, Dynamic, RowMajor> eig_vectors = eig_solver.eigenvectors().real();

            Map<Matrix<double, Dynamic, Dynamic, RowMajor>>(matrix2double, size,
                    size) = eig_vectors; // convert from Matrix to double*

            uint64_t tmp_share;
            for (int i = 0; i < size; i++) {
//                M1[i] = new uint64_t[size];
//                M2[i] = new uint64_t[size];
                for (int j = 0; j < size; j++) {
                    tmp_share = proxy->generateRandom();
                    M1[i][j] = convert2uint64(matrix2double[i * size + j]) -
                               tmp_share; // make sure that this conversion preserves the same orientation - row-major vs column-major!!!
                    M2[i][j] = tmp_share;
                    addVal2CharArray(M1[i][j], &ptr);
                    addVal2CharArray(M2[i][j], &ptr2);
                }
            }

            // eigenvalues - size-by-1
            Matrix<double, Dynamic, 1> vec_eig_values = eig_solver.eigenvalues().real();
            Map<Matrix<double, Dynamic, 1>>(double_eig_values, size,
                    1) = vec_eig_values; // Since it is a vector, do we need to specify that it is row-major?

//        double* delta = new double[size]; // a vector of random values to mask eigenvalues
            double alpha = MIN_DELTA + ((double) rand() / RAND_MAX) * (MAX_DELTA - MIN_DELTA); // a scalar random value to mask eigenvalues
            for (uint32_t i = 0; i < size; i++) {
                delta[i] = MIN_DELTA + ((double) rand() / RAND_MAX) * (MAX_DELTA - MIN_DELTA);
                masked_eig_values[i] = convert2uint64((double_eig_values[i] + epsilon) * delta[i] + alpha);
                addVal2CharArray(masked_eig_values[i], &ptr);
                addVal2CharArray(convert2uint64(delta[i]), &ptr2);
            }

            dbg_masked_eig_values[g] = masked_eig_values;

            addVal2CharArray(convert2uint64(alpha), &ptr2); // add scalar alpha to buffer2
        }

        // 4. send
        //  - the share of eigenvectors (size * size * 8 bits) and masked eigenvalues (size * 8 bits) to P1
        //  - the share of eigenvectors (size * size * 8 bits), the mask delta (size * 8 bits) and scalar alpha (8 bits) to P2
        thr1 = thread(Send, socket_p1, proxy->getBuffer1(), n_gms * (size * size * 8 + size * 8));
        thr2 = thread(Send, socket_p2, proxy->getBuffer2(), n_gms * (size * size * 8 + size * 8 + 8));
        thr1.join();
        thr2.join();

        for (uint32_t i = 0; i < size; i++) {
            delete[] M1[i];
            delete[] M2[i];
        }
        delete[] M1;
        delete[] M2;
        delete[] double_eig_values;
        delete[] masked_eig_values;
        delete[] matrix2double;
        delete[] masked_G;
        delete[] delta;
    }*/
}

uint64_t*** GM2KM(Party* proxy, uint64_t ***G, uint64_t alpha, uint32_t n_gms, uint32_t size) {
    /* This function computes the Gaussian kernel matrices based on the Gram matrices of the samples in bulk form.
     *
     * Input(s)
     * G: n_gms-many size-by-size Gram matrices
     * alpha: the similarity adjustment parameter of the Gaussian kernel function
     * n_gms: number of Gram matrices
     * size: number of samples in the Gram matrices
     *
     * Output(s)
     * Returns n_gms-many Gaussian kernel matrices of size size-by-size
     */
    int p_role = proxy->getPRole();
    if(p_role == P1 || p_role == P2) {
        uint32_t step_size = (size * (size + 1)) / 2;
        cout << "Step size in GM2KM_v2: " << step_size << endl;

        // computes the initial form of the kernel matrices
        uint64_t ***km = new uint64_t**[n_gms];
        // initialization
        for(uint32_t g = 0; g < n_gms; g++) {
            km[g] = new uint64_t*[size];
            for(uint32_t k = 0; k < size; k++) {
                km[g][k] = new uint64_t[size];
            }
        }

        // straighten the initial kernel matrix
        uint64_t *str_km = new uint64_t[n_gms * step_size];
        uint32_t ind = 0;
        for(uint32_t g = 0; g < n_gms; g++) {
            for(uint32_t i = 0; i < size; i++) {
                for(uint32_t j = i; j < size; j++) {
                    str_km[ind] = local_MUL(alpha, (G[g][i][j] - ((uint64_t) 1 << FRAC) * p_role));
                    ind++;
                }
            }
        }

        // compute the exponential of the values in the kernel matrix
        uint64_t *tmp_exp = EXP(proxy, str_km, n_gms * step_size);

        // for debugging -- to see if we compute the exponentials correctly
//        uint64_t*** atkafa = new uint64_t**[n_gms];
//        for(int i = 0; i < n_gms; i++) {
//            atkafa[i] = new uint64_t*[2];
//            atkafa[i][0] = new uint64_t[step_size];
//            atkafa[i][1] = new uint64_t[step_size];
//            for(int j = 0; j < step_size; j++) {
//                atkafa[i][0][j] = str_km[i * step_size + j];
//                atkafa[i][1][j] = tmp_exp[i * step_size + j];
//            }
//            print2DArrayRecAndConv("" + to_string(i), atkafa[i], 2, step_size, false);
//        }

        cout << "GM2KM: EXP computation is done!" << endl;

        // the first k layer is the first part of the resulting tmp_exp
        ind = 0;
        for(int i = 0; i < size; i++) {
            for(int j = i; j < size; j++) {
                km[0][i][j] = tmp_exp[ind];
                km[0][j][i] = km[0][i][j];
                ind++;
            }
        }

        // multiplication of the individually computed kernel matrices and the restoration of these multiplications
        uint64_t *tmp_mul = tmp_exp;
        for (int g = 1; g < n_gms; g++) {
            tmp_mul = MUL(proxy, tmp_mul, &tmp_exp[g * step_size], step_size);
            ind = 0;
            for(int i = 0; i < size; i++) {
                for(int j = i; j < size; j++) {
                    km[g][i][j] = tmp_mul[ind];
                    km[g][j][i] = km[g][i][j];
                    ind++;
                }
            }
        }
        return km;
    }
    else if( p_role == HELPER) {
        uint32_t step_size = (size * (size + 1)) / 2;

        EXP(proxy, 0, n_gms * step_size);
        for (int g = 1; g < n_gms; g++) {
            MUL(proxy, 0, 0, step_size);
        }
        return NULL;
    }
    return nullptr;
}

uint64_t** INVSQRT(Party* proxy, uint64_t **G, uint32_t size, double epsilon = 0.01) {
    /*
     * Computes the inverse square root of the given size-by-size gram matrix G by employing eigenvalue decomposition.
     * We based our solution ma
     * Reference: Zhou, Lifeng, and Chunguang Li. "Outsourcing eigen-decomposition and singular value decomposition of
     * large matrix to a public cloud." IEEE Access 4 (2016): 869-879.
     *
     * Input(s)
     * G: size-by-size gram matrix
     * size: the size of the gram matrix
     *
     * Output(s)
     *
     */
    int p_role = proxy->getPRole();
    if (p_role == P1 || p_role == P2) {
        // generate a scalar value whose max is MAXSCALAR
        uint64_t scalar_s = proxy->generateCommonRandom() & MAXSCALAR;
        printValue("s", convert2double(scalar_s));
        uint64_t scalar_a = proxy->generateCommonRandom() & MAXA; // multiplier
        double d_scalar_a = convert2double(scalar_a);
        printValue("a", d_scalar_a);

//        if (p_role == P1) {
//            // v1: compute A1 + sI
//            for (int i = 0; i < size; i++) {
//                G[i][i] += scalar_s;
//            }
        // v2: compute A1 + sI
        for( int j = 0; j < size; j++) {
            for (int i = 0; i < size; i++) {
                G[i][j] = local_MUL(G[i][j], scalar_a);
//                if(i == 0 && j == 0)
//                    printValueRecAndConv("G[0][0]", G[i][j]);
                if(i == j && p_role == P1) {
                    G[i][i] += scalar_s;
                }
            }
        }
//        }

//        print2DArrayRecAndConv("a * G + s * I", G, size, size);
        cout << "C1" << endl;
        // generate a random orthogonal matrix
        uint64_t **M = randomOrthogonalMatrix(proxy, size);
        cout << "C2" << endl;
        // compute M * (A + sI) * M^T
        uint64_t **masked_gram_matrix = local_MATMATMUL(M, G, size, size, size);
        cout << "C3" << endl;
        uint64_t **trM;
//        if (p_role == P1) {
        trM = new uint64_t *[size];
        for (int i = 0; i < size; i++) {
            trM[i] = new uint64_t[size];
            for (int j = 0; j < size; j++) {
                trM[i][j] = M[j][i];
            }
        }
        masked_gram_matrix = local_MATMATMUL(masked_gram_matrix, trM, size, size, size);
//        } else {
//            trM = M;
//            masked_gram_matrix = local_NF_MATMATMUL(masked_gram_matrix, trM, size, size, size);
//        }

//        print2DArrayRecAndConv("M (a * G + s * I) M.T", masked_gram_matrix, size, size);
        cout << "C4" << endl;
        // send the mask Gram matrix to Helper
//        proxy->SendBytes(RKN_EIG, size);
        EIG(proxy, size, epsilon);
        cout << "C4.5" << endl;
        unsigned char *ptr = proxy->getBuffer1();
        addArray2CharArray(masked_gram_matrix, &ptr, size, size);
        Send(proxy->getSocketHelper(), proxy->getBuffer1(), size * size * 8);

        /* receive the corresponding part of the resulting eigenvalue decomposition from Helper - note that these are
        specific for the computation of the inverse square root of the Gram matrix */
        // First size * size * 8 bits: masked eigenvalues (or unmasker in case P2)
        // Second size * size * 8 bits: masked eigenvectors
        // Additional 8 bits to get the scalar alpha in P2
        Receive(proxy->getSocketHelper(), proxy->getBuffer1(), size * size * 8 + size * 8 + (8 * p_role));
        if(DEBUG_FLAG >= 2)
            cout << "Received the components of eigenvalue decomposition" << endl;
        ptr = proxy->getBuffer1();
        uint64_t alpha;
        uint64_t **masked_eig_vecs;
        uint64_t *eig_vals; // = new uint64_t[size];
        cout << "C4.6" << endl;
        convert22DArray(&ptr, masked_eig_vecs, size, size); // the share of the masked eigenvectors
        cout << "C4.7" << endl;
        convert2Array(&ptr, eig_vals, size); // eigenvalue related things
        cout << "C4.8" << endl;

//        if(DEBUG_FLAG >= 3)
        print1DArray("Received eigenvalue related things", convert2double(eig_vals, size), size);

        if (p_role == P2) {
            alpha = convert2Long(&ptr); // scalar alpha
            if(DEBUG_FLAG >= 3)
                printValue("Alpha", convert2double(alpha));
        }

        // unmasking the eigenvectors
        uint64_t **eig_vecs = local_MATMATMUL(trM, masked_eig_vecs, size, size, size);

        if(DEBUG_FLAG >= 3)
            print2DArray("Reconstructed eigenvectors", convert2double(REC(proxy, eig_vecs, size, size), size, size), size, size);

        // unmasking the inverse square root of the eigenvalues
        uint64_t *unmasker = new uint64_t[size];;
        if (p_role == P1) {
            if(DEBUG_FLAG >= 2)
                cout << "Receiving unmasker from P2..." << endl;
            Receive(proxy->getSocketP2(), proxy->getBuffer1(), size * 8);
            if(DEBUG_FLAG >= 2)
                cout << "Done!" << endl;
            ptr = proxy->getBuffer1();
            convert2Array(&ptr, unmasker, size);
            for (uint32_t i = 0; i < size; i++) {
//                // v1
//                eig_vals[i] = eig_vals[i] - unmasker[i];
                // v2
                eig_vals[i] = convert2uint64(convert2double(eig_vals[i]) / d_scalar_a) - unmasker[i];
            }
        } else {
            // s * delta + alpha
            for (uint32_t i = 0; i < size; i++) {
//                // v1
//                unmasker[i] = local_NF_MUL(scalar_s, eig_vals[i]) + alpha;
                // v2
                unmasker[i] = convert2uint64(convert2double(local_MUL(scalar_s, eig_vals[i]) + alpha) / d_scalar_a);
            }

            if(DEBUG_FLAG >= 3)
                print1DArray("Unmasker", convert2double(unmasker, size), size);

            ptr = proxy->getBuffer1();
            addVal2CharArray(unmasker, &ptr, size);
            if(DEBUG_FLAG >= 2)
                cout << "Sending the unmasker to P1..." << endl;
            Send(proxy->getSocketP1(), proxy->getBuffer1(), size * 8);
            if(DEBUG_FLAG >= 2)
                cout << "Done!" << endl;
        }

        uint64_t *zero_vec = new uint64_t[size];
        for (uint32_t i = 0; i < size; i++) {
            zero_vec[i] = 0;
            if(p_role == P1) {
                eig_vals[i] = convert2uint64(1.0 / pow(convert2double(eig_vals[i]), 0.5));
            }
            else {
                eig_vals[i] = convert2uint64(pow(convert2double(eig_vals[i]), 0.5));
            }
        }

        if(DEBUG_FLAG >= 3)
            print1DArray("Revised eigenvalue related things", convert2double(eig_vals, size), size);

        if(DEBUG_FLAG >= 2)
            cout << "Computing the inverse square root of the eigenvalues..." << endl;
        uint64_t *sqrt_eig_vals;
        if (p_role == P1) {
            sqrt_eig_vals = MUL(proxy, eig_vals, zero_vec, size);
        } else {
            sqrt_eig_vals = MUL(proxy, zero_vec, eig_vals, size);
        }

        // if(DEBUG_FLAG >= 3)
        print1DArray("Inverse square root of the eigenvalues", convert2double(REC( proxy,sqrt_eig_vals, size), size), size);

        if(DEBUG_FLAG >= 2)
            cout << "The inverse square root of the eigenvalues are computed." << endl;

        // construct the inverse square root of the Gram matrix
        uint64_t **tr_eig_vecs = new uint64_t *[size];
        uint64_t **eig_vals_mat = new uint64_t *[size];
        for (uint32_t i = 0; i < size; i++) {
            tr_eig_vecs[i] = new uint64_t[size];
            eig_vals_mat[i] = new uint64_t[size];
            eig_vals_mat[i][i] = sqrt_eig_vals[i];
            for (uint32_t j = 0; j < size; j++) {
                tr_eig_vecs[i][j] = eig_vecs[j][i];

                // if the elements of the two dimensional dynamic array are zero after the initialization, no need for this part
                if (j != i) {
                    eig_vals_mat[i][j] = 0;
                }
            }
        }
        cout << "C5" << endl;
        uint64_t **invsqrt_G = MATMATMUL(proxy, MATMATMUL(proxy, eig_vecs, eig_vals_mat, size, size, size),
                                            tr_eig_vecs, size, size, size);
        if(DEBUG_FLAG >= 1)
            cout << "Returning from Party::INVSQRT...\n************************************************************" << endl;
        return invsqrt_G;
    }
    else if( p_role == HELPER) {
        EIG(proxy, size);
        MUL(proxy, 0, 0, size);
        MATMATMUL(proxy, 0, 0, size * size * size, 0, 0);
        MATMATMUL(proxy, 0, 0, size * size * size, 0, 0);
        return NULL;
    }
    return NULL;
}

uint64_t*** INVSQRT(Party* proxy, uint64_t ***G, uint32_t n_gms, uint32_t size, double epsilon = 0.01) {
    /*
     * Computes the inverse square root of the given size-by-size n_gms gram matrices in G by employing eigenvalue
     * decomposition. We based our solution the following:
     * Reference: Zhou, Lifeng, and Chunguang Li. "Outsourcing eigen-decomposition and singular value decomposition of
     * large matrix to a public cloud." IEEE Access 4 (2016): 869-879.
     *
     * Input(s)
     * G: n_gms-by-size-by-size gram matrix
     * n_gms: the number of gram matrices in G
     * size: the size of the gram matrix
     *
     * Output(s)
     *
     */
    int p_role = proxy->getPRole();
    if (p_role == P1 || p_role == P2) {
        // generate a scalar value whose max is MAXSCALAR
        uint64_t scalar_s = proxy->generateCommonRandom() & MAXSCALAR;
        double d_scalar_s = convert2double(scalar_s); // debugging purposes
//        cout << "INVSQRT: Scalar added to the gram matrices in double: " << d_scalar_s << endl; // debugging purposes
//        cout << "INVSQRT: Scalar added to the gram matrices in uint64_t: " << scalar_s << endl; // debugging purposes

        if (p_role == P1) {
            // compute A1 + sI
            for(int g = 0; g < n_gms; g++) {
                for (int i = 0; i < size; i++) {
                    G[g][i][i] += scalar_s;
                }
            }
        }

        // generate random orthogonal matrices
        uint64_t ***M = new uint64_t**[n_gms];
        for(int i = 0; i < n_gms; i++) {
            M[i] = randomOrthogonalMatrix(proxy, size);
            double** tmp_M = convert2double(M[i], size, size);
            double** tr_tmp_M = new double*[size];
            for(int j = 0; j < size; j++) {
                tr_tmp_M[j] = new double[size];
                for(int k = 0; k < size; k++) {
                    tr_tmp_M[j][k] = tmp_M[k][j];
                }
            }
        }

        // compute M * (A + sI) * M^T
//        uint64_t ***masked_gram_matrix = MNF_MATMATMUL(M, G, n_gms, size, size, size);
        uint64_t ***masked_gram_matrix = local_MATMATMUL(M, G, n_gms, size, size, size);

        uint64_t ***trM;
//        if (p_role == P1) {
        trM = new uint64_t**[n_gms];
        for(int g = 0; g < n_gms; g++) {
            trM[g] = new uint64_t *[size];
            for (int i = 0; i < size; i++) {
                trM[g][i] = new uint64_t[size];
                for (int j = 0; j < size; j++) {
                    trM[g][i][j] = M[g][j][i];
                }
            }
        }
//            masked_gram_matrix = MNF_MATMATMUL(masked_gram_matrix, trM, n_gms, size, size, size);
        masked_gram_matrix = local_MATMATMUL(masked_gram_matrix, trM, n_gms, size, size, size);
//        } else {
//            trM = M;
////            masked_gram_matrix = MNF_MATMATMUL(masked_gram_matrix, trM, n_gms, size, size, size);
//            masked_gram_matrix = local_MNF_MATMATMUL(masked_gram_matrix, trM, n_gms, size, size, size);
//        }

        // debugging purposes
//        double*** rec_kmer_kms = new double**[n_gms];
//        double*** rec_M = new double**[n_gms];
//        double*** rec_tr_M = new double**[n_gms];
//        for(int g = 0; g < n_gms; g++) {
//            rec_kmer_kms[g] = convert2double(REC(proxy, G[g], size, size), size, size);
//            rec_M[g] = convert2double(M[g], size, size);
//            rec_tr_M[g] = convert2double(trM[g], size, size);
//        }
//        double*** MGMT = double_MATMATMUL(double_MATMATMUL(rec_M, rec_kmer_kms, n_gms, size, size, size), rec_tr_M,
//                                          n_gms, size, size, size);
//
//        double** eig_MGMT = new double*[n_gms];
//        for(int i = 0; i < n_gms; i++) {
//            EigenSolver<Matrix<double, Dynamic, Dynamic, RowMajor>> AT_ges;
//            Map<Matrix<double, Dynamic, Dynamic, RowMajor>> AT_matrix_G(straighten2DArray(MGMT[i], size, size), size, size);
//            AT_ges.compute(AT_matrix_G);
//            Matrix<double, Dynamic, 1> AT_eig_vals = AT_ges.eigenvalues().real();
//            eig_MGMT[i] = new double[size];
//            Map<Matrix<double, Dynamic, 1>>(eig_MGMT[i], size) = AT_eig_vals;
//            for(int j = 0; j < size; j++) {
//                eig_MGMT[i][j] -= d_scalar_s;
//            }
//            bubbleSort(eig_MGMT[i], size);
//        }
////        print2DArray("Eigenvalues of MGMT", eig_MGMT, n_gms, size, false);
//
//        double*** rec_masked_gram_matrix = new double**[n_gms];
//        for(int i = 0; i < n_gms; i++) {
//            rec_masked_gram_matrix[i] = convert2double(REC(proxy, masked_gram_matrix[i], size, size), size, size);
//        }
//
//        double*** diff_MGMT = new double**[n_gms];
//        double* total_diffs = new double[n_gms];
//
//        for(int i = 0; i < n_gms; i++) {
////            saveAsCsv("output/MGMT_" + to_string(i) + ".csv", MGMT[i], size, size);
////            saveAsCsv("output/rec_masked_gram_matrix_" + to_string(i) + ".csv", rec_masked_gram_matrix[i], size, size);
//            double tmp = 0;
//            int cntr_invalid = 0;
//            int cntr_negative = 0;
//            diff_MGMT[i] = new double*[size];
//            for(int j = 0; j < size; j++) {
//                diff_MGMT[i][j] = new double[size];
//                for(int k = 0; k < size; k++) {
//                    diff_MGMT[i][j][k] = MGMT[i][j][k] - rec_masked_gram_matrix[i][j][k];
//                    tmp += diff_MGMT[i][j][k];
//                    if(abs(MGMT[i][j][k] - rec_masked_gram_matrix[i][j][k]) > 1) {
//                        cout << "MGMT[" << i << "][" << j << "][" << k << "]: " << MGMT[i][j][k] << "\tComputed: " << rec_masked_gram_matrix[i][j][k] << endl;
//                        cntr_invalid++;
//                    }
//                }
//            }
//            total_diffs[i] = tmp;
//            cout << "In " << i << "th kernel matrix, there are " << cntr_invalid << " invalid entries." << endl;
//        }
//        print1DArray("Differences between GT MGMT and computed one", total_diffs, n_gms);
        // end of debugging

        // send the mask Gram matrix to Helper
        EIG(proxy, n_gms, size, epsilon);
        unsigned char *ptr = proxy->getBuffer1();
        for(uint32_t i = 0; i < n_gms; i++) {
            addArray2CharArray(masked_gram_matrix[i], &ptr, size, size);
        }
        Send(proxy->getSocketHelper(), proxy->getBuffer1(), n_gms * size * size * 8);

        /* receive the corresponding part of the resulting eigenvalue decomposition from Helper - note that these are
        specific for the computation of the inverse square root of the Gram matrix */
        // We have n_gms chunks. For each chunk of (size * size * 8 + size * 8 + p_role * 8) bits block:
        // First size * size * 8 bits: masked eigenvalues (or unmasker in case P2)
        // Second size * size * 8 bits: masked eigenvectors
        // Additional 8 bits to get the scalar alpha in P2
        Receive(proxy->getSocketHelper(), proxy->getBuffer1(), n_gms * (size * size * 8 + size * 8 + (8 * p_role)));
        ptr = proxy->getBuffer1();
        uint64_t *alpha = new uint64_t[n_gms];
        uint64_t ***masked_eig_vecs = new uint64_t**[n_gms];
        uint64_t **eig_vals = new uint64_t*[n_gms];

        for(uint64_t g = 0; g < n_gms; g++) {
            convert22DArray(&ptr, masked_eig_vecs[g], size, size); // the share of the masked eigenvectors
            convert2Array(&ptr, eig_vals[g], size); // eigenvalue related things

            if (p_role == P2) {
                alpha[g] = convert2Long(&ptr); // scalar alpha
            }
        }

        // unmasking the eigenvectors
//        uint64_t ***eig_vecs = MNF_MATMATMUL(trM, masked_eig_vecs, n_gms, size, size, size);
        uint64_t ***eig_vecs = local_MATMATMUL(trM, masked_eig_vecs, n_gms, size, size, size);

        // unmasking the inverse square root of the eigenvalues
        uint64_t **unmasker;
        if (p_role == P1) {
            Receive(proxy->getSocketP2(), proxy->getBuffer1(), n_gms * size * 8);
            ptr = proxy->getBuffer1();
            convert22DArray(&ptr, unmasker, n_gms, size);
            for(uint32_t g = 0; g < n_gms; g++) {
                for (uint32_t i = 0; i < size; i++) {
                    eig_vals[g][i] = eig_vals[g][i] - unmasker[g][i];
                }
            }
        } else {
            ptr = proxy->getBuffer1();
            unmasker = new uint64_t*[n_gms];
            for(uint32_t g = 0; g < n_gms; g++) {
                // s * delta + alpha
                unmasker[g] = new uint64_t[size];
                for (uint32_t i = 0; i < size; i++) {
                    unmasker[g][i] = local_MUL(scalar_s, eig_vals[g][i]) + alpha[g];
                }
                addVal2CharArray(unmasker[g], &ptr, size);
            }
            Send(proxy->getSocketP1(), proxy->getBuffer1(), n_gms * size * 8);
        }

//        uint64_t **zero_vec = new uint64_t*[n_gms];
        uint64_t *str_zero_vec = new uint64_t[n_gms * size];
        uint64_t *str_untouched_eigvals = new uint64_t[n_gms * size]; // for debugging purposes
        uint64_t *str_processed_eigvals = new uint64_t[n_gms * size];
        for(uint32_t g = 0; g < n_gms; g++) {
//            zero_vec[g] = new uint64_t[size];
            for (uint32_t i = 0; i < size; i++) {
//                zero_vec[g][i] = 0;
                str_zero_vec[g * size + i] = 0;
                if(p_role == P1) {
                    str_untouched_eigvals[g * size + i] = eig_vals[g][i];
                    str_processed_eigvals[g * size + i] = convert2uint64(1.0 / pow(convert2double(eig_vals[g][i]), 0.5));
                }
                else {
                    str_untouched_eigvals[g * size + i] = convert2uint64(1.0 / convert2double(eig_vals[g][i]));
                    str_processed_eigvals[g * size + i] = convert2uint64(pow(convert2double(eig_vals[g][i]), 0.5));
                }
//                str_processed_eigvals[g * size + i] = eig_vals[g][i];
            }
        }

//        cout << "C3: MUL: Size: " << (n_gms * size) << endl;
        uint64_t *tmp_res;
        if (p_role == P1) {
//            for(uint32_t g = 0; g < n_gms; g++) {
//                for(uint32_t i = 0; i < size; i++) {
//                    tmp_all_eig_vals[g * size + i] = eig_vals[g][i];
//                }
//            }
            tmp_res = MUL(proxy, str_processed_eigvals, str_zero_vec, n_gms * size);
        } else {
//            for(uint32_t g = 0; g < n_gms; g++) {
//                for(uint32_t i = 0; i < size; i++) {
//                    tmp_all_eig_vals[g * size + i] = eig_vals[g][i];
//                }
//            }
            tmp_res = MUL(proxy, str_zero_vec, str_processed_eigvals, n_gms * size);
        }

        uint64_t **sqrt_eig_vals = new uint64_t*[n_gms];
        for(uint32_t g = 0; g < n_gms; g++) {
            sqrt_eig_vals[g] = new uint64_t[size];
            for(uint32_t i = 0; i < size; i++) {
                sqrt_eig_vals[g][i] = tmp_res[g * size + i];
            }
        }

//        print2DArray("Square root of eigenvalues", convert2double(REC(proxy, sqrt_eig_vals, n_gms, size), n_gms, size), n_gms, size);

        // construct the inverse square root of the Gram matrix
        uint64_t ***tr_eig_vecs = new uint64_t**[n_gms];
        uint64_t ***eig_vals_mat = new uint64_t**[n_gms];
        for(uint32_t g = 0; g < n_gms; g++) {
            tr_eig_vecs[g] = new uint64_t*[size];
            eig_vals_mat[g] = new uint64_t*[size];
            for (uint32_t i = 0; i < size; i++) {
                tr_eig_vecs[g][i] = new uint64_t[size];
                eig_vals_mat[g][i] = new uint64_t[size];
                eig_vals_mat[g][i][i] = sqrt_eig_vals[g][i];
                for (uint32_t j = 0; j < size; j++) {
                    tr_eig_vecs[g][i][j] = eig_vecs[g][j][i];

                    // if the elements of the two dimensional dynamic array are zero after the initialization, no need for this part
                    if (j != i) {
                        eig_vals_mat[g][i][j] = 0;
                    }
                }
            }
        }

        uint64_t*** tmp_G = MATMATMUL(proxy, eig_vecs, eig_vals_mat, n_gms, size, size, size);
        uint64_t ***invsqrt_G = MATMATMUL(proxy, tmp_G,tr_eig_vecs, n_gms, size, size, size);
//        uint64_t ***invsqrt_G = local_MNF_MATMATMUL(local_MNF_MATMATMUL(eig_vecs, eig_vals_mat, n_gms, size, size, size),
//                                            tr_eig_vecs, n_gms, size, size, size);

//        cout << "C7" << endl;
        for(uint32_t g = 0; g < n_gms; g++) {
            delete [] sqrt_eig_vals[g];
            for(uint32_t i = 0; i < size; i++) {
                delete [] tr_eig_vecs[g][i];
                delete [] eig_vals_mat[g][i];
            }
            delete [] tr_eig_vecs[g];
            delete [] eig_vals_mat[g];
        }
        delete [] sqrt_eig_vals;
        delete [] tr_eig_vecs;
        delete [] eig_vals_mat;
        delete [] str_processed_eigvals;
        delete [] str_zero_vec;
        delete [] tmp_res;

//        cout << "Returning from INVSQRT...\n************************************************************" << endl;
        return invsqrt_G;
    }
    else if( p_role == HELPER) {
        EIG(proxy, n_gms, size);
        MUL(proxy, 0, 0, n_gms * size);
        MATMATMUL(proxy, 0, 0, 0, n_gms * size * size * size, 0, 0);
        MATMATMUL(proxy, 0, 0, 0, n_gms * size * size * size, 0, 0);
        return NULL;
    }
    return NULL;
}

uint64_t* RKN_ITERATION(Party* proxy, uint64_t* x, uint64_t* z, uint64_t* ct1, uint32_t n_dim, uint32_t n_anc, uint32_t k_mer, double lambda, double alpha) {
    /* This function performs a single time point iteration of ppRKN inference.
     *
     * Input(s)
     * proxy: Party instance
     * x: data vector of size n_dim at time t
     * z: anchor points of size (k_mer * n_anc * n_dim)
     * ct1: the mapping from time t-1 of size ((k_mer + 1) * n_anc) -- the first n_anc elements are 1
     * n_dim: the number of elements to represent characters in the sequence
     * n_anc: the number of anchor points
     * k_mer: the k value of k-mers
     * lambda: the downscaling factor to decrease the effect of the previous time point
     *
     * Return(s)
     * ckt: the output of the mapping of the sequence after time point t
     */

    int p_role = proxy->getPRole();

    if( p_role == P1 || p_role == P2) {
        uint32_t size = k_mer * n_anc * n_dim;
        uint32_t size2 = k_mer * n_anc;
        uint64_t* rep_x = new uint64_t[size];

        for(int i = 0; i < k_mer; i++) {
            for(int j = 0; j < n_anc; j++) {
                for(int k = 0; k < n_dim; k++) {
                    rep_x[(i * n_anc * n_dim) + (j * n_dim) + k] = x[k];
                }
            }
        }
    //    print1DArray("rep_x in pprkn_iteration", Mconvert2double(MReconstruct(rep_x, size), size), size);

        // computation of b_{k}[t]
        uint64_t* dp = DP(proxy, rep_x, z, size, n_dim);
        uint64_t tmp_alpha = convert2uint64(alpha);
        uint64_t tmp_minus_one = convert2uint64(-1);
        for(uint32_t i = 0; i < size2; i++) {
    //        uint64_t tmp = dp[i] + p_role * tmp_minus_one;
            uint64_t tmp = dp[i] - (((uint64_t) 1 << FRAC) * p_role);
    //        printValue("tmp " + to_string(i), convert2double(Reconstruct(tmp)));
            dp[i] = local_MUL(tmp_alpha, tmp);
    //        printValue("dp[" + to_string(i) + "]", convert2double(Reconstruct(dp[i])));
        }

    //    print1DArray("dp in pprkn_iteration", Mconvert2double(MReconstruct(dp, size2), size2), size2);
        uint64_t* res_b = EXP(proxy, dp, size2);
//        print1DArray("res_b in pprkn_iteration", convert2double(REC(proxy, res_b, size2), size2), size2);

        // computation of c_{k-1}[t-1] * b_{k}[t]
        uint64_t* skt = MUL(proxy, res_b, ct1, size2);
//        print1DArray("skt in pprkn_iteration", convert2double(REC(proxy, skt, size2), size2), size2);

        // computation of c_{k-1}[t-1] * b_{k}[t]
    //    cout << "********* lambda: " << convert2uint64(lambda) << endl;
        uint64_t* ckt = new uint64_t[size2];
        for(uint32_t i = 0; i < size2; i++) {
    //        ckt[i] = skt[i] + convert2uint64(lambda * convert2double(ct1[i + n_anc]));
            ckt[i] = local_MUL(convert2uint64(lambda), skt[i]) + local_MUL(convert2uint64(1 - lambda), ct1[i + n_anc]); // this provides more precision
    //        ckt[i] = skt[i] + local_NF_MUL(convert2uint64(lambda), ct1[i + n_anc]); // this provides more precision
        }

        delete[] dp;
        delete[] res_b;
        delete[] skt;
        delete[] rep_x;

        return ckt;
    }
    else if( p_role == HELPER) {
        // n_dim contains "size"
        // n_anc contains "size2"
        uint32_t size = n_dim;
        uint32_t size2 = n_anc;

        DP(proxy, 0, 0, size, 0);
        EXP(proxy, 0, size2);
        MUL(proxy, 0, 0, size2);

        return NULL;
    }
    return NULL;
}


#endif //CECILIA_RKN_H
