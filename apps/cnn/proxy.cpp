#include <cstdlib>
#include <iostream>
#include <deque>
#include <chrono>
#include <iomanip>
#include <assert.h>
#include "../../core/cnn.h"
using namespace std;


void print1DArray(string const &str1, double* x, uint32_t size) {
    cout << "======================= " << str1 << " =======================" << endl;
    for(uint32_t i = 0; i < size; i++) {
        cout << x[i] << "\t";
    }
    cout << endl;
    cout << "==============================================================" << endl;
}

/**
 * Prints the given matrix in a way to show the values in a column and row based format while alsow different windows of
 * the matrix are seperated from each other visually.
 * @param str1 - String to be printed out before the visualization of the matrix.
 * @param matrix - the matrix to be shown represented by a vector where each row of the matrix contains m_col / w_col values per window.
 * @param m_row - number of rows in matrix
 * @param m_col - number of columns in matrix
 * @param w_row - number of rows per window
 * @param w_col - number of columns per window
 */
void print1DMatrixByWindows(string const &str1, uint64_t* matrix, uint32_t m_row, uint32_t m_col, uint32_t w_row, uint32_t w_col) {
    uint32_t rowNumberOfWin = floor(m_col/w_col);
    uint64_t mSize = m_col*m_row;
    double *mConverted = convert2double(matrix, mSize);
    cout << "======================= " << str1 << " =======================" << endl << endl;
    for(uint32_t i = 0; i < m_row; i++) {
        //delimiter between windows in horizontal direction
        if(i % w_row == 0){
            for(uint8_t d = 0; d < m_col; d++){
                cout << " _ _ _ _ _ _ _ _";
            }
            cout << endl;
        }
        for(uint32_t j = 0; j < m_col; j++){
            //delimiter between windows in vertical direction
            if(j % w_col == 0){
                cout << "|\t";
            }
            cout << mConverted[i*m_col + j] << " \t " ;
        }
        cout << "|" << endl;
    }
    for(uint8_t d = 0; d < m_col; d++){
        cout << " _ _ _ _ _ _ _ _";
    }
    cout << endl << "==============================================================" << endl;
}

int main(int argc, char* argv[]) {
    uint8_t role = atoi(argv[1]);
    uint16_t cport = atoi(argv[2]);
    string caddress(argv[3]);
    uint16_t hport = atoi(argv[4]);
    string haddress(argv[5]);

    uint16_t window_size = atoi(argv[6]);
    string matrix_string = argv[7];

    uint32_t matrix_col_size = atoi(argv[8]);
    uint32_t matrix_row_size = atoi(argv[9]);

    //ensure ports are not 0 for helper or client
    cout << "Setting ports for helper/client...";
    if (cport != 0) {
        assert(cport < 1 << (sizeof(uint16_t) * 8));
    }
    if (hport != 0) {
        assert(hport < 1 << (sizeof(uint16_t) * 8));
    }

    Party *proxy;
    cout << "Creating Proxy...";
    if (role==0)
        proxy = new Party(P1,hport, haddress, cport, caddress);
    else
        proxy = new Party(P2,hport, haddress, cport, caddress);
    //cout << "Secret matrix (" << matrix_row_size << " x " << matrix_col_size << ")" << endl;

    //init
    uint16_t matrix_size = matrix_row_size * matrix_col_size;
    double *doubleMatrix = new double [matrix_size];

    //fill matrix (as vector) with params from input
    stringstream mss(matrix_string);
    uint32_t i = 0;
    while(mss.good() && i<matrix_size){
        string elements;
        getline(mss, elements, ',');
        double element = std::atof(elements.c_str());
        doubleMatrix[i] = element;
        //cout << element << " ";
        i++;
    }
    //convert matrix elements
    uint64_t *matrix = convert2uint64(doubleMatrix, matrix_size, 15);
    uint64_t *resortedMatrix = new uint64_t [matrix_size];
    //                                                     here same param twice but could be different sizes
    resortedMatrix = resortMatrix(matrix, matrix_col_size, matrix_row_size, window_size, window_size, resortedMatrix);

    //helper variable to create shares of v1 and v2
    uint64_t *mTmp = new uint64_t [matrix_size]; //pick random value in ring --> shareOfMatrix is calculated from that
    uint64_t *shareOfMatrix;

    //fill helper variables with random values
    for(uint32_t i = 0; i < matrix_size; i++){
        mTmp[i] = proxy->generateCommonRandom();
    }
    cout << "Start creating share for proxy " << to_string(role) << endl;
    //create shares
    if(role == P1) {
        for(uint32_t i = 0; i < matrix_size; i++){
            shareOfMatrix[i] = resortedMatrix[i] - mTmp[i];
        }
    }
    else {
        shareOfMatrix = mTmp;
    }

    cout << "___MMAXPOOL___: " << endl;
    cout << "Ground Truth: " << endl;
    print1DMatrixByWindows("matrix", matrix, matrix_row_size, matrix_col_size, window_size, window_size);
    delete []matrix;
    print1DMatrixByWindows("resorted matrix", resortedMatrix, matrix_row_size, matrix_col_size, window_size, window_size);

    uint64_t *maxElements = nullptr;
    maxElements = MAX(proxy, shareOfMatrix, matrix_size, window_size, 40000000);

    uint16_t window_length = window_size*window_size;
    uint16_t numberOfWins = static_cast<uint16_t >(floor(matrix_size / window_length));
    print1DArray("Final result", convert2double(REC(proxy, maxElements, numberOfWins), numberOfWins), numberOfWins);
    //proxy->print1DArray("Final result", proxy->Mconvert2double(proxy->REC(maxElements, numberOfWins, N), numberOfWins), numberOfWins);


    // CNN inference pipeline, call functions sequantially for inference (Matrix MUL, MAXPOOL, RELU, vs...)


    proxy->SendBytes(CORE_END);
    proxy->PrintBytes();
    return 0;
}

