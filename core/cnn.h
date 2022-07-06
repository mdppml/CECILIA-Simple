//
// Created by Mete Akgun on 17.01.22.
//

#ifndef PPAUC_CNN_H
#define PPAUC_CNN_H

uint64_t * FLT(uint64_t*** images, uint32_t i_dim, uint32_t i_number);
uint64_t ***INC(uint64_t ***input, uint32_t channel, uint32_t conv_height, uint32_t conv_width, uint32_t kernel_dim,
                uint32_t stride);
uint64_t ** transpose(uint64_t** matrix, uint32_t rows, uint32_t cols);


#include "core.h"
#include "../utils/flib.h"
#include "bitset"


void print1DMatrixByWindows(string const &str1, double *matrix, uint32_t m_row, uint32_t m_col, uint32_t w_row,
                             uint32_t w_col) {
    cout << "======================= " << str1 << " =======================" << endl << endl;
    for(uint32_t i = 0; i < m_row; i++) {
        //delimiter between windows in horizontal direction
        if(i % w_row == 0){
            for(uint32_t d = 0; d < m_col; d++){
                cout << " _ _ _ _ _ _ _ _";
            }
            cout << endl;
        }
        for(uint32_t j = 0; j < m_col; j++){
            //delimiter between windows in vertical direction
            if(j % w_col == 0){
                cout << "|\t";
            }
            cout << matrix[i*m_col + j] << " \t " ;
        }
        cout << "|" << endl;
    }
    for(uint32_t d = 0; d < m_col; d++){
        cout << " _ _ _ _ _ _ _ _";
    }
    cout << endl << "==============================================================" << endl;
}


/**
 * Vectorized subtraction: subtract elementwise b from a (a-b)
 * @param a from values of this vector the values of b will be subtracted
 * @param b the vector whose values are going to be subtracted from a
 * @param length of a and b (must have same length)
 * @return vector d where each element is the result of the according elements a-b.
 */
uint64_t* SUB(const uint64_t *a, const uint64_t *b, uint32_t length){
    uint64_t *subtractedValues = new uint64_t [length];
    for(uint32_t i = 0; i<length; i++){
        subtractedValues[i] = *(a+i) - *(b+i);
        //cout << i << ": (a-b = c) " << to_string(a[i]) << " - " << to_string(b[i]) << " = " << to_string(subtractedValues[i]) << endl;
    }
    return subtractedValues;
}


/**
 * Resort (RST) the given matrix so that elements of one window are found as a sequence of w_rows*w_cols.
 * Elements of a window can be found in matrix starting from index
 * i up to i + w_cols and
 * i + m_cols * win_row up to i + m_cols * win_row for each row of window  win_row.
 * compute window_size vectors:
 * matrix vector M:  a, b, c, d,
 *                   e, f, g, h,
 *                   i, j, k, l,
 *                   m, n, o, p              shall become
 * 4 window_size vectors  a, b, e, f,
 *                        c, d, g, h,
 *                        i, j, m, n,
 *                        k, l, o, p
 * Note: this method can handle and resort only values from windows which are subsequent. Assuming the true matrix M but
 *       the method receives only
 *       vector m:  a, b,
 *                  e, f
 *       This will not correctly be resorted as matrix contains windows c and d between b and e which are
 *       unknown to this method.
 * @param matrix the matrix to be resorted
 * @param m_cols number of columns in matrix
 * @param m_rows number of rows in matrix
 * @param w_cols number of columns in window
 * @param w_rows number of rows in window
 * @param resortedMatrix pointer to the resulting resorted Matrix, which will have length of m_cols * m_rows
 *
 * CAUTION: can not process matrices greater than a size of xx by now.  TODO find out the max size
 */
void RST(const uint64_t* matrix, uint32_t m_cols, uint32_t m_rows, uint32_t w_cols, uint32_t w_rows, uint64_t* resortedMatrix){
    uint32_t winSize = w_cols * w_rows;
    uint32_t numberOfWins = (m_cols * m_rows) / winSize;
    uint32_t winsPerRow = m_cols / w_cols;
    uint32_t w_count = 0;

    while(w_count < numberOfWins){
        uint32_t windowStart = static_cast<uint32_t>(w_rows * m_cols * floor(w_count / winsPerRow) + (w_cols * (w_count & (winsPerRow-1))));
        for(uint32_t i = 0; i<w_rows; i++){
            uint32_t m_index = windowStart + i*m_cols;
            for(uint32_t  j = 0; j<w_cols; j++){
                uint32_t w_index = w_count * winSize + (i*w_cols + j);
                resortedMatrix[w_index] = matrix[m_index + j];
            }
        }
        w_count++;
    }
}

/**
 * Selects the maximum element from the given matrix in secret shared form.
 * @param mShare - secret share of the matrix from which the maximal element shall be computed.
 * @param matrix_size - size of mShare.
 * @return The maximum element which was found in mShare.
 */
uint64_t MAX(Party* proxy, uint64_t *mShare, uint32_t matrix_size){
    /** MAIN IDEA:
     * Compare values by splitting the matrix in two halves and
     * comparing each value to its counterpart at the same position in the other half.
     * If size of the given matrix is odd, there will be a residue, which is stored in residue.
     */
    uint32_t cmpVectorSize = matrix_size; //size of resulting vector after cmp, MUX and its divided by 2 is size of each half.
    bool isResidueStored = false;
    if ( proxy->getPRole() == P1 ||  proxy->getPRole() == P2) {
        uint64_t *maxElements = mShare;
        uint32_t maxHalfSizes = static_cast<uint32_t>(ceil(matrix_size / 2)); // ceil because residue might be added.
        uint64_t firstHalf [maxHalfSizes];
        uint64_t secondHalf [maxHalfSizes];
        uint64_t residue; //there is at most one residual element.

        while (cmpVectorSize > 1) {
            uint32_t halfSize = static_cast<uint32_t>(floor(cmpVectorSize / 2));

            for (uint32_t i = 0; i < halfSize; i++){
                firstHalf[i] = maxElements[i];
                secondHalf[i] = maxElements[i + halfSize];
            }
            if (cmpVectorSize & 0x1) {                          //there is a residue remaining
                if (isResidueStored) {                          //second residue found --> add stored and current residue to one half each.
                    halfSize++;                                 //each half of window increases size by 1 because of residues
                    *(firstHalf + halfSize) = residue;
                    *(secondHalf + halfSize) = *(firstHalf + cmpVectorSize - 1); //last element of cmpVector
                } else {                                        //no residue stored up to now:
                    isResidueStored = true;
                    residue = *(firstHalf + cmpVectorSize - 1); // store last element in residue
                }
            }
            if (halfSize > 0) {                                 // maximums are not yet found
                //compare: a-b =c and then MSB(c) =d
                uint64_t *c = SUB(firstHalf, secondHalf, halfSize);
                uint64_t *d = MSB(proxy, c, halfSize);

                //MUX:
                maxElements = MUX(proxy, firstHalf, secondHalf, d, halfSize);
            }
            //prepare next round:
            cmpVectorSize = halfSize;
        }
        uint64_t max = maxElements[0];                          // should only contain one element at the end.
        delete [] maxElements;
        return max;
    }
    else if ( proxy->getPRole() == HELPER) {
        /**Compares values in a given window by splitting the window in two halves and comparing each value to its counterpart at the same position in the other half.
        If size of the given windowVector is odd, there will be a residue, which is stored in residue. */
        while (cmpVectorSize > 1) {
            uint32_t halfSize = static_cast<uint32_t>(floor(cmpVectorSize / 2));
            if (cmpVectorSize % 2 == 1) {                   //there is an residue remaining
                if (isResidueStored) {                            //second residue found --> add stored and current residue each to one half.
                    halfSize++;                                //each half of window increases size by 1 because of residues
                } else {                                       //no residue stored up to now:
                    isResidueStored = true;
                }
            }
            //if cmpVectorSize is odd, store the last element as residue
            if (halfSize > 0) {          // maximums are not yet found
                //compare: a-b =c and then MSB(c) =d
                MSB(proxy, nullptr, halfSize);

                //MUX:
                MUX(proxy, nullptr, nullptr, nullptr, halfSize);
            }
            //prepare next round:
            cmpVectorSize = halfSize;
        }
        return 0;
    }
    return -1;
}


/**
 * Maxpooling of a matrix of size matrix_size and non overlapping windows of size window_size. Address of the first
 * address in the matrix is mShare. as the matrix is processed as secret shares.
 * @param mShare - share of the matrix where values of a window are to be expected subsequently and then values of the
 * next value will follow (so not rows are concatenated containing values of all windows over all columns but windows
 * are concatenated).
 * @param m_rows - number of rows in the matrix mShare.
 * @param m_cols - number of columns in the matrix mShare.
 * @param window_size - size of the window in both dimensions (symmetric windows possible only)
 * @return the maximum element per window, therefore floor(matrix_size / (window_size*window_size)) elements in a vector
 */
uint64_t* MAX(Party* proxy, uint64_t *mShare, uint32_t m_rows, uint32_t m_cols, uint32_t window_size){
    uint32_t matrix_size = m_rows * m_cols;

    //compare within one window for each window:
    uint32_t window_length = window_size * window_size;
    uint32_t cmpWindowVectorSize = window_length; //size of resulting vector after cmp, MUX and its divided by 2 is size of each half.
    bool isResidueStored = false;
    bool isResidueInBuffer = false;

    uint32_t numberOfWins;
    if ( proxy->getPRole() == P1 ||  proxy->getPRole() == P2) {
        uint64_t *resorted = new uint64_t [matrix_size];                // RESORT matrix to have all values of a window subsequently
        // only rows calculated because always complete rows of windows are processed
        RST(mShare, m_cols, m_rows, window_size, window_size, resorted);
        numberOfWins = matrix_size / window_length;
        uint64_t *maxElements = new uint64_t [numberOfWins];
        uint32_t comparisons;

        while (cmpWindowVectorSize > 1) {
            /**Compares values in a given window by splitting the window in two halves and comparing each value to
             * its counterpart at the same position in the other half. If size of the given windowVector is odd,
             * there will be a residue, which is stored in residue.
             * All three vectors are initialized within the loop to free space as less and less space is required.*/
            uint32_t halfSize = cmpWindowVectorSize / 2; // elements per win in halves.

            comparisons = halfSize * numberOfWins;
            uint64_t firstHalf [comparisons];
            uint64_t secondHalf [comparisons];
            uint64_t residue [numberOfWins]; //if there are residues, for each window there is one element --> size is number of windows.


            for (uint32_t i = 0; i < numberOfWins; i++) {
                uint64_t *currWindowStart = (resorted + i * cmpWindowVectorSize);
                uint64_t *currWindowMiddle = currWindowStart + halfSize;
                //if cmpWindowVectorSize is odd, store the last element as residue
                if (cmpWindowVectorSize & 1) {                        //there is a residue remaining
                    if (isResidueStored) {                            //second residue found --> add stored and current residue to one half each.
                        if (i == 0){                                  // only once for all windows:
                            isResidueInBuffer = false;                //after processing all windows, buffer is remembered to be empty; dont set isResidueStored directly otherwise residues in one loop iteration are treated differently
                            halfSize++;                               //one half of window increases size by 1 because of residue
                            comparisons += numberOfWins;              // for each window, 1 element more
                            cmpWindowVectorSize++;
                        }
                        uint32_t posInHalfes = (i + 1) * halfSize - 1;

                        firstHalf[posInHalfes] = residue[i];
                        secondHalf[posInHalfes] = *(currWindowStart + cmpWindowVectorSize - 1); //last element of current window vector
                    }
                    else {                                          //no residue stored up to now:
                        isResidueInBuffer = true;                   //dont set isResidueStored directly, as then residues in one loop iteration would be treated differently
                        residue[i] = *(currWindowStart + cmpWindowVectorSize - 1);
                    }
                }
                uint64_t vHalfIndex = i * halfSize;         // index at which values of window i are to be stored.
                for (uint32_t j = 0; j < halfSize; j++) {   // splitting the windows on firstHalf and secondHalf
                    firstHalf[vHalfIndex + j] = *(currWindowStart + j);
                    secondHalf[vHalfIndex + j] = *(currWindowMiddle + j);
                }
            }
            if (comparisons > 0) {          // maximums are not yet found
                //compare: a-b =c and then MSB(c) =d
                uint64_t *c = SUB(firstHalf, secondHalf, comparisons);
                uint64_t *d = MSB(proxy, c, comparisons);
                //MUX: returns for each position i: firstHalf[i] if d[i] = 0; secondHalf[i] if d[i] = 1
                resorted = MUX(proxy, firstHalf, secondHalf, d, comparisons);
            }
            //prepare next round:
            cmpWindowVectorSize = halfSize;
            isResidueStored = isResidueInBuffer;
        }
        for (uint32_t m = 0; m < numberOfWins; m++){
            maxElements[m] = resorted[m];
        }
        return maxElements;
    }
    else if ( proxy->getPRole() == HELPER) {
        numberOfWins = matrix_size / window_length;
        /**Compares values in a given window by splitting the window in two halves and comparing each value to its counterpart at the same position in the other half.
        If size of the given windowVector is odd, there will be a residue, which is stored in residue. */
        while (cmpWindowVectorSize > 1) {
            uint32_t halfSize = cmpWindowVectorSize / 2;
            if (cmpWindowVectorSize & 1) {                   //there is an residue remaining
                if (isResidueStored) {                            //second residue found --> add stored and current residue each to one half.
                    isResidueInBuffer = false;                    //after processing all windows, buffer is remembered to be empty; dont set isResidueStored directly otherwise residues in one loop iteration are treated differently
                    halfSize++;                                   //one half of window increases size by 1 because of residue
                }
                else {
                    isResidueInBuffer = true;
                }
            }

            uint32_t vectorLength = static_cast<uint32_t>(floor(halfSize * numberOfWins));
            if (vectorLength > 0) {          // maximums are not yet found
                //compare: a-b =c and then MSB(c) =d
                MSB(proxy, nullptr, vectorLength);
                //MUX:
                MUX(proxy, nullptr, nullptr, nullptr, vectorLength);
            }
            //prepare next round:
            cmpWindowVectorSize = halfSize;
            isResidueStored = isResidueInBuffer;
        }
        return nullptr;
    }
    return nullptr;
}

/**
 * Method for private computation of the RELU function.
 * @param proxy
 * @param x - variable x for which to compute RELU(x)
 * @return
 */
uint64_t RELU(Party* proxy, uint64_t x){
    uint64_t K = (RING_N>>1); // N is the ring size - 1 = 2^64 -1

    if (proxy->getPRole() == P1 ||  proxy->getPRole() == P2) {
        uint64_t commonValues[4];
        for (unsigned long & commonValue : commonValues) {
            commonValue = proxy->generateCommonRandom() | 0x1; // common values must be odd
        }
        // create even random shares
        uint64_t e[] = {proxy->generateRandom() & EVEN_MASK, proxy->generateRandom() & EVEN_MASK};

        // init
        int f = proxy->generateCommonRandom() & 0x1;
        int g = proxy->generateCommonRandom() & 0x1;
        int h = proxy->generateCommonRandom() & 0x1;

        // make the shares more random by adding i to e_i for S_i:
        e[f] += proxy->getPRole();

        uint64_t t = x & K; // get first L-1 bit of the share

        uint64_t d = MOC(proxy, t);
        uint64_t z = x - d;

        // compute parts of a, b and c:
        uint64_t values[6];
        values[0] = proxy->getPRole() * f * (K+1) - z;                  // a_0
        values[1] = proxy->getPRole() * (1 - f) * (K+1) - z;            // a_1
        values[2] = (x + e[g]) * commonValues[0];                       // b_0
        values[3] = (x + e[1 - g]) * commonValues[1];                   // b_1
        values[4] = (e[h]) * commonValues[2];                           // c_0
        values[5] = (e[1 - h]) * commonValues[3];                       // c_1

        // proxy sends a,b and c to HELPER:
        unsigned char *ptr_out = proxy->getBuffer1();
        addVal2CharArray(values, &ptr_out, 6);
        Send(proxy->getSocketHelper(), proxy->getBuffer1(), 6 * 8);

        // receive fresh share from helper
        Receive(proxy->getSocketHelper(), proxy->getBuffer1(),
                8 * 8); // order is ab[0], ab[1], ab[2], ab[3], ac[0], ac[1], ac[2], ac[3]

        ptr_out = proxy->getBuffer1();
        uint64_t ab[4];
        uint64_t ac[4];
        convert2Array(&ptr_out, &ab[0], 4);
        convert2Array(&ptr_out, &ac[0], 4);

        uint64_t em;
        if (g == h){
            uint64_t r2_inverse = getModularInverse(commonValues[2]);
            em = ac[2*f] * r2_inverse;
        }else{
            uint64_t r3_inverse = getModularInverse(commonValues[3]);
            em = ac[2*f+1] * r3_inverse;
        }

        uint64_t r0_inverse = getModularInverse(commonValues[0]);
        uint64_t xm = ab[2*f] * r0_inverse - em;
        z = x - xm;
        return z;
    }
    else if (proxy->getPRole() == HELPER) {


        MOC(proxy, 0);

        Receive(proxy->getSocketP1(), proxy->getBuffer1(), 6 * 8);
        Receive(proxy->getSocketP2(), proxy->getBuffer2(), 6 * 8);
        unsigned char *ptr = proxy->getBuffer1();
        unsigned char *ptr2 = proxy->getBuffer2();

        uint64_t reconstructedVals[6];
        for (uint i = 0; i < 6; i++) {
            reconstructedVals[i] = (convert2Long(&ptr) + convert2Long(&ptr2));
            // 6 values expected per party --> first two (a values) are reconstructed and divided by K
            if (i < 2){
                reconstructedVals[i] /= (K+1);
            }

        }
        // reconstruct values a, b and c obtained from proxy 1 and 2
        uint64_t ab [4];
        uint64_t ac [4];
        ab[0] = reconstructedVals[0] * reconstructedVals[2]; // a_0 * b_0
        ab[1] = reconstructedVals[0] * reconstructedVals[3]; // a_0 * b_1
        ab[2] = reconstructedVals[1] * reconstructedVals[2]; // a_1 * b_0
        ab[3] = reconstructedVals[1] * reconstructedVals[3]; // a_1 * b_1

        ac[0] = reconstructedVals[0] * reconstructedVals[4]; // a_0 * c_0
        ac[1] = reconstructedVals[0] * reconstructedVals[5]; // a_0 * c_1
        ac[2] = reconstructedVals[1] * reconstructedVals[4]; // a_1 * c_0
        ac[3] = reconstructedVals[1] * reconstructedVals[5]; // a_1 * c_1
        //cout << "ab and ac are calculated" << endl;

        // create fresh shares
        uint64_t share1 [8];
        uint64_t share2 [8];
        for (uint i = 0; i < 4; i++) {
            uint64_t tmp = proxy->generateRandom();
            share1[i] = tmp;
            share2[i] = ab[i] - tmp;

            tmp = proxy->generateRandom();
            share1[i+4] = tmp;
            share2[i+4] = ac[i] - tmp;
        }

        // send shares to proxy 1 and 2
        unsigned char *ptr_back = proxy->getBuffer1();
        unsigned char *ptr_back2 = proxy->getBuffer2();
        addVal2CharArray(share1, &ptr_back, 8);
        addVal2CharArray(share2, &ptr_back2, 8);

        thread thr1 = thread(Send, proxy->getSocketP1(), proxy->getBuffer1(), 8 * 8);
        thread thr2 = thread(Send, proxy->getSocketP2(), proxy->getBuffer2(), 8 * 8);

        thr1.join();
        thr2.join();
        return 0;
    }
    return -1;
}

/**
 * Method for private computation of the RELU function.
 * @param proxy
 * @param x - vector of variables for which to compute RELU
 * @param size - size of vector x
 * @return vector of resulting values for each position in x
 */
uint64_t* RELU(Party* proxy, uint64_t* x, uint64_t size){
    uint64_t K = (RING_N>>1); // N is the ring size - 1 = 2^64 -1

    if (proxy->getPRole() == P1 ||  proxy->getPRole() == P2) {
        uint64_t commonValues[4*size];
        for (unsigned long & commonValue : commonValues) {
            commonValue = proxy->generateCommonRandom() | 0x1; // common values must be odd
        }

        uint64_t* e = new uint64_t [2*size]; // 2 shares per value to compute Relu for
        int* f = new int[size];
        int* g = new int[size];
        int* h = new int[size];
        uint64_t* t = new uint64_t [size];
        for (uint64_t i = 0; i < size; i++){
            // create even random shares: store first all 'e_0', then all 'e_1'
            e[i] = proxy->generateRandom() & EVEN_MASK;
            e[i+size] = proxy->generateRandom() & EVEN_MASK;

            // init f, g, and h per value in x
            f[i] = proxy->generateCommonRandom() & 0x1;
            g[i] = proxy->generateCommonRandom() & 0x1;
            h[i] = proxy->generateCommonRandom() & 0x1;

            // make the shares more random by adding i to e_i for S_i:
            // is 0*size or 1*size --> add pRole to either according e_0 or e_1
            e[i + f[i]*size] += proxy->getPRole();

            t[i] = x[i] & K; // get first L-1 bit of the share
        }
        uint64_t* d = MOC(proxy, t, size);

        uint64_t* z = new uint64_t [size];
        uint64_t values[6*size];
        for(uint64_t j = 0; j<size; j++){
            z[j] = x[j] - d[j];

            // compute parts of a, b and c:
            values[j*6] = proxy->getPRole() * f[j] * (K+1) - z[j];                  // a_0
            values[j*6 + 1] = proxy->getPRole() * (1 - f[j]) * (K+1) - z[j];        // a_1
            values[j*6 + 2] = (x[j] + e[j + g[j]*size]) * commonValues[j*4];        // b_0
            values[j*6 + 3] = (x[j] + e[j + (1-g[j])*size]) * commonValues[j*4+1];  // b_1
            values[j*6 + 4] = (e[j + h[j]*size]) * commonValues[j*4+2];             // c_0
            values[j*6 + 5] = (e[j + (1-h[j])*size]) * commonValues[j*4+3];         // c_1
        }
        // proxy sends a,b and c to HELPER:
        unsigned char *ptr_out = proxy->getBuffer1();
        addVal2CharArray(values, &ptr_out, 6*size);
        Send(proxy->getSocketHelper(), proxy->getBuffer1(), 6 * size * 8);

        // receive fresh share from helper
        Receive(proxy->getSocketHelper(), proxy->getBuffer1(),
                8 * size * 8); // order is all (ab[0], ab[1], ab[2], ab[3]), then all (ac[0], ac[1], ac[2], ac[3])

        ptr_out = proxy->getBuffer1();
        uint64_t ab[4*size];
        uint64_t ac[4*size];
        convert2Array(&ptr_out, &ab[0], 4*size);
        convert2Array(&ptr_out, &ac[0], 4*size);

        for(uint64_t i = 0; i<size; i++) {
            uint64_t em;
            if (g[i] == h[i]) {
                uint64_t r2_inverse = getModularInverse(commonValues[i*4 + 2]);
                em = ac[i*4 + 2 * f[i]] * r2_inverse;
            } else {
                uint64_t r3_inverse = getModularInverse(commonValues[i*4 + 3]);
                em = ac[i*4 + 2 * f[i] + 1] * r3_inverse;
            }

            uint64_t r0_inverse = getModularInverse(commonValues[i*4]);
            uint64_t xm = ab[i*4 + 2 * f[i]] * r0_inverse - em;
            z[i] = x[i] - xm;
        }
        return z;
    }
    else if (proxy->getPRole() == HELPER) {
        MOC(proxy, 0, size);

        Receive(proxy->getSocketP1(), proxy->getBuffer1(), 6 * size * 8);
        Receive(proxy->getSocketP2(), proxy->getBuffer2(), 6 * size * 8);
        unsigned char *ptr = proxy->getBuffer1();
        unsigned char *ptr2 = proxy->getBuffer2();

        uint64_t reconstructedVals[6*size];
        for (uint i = 0; i < 6*size; i++) {
            reconstructedVals[i] = (convert2Long(&ptr) + convert2Long(&ptr2));
            // 6 values expected per party --> first two (a values) are reconstructed and divided by K
            if ((i % 6) < 2){
                reconstructedVals[i] /= (K+1);
            }

        }
        // reconstruct values a, b and c obtained from proxy 1 and 2
        uint64_t ab [4*size];
        uint64_t ac [4*size];
        for (uint i = 0; i < size; i++) {
            ab[i*4] = reconstructedVals[i*6] * reconstructedVals[i*6 + 2]; // a_0 * b_0
            ab[i*4 + 1] = reconstructedVals[i*6] * reconstructedVals[i*6 + 3]; // a_0 * b_1
            ab[i*4 + 2] = reconstructedVals[i*6 + 1] * reconstructedVals[i*6 + 2]; // a_1 * b_0
            ab[i*4 + 3] = reconstructedVals[i*6 + 1] * reconstructedVals[i*6 + 3]; // a_1 * b_1

            ac[i*4] = reconstructedVals[i*6] * reconstructedVals[i*6 + 4]; // a_0 * c_0
            ac[i*4 + 1] = reconstructedVals[i*6] * reconstructedVals[i*6 + 5]; // a_0 * c_1
            ac[i*4 + 2] = reconstructedVals[i*6 + 1] * reconstructedVals[i*6 + 4]; // a_1 * c_0
            ac[i*4 + 3] = reconstructedVals[i*6 + 1] * reconstructedVals[i*6 + 5]; // a_1 * c_1
        }
        // create fresh shares
        uint64_t share1 [8*size];
        uint64_t share2 [8*size];
        uint64_t* tmp = convert2uint64(random_1D_data(proxy, 8*size), size);
        for (uint i = 0; i < 4*size; i++) {
            share1[i] = tmp[i];
            share2[i] = ab[i] - tmp[i];

            share1[i+4*size] = tmp[i+4*size];
            share2[i+4*size] = ac[i+4*size] - tmp[i+4*size];
        }

        // send shares to proxy 1 and 2
        unsigned char *ptr_back = proxy->getBuffer1();
        unsigned char *ptr_back2 = proxy->getBuffer2();
        addVal2CharArray(share1, &ptr_back, 8*size);
        addVal2CharArray(share2, &ptr_back2, 8*size);

        thread thr1 = thread(Send, proxy->getSocketP1(), proxy->getBuffer1(), 8 * size * 8);
        thread thr2 = thread(Send, proxy->getSocketP2(), proxy->getBuffer2(), 8 * size * 8);

        thr1.join();
        thr2.join();
        return nullptr;
    }
    return nullptr;
}

/**
 * Method for private computation of the derivative of the RELU function.
 * @param proxy
 * @param x - variable x for which to compute RELU'(x), the derivative of the RELU function.
 * @return
 */
uint64_t DRELU(Party* proxy, uint64_t x){
    uint64_t K = (RING_N>>1); // N is the ring size - 1 = 2^64 -1
    // K is 2^63 - 1
    uint8_t exchangingBit = 2;
    if (proxy->getPRole() == P1 ||  proxy->getPRole() == P2) {

        // init
        uint64_t f = proxy->generateCommonRandom() & 0x1; // generate common random bit

        uint64_t t = x & K; // get first L-1 bit of the share
        K += 1; // increase K by 1 K is 2^63
        uint64_t d = MOC(proxy, t);
        uint64_t z = x - d;


        // compute parts of a:
        uint64_t values[exchangingBit];
        values[0] = proxy->getPRole() * f * K - z;       // a_0
        values[1] = proxy->getPRole() * (1 - f) * K - z; // a_1

        // proxy sends a to HELPER:
        unsigned char *ptr = proxy->getBuffer1();
        addVal2CharArray(values, &ptr, exchangingBit);
        Send(proxy->getSocketHelper(), proxy->getBuffer1(), exchangingBit * 8);

        // receive fresh share from helper
        Receive(proxy->getSocketHelper(), proxy->getBuffer1(), exchangingBit * 8); // order is share of a[0], a[1]
        ptr = proxy->getBuffer1();
        z = proxy->getPRole() -  convert2Long(&ptr);
        if (f)  // if f is 1  we get the next long value in the buffer.
            z = proxy->getPRole() - convert2Long(&ptr);

        return z;
    }
    else if (proxy->getPRole() == HELPER) {
        K += 1;
        MOC(proxy, 0);

        Receive(proxy->getSocketP1(), proxy->getBuffer1(), exchangingBit * 8);
        Receive(proxy->getSocketP2(), proxy->getBuffer2(), exchangingBit * 8);
        unsigned char *ptr1 = proxy->getBuffer1();
        unsigned char *ptr2 = proxy->getBuffer2();

        uint64_t reconstructedVals[exchangingBit];
        reconstructedVals[0] = (convert2Long(&ptr1) + convert2Long(&ptr2)) / K;
        reconstructedVals[1] = (convert2Long(&ptr1) + convert2Long(&ptr2)) / K;

        //reassign buffer because ptr1 and ptr2 were incremented by convert2Long calls.
        ptr1 = proxy->getBuffer1();
        ptr2 = proxy->getBuffer2();

        uint64_t tmp = proxy->generateRandom();
        addVal2CharArray(tmp,&ptr1);
        addVal2CharArray(reconstructedVals[0]-tmp,&ptr2);
        tmp = proxy->generateRandom();
        addVal2CharArray(tmp,&ptr1);
        addVal2CharArray(reconstructedVals[1]-tmp,&ptr2);

        thread thr1 = thread(Send, proxy->getSocketP1(), proxy->getBuffer1(), exchangingBit * 8);
        thread thr2 = thread(Send, proxy->getSocketP2(), proxy->getBuffer2(), exchangingBit * 8);
        thr1.join();
        thr2.join();

        return 0;
    }
    return -1;
}

/**
 * Method for private computation of the derivative of the RELU function.
 * @param proxy
 * @param x - vector of variables for which to compute RELU' (the derivative of the RELU function).
 * @param size - size of vector x
 * @return
 */
uint64_t* DRELU(Party* proxy, uint64_t* x, uint64_t size){
    uint64_t K = (RING_N>>1); // N is the ring size - 1 = 2^64 -1
    // K is 2^63 - 1
    uint8_t exchangingBit = 2 * size;
    role pRole = proxy->getPRole();
    if (pRole == P1 || pRole == P2) {
        // init
        uint64_t* f = new uint64_t [size];
        uint64_t* t = new uint64_t [size];
        for (uint64_t counter = 0; counter < size; counter++){
            f[counter] = proxy->generateCommonRandom() & 0x1;
            t[counter] = x[counter] & K; // get first L-1 bit of the share
        }

        uint64_t* d = MOC(proxy, t, size);

        uint64_t* z = new uint64_t [size];
        for (uint64_t counter = 0; counter < size; counter++){
            z[counter] = x[counter] - d[counter];
        }
        K += 1; // increase K by 1 K is 2^63

        // compute parts of a:
        uint64_t values[exchangingBit];
        for (uint64_t exBit = 0; exBit < exchangingBit; exBit++){
            uint64_t c = floor(exBit / 2);
            values[exBit] = pRole * f[c] * K - z[c];         // a_0 of thee exBit/2-th variable
            exBit++;
            c = floor(exBit / 2);
            values[exBit] = pRole * (1 - f[c]) * K - z[c];   // a_1
        }

        // proxy sends a to HELPER:
        unsigned char *ptr = proxy->getBuffer1();
        addVal2CharArray(values, &ptr, exchangingBit);
        Send(proxy->getSocketHelper(), proxy->getBuffer1(), exchangingBit * 8);

        // receive fresh share from helper
        Receive(proxy->getSocketHelper(), proxy->getBuffer1(), exchangingBit * 8); // order is share of a[0], a[1]
        ptr = proxy->getBuffer1();

        for (uint64_t c = 0; c < size; c++){
            z[c] = pRole - convert2Long(&ptr);
            if (f[c])  // if f is 1  we get the next long value in the buffer.
                z[c] = pRole - convert2Long(&ptr);
        }
        return z;
    }
    else if (proxy->getPRole() == HELPER) {
        MOC(proxy, 0, size);

        K += 1;
        Receive(proxy->getSocketP1(), proxy->getBuffer1(), exchangingBit * 8);
        Receive(proxy->getSocketP2(), proxy->getBuffer2(), exchangingBit * 8);
        unsigned char *ptr1 = proxy->getBuffer1();
        unsigned char *ptr2 = proxy->getBuffer2();

        uint64_t reconstructedVals[exchangingBit];
        for (uint64_t v = 0; v < exchangingBit; v++) {
            reconstructedVals[v] = (convert2Long(&ptr1) + convert2Long(&ptr2)) / K;
            v++;
            reconstructedVals[v] = (convert2Long(&ptr1) + convert2Long(&ptr2)) / K;
        }
        //reassign buffer because ptr1 and ptr2 were incremented
        ptr1 = proxy->getBuffer1();
        ptr2 = proxy->getBuffer2();

        uint64_t* tmp = convert2uint64(random_1D_data(proxy, exchangingBit), exchangingBit);  // values for P1
        addVal2CharArray(tmp, &ptr1, exchangingBit);

        uint64_t * share = new uint64_t [exchangingBit];                                // values for P2
        for (uint64_t v = 0; v < exchangingBit; v++) {
            share[v] = reconstructedVals[v] - tmp[v];
        }
        addVal2CharArray(share, &ptr2, exchangingBit);

        thread thr1 = thread(Send, proxy->getSocketP1(), proxy->getBuffer1(), exchangingBit * 8);
        thread thr2 = thread(Send, proxy->getSocketP2(), proxy->getBuffer2(), exchangingBit * 8);
        thr1.join();
        thr2.join();

        return nullptr;
    }
    return nullptr;
}


uint64_t DIV(Party* proxy, uint64_t a, uint64_t b) {
    //compare with Algorithm 8 Division by SecureNN
    if (proxy->getPRole() == P1 || proxy->getPRole() == P2) {

        uint64_t result = 0;
        uint64_t u = proxy->createShare(0);
        for (int16_t i = FRAC - 1; i >= 0; i--) {
            //___3. step
            uint64_t y = a << i; // equals a * 2^i not as pseudocode in paper
            uint64_t z = y - u - b;

            //___4. step (DRELU)
            uint64_t drelu = DRELU(proxy, z);

            //___5. step (MUL)
            uint64_t mul = MUL(proxy, drelu, b);

            //___6. step and 8. step (combined)
            result += (drelu << i);
            //___7. step
            u += mul;
            cout << "current result: " << bitset<L_BIT>(REC(proxy, result)) << endl;
        }
        return result;
    } else if (proxy->getPRole() == HELPER) {
        for (int16_t i = FRAC-1; i >= 0; i--) {
            //___4. step (DRELU)
            DRELU(proxy, 0);

            //___5. step (MUL)
            MUL(proxy, 0, 0);
        }
        return 0;
    }
    return -1;
}

uint64_t DIV_BIT(Party* proxy, uint64_t a, uint64_t b) {
    cout << "a = " << bitset<L_BIT>(REC(proxy, a)) << "; b = " << bitset<L_BIT>(REC(proxy, b)) << endl;
    if (proxy->getPRole()  == P1 || proxy->getPRole()  == P2){
        uint64_t remainder = 0;
        cout << "remainder = " << bitset<L_BIT>(REC(proxy, remainder)) << endl;
        uint64_t result = 0;
        for (int16_t i = L_BIT-1; i >= 0; i--){
            remainder = (remainder << i);
            remainder |= (a & (1 << i));
            cout << "remainder = " << bitset<L_BIT>(remainder) << endl;
            uint64_t z = remainder - b;
            cout << "z = " << bitset<L_BIT>(z) << endl;
            uint64_t quotient = MSB(proxy, z); // MSB of z is intermediate result/quotient
            cout << "MSB(z) = " << bitset<L_BIT>(REC(proxy, quotient)) << endl;

            // subtract to calc remainder
            if(quotient){
                remainder -= b;
                result |= (1 << i);
            }
        }
        return result;
    }
    else if (proxy->getPRole() == HELPER) {
        for (int16_t i = L_BIT-1; i >= 0; i--) {
            MSB(proxy, 0);
        }
        return 0;
    }
    return 0;
}

// LAYER FUNCTIONS:

/**
 * Implements the function of a convolutional layer (CL) using ReLU and then Maxpool with a 2x2 filter as activation function.
 * @param proxy
 * @param input data on which convolution is performed using provided kernels. TODO make method to be called for batch of inputs
 * @param i_dim size of input in symmetric shape. So input has size of i_size x i_size TODO what about third dimension (channel param)?
 * @param kernel all kernel to be used for convolution. Those are all expected to be vectors of length k_size x k_size.
 * @param k_dim dimension of the kernel TODO first only kernel_size * kernel_size * 1 but later also for 3 channel
 * @param k_number number of kernels to be used, which will define the third dimension of the outcome.
 * @param stride step size of each kernel to shift per iteration
 * @return Output of the input convoluted by the given kernels.
 *         Output size will be: ZxZxk_number with Z = floor((28 - 5)/stride) + 1
 */
 /*
uint64_t*** CL(Party* proxy, uint64_t** input, uint64_t i_dim, uint64_t** kernel, uint32_t k_dim, uint32_t k_number, uint8_t stride){
    uint32_t k_size = k_dim * k_dim;
    uint32_t conv_dim = static_cast<uint32_t>(floor((i_dim - k_dim) / stride) + 1);
    uint32_t conv_len = conv_dim * conv_dim;
    uint32_t out_size = conv_dim / 2; // divide by 2 because maxpool has window size of 2 --> reduces by half.
    uint32_t last_start = i_dim - k_dim + 1; // TODO adapt
    cout << "CL1: started init: conv_dim= " << conv_dim << " , outSize= " << out_size << ", lastStart = " << last_start << endl;

    // stretch the input for vectorized MATVECMUL
    if (proxy->getPRole() == P1 || proxy->getPRole() == P2) {
        uint64_t ** stretched_input = INC(input, conv_dim, last_start, k_dim, stride);
        //print2DArray("Stretched matrix for MATVECMUL", convert2double(REC(proxy, stretched_input, conv_len, k_size), conv_len, k_size), conv_len, k_size);

        // convolution:
        uint64_t ***conv_layer = new uint64_t **[k_number];
        for (uint32_t k = 0; k < k_number; k++) {                          // for each kernel
            // multiply stretched_input with kernel
            uint64_t *conv_result = MATVECMUL(proxy, stretched_input, kernel[k], conv_len, k_size);
            //print1DArray("Result of MATVECMUL", convert2double(REC(proxy, conv_result, conv_len), conv_len),conv_len);

            uint64_t* conv_activated = new uint64_t [conv_len];
            // ACTIVATION:
            // ReLU
            for (uint32_t i = 0; i < conv_len; i++) {
                conv_activated[i] = RELU(proxy, conv_result[i]);
                //cout << "RELU of " << convert2double(REC(proxy, conv_result[i])) << " is " << convert2double(REC(proxy, conv_activated[i])) << endl;
            }

            // Maxpool:
            uint64_t *reduced = MAX(proxy, conv_activated, conv_dim, conv_dim, 2);
            // bring result in matrix shape
            conv_layer[k] = new uint64_t *[out_size];
            for (uint32_t row = 0; row < out_size; row++) {
                conv_layer[k][row] = new uint64_t[out_size];
                for (uint32_t col = 0; col < out_size; col++) {
                    conv_layer[k][row][col] = reduced[row * out_size + col];
                }
            }
        }
        return conv_layer;
    }
    else if (proxy->getPRole() == HELPER){
        // convolution:
        for (uint32_t k = 0; k < k_number; k++) {
            MATVECMUL(proxy, nullptr, nullptr, conv_dim * conv_dim * k_size, 0);

            // ACTIVATION:
            // ReLU
            for (uint32_t i = 0; i < conv_dim * conv_dim; i++) {
                RELU(proxy, 0);
            }
            // Maxpool:
            MAX(proxy, nullptr, conv_dim, conv_dim, 2);
        }
        return nullptr;
    }
    return nullptr;

}*/

/**
 * Implements the function of a convolutional layer (CL) using ReLU and then Maxpool with a 2x2 filter as activation function.
 * The input is supposed to be in 3D shape such as images/matrices with one or several channel.
 * @param proxy
 * @param input data on which convolution is performed using the provided kernels.
 * Input has a shape of i_width x i_height x i_channel.
 * @param i_channel number of channels of the input; must be > 0
 * @param i_height size of input in the other dimension; must be > 0
 * @param i_width size of input in one dimension; must be > 0
 * @param kernel vector containing all kernel to be used for convolution.
 * The parameter kernel has shape output_channel * i_channel x k_dim * k_dim.
 * So each channel of a kernel is represented by a vector of length k_dim * k_dim.
 * The kernel's number of channels is defined by i_channel and all channel of a kernel are stored in 'kernel' sequentially.
 * Then all the channel of the next kernel follow.
 * @param k_dim dimension of the symmetric kernel; must be > 0
 * @param output_channel number of kernels with shape i_channel x k_dim x k_dim; must be > 0
 * For each kernel there will be one output channel in the result, where the value at the according location is
 * the sum of the single multiplications between input and kernel of same channel number.
 * @param stride step size of each kernel to shift per iteration; must be > 0
 * @return Output of the input convoluted by the given kernels.
 *         Output size will be: number with Z = floor((28 - 5)/stride) + 1
 */
uint64_t*** CL(Party* proxy, uint64_t*** input, uint32_t i_channel, uint32_t i_height, uint32_t i_width, uint64_t*** kernel, uint32_t k_dim, uint32_t output_channel, uint32_t stride){
    uint32_t k_size = k_dim * k_dim;
    uint32_t conv_width = static_cast<uint32_t>(floor((i_width - k_dim) / stride) + 1);
    uint32_t conv_height = static_cast<uint32_t>(floor((i_height - k_dim) / stride) + 1);
    uint32_t conv_len = conv_width * conv_height;

    uint32_t out_width = conv_width / 2; // divide by 2 because maxpool has window size of 2 --> reduces by half.
    uint32_t out_height = conv_height / 2;
    // stretch the input for vectorized MATVECMUL
    if (proxy->getPRole() == P1 || proxy->getPRole() == P2) {
        uint64_t *** stretched_input = INC(input, i_channel, i_width, i_height, k_dim, stride);
        //print2DArray("Stretched matrix for MATVECMUL", convert2double(REC(proxy, stretched_input[0], conv_len, k_size), conv_len, k_size), conv_len, k_size);

        // convolution:
        uint64_t ***conv_layer = new uint64_t **[output_channel];
        for (uint32_t k = 0; k < output_channel; k++) {                          // for each kernel
            // multiply stretched_input with kernel (for all channel)
            cout << "calling MATVECMUL" << endl;
            uint64_t **conv_result = MATVECMUL(proxy, stretched_input, kernel[k], i_channel, conv_len, k_size);
            cout << "returned from MATVECMUL" << endl;
            //print2DArray("Result of MATVECMUL", convert2double(REC(proxy, conv_result, i_channel, conv_len), i_channel, conv_len), i_channel, conv_len);

            // sum up the channel results to obtain one output_channel for all input channel
            uint64_t* summed_channel_conv = ADD(proxy, conv_result, i_channel, conv_len);
            // ACTIVATION:
            uint64_t* conv_activated = new uint64_t [conv_len];
            // ReLU
            for (uint32_t i = 0; i < conv_len; i++) {
                conv_activated[i] = RELU(proxy, summed_channel_conv[i]);
            }

            // Maxpool:
            uint64_t *reduced = MAX(proxy, conv_activated, conv_width, conv_height, 2);
            // bring result in matrix shape
            conv_layer[k] = new uint64_t *[out_height];
            for (uint32_t row = 0; row < out_width; row++) {
                conv_layer[k][row] = new uint64_t[out_height];
                for (uint32_t col = 0; col < out_width; col++) {
                    conv_layer[k][row][col] = reduced[row * out_height + col];
                }
            }
        }
        return conv_layer;
    }
    else if (proxy->getPRole() == HELPER){
        // convolution:
        for (uint32_t k = 0; k < output_channel; k++) {
            MATVECMUL(proxy, nullptr, nullptr, 0, i_channel * conv_len * k_size, 0);

            // ACTIVATION:
            // ReLU
            for (uint32_t i = 0; i < conv_len; i++) {
                RELU(proxy, 0);
            }
            // Maxpool:
            MAX(proxy, nullptr, conv_width,  conv_height, 2);
        }
        return nullptr;
    }
    return nullptr;

}
/**
 *
 * Increase the size of input by stretching its values per channel so that
 * in each row all values by a kernels channel with dimension k_dim * k_dim are stored.
 * @param input input matrix to be stretched accordingly. Its size is channel x conv_width x conv_height
 * @param channel number of channels (size of first dimension) for input
 * @param height height of the matrix after convolution
 * @param width width of the matrix after convolution
 * @param k_dim kernel dimension. Each kernel's channel holds k_dim * k_dim values in total.
 * @param stride step size of the kernels
 * @return the stretched input image with shape channel x conv_width*conv_height x k_dim*k_dim
 *         with conv_width = (width - k_dim)/stride + 1 and conv_height = (height - k_dim)/stride + 1
 */
uint64_t ***INC(uint64_t ***input, uint32_t channel, uint32_t height, uint32_t width, uint32_t k_dim,
                uint32_t stride) {
    uint32_t k_size = k_dim * k_dim;
    uint32_t conv_height = (height - k_dim)/stride + 1;
    uint32_t last_row_start = height - k_dim + 1;
    uint32_t last_col_start = width - k_dim + 1;
    // stretch the input for vectorized MATVECMUL
    uint64_t ***stretched_input = new uint64_t **[channel];
    for(uint32_t c = 0; c < channel; c++){
        stretched_input[c] = new uint64_t *[height * width];
        for(uint32_t row = 0; row < last_row_start; row += stride) {
            for(uint32_t col = 0; col < last_col_start; col += stride) {
                uint32_t index = (row * conv_height + col) / stride;
                stretched_input[c][index] = new uint64_t[k_size];
                // shift kernel over all cols in according row
                for (uint32_t l = 0; l < k_size; l++) {
                    //                    row           columns have same length as kernel size
                    stretched_input[c][index][l] = input[c][row + l / k_dim][col + l % k_dim];
                }
            }
        }
    }
    return stretched_input;
}


/**
 * Implements functionality of a fully connected layer (FCL) TODO which activation?
 * @param proxy
 * @param input the input matrices to be fully connected to the output nodes of shape i_dim * i_dim * i_number
 * @param i_dim dimension of the input matrices
 * @param i_number number of the input matrices
 * @param weights to be used must be of shape node_number * (i_dim * i_dim * i_number) as provided by Chameleon files
 * @param node_number number of output nodes of this layer
 * @return the output layer that have been computed by using the dot product between input values and weights. It will be of length node_number.
 */
uint64_t* FCL(Party* proxy, uint64_t*** input, uint64_t i_dim, uint32_t i_number, uint64_t** weights, uint32_t node_number){
    uint64_t out_size = i_dim * i_dim * i_number; //first dimension of output
    if (proxy->getPRole() == P1 || proxy->getPRole() == P2){
        uint64_t *flattened = FLT(input, i_dim, i_number);
        uint64_t ** w = transpose(weights, out_size, node_number);
        uint64_t *output = new uint64_t [node_number];
       // uint64_t *stretched_flattened = new uint64_t [out_size*node_number];
        //memcpy(stretched_flattened, flattened, node_number);
        for (uint32_t node = 0; node < node_number; node++){ //TODO use metes stretching method to make for loop obsolete.
            output[node] = DP(proxy, flattened, w[node], out_size);
        }
        return output;
    }
    else if (proxy->getPRole() == HELPER){
        for (uint32_t node = 0; node < node_number; node++){ //TODO use metes stretching method to make for loop obsolete.
            DP(proxy, nullptr, nullptr, out_size);
        }
        return nullptr;
    }
    else{
        cout << "Error: unknown proxy role for fully connected layer." << endl;
        return nullptr;
    }
}


/**
 * Flatten the values of a matrix concatenating its values so that a single vector is the result
 * @param images matrices to be flattened to one single vactor
 * @param i_dim dimension of a single matrix in either direction (symmetric shape)
 * @param i_number number of matrices in images
 * @return the flattened vector of length i_dim * i_dim * i_number
 */
uint64_t * FLT(uint64_t*** images, uint32_t i_dim, uint32_t i_number){
    uint64_t i_size = i_dim * i_dim;
    // stretch the input for vectorized MATVECMUL
    uint64_t * flattened = new uint64_t [i_size * i_number];
    for (uint32_t i = 0; i < i_number; i++){
        for (uint32_t el = 0; el < i_size; el ++) {
            flattened[el + i * i_size] = images[i][el / i_dim][el % i_dim];
        }
    }
    return flattened;
}

uint64_t ** transpose(uint64_t** matrix, uint32_t rows, uint32_t cols){
    uint64_t ** res = new uint64_t *[cols];
    for (int c = 0; c < cols; ++c) {
        res[c] = new uint64_t [rows];
        for(int r = 0; r<rows; r++) {
            res[c][r] = matrix[r][c];
        }
    }
    return res;
}

#endif //PPAUC_CNN_H

