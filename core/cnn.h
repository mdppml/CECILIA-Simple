//
// Created by Mete Akgun on 17.01.22.
//

#ifndef PPAUC_CNN_H
#define PPAUC_CNN_H

#include "core.h"
#include "../utils/flib.h"
#include "bitset"

/**
 * Vectorized subtraction: subtract elementwise b from a (a-b)
 * @param a from values of this vector the values of b will be subtracted
 * @param b the vector whose values are going to be subtracted from a
 * @param length of a and b (must have same length)
 * @return vector d where each element is the result of the according elements a-b.
 */
uint64_t* SUB(const uint64_t *a, const uint64_t *b, uint32_t length){
    auto *subtractedValues = new uint64_t [length];
    for(uint32_t i = 0; i<length; i++){
        subtractedValues[i] = *(a+i) - *(b+i);
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
 * CAUTION: only matrices up to a size of 9000 x 9000 can be granted to be processed.
 */
void RST(const uint64_t* matrix, uint32_t m_cols, uint32_t m_rows, uint32_t w_cols, uint32_t w_rows, uint64_t* resortedMatrix){
    uint32_t winSize = w_cols * w_rows;
    uint32_t numberOfWins = (m_cols * m_rows) / winSize;
    uint32_t winsPerRow = m_cols / w_cols;
    uint32_t w_count = 0;

    while(w_count < numberOfWins){
        auto windowStart = static_cast<uint32_t>(w_rows * m_cols * floor(w_count / winsPerRow) + (w_cols * (w_count & (winsPerRow-1))));
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
 * Finds the index of the maximum element from the given matrix in secret shared form.
 * @param matrix - secret share of the matrix from which the index of the maximal element shall be computed.
 * @param size - size of matrix
 * @return The index of the maximum element in matrix.
 */
uint64_t ARGMAX(Party* proxy, uint64_t *mShare, uint32_t matrix_size){
    /** MAIN IDEA:
     * As the MAX is performed, a second matrix of same length is created containing the indices only.
     * However values are selected from matrix, is also done for the indices-matrix.
     */
    int L1 = L_BIT - 1;
    uint32_t cmpVectorSize = matrix_size; //size of resulting vector after cmp, MUX.
    bool isResidueStored = false;
    bool isSecondHalfFilled = false;
    if ( proxy->getPRole() == P1 ||  proxy->getPRole() == P2) {
        uint64_t *maxElements = mShare;
        // generate matrix with indices
        auto *args = new uint64_t [matrix_size];
        uint64_t * tmp = convert2uint64(random_1D_data(proxy, matrix_size), matrix_size);
        for (uint32_t i = 0; i < matrix_size; ++i) {
            if(proxy->getPRole() == P1){
                double index = i;
                args[i] = convert2uint64(index) - tmp[i];
            }
            else{
                args[i] = tmp[i];
            }
        }
        auto maxHalfSizes = static_cast<uint32_t>(ceil(matrix_size / 2)); // ceil because residue might be added.
        uint64_t elements1 [maxHalfSizes];
        uint64_t elements2 [maxHalfSizes];
        uint64_t args1 [maxHalfSizes];
        uint64_t args2 [maxHalfSizes];
        uint64_t residue[2]; //there is at most one residual element and its index

        while (cmpVectorSize > 1) {
            auto halfSize = static_cast<uint32_t>(floor(cmpVectorSize / 2));

            if (cmpVectorSize & 0x1) {                          //there is a residue remaining
                if (isResidueStored) {                          //second residue found --> add stored and current residue to one half each.
                    halfSize++;                                 //each half of window increases size by 1 because of residues
                    memcpy(elements2, &maxElements[halfSize], (halfSize-1)*sizeof(maxElements[0]));
                    elements2[halfSize-1] = residue[0];
                    memcpy(args2, &args[halfSize], (halfSize-1)*sizeof(maxElements[0]));
                    args2[halfSize-1] = residue[1];
                    isSecondHalfFilled = true;
                } else {                                        //no residue stored up to now:
                    isResidueStored = true;
                    residue[0] = maxElements[cmpVectorSize-1]; // store last element in residue
                    residue[1] = args[cmpVectorSize-1];
                }
            }
            memcpy(elements1, maxElements, halfSize * sizeof(maxElements[0]));
            memcpy(args1, args, halfSize * sizeof(maxElements[0]));
            if(!isSecondHalfFilled){
                memcpy(elements2, &maxElements[halfSize], halfSize * sizeof(maxElements[0]));
                memcpy(args2, &args[halfSize], halfSize * sizeof(maxElements[0]));
            }

            if (halfSize > 0) {                                 // maximums not yet found
                //compare: a-b =c and then MSB(c) =d
                uint64_t *c = SUB(elements1, elements2, halfSize);
                uint64_t *d = MSB(proxy, c, halfSize);

                //MUX:
                maxElements = MUX(proxy, elements1, elements2, d, halfSize);
                args = MUX(proxy, args1, args2, d, halfSize);
                delete[] c;
                delete[] d;
            }
            //prepare next round:
            cmpVectorSize = halfSize;
            isSecondHalfFilled = false;
        }
        delete [] maxElements;
        uint64_t argmax = args[0];                          // should only contain one element at the end.
        delete [] args;
        return argmax;
    }
    else if ( proxy->getPRole() == HELPER) {
        /**Compares values in a given window by splitting the window in two halves and comparing each value to its counterpart at the same position in the other half.
        If size of the given windowVector is odd, there will be a residue, which is stored in residue. */
        while (cmpVectorSize > 1) {
            auto halfSize = static_cast<uint32_t>(floor(cmpVectorSize / 2));
            if (cmpVectorSize % 2 == 1) {                   //there is a residue remaining
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
    uint32_t cmpVectorSize = matrix_size; //size of resulting vector after cmp, MUX.
    bool isResidueStored = false;
    bool isSecondHalfFilled = false;
    if ( proxy->getPRole() == P1 ||  proxy->getPRole() == P2) {
        uint64_t *maxElements = mShare;
        auto maxHalfSizes = static_cast<uint32_t>(ceil(matrix_size / 2)); // ceil because residue might be added.
        uint64_t firstHalf [maxHalfSizes];
        uint64_t secondHalf [maxHalfSizes];
        uint64_t residue; //there is at most one residual element.

        while (cmpVectorSize > 1) {
            auto halfSize = static_cast<uint32_t>(floor(cmpVectorSize / 2));

            if (cmpVectorSize & 0x1) {                          //there is a residue remaining
                if (isResidueStored) {                          //second residue found --> add stored and current residue to one half each.
                    halfSize++;                                 //each half of window increases size by 1 because of residues
                    memcpy(secondHalf, &maxElements[halfSize], (halfSize-1)*sizeof(maxElements[0]));
                    secondHalf[halfSize-1] = residue;
                    isSecondHalfFilled = true;
                } else {                                        //no residue stored up to now:
                    isResidueStored = true;
                    residue = maxElements[cmpVectorSize-1]; // store last element in residue
                }
            }
            memcpy(firstHalf, maxElements, halfSize*sizeof(maxElements[0]));
            if(!isSecondHalfFilled){
                memcpy(secondHalf, &maxElements[halfSize], halfSize*sizeof(maxElements[0]));
            }

            if (halfSize > 0) {                                 // maximums not yet found
                //compare: a-b =c and then MSB(c) =d
                uint64_t *c = SUB(firstHalf, secondHalf, halfSize);
                uint64_t *d = MSB(proxy, c, halfSize);

                //MUX:
                maxElements = MUX(proxy, firstHalf, secondHalf, d, halfSize);
                delete[] c;
                delete[] d;
            }
            //prepare next round:
            cmpVectorSize = halfSize;
            isSecondHalfFilled = false;
        }

        uint64_t max = maxElements[0];                          // should only contain one element at the end.
        delete [] maxElements;
        return max;
    }
    else if ( proxy->getPRole() == HELPER) {
        /**Compares values in a given window by splitting the window in two halves and comparing each value to its counterpart at the same position in the other half.
        If size of the given windowVector is odd, there will be a residue, which is stored in residue. */
        while (cmpVectorSize > 1) {
            auto halfSize = static_cast<uint32_t>(floor(cmpVectorSize / 2));
            if (cmpVectorSize % 2 == 1) {                   //there is a residue remaining
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
 * @param matrix - share of the matrix where values of a window are to be expected subsequently and then values of the
 * next value will follow (so not rows are concatenated containing values of all windows over all columns but windows
 * are concatenated).
 * @param m_rows - number of rows in the matrix mShare.
 * @param m_cols - number of columns in the matrix mShare.
 * @param win_rows - number of rows of the sliding window for which maxpool is performed.
 * @param win_cols - number of columns of the sliding window for which maxpool is performed.
 * @return the maximum element per window, therefore floor(matrix_size / (win_rows*win_cols)) elements in a vector
 * return vector must be deleted if not needed anymore.
 */
uint64_t* MAX(Party* proxy, uint64_t *matrix, uint32_t m_rows, uint32_t m_cols, uint32_t win_rows, uint32_t win_cols){
    uint32_t matrix_size = m_rows * m_cols;
    uint32_t window_length = win_cols * win_rows;
    uint32_t cmpWindowVectorSize = window_length; //size of resulting vector after cmp, MUX, and it's divided by 2 is size of each half.
    bool isResidueStored = false;
    bool isResidueInBuffer = false;

    uint32_t numberOfWins;
    if ( proxy->getPRole() == P1 ||  proxy->getPRole() == P2) {
        auto *resorted = new uint64_t [matrix_size];                // RESORT matrix to have all values of a window subsequently
        RST(matrix, m_cols, m_rows, win_cols, win_rows, resorted);
        numberOfWins = matrix_size / window_length;

        auto *maxElements = new uint64_t [numberOfWins];
        uint64_t residue [numberOfWins];
        uint32_t comparisons;
        while (cmpWindowVectorSize > 1) {
            /**Compares values in a given window by splitting the window in two halves and comparing each value to
             * its counterpart at the same position in the other half. If size of the given windowVector is odd,
             * there will be a residue, which is stored in residue.*/
            uint32_t halfSize = cmpWindowVectorSize / 2;

            comparisons = halfSize * numberOfWins;
            uint64_t firstHalf [comparisons+numberOfWins]; // allocate space for 1 value more in case residue from store is used
            uint64_t secondHalf [comparisons+numberOfWins];

            for (uint32_t i = 0; i < numberOfWins; i++) {
                uint64_t *currWindowStart = (resorted + i * cmpWindowVectorSize);
                uint64_t *currWindowMiddle = currWindowStart + halfSize;
                if (cmpWindowVectorSize & 1) {                        //there is a residue remaining
                    if (isResidueStored) {                            //second residue found --> add stored and current residue to one half each.
                        if(i == 0){
                            isResidueInBuffer = false;                    //after processing all windows, buffer is remembered to be empty; don't set isResidueStored directly otherwise residues in one loop iteration are treated differently
                            halfSize++;
                            currWindowMiddle++; //halfSize was too small in line 423
                            comparisons += numberOfWins;
                        }
                        //fill secondHalf directly, because inserting with pointer is not possible without loss
                        memcpy(&secondHalf[i * halfSize], currWindowMiddle, (halfSize-1)*sizeof(resorted[0]));
                        secondHalf[i * halfSize + halfSize - 1] = residue[i];
                    }
                    else {                                          //no residue stored up to now:
                        isResidueInBuffer = true;                   //don't set isResidueStored directly, as then residues in one loop iteration would be treated differently
                        residue[i] = *(currWindowStart + cmpWindowVectorSize - 1);
                    }
                }
                memcpy(&firstHalf[i * halfSize], currWindowStart, halfSize*sizeof(resorted[0]));
                if(isResidueInBuffer or !isResidueStored){
                    //otherwise, memcpy was already called for secondHalf
                    memcpy(&secondHalf[i * halfSize], currWindowMiddle, halfSize*sizeof(resorted[0]));
                }
            }
            if (comparisons > 0) {          // maximums are not yet found
                //compare: a-b =c and then MSB(c) =d
                uint64_t *c = SUB(firstHalf, secondHalf, comparisons);
                uint64_t *d = MSB(proxy, c, comparisons);
                //MUX: returns for each position i: firstHalf[i] if d[i] = 0; secondHalf[i] if d[i] = 1
                resorted = MUX(proxy, firstHalf, secondHalf, d, comparisons);
                delete[] c;
                delete[] d;
            }
            //prepare next round:
            cmpWindowVectorSize = halfSize;
            isResidueStored = isResidueInBuffer;
        }
        if (isResidueStored) {
            uint64_t *c = SUB(resorted, residue, numberOfWins);
            uint64_t *d = MSB(proxy, c, numberOfWins);
            resorted = MUX(proxy, resorted, residue, d, numberOfWins);
            delete[] c;
            delete[] d;
        }
        for (uint32_t m = 0; m < numberOfWins; m++){
            maxElements[m] = resorted[m];
        }
        delete [] resorted;
        return maxElements;
    }
    else if ( proxy->getPRole() == HELPER) {
        numberOfWins = matrix_size / window_length;
        /**Compares values in a given window by splitting the window in two halves and comparing each value to its counterpart at the same position in the other half.
        If size of the given windowVector is odd, there will be a residue, which is stored in residue. */
        while (cmpWindowVectorSize > 1) {
            uint32_t halfSize = cmpWindowVectorSize / 2;

            if (cmpWindowVectorSize & 1) {
                if (isResidueStored) {                            //second residue found --> add stored and current residue to one half each.
                    isResidueInBuffer = false;                    //after processing all windows, buffer is remembered to be empty; don't set isResidueStored directly otherwise residues in one loop iteration are treated differently
                    halfSize++;
                }
                else {                                          //no residue stored up to now:
                    isResidueInBuffer = true;                   //don't set isResidueStored directly, as then residues in one loop iteration would be treated differently
                }
            }

            auto vectorLength = static_cast<uint32_t>(floor(halfSize * numberOfWins));
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
        if (isResidueStored) {
            MSB(proxy, nullptr, numberOfWins);
            MUX(proxy, nullptr, nullptr, nullptr, numberOfWins);
        }
        return nullptr;
    }
    return nullptr;
}

/**
 * Method for private computation of the RELU function.
 * @param proxy
 * @param x - secret share of variable x for which to compute RELU(x)
 * @return
 */
uint64_t RELU(Party* proxy, uint64_t x){
    uint64_t K = (RING_N>>1); // N is the ring size - 1 = 2^64 -1

    if (proxy->getPRole() == P1 ||  proxy->getPRole() == P2) {
        uint64_t commonValues[3];
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
        uint64_t values[5];
        values[0] = proxy->getPRole() * f * (K+1) - z;                  // a_0
        values[1] = proxy->getPRole() * (1 - f) * (K+1) - z;            // a_1
        values[2] = (x + e[g]) * commonValues[0];                       // b
        values[3] = (e[h]) * commonValues[1];                           // c_0
        values[4] = (e[1 - h]) * commonValues[2];                       // c_1

        // proxy sends a,b and c to HELPER:
        unsigned char *ptr_out = proxy->getBuffer1();
        addVal2CharArray(values, &ptr_out, 5);
        Send(proxy->getSocketHelper(), proxy->getBuffer1(), 5 * 8);

        // receive fresh share from helper
        Receive(proxy->getSocketHelper(), proxy->getBuffer1(),
                6 * 8); // order is ab[0], ab[1], ac[0], ac[1], ac[2], ac[3]

        ptr_out = proxy->getBuffer1();
        uint64_t ab[2];
        uint64_t ac[4];
        convert2Array(&ptr_out, &ab[0], 2);
        convert2Array(&ptr_out, &ac[0], 4);

        uint64_t em;
        if (g == h){
            uint64_t r2_inverse = getModularInverse(commonValues[1]);
            em = ac[2*f] * r2_inverse;
        }else{
            uint64_t r3_inverse = getModularInverse(commonValues[2]);
            em = ac[2*f+1] * r3_inverse;
        }

        uint64_t r0_inverse = getModularInverse(commonValues[0]);
        uint64_t xm = ab[f] * r0_inverse - em;
        z = x - xm;
        return z;
    }
    else if (proxy->getPRole() == HELPER) {
        MOC(proxy, 0);

        Receive(proxy->getSocketP1(), proxy->getBuffer1(), 5 * 8);
        Receive(proxy->getSocketP2(), proxy->getBuffer2(), 5 * 8);
        unsigned char *ptr = proxy->getBuffer1();
        unsigned char *ptr2 = proxy->getBuffer2();

        uint64_t reconstructedVals[5];
        for (uint i = 0; i < 5; i++) {
            reconstructedVals[i] = (convert2Long(&ptr) + convert2Long(&ptr2));
            // 6 values expected per party --> first two (a values) are reconstructed and divided by K
            if (i < 2){
                reconstructedVals[i] /= (K+1);
            }
        }
        // reconstruct values a, b and c obtained from proxy 1 and 2
        uint64_t ab [2];
        uint64_t ac [4];
        ab[0] = reconstructedVals[0] * reconstructedVals[2]; // a_0 * b
        ab[1] = reconstructedVals[1] * reconstructedVals[2]; // a_1 * b

        ac[0] = reconstructedVals[0] * reconstructedVals[3]; // a_0 * c_0
        ac[1] = reconstructedVals[0] * reconstructedVals[4]; // a_0 * c_1
        ac[2] = reconstructedVals[1] * reconstructedVals[3]; // a_1 * c_0
        ac[3] = reconstructedVals[1] * reconstructedVals[4]; // a_1 * c_1
        //cout << "ab and ac are calculated" << endl;

        // create fresh shares
        uint64_t share1 [6];
        uint64_t share2 [6];

        uint64_t tmp = proxy->generateRandom();
        share1[0] = tmp;
        share2[0] = ab[0] - tmp;
        tmp = proxy->generateRandom();
        share1[1] = tmp;
        share2[1] = ab[1] - tmp;
        for (uint i = 0; i < 4; i++) {
            tmp = proxy->generateRandom();
            share1[i+2] = tmp;
            share2[i+2] = ac[i] - tmp;
        }

        // send shares to proxy 1 and 2
        unsigned char *ptr_back = proxy->getBuffer1();
        unsigned char *ptr_back2 = proxy->getBuffer2();
        addVal2CharArray(share1, &ptr_back, 6);
        addVal2CharArray(share2, &ptr_back2, 6);

        thread thr1 = thread(Send, proxy->getSocketP1(), proxy->getBuffer1(), 6 * 8);
        thread thr2 = thread(Send, proxy->getSocketP2(), proxy->getBuffer2(), 6 * 8);

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
 * @return vector of resulting values for each position in x, must be deleted if not needed anymore.
 */
uint64_t* RELU(Party* proxy, const uint64_t* x, uint64_t size){
    uint64_t K = (RING_N>>1); // N is the ring size - 1 = 2^64 -1
    if (proxy->getPRole() == P1 ||  proxy->getPRole() == P2) {
        uint64_t commonValues[3*size];
        for (unsigned long & commonValue : commonValues) {
            commonValue = proxy->generateCommonRandom() | 0x1; // common values must be odd
        }

        auto* e = new uint64_t [2*size]; // 2 shares per value to compute Relu for
        auto* f = new int[size];
        auto* g = new int[size];
        auto* h = new int[size];
        auto* t = new uint64_t [size];
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
        delete [] t;

        auto* z = new uint64_t [size];
        uint64_t values[5*size];
        for(uint64_t j = 0; j<size; j++){
            z[j] = x[j] - d[j];

            // compute parts of a, b and c:
            values[j*5] = proxy->getPRole() * f[j] * (K+1) - z[j];                  // a_0
            values[j*5 + 1] = proxy->getPRole() * (1 - f[j]) * (K+1) - z[j];        // a_1
            values[j*5 + 2] = (x[j] + e[j + g[j]*size]) * commonValues[j*3];        // b
            values[j*5 + 3] = (e[j + h[j]*size]) * commonValues[j*3+1];             // c_0
            values[j*5 + 4] = (e[j + (1-h[j])*size]) * commonValues[j*3+2];         // c_1
        }
        delete[] e;
        delete[] d;

        // proxy sends a,b and c to HELPER:
        unsigned char *ptr_out = proxy->getBuffer1();
        addVal2CharArray(values, &ptr_out, 5*size);
        Send(proxy->getSocketHelper(), proxy->getBuffer1(), 5 * size * 8);

        // receive fresh share from helper
        Receive(proxy->getSocketHelper(), proxy->getBuffer1(),
                6 * size * 8); // order is (ab[0], ab[1], then all (ac[0], ac[1], ac[2], ac[3])

        ptr_out = proxy->getBuffer1();
        uint64_t ab[2*size];
        uint64_t ac[4*size];
        convert2Array(&ptr_out, &ab[0], 2*size);
        convert2Array(&ptr_out, &ac[0], 4*size);

        for(uint64_t i = 0; i<size; i++) {
            uint64_t em;
            if (g[i] == h[i]) {
                uint64_t r2_inverse = getModularInverse(commonValues[i*3 + 1]);
                em = ac[i*4 + 2 * f[i]] * r2_inverse;
            } else {
                uint64_t r3_inverse = getModularInverse(commonValues[i*3 + 2]);
                em = ac[i*4 + 2 * f[i] + 1] * r3_inverse;
            }

            uint64_t r0_inverse = getModularInverse(commonValues[i*3]);
            uint64_t xm = ab[i*2 + f[i]] * r0_inverse - em;
            z[i] = x[i] - xm;
        }
        delete [] f;
        delete [] g;
        delete [] h;
        return z;
    }
    else if (proxy->getPRole() == HELPER) {
        MOC(proxy, nullptr, size);

        Receive(proxy->getSocketP1(), proxy->getBuffer1(), 5 * size * 8);
        Receive(proxy->getSocketP2(), proxy->getBuffer2(), 5 * size * 8);
        unsigned char *ptr = proxy->getBuffer1();
        unsigned char *ptr2 = proxy->getBuffer2();

        uint64_t reconstructedVals[5*size];
        for (uint i = 0; i < 5*size; i++) {
            reconstructedVals[i] = (convert2Long(&ptr) + convert2Long(&ptr2));
            if ((i % 5) < 2){
                reconstructedVals[i] /= (K+1);
            }
        }
        // reconstruct values a, b and c obtained from proxy 1 and 2
        uint64_t ab [2*size];
        uint64_t ac [4*size];
        for (uint i = 0; i < size; i++) {
            ab[i*2] = reconstructedVals[i*5] * reconstructedVals[i*5 + 2]; // a_0 * b
            ab[i*2 + 1] = reconstructedVals[i*5 + 1] * reconstructedVals[i*5 + 2]; // a_1 * b

            ac[i*4] = reconstructedVals[i*5] * reconstructedVals[i*5 + 3]; // a_0 * c_0
            ac[i*4 + 1] = reconstructedVals[i*5] * reconstructedVals[i*5 + 4]; // a_0 * c_1
            ac[i*4 + 2] = reconstructedVals[i*5 + 1] * reconstructedVals[i*5 + 3]; // a_1 * c_0
            ac[i*4 + 3] = reconstructedVals[i*5 + 1] * reconstructedVals[i*5 + 4]; // a_1 * c_1
        }
        // create fresh shares
        uint64_t share1 [6*size];
        uint64_t share2 [6*size];
        uint64_t* tmp = convert2uint64(random_1D_data(proxy, 6*size), 6*size);
        for (uint i = 0; i < 2*size; i++) {
            share1[i] = tmp[i];
            share2[i] = ab[i] - tmp[i];

            share1[i+2*size] = tmp[i+2*size];
            share2[i+2*size] = ac[i] - tmp[i+2*size];

            share1[i+4*size] = tmp[i+4*size];
            share2[i+4*size] = ac[i + 2*size] - tmp[i+4*size];
        }

        // send shares to proxy 1 and 2
        unsigned char *ptr_back = proxy->getBuffer1();
        unsigned char *ptr_back2 = proxy->getBuffer2();
        addVal2CharArray(share1, &ptr_back, 6*size);
        addVal2CharArray(share2, &ptr_back2, 6*size);

        thread thr1 = thread(Send, proxy->getSocketP1(), proxy->getBuffer1(), 6 * size * 8);
        thread thr2 = thread(Send, proxy->getSocketP2(), proxy->getBuffer2(), 6 * size * 8);

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

        z = proxy->getPRole() - convert2Long(&ptr);
        if (f) { // if f is 1  we get the next long value in the buffer.
            z = proxy->getPRole() - convert2Long(&ptr);
            cout << convert2double(REC(proxy, z)) << endl;
        }
        else{
            cout << convert2double(REC(proxy, convert2Long(&ptr))) << endl;
        }
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
 * @return vector containing the resulting DRELU(x), must be deleted if not needed anymore.
 */
uint64_t* DRELU(Party* proxy, uint64_t* x, uint32_t size){
    uint64_t K = (RING_N>>1); // N is the ring size - 1 = 2^64 -1
    // K is 2^63 - 1
    uint32_t exchangingBit = 2 * size;
    if (proxy->getPRole() == P1 || proxy->getPRole() == P2) {
        // init
        auto* f = new uint64_t [size];
        auto* t = new uint64_t [size];
        for (uint32_t i = 0; i < size; i++){
            f[i] = proxy->generateCommonRandom() & 0x1;
            t[i] = x[i] & K; // get first L-1 bit of the share
        }
        K += 1;
        uint64_t* d = MOC(proxy, t, size);
        delete [] t;
        uint64_t* z = SUB(x, d, size);
        delete[] d;

        // compute parts of a:
        uint64_t values[exchangingBit];
        for (uint32_t i = 0; i < size; i++){
            values[i * 2] = proxy->getPRole() * f[i] * K - z[i];         // a_0 of the i-th variable
            values[i * 2 + 1] = proxy->getPRole() * (1 - f[i]) * K - z[i];   // a_1
        }

        // proxy sends a to HELPER:
        unsigned char *ptr = proxy->getBuffer1();
        addVal2CharArray(values, &ptr, exchangingBit);
        Send(proxy->getSocketHelper(), proxy->getBuffer1(), exchangingBit * 8);

        // receive fresh share from helper
        Receive(proxy->getSocketHelper(), proxy->getBuffer1(), exchangingBit * 8); // order is share of a[0], a[1]
        ptr = proxy->getBuffer1();

        for (uint64_t c = 0; c < size; c++){
            z[c] = proxy->getPRole() - convert2Long(&ptr);
            if (f[c]) {  // if f is 1  we get a1.
                z[c] = proxy->getPRole() - convert2Long(&ptr);
            }
            else{// increase ptr so that in next iteration not a1 is read:
                (*ptr)+=7;
            }
        }
        delete[] f;
        return z;
    }
    else if (proxy->getPRole() == HELPER) {
        K += 1;
        MOC(proxy, nullptr, size);

        Receive(proxy->getSocketP1(), proxy->getBuffer1(), exchangingBit * 8);
        Receive(proxy->getSocketP2(), proxy->getBuffer2(), exchangingBit * 8);
        unsigned char *ptr1 = proxy->getBuffer1();
        unsigned char *ptr2 = proxy->getBuffer2();

        uint64_t reconstructedVals[exchangingBit];
        for (uint64_t v = 0; v < exchangingBit; v++) {
            reconstructedVals[v] = (convert2Long(&ptr1) + convert2Long(&ptr2)) / K;
        }
        //reassign buffer because ptr1 and ptr2 were incremented
        ptr1 = proxy->getBuffer1();
        ptr2 = proxy->getBuffer2();

        uint64_t* tmp = convert2uint64(random_1D_data(proxy, exchangingBit), exchangingBit);  // values for P1
        addVal2CharArray(tmp, &ptr1, exchangingBit);
        uint64_t * share = SUB(reconstructedVals, tmp, exchangingBit);
        addVal2CharArray(share, &ptr2, exchangingBit);

        thread thr1 = thread(Send, proxy->getSocketP1(), proxy->getBuffer1(), exchangingBit * 8);
        thread thr2 = thread(Send, proxy->getSocketP2(), proxy->getBuffer2(), exchangingBit * 8);
        thr1.join();
        thr2.join();
        delete [] share;

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

// LAYER HELPER FUNCTIONS:

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
 *         with conv_width = (width - k_dim + 1)/stride and conv_height = (height - k_dim + 1)/stride
 *         return vector must be deleted if not needed anymore.
 */
uint64_t ***INC(uint64_t ***input, uint32_t channel, uint32_t height, uint32_t width, uint32_t k_dim,
                uint32_t stride) {
    uint32_t k_size = k_dim * k_dim;
    uint32_t last_row_start = height - k_dim + 1;
    uint32_t last_col_start = width - k_dim + 1;
    uint32_t conv_height = last_row_start/stride;
    uint32_t conv_width = last_col_start/stride;
    // stretch the input for vectorized MATVECMUL
    auto ***stretched_input = new uint64_t **[channel];
    for(uint32_t c = 0; c < channel; c++){
        stretched_input[c] = new uint64_t *[conv_width * conv_height];
        for(uint32_t row = 0; row < last_row_start; row += stride) {
            for(uint32_t col = 0; col < last_col_start; col += stride) {
                uint32_t index = (row * conv_height + col)/stride;
                stretched_input[c][index] = new uint64_t[k_size];
                // shift kernel over all cols in according row
                for (uint32_t l = 0; l < k_size; l++) {
                    stretched_input[c][index][l] = input[c][row + l / k_dim][col + l % k_dim];
                }
            }
        }
    }
    return stretched_input;
}

/**
 * Symmetric padding of the input matrix.
 * @param input matrix to be padded
 * @param rows number of rows in input
 * @param cols number of cols in input
 * @param padding_value value to be inserted for padding
 * @param padding_size number of padding_values to be inserted in each direction: top, bottom, right and left
 * @return the padded input matrix, return vector must be deleted if not needed anymore.
 */
uint64_t **PAD(uint64_t** input, uint32_t rows, uint64_t cols, uint64_t padding_value, uint32_t padding_size){
    uint32_t padded_row_length = 2*padding_size+cols;
    uint32_t padded_col_length = 2*padding_size+rows;
    auto** padded_input = new uint64_t *[padded_col_length];

    for (uint32_t i = 0; i<padded_col_length; i++){
        padded_input[i] = new uint64_t[padded_row_length];
        // init whole matrix with padding value
        memset(padded_input[i], padding_value, padded_row_length*sizeof (uint64_t));
        if(i >= padding_size && i < (padding_size + rows)){
            memcpy(&padded_input[i][padding_size], input[i - padding_size], cols * sizeof(uint64_t));        // copy values
        }
    }
    return padded_input;
}

/**
 * Flatten the values of a matrix concatenating its values so that a single vector is the result
 * @param images matrices to be flattened to one single vactor
 * @param i_height height of a single image in images
 * @param i_width width of a single image in images
 * @param i_number number of matrices in images
 * @return the flattened vector of length i_dim * i_dim * i_number
 * return vector must be deleted if not needed anymore.
 */
uint64_t * FLT(uint64_t*** images, uint32_t i_height, uint32_t i_width, uint32_t i_number){
    uint64_t i_size = i_height * i_width;
    auto * flattened = new uint64_t [i_size * i_number];
    for (uint32_t i = 0; i < i_number; i++){
        for (uint32_t el = 0; el < i_size; el ++) {
            flattened[el + i * i_size] = images[i][el / i_width][el % i_width];
        }
    }
    return flattened;
}


// LAYER FUNCTIONS:

/**
 * Implements the function of a convolutional layer (CL) using ReLU as activation function
 * and then Maxpool with a 2x2 filter if according parameter is set.
 *
 * @param proxy
 * @param input data on which convolution is performed using the provided kernels.
 * The input is supposed to be in 2D shape where each channel is appended in the second dimension.
 * So a input of shape WxHxC is supposed to be represented in shape WxH*C.
 * @param i_channel number of channels of the input; must be > 0
 * @param i_height size of input in one dimension, divided by i_channel; must be > 0
 * @param i_width size of input in the other dimension; must be > 0
 * @param kernel vector containing all kernel to be used for convolution.
 * The parameter kernel has shape output_channel x i_channel * k_dim * k_dim.
 * So each channel of a kernel is represented by a vector of length k_dim * k_dim.
 * The kernel's number of channels is defined by i_channel.
 * @param k_dim dimension of the symmetric kernel; must be > 0
 * @param output_channel number of kernels with length i_channel * k_dim * k_dim; must be > 0
 * For each kernel there will be one output channel in the result, where the value at the according location is
 * the sum of the single multiplications between input and kernel of same channel number.
 * @param stride step size of each kernel to shift per iteration; must be > 0
 * @param max_win_height defines window heigth of maxpool, which will be realized after Relu activation.
 * A max_win_height <= 0 indicates that no maxpool shall be performed after Relu activation.
 * @param max_win_width defines window width of maxpool, which will be realized after Relu activation.
 * A max_win_width <= 0 indicates that no maxpool shall be performed after Relu activation.
 * In order for maxpool to be performed it must hold: (max_win_width > 0) and (max_win_height > 0)
 * @param maxpool_window_dim defines the window size of maxpooling operation, which will be realized after Relu activation.
 * A maxpool_window_dim <= 0 inficates that no maxpooling shall be performed after Relu activation.
 * @param bias vector of length output_channel. For each kernel there is one bias value which is added to every value of the according output_channel.
 * @param last_conv defines if this convolutional layer will be the last one of the cnn.
 * If so, the output will be the flattened matrix (as required for fully connected layer), so that output[0][0] = flattened matrix.
 * @return Output of the input convoluted by the given kernels.
 *         Shape of output will be:  c x h x w if last_conv == false
 *         otherwise the shape of the output will be 1 x 1 x c*h*w with
 *         c = output_channel
 *         h = floor((i_height - k_dim + 1)/stride) if doMaxpooling is false; otherwise h = floor((i_height - k_dim + 1)/(2*stride))
 *         w = floor((i_weight - k_dim + 1)/stride) if doMaxpooling is false; otherwise w = floor((i_weight - k_dim + 1)/(2*stride))
 *         return vector must be deleted if not needed anymore.
 */
uint64_t*** CL(Party* proxy, uint64_t*** input, uint32_t i_channel, uint32_t i_height, uint32_t i_width, uint64_t*** kernel, uint32_t k_dim, uint32_t output_channel, uint32_t stride, uint32_t max_win_height, uint32_t max_win_width, uint64_t* bias, bool last_conv = false){
    uint32_t k_size = k_dim * k_dim;
    auto conv_width = static_cast<uint32_t>(floor((i_width - k_dim + 1) / stride));
    auto conv_height = static_cast<uint32_t>(floor((i_height - k_dim + 1) / stride));
    uint32_t conv_len = conv_width * conv_height;
    bool doMaxpooling = ((max_win_width > 0) && (max_win_height > 0));
    // stretch the input for vectorized MATVECMUL
    if (proxy->getPRole() == P1 || proxy->getPRole() == P2) {
        //print1DArray("bias in CL: ", convert2double(REC(proxy, bias, output_channel), output_channel), output_channel);
        uint64_t *** stretched_input = INC(input, i_channel, i_height, i_width, k_dim, stride);
        //print2DArray("INC ", convert2double(REC(proxy,stretched_input[0], conv_len, k_size), conv_len, k_size), conv_len, k_size);
        uint32_t out_width = conv_width;
        uint32_t out_height = conv_height;
        if (doMaxpooling) {
            out_height /= max_win_height;
            out_width /= max_win_width;
        }
        // convolution:
        auto ***conv_layer = new uint64_t **[output_channel];
        for (uint32_t k = 0; k < output_channel; k++) {
            // multiply stretched_input with kernel (for all channel)
            uint64_t **conv_result = MATVECMUL(proxy, stretched_input, kernel[k], i_channel, conv_len, k_size);
            //print2DArray("conv Result ", convert2double(REC(proxy,conv_result, i_channel, conv_len), i_channel, conv_len), i_channel, conv_len);
            uint64_t* summed_channel_conv;
            if(i_channel == 1){
                summed_channel_conv = conv_result[0];
            }
            else{
                summed_channel_conv = ADD(proxy, conv_result, i_channel, conv_len);
            }
            // ACTIVATION:
            uint64_t* conv_activated = RELU(proxy, summed_channel_conv, conv_len);
            delete[] conv_result;
            //print1DArray("ReLU result", convert2double(REC(proxy,conv_activated, conv_len), conv_len), conv_len);
            delete[] summed_channel_conv;
            if (doMaxpooling){
                // Maxpool:
                conv_activated = MAX(proxy, conv_activated, conv_height, conv_width, max_win_height, max_win_width);
            }
            if(last_conv){
                if(k == 0){
                    conv_layer = new uint64_t **[1];
                    conv_layer[0] = new uint64_t*[1];
                    conv_layer[0][0] = new uint64_t[output_channel*out_width*out_height];
                }
                for (uint32_t row = 0; row < out_height; row++) {
                    for (uint32_t col = 0; col < out_width; col++) {
                        conv_layer[0][0][k*out_width*out_height + row*out_width + col] = conv_activated[row * out_width + col] + bias[k];
                    }
                }
            }
            else{
                // bring result in matrix shape
                conv_layer[k] = new uint64_t *[out_height];
                for (uint32_t row = 0; row < out_height; row++) {
                    conv_layer[k][row] = new uint64_t[out_width];
                    for (uint32_t col = 0; col < out_width; col++) {
                        conv_layer[k][row][col] = conv_activated[row * out_width + col] + bias[k];
                    }
                }
            }
            delete[] conv_activated;
        }
        for (int i = 0; i < i_channel; ++i) {
            for (int row = 0; row < conv_len; ++row) {
                delete[] stretched_input[i][row];
            }
            delete[] stretched_input[i];
        }
        delete [] stretched_input;
        return conv_layer;
    }
    else if (proxy->getPRole() == HELPER){
        // convolution:
        for (uint32_t k = 0; k < output_channel; k++) {
            MATVECMUL(proxy, nullptr, nullptr, 0, i_channel * conv_len * k_size, 0);

            // ACTIVATION:
            // ReLU
            RELU(proxy, nullptr, conv_len);
            if(doMaxpooling){
                // Maxpool:
                MAX(proxy, nullptr, conv_height, conv_width, max_win_height, max_win_width);
            }
        }
        return nullptr;
    }
    return nullptr;
}

/**
 * Implements functionality of a fully connected layer (FCL) with ReLU as activation function.
 * @param proxy
 * @param input the input vector of length in_size to be fully connected to the output nodes.
 * @param in_size length of the input vector (when previous layer is a convolutional layer, use the Flattening method FLT before)
 * @param weights to be used must be of shape in_size x node_number as provided by Chameleon files
 * @param node_number number of output nodes of this layer
 * @param bias vector of length node_number. For each output node there is one bias value which is added.
 * @return the output layer that have been computed by using the dot product between input values and weights,
 * activated with ReLu and the according bias added. It will be of length node_number.
 * return vector must be deleted if not needed anymore.
 */
uint64_t* FCL(Party* proxy, uint64_t* input, uint32_t in_size, uint64_t** weights, uint32_t node_number, uint64_t* bias){
    if (proxy->getPRole() == P1 || proxy->getPRole() == P2){
        uint64_t *output = MATVECMUL(proxy, weights, input, node_number, in_size);

        uint64_t* relu = RELU(proxy, output, node_number);
        output = ADD(proxy, relu, bias, node_number);
        delete[] relu;
        return output;
    }
    else if (proxy->getPRole() == HELPER){
        MATVECMUL(proxy, nullptr, nullptr, node_number*in_size, 0);
        RELU(proxy, nullptr, node_number);
        return nullptr;
    }
    else{
        cout << "Error: unknown proxy role for fully connected layer." << endl;
        return nullptr;
    }
}

#endif //PPAUC_CNN_H

