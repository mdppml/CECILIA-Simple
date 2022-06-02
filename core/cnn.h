//
// Created by Mete Akgun on 17.01.22.
//

#ifndef PPAUC_CNN_H
#define PPAUC_CNN_H

#include "core.h"
#include "../utils/flib.h"
#include "bitset"

/**
 * Prints the given matrix in a way to show the values in a column and row based format while also different windows of
 * the matrix are seperated from each other visually. The matrix will be reconstructed and converted to double.
 * @param str1 - String to be printed out before the visualization of the matrix.
 * @param matrix - the matrix to be shown represented by a vector where each row of the matrix contains m_col / w_col values per window.
 * @param m_row - number of rows in matrix
 * @param m_col - number of columns in matrix
 * @param w_row - number of rows per window
 * @param w_col - number of columns per window
 */
void print1DMatrixByWindows(string const &str1, uint64_t* matrix, uint32_t m_row, uint32_t m_col, uint32_t w_row, uint32_t w_col) {
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

void print1DMatrixByWindows(string const &str1, double *matrix, uint32_t m_row, uint32_t m_col, uint32_t w_row,
                             uint32_t w_col) {
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
            cout << matrix[i*m_col + j] << " \t " ;
        }
        cout << "|" << endl;
    }
    for(uint8_t d = 0; d < m_col; d++){
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
uint64_t* SUB(uint64_t *a, uint64_t *b, uint32_t length){
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
    //cout << "entered RST..." << endl;
    uint32_t winSize = w_cols * w_rows;
    uint32_t numberOfWins = (m_cols * m_rows) / winSize;
    uint32_t winsPerRow = m_cols / w_cols;
    uint32_t w_count = 0;

    while(w_count < numberOfWins){
        auto windowStart = static_cast<uint32_t>(w_rows * m_cols * floor(w_count / winsPerRow) + (w_cols * (w_count & (winsPerRow-1))));
        // cout << w_count << " out of " << numberOfWins << " windows: "<< windowStart << endl;
        for(uint32_t i = 0; i<w_rows; i++){
            uint32_t m_index = windowStart + i*m_cols;
            for(uint32_t  j = 0; j<w_cols; j++){
                uint32_t w_index = w_count * winSize + (i*w_cols + j);
                //cout << "Matrix: " << m_index+j << " = (start: " << windowStart << " + wrow: " << i << " * cols: " << m_cols << " + wcol: " << j << ");  resorted: " << w_index << "(value = " << matrix[m_index+j] << endl;
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
            if (cmpVectorSize % 2 == 1) {                   //there is an residue remaining
                if (isResidueStored) {                            //second residue found --> add stored and current residue each to one half.
                    halfSize++;                                //each half of window increases size by 1 because of residues
                    *(firstHalf + halfSize) = residue;
                    *(secondHalf + halfSize) = *(firstHalf + cmpVectorSize - 1); //last element of cmpVector
                } else {                                       //no residue stored up to now:
                    isResidueStored = true;
                    residue = *(firstHalf + cmpVectorSize - 1); // store last element in residue
                }
            }
            //if cmpVectorSize is odd, store the last element as residue
            if (halfSize > 0) {          // maximums are not yet found
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
        auto *resorted = new uint64_t [matrix_size];                // RESORT matrix to have all values of a window subsequently
        // only rows calculated because always complete rows of windows are processed
        RST(mShare, m_cols, m_rows, window_size, window_size, resorted);
        numberOfWins = matrix_size / window_length;
        auto *maxElements = new uint64_t [numberOfWins];
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
                if (cmpWindowVectorSize & 1) {                   //there is an residue remaining
                    if (isResidueStored) {                            //second residue found --> add stored and current residue each to one half.
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
    return NULL;
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
        for (uint8_t i = 0; i < 4; i++) {
            commonValues[i] = proxy->generateCommonRandom() | 0x1; // common values must be odd
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
        values[2] = (x + e[g]) * commonValues[0];        // b_0
        values[3] = (x + e[1 - g]) * commonValues[1];    // b_1
        values[4] = (e[h]) * commonValues[2];        // c_0
        values[5] = (e[1 - h]) * commonValues[3];    // c_1

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
        uint64_t f = proxy->generateCommonRandom() & 0x1;

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


uint64_t DIV(Party* proxy, uint64_t a, uint64_t b) {
    uint64_t K = (RING_N>>1);
    cout << "a = " << bitset<L_BIT>(REC(proxy, a)) << "; b = " << bitset<L_BIT>(REC(proxy, b)) << endl;
    if (proxy->getPRole()  == P1 || proxy->getPRole()  == P2){
        uint64_t dividend = MSB(proxy, a);
        uint64_t result = 0;
        for (int16_t i = L_BIT-1; i >= 0; i--){
            //cout << "______________ i : " << i << "; dividend = " << bitset<L_BIT>(REC(proxy, dividend)) << endl;
            //cout << "result: " << bitset<L_BIT>(convert2double(REC(proxy, result))) << endl;

            uint64_t z = ADD(proxy, b, -dividend);
            uint64_t quotient = MSB(proxy, z); // MSB of z is intermediate result/quotient //TODO what if b == dividend? quotient should be 1 then.
//            cout << "z: " << bitset<L_BIT>(convert2double(REC(proxy, z))) << "; quotient: " << bitset<L_BIT>(REC(proxy, quotient)) << endl;

            // set the i-th bit of result to the value of v
            result |= quotient & 1;
           /* if(quotient > 0 ){
                result |= (1UL << i);
            }
            else{
                result &= ~(1UL << i);
            }*/

          //  cout << "dividend before: " << bitset<L_BIT>(convert2double(REC(proxy, dividend))) << endl;
            dividend -= MUL(proxy, quotient, b); // next dividend is remainder (= dividend - (quotient * b))
            dividend |= a & (1 << i);
        }
        //result += zeroShare;
       // cout << "result: " << REC(proxy, result) << endl;
        return result;
    }
    else if (proxy->getPRole() == HELPER) {
        MSB(proxy, 0);
        for (int16_t i = L_BIT-1; i >= 0; i--) {
            ADD(proxy, 0, 0);
            MSB(proxy, 0);
            MUL(proxy, 0, 0);
        }
        return 0;
    }
}



#endif //PPAUC_CNN_H

