//
// Created by Mete Akgun on 17.01.22.
//

#ifndef PPAUC_CNN_H
#define PPAUC_CNN_H


#include "core.h"
#include "../utils/flib.h"

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
    std::unique_ptr<uint64_t[]> subtractedValues(new uint64_t (length));
    if( subtractedValues.get() == nullptr ) {
        cout << "Insufficient memory" << endl;
        return nullptr;
    }
    for(uint32_t i = 0; i<length; i++){
        subtractedValues[i] = *(a+i) - *(b+i);
        //cout << "(a-b = c) " << to_string(a[i]) << to_string(b[i]) << " = " << to_string(subtractedValues[i]) << endl;
    }
    return subtractedValues.get();
}


/**
 * Resort (RST) the given matrix so that elements of one window are found as a sequence of w_rows*w_cols.
 * Elements of a window can be found in matrix starting from index
 * i up to i + w_cols and
 * i + m_cols * win_row up to i + m_cols * win_row for each row of window  win_row.
 * compute window_size vectors:
 * matrix vector     a, b, c, d,
 *                   e, f, g, h,
 *                   i, j, k, l,
 *                   m, n, o, p              shall become
 * 4 window_size vectors  a, b, e, f,
 *                        c, d, g, h,
 *                        i, j, m, n,
 *                        k, l, o, p
 * @param matrix the matrix to be resorted
 * @param m_cols number of columns in matrix
 * @param m_rows number of rows in matrix
 * @param w_cols number of columns in window
 * @param w_rows number of rows in window
 * @param resortedMatrix pointer to the resulting resorted Matrix, which will have length of m_cols * m_rows
 *
 * CAUTION: can not process matrices greater than a size of xx by now.
 */
void RST(const uint64_t* matrix, uint32_t m_cols, uint32_t m_rows, uint32_t w_cols, uint32_t w_rows, uint64_t* resortedMatrix){
    //cout << "entered RST..." << endl;
    uint32_t winSize = w_cols * w_rows;
    uint32_t numberOfWins = (m_cols * m_rows) / winSize;
    uint32_t winsPerRow = m_cols / w_cols;
    uint32_t w_count = 0;

    while(w_count < numberOfWins){
        auto  windowStart = static_cast<uint32_t>(w_rows * m_cols * floor(w_count / winsPerRow) + (w_cols * w_count));
        //cout << w_count << " out of " << numberOfWins << " windows: "<< windowStart << endl;
        for(uint32_t  j = 0; j<w_cols; j++){
            for(uint32_t i = 0; i<w_rows; i++){
                //start of a window     shift to right row of windows                 shift to right column    shift to right row
                uint32_t m_index = windowStart + (i*m_cols + j);
                uint32_t w_index = w_count * winSize + (i*w_cols + j);
                //cout << "Matrix: " << m_index << " = (" << windowStart << " + " << i*m_cols + j << ");  resorted: " << w_index << endl;
                resortedMatrix[w_index] = matrix[m_index];
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
        auto *firstHalf = new uint64_t[matrix_size / 2];
        uint64_t *secondHalf;
        uint64_t residue; //there is at most one residual element.

        while (cmpVectorSize > 1) {
            uint32_t halfSize = static_cast<uint8_t>(floor(cmpVectorSize / 2));

            firstHalf = maxElements;
            secondHalf = firstHalf + halfSize;
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
        delete[]firstHalf;
        cout << "max = " << convert2double(REC(proxy, maxElements[0])) << endl;
        return maxElements[0]; // should only contain one element at the end.
    }
    else if ( proxy->getPRole() == HELPER) {
        /**Compares values in a given window by splitting the window in two halves and comparing each value to its counterpart at the same position in the other half.
        If size of the given windowVector is odd, there will be a residue, which is stored in residue. */
        uint64_t maxElement;
        while (cmpVectorSize > 1) {
            uint32_t halfSize = static_cast<uint8_t>(floor(cmpVectorSize / 2));
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

    uint64_t mProcessing_start = 0;                     // variables needed if matrix_size increases MAX_MATRIX_SIZE
    uint32_t mProcessing_length = matrix_size;          // first compute maxpool for max number of elements possible:

    if(matrix_size > MAX_MATRIX_SIZE){
        mProcessing_length = MAX_MATRIX_SIZE - (MAX_MATRIX_SIZE % window_length); // ensure that full windows are processed
        cout << "matrix contains to many elements. " << mProcessing_length << " will be processed which are " << to_string(mProcessing_length/window_length) << " windows of length " << window_length << endl;
    }
    uint32_t numberOfWins;
    auto maxElements = new uint64_t[matrix_size / window_length];
    cout << "Individual parts are starting...: processingLength= " << mProcessing_length << ", windowLength = " << window_length << endl;
    if ( proxy->getPRole() == P1 ||  proxy->getPRole() == P2) {
        uint32_t storedBufferVals = 0; // number of elements stored in the buffer --> needed to check if there is still capacity
        unsigned char *buffer = proxy->getBuffer1();

        do { // loop through parts of the matrix if needed to be particioned.
            numberOfWins = static_cast<uint32_t >(floor(mProcessing_length / window_length));
            cout << "numberOfWins= " << to_string(numberOfWins) << endl;

            while (cmpWindowVectorSize > 1) {
                /**Compares values in a given window by splitting the window in two halves and comparing each value to
                 * its counterpart at the same position in the other half. If size of the given windowVector is odd,
                 * there will be a residue, which is stored in residue.
                 * All three vectors are initialized within the loop to free space as less and less space is required.*/
                 cout << "started iterating through window...: cmpWindowVectorSize= " << to_string(cmpWindowVectorSize) << endl;
                auto halfSize = static_cast<uint32_t>(floor(cmpWindowVectorSize / 2)); // elements per win in halves.

                uint64_t comparisons = (halfSize * numberOfWins);
                auto *firstHalf = new uint64_t[comparisons];
                auto *secondHalf = new uint64_t[comparisons];
                auto *residue = new uint64_t[numberOfWins]; //if there are residues, for each window there is one element --> size is number of windows.
                cout << "first, second half and residue are created. resort..." << endl;

                // RESORT matrix to have all values of a window subsequently
                auto *resorted = new uint64_t[mProcessing_length];
                //TODO recalculate mCols and mRows so that they match processing legnth.
                RST(mShare, window_size*numberOfWins, window_size*numberOfWins, window_size, window_size, resorted);
                cout << "matrix was resorted. vector size to be compared per win: " << cmpWindowVectorSize << "; halfsize= " << halfSize << endl;
                for (uint32_t i = 0; i < numberOfWins; i++) {
                    uint64_t *currWindowStart = mShare + mProcessing_start + i * cmpWindowVectorSize;
                    uint64_t *currWindowMiddle = currWindowStart + halfSize;

                    //if cmpWindowVectorSize is odd, store the last element as residue
                    if (cmpWindowVectorSize % 2 == 1) {                   //there is an residue remaining
                        if (isResidueStored) {                            //second residue found --> add stored and current residue each to one half.
                            isResidueInBuffer = false;                    //after processing all windows, buffer is remembered to be empty; dont set isResidueStored directly otherwise residues in one loop iteration are treated differently
                            currWindowMiddle += 1;                        //one more for residues
                            halfSize++;                                   //one half of window increases size by 1 because of residue
                            uint32_t posInHalfes = (i + 1) * halfSize - 1;

                            firstHalf[posInHalfes] = residue[i];
                            secondHalf[posInHalfes] = *(currWindowStart + cmpWindowVectorSize -
                                                            1); //last element of current window vector
                        } else {                           //no residue stored up to now:
                            cout << "store " << i << "th value in buffer" << endl;
                            isResidueInBuffer = true;                   //dont set isResidueStored directly, as then residues in one loop iteration would be treated differently
                            residue[i] = *(currWindowStart + cmpWindowVectorSize - 1);
                        }
                    }
                    uint32_t vHalfIndex = i * halfSize;         // index at which values of window i are to be stored.
                    for (uint32_t j = 0; j < halfSize; j++) {   // splitting the windows on firstHalf and secondHalf
                        //cout << to_string(j) << "currWinStart = " << currWindowStart + j << "of " <<  " was store? ";
                        firstHalf[vHalfIndex + j] = *(currWindowStart + j);
                        //cout << to_string(firstHalf[vHalfIndex + j]) << endl << "currWinMiddle = " << to_string(*(currWindowMiddle+ j)) << " was store? ";
                        secondHalf[vHalfIndex + j] = *(currWindowMiddle + j);
                        //cout << to_string(secondHalf[vHalfIndex + j]) << endl;
                    }/*
                    firstHalf[vHalfIndex] = *(currWindowStart);
                    secondHalf[vHalfIndex] = *(currWindowMiddle);
                    cout << "this is correct: " << *(currWindowStart + 3) << " = " << firstHalf[vHalfIndex+3] << endl;*/
                }
                uint32_t vectorLength = halfSize * numberOfWins;
                if (vectorLength > 0) {          // maximums are not yet found
                    //compare: a-b =c and then MSB(c) =d
                    cout << "start SUB... vectorLength= " << vectorLength << endl;
                    uint64_t *c = SUB(firstHalf, secondHalf, vectorLength);
                    cout << "finished SUB: " << &c << "; start MSB..." << endl;
                    uint64_t *d = MSB(proxy, c, vectorLength);
                    cout << "finished MSB: " << d << "; start MUX..." << endl;
                    //MUX:
                    mShare = MUX(proxy, firstHalf, secondHalf, d, vectorLength);
                    cout << "MUX = " << maxElements << endl;
                }
                //prepare next round:
                cmpWindowVectorSize = halfSize;
                isResidueStored = isResidueInBuffer;
            }
            mProcessing_start += mProcessing_length + 1; // next start where current part ended.
            mProcessing_length = matrix_size - mProcessing_start; // remaining length to process
            if (mProcessing_length > MAX_MATRIX_SIZE) {
                // remaining data is still too big
                mProcessing_length = MAX_MATRIX_SIZE - (MAX_MATRIX_SIZE % window_length);
            }
            if ((numberOfWins + storedBufferVals) <= BUFFER_SIZE){
                // there is still capacity to store more elements:
                for (uint32_t win = 0; win < numberOfWins; win++) {
                    uint64_t maxVal = *(maxElements + win);
                    addVal2CharArray(maxVal, &buffer);
                }
                storedBufferVals += numberOfWins;
            }
            else{
                // empty buffer TODO
            }
        } while (mProcessing_start + mProcessing_length <= matrix_size);
        // all windows were processed --> put data in maxElements and empty buffer
        for (uint32_t win = 0; win < storedBufferVals; win++) {
            maxElements[win] = buffer[win];
        }
        buffer = nullptr;
        return maxElements;
    }
    else if ( proxy->getPRole() == HELPER) {
        do {
            numberOfWins = static_cast<uint32_t >(floor(mProcessing_length / window_length));
            cout << "numberOfWins= " << to_string(numberOfWins) << endl;
            /**Compares values in a given window by splitting the window in two halves and comparing each value to its counterpart at the same position in the other half.
            If size of the given windowVector is odd, there will be a residue, which is stored in residue. */
            while (cmpWindowVectorSize > 1) {
                cout << "start a round of the while loop..." << endl;
                auto halfSize = static_cast<uint32_t>(floor(cmpWindowVectorSize / 2));

                // doesn not need to be in for loop (see P0/P1) as there, bool values are set to same in each iteration
                if (cmpWindowVectorSize % 2 == 1) {                   //there is an residue remaining
                    if (isResidueStored) {                            //second residue found --> add stored and current residue each to one half.
                        isResidueInBuffer = false;
                        halfSize++;                                   //one half of window increases size by 1 because of residue
                    }
                    else {
                        isResidueInBuffer = true;
                    }
                }

                uint32_t vectorLength = halfSize * numberOfWins;
                if (vectorLength > 0) {          // maximums are not yet found
                    //compare: a-b =c and then MSB(c) =d
                    SUB(0, 0, vectorLength);
                    cout <<"vectorLength = " << vectorLength << endl;
                    MSB(proxy, 0, vectorLength);
                    cout << "finished MSB; start MUX..." << endl;
                    //MUX:
                    MUX(proxy, nullptr, nullptr, nullptr, vectorLength);
                    cout << "finished MUX." << endl;
                }
                //prepare next round:
                cmpWindowVectorSize = halfSize;
                isResidueStored = isResidueInBuffer;
            }
            mProcessing_start += mProcessing_length + 1; // next start where current part ended.
            mProcessing_length = matrix_size - mProcessing_start; // remaining length to process
            if (mProcessing_length > MAX_MATRIX_SIZE) {
                // remaining data is still too big
                mProcessing_length = MAX_MATRIX_SIZE;
            }
        } while (mProcessing_start + mProcessing_length <= matrix_size);
        return NULL;
    }
}



/**
 * Method for private computation of the RELU function.
 * @param proxy
 * @param x - variable x for which to compute RELU(x)
 * @return
 */
uint64_t RELU(Party* proxy, uint64_t x){
    uint64_t K = (N>>1); // N is the ring size - 1 = 2^64 -1

    role p_role = proxy->getPRole();
    if (proxy->getPRole() == P1 ||  proxy->getPRole() == P2) {
        auto *commonValues = new uint64_t [3];
        for (uint8_t i=0; i<3; i++){
            commonValues[i] = proxy->generateCommonRandom();
        }
        uint8_t *buffer = proxy->getBuffer1();
        int socket_helper = proxy->getSocketHelper();

        // create even random shares
        uint64_t e_0 = proxy->generateRandom() & EVEN_MASK; //ensure e_0 is even
        uint64_t e_1 = proxy->generateRandom() & EVEN_MASK; //ensure e_1 is even
        uint64_t e[] = {e_0, e_1};

        // init
        bool f, g, h;
        f = proxy->generateCommonRandom() & 0x1;
        g = proxy->generateCommonRandom() & 0x1;
        h = proxy->generateCommonRandom() & 0x1;
        //cout << "bools: " << f << " " << g << " " << h << " " << endl;

        // make the shares more random by adding i to e_i for S_i:
        e[f] += p_role;

        uint64_t t = x & K; // get first L-1 bit of the share

        uint64_t d = MOC(proxy, t);
        uint64_t z = x - d;

        // compute parts of a, b and c:
        auto *values = new uint64_t [6];
        values[0] = p_role * f * K - z;                  // a_0
        values[1] = p_role * (1 - f) * K - z;            // a_1
        values[2] = (x + e[g]) * commonValues[0];        // b_0
        values[3] = (x + e[1 - g]) * commonValues[1];    // b_1
        values[4] = (x + e[h]) * commonValues[2];        // c_0
        values[5] = (x + e[1 - h]) * commonValues[3];    // c_1

        // proxy sends a,b and c to HELPER:
        unsigned char *ptr_out = &buffer[0];
        addVal2CharArray(values, &ptr_out, 6);
        Send(socket_helper, buffer, 6 * 8);

        // receive fresh share from helper
        Receive(socket_helper, buffer, 8 * 8); // order is ab[0], ab[1], ab[2], ab[3], ac[0], ac[1], ac[2], ac[3]

        uint64_t r2_inverse = getModularInverse(commonValues[2]);
        // ac[2f]
        uint64_t em = buffer[2*f+4] * r2_inverse;

        uint64_t r0_inverse = getModularInverse(commonValues[0]);
        // ab[2f]
        uint64_t xm = buffer[2*f] * r0_inverse - em;
        z = x - xm;
        return z;
    }
    else if (p_role == HELPER) {
        K += 1; // increase K by 1 K is 2^63
        int socket_p1 = proxy->getSocketP1();
        int socket_p2 = proxy->getSocketP2();

        uint8_t * buffer = proxy->getBuffer1();
        uint8_t * buffer2 = proxy->getBuffer2();

        MOC(proxy, NULL);

        Receive(socket_p1, buffer, 6 * 8);
        Receive(socket_p2, buffer2, 6 * 8);
        unsigned char *ptr = &buffer[0];
        unsigned char *ptr2 = &buffer2[0];

        auto *reconstructedVals = new uint64_t [6];
        for (uint i = 0; i < 6; i++) {
            reconstructedVals[i] = (*(ptr + i) + *(ptr2 + i));
            // 6 values expected per party --> first two (a values) are reconstructed and divided by K
            if (i < 2){
                reconstructedVals[i] /= K;
            }
        }
        // reconstruct values a, b and c obtained from proxy 1 and 2
        auto *ab = new uint64_t [4];
        auto *ac = new uint64_t [4];
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
        auto *share1 = new uint64_t [8];
        auto *share2 = new uint64_t [8];
        for (uint i = 0; i < 4; i++) {
            uint64_t tmp = proxy->generateRandom();
            share1[i] = tmp;
            share2[i] = ab[i] - tmp;

            tmp = proxy->generateRandom();
            share1[i+4] = tmp;
            share2[i+4] = ac[i] - tmp;
        }

        // send shares to proxy 1 and 2
        unsigned char *ptr_back = &buffer[0];
        unsigned char *ptr_back2 = &buffer2[0];
        addVal2CharArray(share1, &ptr_back, 8);
        addVal2CharArray(share2, &ptr_back2, 8);

        thread thr1 = thread(Send, socket_p1, buffer, 8 * 8);
        thread thr2 = thread(Send, socket_p2, buffer2, 8 * 8);

        thr1.join();
        thr2.join();
        return 0;
    }
}


/**
 * Method for private computation of the derivative of the RELU function.
 * @param proxy
 * @param x - variable x for which to compute RELU'(x), the derivative of the RELU function.
 * @return
 */
uint64_t DRELU(Party* proxy, uint64_t x){
    uint64_t K = (N>>1); // N is the ring size - 1 = 2^64 -1
    // K is 2^63 - 1
    uint8_t exchangingBit = 2;
    if (proxy->getPRole() == P1 ||  proxy->getPRole() == P2) {
        //int socket_helper = proxy->getSocketHelper(); // you dont need to assign proxy->getSocketHelper() to a variable.  proxy->getSocketHelper() is a getter for the socket of the helper. delete this line

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

        // proxy sends a,b and c to HELPER:
        unsigned char *ptr = proxy->getBuffer1();
        addVal2CharArray(values, &ptr, exchangingBit);
        Send(proxy->getSocketHelper(), proxy->getBuffer1(), exchangingBit * 8);

        // receive fresh share from helper
        Receive(proxy->getSocketHelper(), proxy->getBuffer1(), exchangingBit * 8); // order is share of a[0], a[1]
        ptr = proxy->getBuffer1();
        z = proxy->getPRole() -  convert2Long(&ptr);
        if (f)  // if f is 1  we get the next long value in the buffer.
            z = proxy->getPRole() -  convert2Long(&ptr);

        //z = p_role - buffer[commonRandom]; we can not do this. buffer[0] is a single byte, we need 8 bytes for the output. delete this line the correct one is above

        return z;
    }
    else if (proxy->getPRole() == HELPER) {
        K += 1;
        MOC(proxy, NULL);

        Receive(proxy->getSocketP1(), proxy->getBuffer1(), exchangingBit * 8);
        Receive(proxy->getSocketP2(), proxy->getBuffer2(), exchangingBit * 8);
        unsigned char *ptr1 = proxy->getBuffer1();
        unsigned char *ptr2 = proxy->getBuffer2();

        uint64_t reconstructedVals[exchangingBit];
        reconstructedVals[0] = (convert2Long(&ptr1) + convert2Long(&ptr2)) / K;
        reconstructedVals[1] = (convert2Long(&ptr1) + convert2Long(&ptr2)) / K;

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
    else{
        cout << "No proxy role recognized." << endl;
        return -1;
    }
}


#endif //PPAUC_CNN_H

