//
// Created by Mete Akgun on 17.01.22.
//

#ifndef PPAUC_CNN_H
#define PPAUC_CNN_H

#include "core.h"

// please write your functions here. (RELU, MAXPOOL, vs...)


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
    }
    return subtractedValues;
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
 * @return the resorted matrix
 */
uint64_t * RST(uint64_t* matrix, uint32_t m_cols, uint32_t m_rows, uint32_t w_cols, uint32_t w_rows, uint64_t* resortedMatrix){
    uint32_t winSize = w_cols * w_rows;
    uint32_t numberOfWins = (m_cols * m_rows) / winSize;
    uint32_t winsPerRow = m_cols / w_cols;
    uint32_t w_count = 0;

    while(w_count < numberOfWins){
        //cout << w_count << " out of " << numberOfWins << " windows."<< endl;
        uint32_t  windowStart = w_rows * m_cols * floor(w_count / winsPerRow);
        uint32_t winsRow = (w_count % winsPerRow) * w_cols;
        for(uint32_t  j = 0; j<w_cols; j++){
            for(uint32_t i = 0; i<w_rows; i++){
                //start of a window     shift to right row of windows                 shift to right column    shift to right row
                uint32_t m_index = windowStart +  winsRow + (i*m_cols + j);
                uint32_t w_index = w_count * winSize + (i*w_cols + j);
                //cout << "Matrix: " << m_index << " = (" << windowStart << " + " << winsRow << " + " << i*m_cols + j << ");  resorted: " << w_index << endl;
                resortedMatrix[w_index] = matrix[m_index];
            }
        }
        w_count++;
    }
    return resortedMatrix;
}

/**
 * Selects the maximum element from the given matrix in secret shared form.
 * @param mShare - secret share of the matrix from which the maximal element shall be computed.
 * @param matrix_size - size of mShare.
 * @return The maximum element which was found in mShare.
 */
uint64_t MAX(Party* proxy, uint64_t *mShare, uint32_t matrix_size){
    /**Compares values by splitting the matrix in two halves and comparing each value to its counterpart at the same position in the other half.
    If size of the given matrix is odd, there will be a residue, which is stored in residue. */

    uint32_t cmpVectorSize =
            matrix_size / 2; //size of resulting vector after cmp, MUX and its divided by 2 is size of each half.
    bool isResidueStored = false;

    if ( proxy->getPRole() == P1 ||  proxy->getPRole() == P2) {
        uint64_t *maxElements = mShare;
        uint64_t *firstHalf = new uint64_t[matrix_size / 2];
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
                MSB(proxy, 0, halfSize);

                //MUX:
                MUX(proxy, 0, 0, 0, halfSize);
            }
            //prepare next round:
            cmpVectorSize = halfSize;
        }
        return NULL;
    }
}


/**
 * Maxpooling of a matrix of size matrix_size and non overlapping windows of size window_size. Address of the first
 * address in the matrix is mShare. as the matrix is processed as secret shares.
 * @param mShare - share of the matrix
 * @param matrix_size - number of elements in the matrix mShare
 * @param window_size - size of the window in both dimensions (symmetric windows possible only)
 * @return the maximum element per window, therefore floor(matrix_size / (window_size*window_size)) elements in a vector
 */
uint64_t* MAX(Party* proxy, uint64_t *mShare, uint32_t matrix_size, uint32_t window_size){
    //compare within one window for each window:
    uint32_t window_length = window_size * window_size;
    uint32_t cmpWindowVectorSize = window_length; //size of resulting vector after cmp, MUX and its divided by 2 is size of each half.
    bool isResidueStored = false;
    bool isResidueInBuffer = false;

    uint64_t mProcessing_start = 0;                     // variables needed if matrix_size increases MAX_MATRIX_SIZE
    uint32_t mProcessing_length = matrix_size;          // first compute maxpool for max number of elements possible:

    /*if(matrix_size > MAX_MATRIX_SIZE){
        mProcessing_length = MAX_MATRIX_SIZE - (MAX_MATRIX_SIZE % window_length); // ensure that full windows are processed
        cout << "matrix contains to many elements. " << mProcessing_length << " will be prcessed which are " << to_string(mProcessing_length/window_length) << " windows of length " << window_length << endl;
    }*/
    uint32_t numberOfWins;
    uint64_t *maxElements = new uint64_t[matrix_size/window_length];

    if ( proxy->getPRole() == P1 ||  proxy->getPRole() == P2) {
        uint32_t storedBufferVals = 0; // number of elements stored in the buffer --> needed to check if there is still capacity
        unsigned char *buffer = proxy->getBuffer1();

        do {
            numberOfWins = static_cast<uint32_t >(floor(mProcessing_length / window_length));
            cout << "number of windows: " << numberOfWins << endl;
            /**Compares values in a given window by splitting the window in two halves and comparing each value to its counterpart at the same position in the other half.
            If size of the given windowVector is odd, there will be a residue, which is stored in residue. */
            uint64_t *firstHalf = new uint64_t[numberOfWins * window_size];
            uint64_t *secondHalf = new uint64_t[numberOfWins * window_size];
            uint64_t *residue = new uint64_t[numberOfWins]; //if there are residues, for each window there is one element --> size is number of windows.
            while (cmpWindowVectorSize > 1) {
                uint32_t halfSize = static_cast<uint8_t>(floor(cmpWindowVectorSize / 2));
                cout << "vector size to be compared per win: " << cmpWindowVectorSize << "; halfsize= " << halfSize << endl;
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

                            *(firstHalf + posInHalfes) = residue[i];
                            *(secondHalf + posInHalfes) = *(currWindowStart + cmpWindowVectorSize -
                                                            1); //last element of current window vector
                        } else {                           //no residue stored up to now:
                            cout << "store " << i << "th value in buffer" << endl;
                            isResidueInBuffer = true;                   //dont set isResidueStored directly, as then residues in one loop iteration would be treated differently
                            *(residue + i) = *(currWindowStart + cmpWindowVectorSize - 1);
                        }
                    }
                    uint32_t vHalfIndex = i * halfSize;
                    for (uint32_t j = 0; j < halfSize; j++) {
                        *(firstHalf + vHalfIndex + j) = *(currWindowStart + j);
                        *(secondHalf + vHalfIndex + j) = *(currWindowMiddle + j);
                    }
                }
                uint32_t vectorLength = halfSize * numberOfWins;
                if (vectorLength > 0) {          // maximums are not yet found
                    //compare: a-b =c and then MSB(c) =d
                    cout << "start SUB... vectorLength= " << vectorLength << endl;
                    uint64_t *c = SUB(firstHalf, secondHalf, vectorLength);
                    cout << "finished SUB: " << c << "; start MSB..." << endl;
                    uint64_t *d = MSB(proxy, c, vectorLength);
                    cout << "finished MSB: " << d << "; start MUX..." << endl;
                    //MUX:
                    maxElements = MUX(proxy, firstHalf, secondHalf, d, vectorLength);
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
            cout << "number of windows: " << numberOfWins << endl;
            /**Compares values in a given window by splitting the window in two halves and comparing each value to its counterpart at the same position in the other half.
            If size of the given windowVector is odd, there will be a residue, which is stored in residue. */
            while (cmpWindowVectorSize > 1) {
                //cout << "start a round of the while loop..." << endl;
                uint32_t halfSize = static_cast<uint8_t>(floor(cmpWindowVectorSize / 2));

                // doesn not need to be in for loop like for parties as bool values are set to same in every iteration
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
                    //SUB(0, 0, vectorLength);
                    cout <<"vectorLength = " << vectorLength << endl;
                    MSB(proxy, 0, vectorLength);
                    //cout << "finished MSB; start MUX..." << endl;
                    //MUX:
                    MUX(proxy, 0, 0, 0, vectorLength);
                    //cout << "finished MUX." << endl;
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

#endif //PPAUC_CNN_H

