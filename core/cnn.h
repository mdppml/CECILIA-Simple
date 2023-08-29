//
// Created by Mete Akgun on 17.01.22.
//

#ifndef PPAUC_CNN_H
#define PPAUC_CNN_H

#include "core.h"
#include "../utils/flib.h"
#include "bitset"
#include "string.h"

using namespace std;

// Node for binary-tree data structure
struct node {
    int index;
    uint64_t multiplier;
    struct node* parent;
    struct node* left_child;
    struct node* right_child;
};

struct node *NewNode(Party *proxy, int index) {
    struct node *n = (struct node *)malloc(sizeof(struct node));

    n->index = index;
    n->multiplier = proxy->GetPRole();
    n->parent = NULL;
    n->left_child = NULL;
    n->right_child = NULL;
    return n;
}

/**
 * Vectorized subtraction: subtract elementwise b from a (a-b)
 * @param a from values of this vector the values of b will be subtracted
 * @param b the vector whose values are going to be subtracted from a
 * @param length of a and b (must have same length)
 * @return vector d where each element is the result of the according elements a-b.
 */
uint64_t* Subtract(const uint64_t *const a, const uint64_t *const b, uint32_t length){
    auto *subtractedValues = new uint64_t [length];
    for(uint32_t i = 0; i<length; i++){
        subtractedValues[i] = *(a+i) - *(b+i);
    }
    return subtractedValues;
}


uint64_t** Transpose(const uint64_t *const *const matrix, uint32_t rows, uint32_t cols){
    uint64_t ** transposed = new uint64_t *[cols];
    for (int c = 0; c < cols; ++c) {
        transposed[c] = new uint64_t [rows];
        for (int r = 0; r < rows; ++r) {
            transposed[c][r] = matrix[r][c];
        }
    }
    return transposed;
}

uint64_t*** Transpose(const uint64_t *const *const *const matrix, uint32_t n_matrices, uint32_t rows, uint32_t cols){
    uint64_t *** transposed = new uint64_t **[n_matrices];
    for (int m = 0; m < n_matrices; ++m) {
        transposed[m] = new uint64_t *[cols];
        for (int c = 0; c < cols; ++c) {
            transposed[m][c] = new uint64_t [rows];
            for (int r = 0; r < rows; ++r) {
                transposed[m][c][r] = matrix[m][r][c];
            }
        }
    }
    return transposed;
}

/**
 * Resort (Resort) the given matrix so that elements of one window are found as a sequence of w_rows*w_cols.
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
void Resort(const uint64_t *const matrix, uint32_t m_cols, uint32_t m_rows, uint32_t w_cols, uint32_t w_rows, uint64_t* resortedMatrix){
    uint32_t winSize = w_cols * w_rows;
    uint32_t numberOfWins = (m_cols * m_rows) / winSize;
    uint32_t winsPerRow = m_cols / w_cols;
    uint32_t w_count = 0;

    while(w_count < numberOfWins){
        auto windowStart = static_cast<uint32_t>(w_rows * m_cols * floor(w_count / winsPerRow) + (w_cols * (w_count % winsPerRow)));
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
uint64_t ArgMax(Party *const proxy, const uint64_t *const mShare, uint32_t matrix_size){
    /** MAIN IDEA:
     * As the Max is performed, a second matrix of same length is created containing the indices only.
     * However values are selected from matrix, is also done for the indices-matrix.
     */
    uint32_t cmpVectorSize = matrix_size; //size of resulting vector after cmp, Multiplex.
    bool isResidueStored = false;
    bool isSecondHalfFilled = false;
    if (proxy->GetPRole() == proxy1 || proxy->GetPRole() == proxy2) {
        uint64_t *maxElements = new uint64_t [2*matrix_size];
        // generate matrix with indices
        uint64_t * tmp = ConvertToUint64(Random1dData(proxy, matrix_size), matrix_size);
        for (int i = 0; i < matrix_size; ++i) {
            maxElements[i] = mShare[i];
            if(proxy->GetPRole() == proxy1){
                maxElements[i + matrix_size] = ConvertToUint64(i) - tmp[i];
            }
            else{
                maxElements[i + matrix_size] = tmp[i];
            }
        }
        auto maxHalfSizes = matrix_size;
        uint64_t elements1 [maxHalfSizes];
        uint64_t elements2 [maxHalfSizes];
        uint64_t residue[2]; //there is at most one residual element and its index

        while (cmpVectorSize > 1) {
            auto halfSize = static_cast<uint32_t>(floor(cmpVectorSize / 2));
            if (cmpVectorSize & 0x1) {                          //there is a residue remaining
                if (isResidueStored) {                          //second residue found --> add stored and current residue to one half each.
                    isResidueStored = false;
                    halfSize++;                                 //each half of window increases size by 1 because of residues
                    memcpy(elements2, maxElements+halfSize, (halfSize-1)*8);
                    elements2[halfSize-1] = residue[0];
                    memcpy(elements2 + halfSize, maxElements + cmpVectorSize + halfSize, (halfSize-1)*8);
                    elements2[cmpVectorSize] = residue[1];
                    isSecondHalfFilled = true;
                } else {                                        //no residue stored up to now:
                    isResidueStored = true;
                    residue[0] = maxElements[cmpVectorSize-1]; // store last element in residue
                    residue[1] = maxElements[2*cmpVectorSize-1];
                }
            }
            memcpy(elements1, maxElements, halfSize * 8);
            memcpy(elements1 + halfSize, maxElements + cmpVectorSize, halfSize*8);
            if(!isSecondHalfFilled){
                memcpy(elements2, maxElements+halfSize, halfSize * 8);
                memcpy(elements2 + halfSize, maxElements + cmpVectorSize + halfSize, halfSize*8);
            }
            if (halfSize > 0) {                                 // maximums not yet found
                //compare: a-b =c and then MostSignificantBit(c) =d
                uint64_t *c = Subtract(elements1, elements2, halfSize);
                uint64_t *d = MostSignificantBit(proxy, c, halfSize);
                uint64_t selection_vector[2*halfSize];
                memcpy(selection_vector, d, halfSize * 8);
                memcpy(selection_vector+halfSize, d, halfSize * 8);
                //Multiplex:
                delete[] maxElements;
                maxElements = Multiplex(proxy, elements1, elements2, selection_vector, 2 * halfSize);
                delete[] c;
                delete[] d;
            }
            //prepare next round:
            cmpVectorSize = halfSize;
            isSecondHalfFilled = false;
        }
        if (isResidueStored) {
            uint64_t c = maxElements[0] - residue[0];
            uint64_t d = MostSignificantBit(proxy, c);
            maxElements[1] = Multiplex(proxy, maxElements[1], residue[1], d);
        }

        uint64_t argmax = maxElements[1];                          // should only contain one element at the end; indices are after values
        delete [] maxElements;
        return argmax;
    }
    else if (proxy->GetPRole() == helper) {
        /**Compares values in a given window by splitting the window in two halves and comparing each value to its counterpart at the same position in the other half.
        If size of the given windowVector is odd, there will be a residue, which is stored in residue. */
        while (cmpVectorSize > 1) {
            auto halfSize = static_cast<uint32_t>(floor(cmpVectorSize / 2));
            if (cmpVectorSize % 2 == 1) {                   //there is a residue remaining
                if (isResidueStored) {                            //second residue found --> add stored and current residue each to one half.
                    isResidueStored = false;
                    halfSize++;                                //each half of window increases size by 1 because of residues
                } else {                                       //no residue stored up to now:
                    isResidueStored = true;
                }
            }
            //if cmpVectorSize is odd, store the last element as residue
            if (halfSize > 0) {          // maximums are not yet found
                //compare: a-b =c and then MostSignificantBit(c) =d
                MostSignificantBit(proxy, nullptr, halfSize);

                //Multiplex:
                Multiplex(proxy, nullptr, nullptr, nullptr, 2 * halfSize);
                //Multiplex(proxy, nullptr, nullptr, nullptr, halfSize);
            }
            //prepare next round:
            cmpVectorSize = halfSize;
        }
        if (isResidueStored) {
            MostSignificantBit(proxy, 0);
            Multiplex(proxy, 0, 0, 0);
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
uint64_t Max(Party *const proxy, const uint64_t *const mShare, uint32_t matrix_size){
    /** MAIN IDEA:
     * Compare values by splitting the matrix in two halves and
     * comparing each value to its counterpart at the same position in the other half.
     * If size of the given matrix is odd, there will be a residue, which is stored in residue.
     */
    uint32_t cmpVectorSize = matrix_size; //size of resulting vector after cmp, Multiplex.
    bool isResidueStored = false;
    bool isSecondHalfFilled = false;
    if (proxy->GetPRole() == proxy1 || proxy->GetPRole() == proxy2) {
        uint64_t *maxElements = new uint64_t[matrix_size];
        std::memcpy(maxElements, mShare, matrix_size*8);
        auto maxHalfSizes = static_cast<uint32_t>(ceil(matrix_size / 2)); // ceil because residue might be added.
        uint64_t firstHalf [maxHalfSizes];
        uint64_t secondHalf [maxHalfSizes];
        uint64_t residue; //there is at most one residual element.

        while (cmpVectorSize > 1) {
            auto halfSize = static_cast<uint32_t>(floor(cmpVectorSize / 2));

            if (cmpVectorSize & 0x1) {                          //there is a residue remaining
                if (isResidueStored) {                          //second residue found --> add stored and current residue to one half each.
                    halfSize++;                                 //each half of window increases size by 1 because of residues
                    memcpy(secondHalf, maxElements+halfSize, (halfSize-1)*8);
                    secondHalf[halfSize-1] = residue;
                    isSecondHalfFilled = true;
                } else {                                        //no residue stored up to now:
                    isResidueStored = true;
                    residue = maxElements[cmpVectorSize-1]; // store last element in residue
                }
            }
            memcpy(firstHalf, maxElements, halfSize*8);
            if(!isSecondHalfFilled){
                memcpy(secondHalf, maxElements+halfSize, halfSize*8);
            }

            if (halfSize > 0) {                                 // maximums not yet found
                //compare: a-b =c and then MostSignificantBit(c) =d
                uint64_t *c = Subtract(firstHalf, secondHalf, halfSize);
                uint64_t *d = MostSignificantBit(proxy, c, halfSize);

                //Multiplex:
                delete[] maxElements;
                maxElements = Multiplex(proxy, firstHalf, secondHalf, d, halfSize);
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
    else if (proxy->GetPRole() == helper) {
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
                //compare: a-b =c and then MostSignificantBit(c) =d
                MostSignificantBit(proxy, nullptr, halfSize);

                //Multiplex:
                Multiplex(proxy, nullptr, nullptr, nullptr, halfSize);
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
uint64_t* Max(Party *const proxy, const uint64_t *const matrix, uint32_t m_rows, uint32_t m_cols, uint32_t win_rows, uint32_t win_cols, bool backprop = false){
    uint32_t matrix_size = m_rows * m_cols;
    uint32_t window_length = win_cols * win_rows;
    uint32_t cmpWindowVectorSize = window_length; //size of resulting vector after cmp, Multiplex, and it's divided by 2 is size of each half.
    bool isResidueStored = false;
    bool isResidueInBuffer = false;

    uint32_t numberOfWins;
    uint64_t **dmax;
    struct node** leaves;
    struct node** current_nodes;
    struct node** previous_nodes;
    struct node** residue_nodes;
    if (proxy->GetPRole() == proxy1 || proxy->GetPRole() == proxy2) {
        auto *resorted = new uint64_t [matrix_size];                // RESORT matrix to have all values of a window subsequently
        Resort(matrix, m_cols, m_rows, win_cols, win_rows, resorted);
        numberOfWins = matrix_size / window_length;

        // for backpropagation - initialize the required variables
        if(backprop) {
            cout << "check 1" << endl;
            // initialize the matrix to hold multiplier of each index
            cout << "depth: " << ceil(log2(window_length)) << endl;
            dmax = new uint64_t*[(int) ceil(log2(window_length))];
            for(int i = 0; i < log2(window_length); i++) {
                dmax[i] = new uint64_t[matrix_size];
                for(int j = 0; j < matrix_size; j++) {
                    dmax[i][j] = ConvertToUint64(proxy->GetPRole());
                    if(i == 0 && j == 0) { cout << "Role:" << proxy->GetPRole() << endl;}
                }
            }
            int tmp_depth = ceil(log2(window_length));
            Print2dArray("dmax - init",
                         ConvertToDouble(Reconstruct(proxy, dmax, tmp_depth, matrix_size), tmp_depth, window_length),
                         tmp_depth, matrix_size);

            Print1dArray("dmax - init", ConvertToDouble(
                    Reconstruct(proxy, Straighten2dArray(dmax, tmp_depth, matrix_size),
                                matrix_size * tmp_depth), matrix_size * tmp_depth), matrix_size * tmp_depth);
            cout << "check 2" << endl;

            // initialize index map
            leaves = new struct node*[matrix_size];
            for(int w = 0; w < numberOfWins; w++) {
                for(int i = 0; i < window_length; i++) {
                    leaves[w * window_length + i] = NewNode(proxy, i);
                }
            }
            cout << "check 3" << endl;

            previous_nodes = leaves;
            residue_nodes = new struct node*[numberOfWins];
            cout << "check 4" << endl;
        }

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

            if(backprop) {
                cout << "check 5 - halfSize: " << halfSize << endl;
                int step_size = halfSize;
                int coreHalfSize = halfSize;
                if (cmpWindowVectorSize & 1) {                        //there is a residue remaining
                    cout << "check 6" << endl;
                    if(isResidueStored) {                        //there is a residue remaining
                        cout << "check 6.1" << endl;
                        current_nodes = new struct node*[comparisons + numberOfWins];
                        step_size++;
                        for(int i = 0; i < numberOfWins; i++) {
                            current_nodes[step_size * (i + 1) - 1] = NewNode(proxy, step_size - 1); // how to index these nodes? individually or thinking them as a whole
                            current_nodes[step_size * (i + 1) - 1]->left_child = previous_nodes[cmpWindowVectorSize * (i + 1) - 1];
                            previous_nodes[cmpWindowVectorSize * (i + 1) - 1]->parent = current_nodes[step_size * (i + 1) - 1];
                            current_nodes[step_size * (i + 1) - 1]->right_child = residue_nodes[i];
                            residue_nodes[i]->parent = current_nodes[step_size * (i + 1) - 1];
                        }
                    }
                    else { //no residue stored up to now:
                        cout << "check 6.2" << endl;
                        current_nodes = new struct node*[comparisons];
                        for(int i = 0; i < numberOfWins; i++) {
                            residue_nodes[i] = previous_nodes[cmpWindowVectorSize * (i + 1) - 1];
                        }
                    }
                }
                else {
                    cout << "check 6.3" << endl;
                    current_nodes = new struct node*[comparisons];
                }
                cout << "check 7" << endl;
                for(int j = 0; j < numberOfWins; j++) {
                     for(int i = 0; i < coreHalfSize; i++) {
                         current_nodes[i + step_size * j] = NewNode(proxy, i);
                         current_nodes[i + step_size * j]->left_child = previous_nodes[i + cmpWindowVectorSize * j];
                         previous_nodes[i + cmpWindowVectorSize * j]->parent = current_nodes[i + step_size * j];
                         current_nodes[i + step_size * j]->right_child = previous_nodes[i + coreHalfSize + cmpWindowVectorSize * j];
                         previous_nodes[i + coreHalfSize + cmpWindowVectorSize * j]->parent = current_nodes[i + step_size * j];
                    }
                }
                cout << "check 8" << endl;
            }

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
                        memcpy(secondHalf+i * halfSize, currWindowMiddle, (halfSize-1)*8);
                        secondHalf[i * halfSize + halfSize - 1] = residue[i];
                    }
                    else {                                          //no residue stored up to now:
                        isResidueInBuffer = true;                   //don't set isResidueStored directly, as then residues in one loop iteration would be treated differently
                        residue[i] = *(currWindowStart + cmpWindowVectorSize - 1);
                    }
                }
                memcpy(firstHalf + i * halfSize, currWindowStart, halfSize*8);
                if(isResidueInBuffer or !isResidueStored){
                    //otherwise, memcpy was already called for secondHalf
                    memcpy(secondHalf + i * halfSize, currWindowMiddle, halfSize*8);
                }
            }
            if (comparisons > 0) {          // maximums are not yet found
                //compare: a-b =c and then MostSignificantBit(c) =d
//                uint64_t *c = Subtract(firstHalf, secondHalf, comparisons);
//                uint64_t *d = MostSignificantBit(proxy, c, comparisons);
                uint64_t *d = Compare(proxy, firstHalf, secondHalf, comparisons);

                //Multiplex: returns for each position i: firstHalf[i] if d[i] = 0; secondHalf[i] if d[i] = 1
//                resorted = Multiplex(proxy, firstHalf, secondHalf, d, comparisons);
                resorted = Multiplex(proxy, secondHalf, firstHalf, d, comparisons);

                if(backprop) {
                    cout << "check 9" << endl;
                    for(int i = 0; i < comparisons; i++) {
                        cout << "check 9." << i << ".1" << endl;
                        current_nodes[i]->left_child->multiplier = d[i];
                        cout << "check 9." << i << ".2" << endl;
                        current_nodes[i]->right_child->multiplier = ConvertToUint64(proxy->GetPRole()) - d[i];
                    }
                    cout << "check 10" << endl;
                    previous_nodes = current_nodes; // update the previous_nodes to current_nodes
                }

//                delete[] c;
                delete[] d;
            }
            //prepare next round:
            cmpWindowVectorSize = halfSize;
            isResidueStored = isResidueInBuffer;
        }
        if (isResidueStored) {
//            uint64_t *c = Subtract(resorted, residue, numberOfWins);
//            uint64_t *d = MostSignificantBit(proxy, c, numberOfWins);
            uint64_t *d = Compare(proxy, resorted, residue, numberOfWins);
//            resorted = Multiplex(proxy, resorted, residue, d, numberOfWins);
            resorted = Multiplex(proxy, residue, resorted, d, numberOfWins);

            if(backprop) {
                cout << "check 11" << endl;
                for(int i = 0; i < comparisons; i++) {
                    cout << "check 11." << i << endl;
                    struct node* tmp_node = NewNode(proxy, 0);
                    tmp_node->left_child = current_nodes[i];
                    current_nodes[i]->parent = tmp_node;
                    tmp_node->left_child->multiplier = d[i];
                    tmp_node->right_child = residue_nodes[i];
                    residue_nodes[i]->parent = tmp_node;
                    tmp_node->right_child->multiplier = ConvertToUint64(proxy->GetPRole()) - d[i];
                }
                cout << "check 12" << endl;
            }

//            delete[] c;
            delete[] d;
        }
        for (uint32_t m = 0; m < numberOfWins; m++){
            maxElements[m] = resorted[m];
        }

        // derivative of maxpool
        if(backprop) {
            cout << "check 13" << endl;
            Print1dArray("dmax - initial", ConvertToDouble(
                                 Reconstruct(proxy, Straighten2dArray(dmax, ceil(log2(window_length)), matrix_size),
                                             matrix_size * ceil(log2(window_length))), matrix_size * ceil(log2(window_length))),
                         matrix_size * ceil(log2(window_length)));
            for(int i = 0; i < matrix_size; i++) {
                cout << "check 13." << i << endl;
                int ind = 0;
                struct node* crnt = leaves[i];
                while(ind < ceil(log2(window_length)) && crnt->parent != NULL) {
                    cout << "check 13." << i << "." << ind << endl;
                    dmax[ind][i] = crnt->multiplier;
                    cout << ind << " -- " << crnt->index << " -> " << crnt->parent->index << endl;
                    cout << "check 13." << i << "." << ind << " - 2" << endl;
                    ind++;
                    crnt = crnt->parent;
                    cout << "check 13." << i << "." << ind << " - 3" << endl;
                }
                cout << "check 13." << i << " - end" << endl;
            }
            cout << "check 14" << endl;

            int n_depth = ceil(log2(window_length));
            int size;
            bool flag = false;
            uint64_t *res = new uint64_t[matrix_size];
            uint64_t *v1;
            uint64_t *v2;
            uint64_t *cur_dmax = Straighten2dArray(dmax, n_depth, matrix_size);
            cout << "check 15" << endl;
            for(int i = 0; i < ceil(log2(ceil(log2(window_length)))); i++) {
//                double *rec_cur_dmax = ConvertToDouble(Reconstruct(proxy, cur_dmax, matrix_size * n_depth), matrix_size * n_depth);
//                Print1dMatrixByWindows("computed max values (test): ", rec_cur_dmax, n_depth, matrix_size, 3, 3);
                Print1dArray("dmax - Step " + to_string(i), ConvertToDouble(
                                     Reconstruct(proxy, cur_dmax, matrix_size * n_depth), matrix_size * n_depth),
                             matrix_size * n_depth);
                size = n_depth / 2;
                cout << "check 16." << i << endl;
                if(n_depth % 2 == 1) {
                    if(flag) {
                        cout << "check 16." << i << " - 1" << endl;
                        size++;
                        v1 = new uint64_t[size * matrix_size];
                        v2 = new uint64_t[size * matrix_size];
                        memcpy(&v1[matrix_size * (size - 1)], &cur_dmax[(n_depth - 1) * matrix_size], matrix_size);
                        memcpy(&v2[matrix_size * (size - 1)], &res, matrix_size);
                        flag = false;
                    }
                    else {
                        cout << "check 16." << i << " - 2" << endl;
                        v1 = new uint64_t[size * matrix_size];
                        v2 = new uint64_t[size * matrix_size];
                        memcpy(res, &cur_dmax[(n_depth - 1) * matrix_size], matrix_size * 8);
                        flag = true;
                    }
                }
                else {
                    cout << "check 16." << i << " - 3" << endl;
                    cout << "Size: " << size << " - Matrix size: " << matrix_size << endl;
                    v1 = new uint64_t[size * matrix_size];
                    v2 = new uint64_t[size * matrix_size];
                }
                cout << "check 17." << i << endl;
                cout << "n_depth: " << n_depth << endl;
                Print1dArray("dmax - Step " + to_string(i) + " before copy", ConvertToDouble(
                                     Reconstruct(proxy, cur_dmax, matrix_size * n_depth), matrix_size * n_depth),
                             matrix_size * n_depth);
                for(int j = 0; j < n_depth / 2; j++) {
                    cout << "check 17." << i << " - 2." << j << endl;
                    memcpy(&v1[j * matrix_size], &cur_dmax[j * matrix_size], matrix_size * 8);
                    cout << " --- start of v2: " << (j + (n_depth / 2)) * matrix_size << endl;
                    memcpy(&v2[j * matrix_size], &cur_dmax[(j + (n_depth / 2)) * matrix_size], matrix_size * 8);
                }
                cout << "check 18." << i << " - mul size: " << size * matrix_size << endl;
                Print1dArray("v1", ConvertToDouble(Reconstruct(proxy, v1, size * matrix_size), size * matrix_size),
                             size * matrix_size);
                Print1dArray("v2", ConvertToDouble(Reconstruct(proxy, v2, size * matrix_size), size * matrix_size),
                             size * matrix_size);
                delete [] cur_dmax;
                cur_dmax = Multiply(proxy, v1, v2, size * matrix_size);
                n_depth = size;

                delete [] v1;
                delete [] v2;
            }
            if(flag) {
                cout << "check 19" << endl;
                cur_dmax = Multiply(proxy, cur_dmax, res, matrix_size); // return this if backprop is true
            }
            Print1dMatrixByWindows("cur_dmax: ",
                                   ConvertToDouble(Reconstruct(proxy, cur_dmax, matrix_size), matrix_size),
                                   numberOfWins, window_length, 1, window_length);
//            Print1dArray("dmax", ConvertToDouble(Reconstruct(proxy, cur_dmax, matrix_size), matrix_size), matrix_size);
        }

        delete [] resorted;
        return maxElements;
    }
    else if (proxy->GetPRole() == helper) {
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
                //compare: a-b =c and then MostSignificantBit(c) =d
//                MostSignificantBit(proxy, nullptr, vectorLength);
                Compare(proxy, nullptr, nullptr, vectorLength);
                //Multiplex:
                Multiplex(proxy, nullptr, nullptr, nullptr, vectorLength);
            }
            //prepare next round:
            cmpWindowVectorSize = halfSize;
            isResidueStored = isResidueInBuffer;
        }
        if (isResidueStored) {
//            MostSignificantBit(proxy, nullptr, numberOfWins);
            Compare(proxy, nullptr, nullptr, numberOfWins);
            Multiplex(proxy, nullptr, nullptr, nullptr, numberOfWins);
        }

        if(backprop) {
            int n_depth = ceil(log2(window_length));
            int size;
            bool flag = false;
            for (int i = 0; i < ceil(log2(ceil(log2(window_length)))); i++) {
                size = n_depth / 2;
                if (n_depth % 2 == 1) {
                    if (flag) {
                        size++;
                        flag = false;
                    } else {
                        flag = true;
                    }
                }
                Multiply(proxy, 0, 0, size * matrix_size);
                n_depth = size;
            }
            if (flag) {
                Multiply(proxy, 0, 0, matrix_size); // return this if backprop is true
            }
        }

        return nullptr;
    }
    return nullptr;
}

/**
 * Method for private computation of the ReLU function.
 * @param proxy
 * @param x - secret share of variable x for which to compute ReLU(x)
 * @return
 */
uint64_t ReLU(Party *const proxy, uint64_t x){
    uint64_t K = (RING_SIZE >> 1); // N is the ring size - 1 = 2^64 -1

    if (proxy->GetPRole() == proxy1 || proxy->GetPRole() == proxy2) {
        uint64_t commonValues[2];
        for (uint64_t & commonValue : commonValues) {
            commonValue = proxy->GenerateCommonRandom() | 0x1; // common values must be odd
        }
        // create even random shares
        uint64_t e[] = {proxy->GenerateRandom() & EVEN_MASK, proxy->GenerateRandom() & EVEN_MASK};

        // init
        int f = proxy->GenerateCommonRandom() & 0x1;
        int g = proxy->GenerateCommonRandom() & 0x1;
        int h = proxy->GenerateCommonRandom() & 0x1;

        // make the shares more random by adding i to e_i for S_i:
        e[f] += proxy->GetPRole();

        uint64_t t = x & K; // get first L-1 bit of the share

        uint64_t d = ModularConversion(proxy, t);
        uint64_t z = x - d;

        // compute parts of a, b and c:
        uint64_t values[5];
        values[0] = proxy->GetPRole() * f * (K + 1) - z;                  // a_0
        values[1] = proxy->GetPRole() * (1 - f) * (K + 1) - z;            // a_1
        values[2] = (x + e[g]);                                         // b
        values[3] = (e[h]) * commonValues[0];                           // c_0
        values[4] = (e[1 - h]) * commonValues[1];                       // c_1

        // proxy sends a,b and c to helper:
        unsigned char *ptr_out = proxy->GetBuffer1();
        AddValueToCharArray(values, &ptr_out, 5);
        Send(proxy->GetSocketHelper(), proxy->GetBuffer1(), 5 * 8);

        // receive fresh share from helper
        Receive(proxy->GetSocketHelper(), proxy->GetBuffer1(),
                6 * 8); // order is ab[0], ab[1], ac[0], ac[1], ac[2], ac[3]

        ptr_out = proxy->GetBuffer1();
        uint64_t ab[2];
        uint64_t ac[4];
        ConvertToArray(&ptr_out, &ab[0], 2);
        ConvertToArray(&ptr_out, &ac[0], 4);

        uint64_t em;
        if (g == h){
            uint64_t r0_inverse = GetModularInverse(commonValues[0]);
            em = ac[2*f] * r0_inverse;
        }else{
            uint64_t r1_inverse = GetModularInverse(commonValues[1]);
            em = ac[2*f+1] * r1_inverse;
        }

        //uint64_t r0_inverse = GetModularInverse(commonValues[0]);
        uint64_t xm = ab[f] - em;
        z = x - xm;
        return z;
    }
    else if (proxy->GetPRole() == helper) {
        ModularConversion(proxy, 0);

        Receive(proxy->GetSocketP1(), proxy->GetBuffer1(), 5 * 8);
        Receive(proxy->GetSocketP2(), proxy->GetBuffer2(), 5 * 8);
        unsigned char *ptr = proxy->GetBuffer1();
        unsigned char *ptr2 = proxy->GetBuffer2();

        uint64_t reconstructedVals[5];
        for (uint i = 0; i < 5; i++) {
            reconstructedVals[i] = (ConvertToLong(&ptr) + ConvertToLong(&ptr2));
            // 6 values expected per party --> first two (a values) are reconstructed and divided by K
            if (i < 2){
                reconstructedVals[i] /= (K+1);
            }
        }
        // BenchmarkReconstruct values a, b and c obtained from proxy 1 and 2
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

        uint64_t tmp = proxy->GenerateRandom();
        share1[0] = tmp;
        share2[0] = ab[0] - tmp;
        tmp = proxy->GenerateRandom();
        share1[1] = tmp;
        share2[1] = ab[1] - tmp;
        for (uint i = 0; i < 4; i++) {
            tmp = proxy->GenerateRandom();
            share1[i+2] = tmp;
            share2[i+2] = ac[i] - tmp;
        }

        // send shares to proxy 1 and 2
        unsigned char *ptr_back = proxy->GetBuffer1();
        unsigned char *ptr_back2 = proxy->GetBuffer2();
        AddValueToCharArray(share1, &ptr_back, 6);
        AddValueToCharArray(share2, &ptr_back2, 6);

        thread thr1 = thread(Send, proxy->GetSocketP1(), proxy->GetBuffer1(), 6 * 8);
        thread thr2 = thread(Send, proxy->GetSocketP2(), proxy->GetBuffer2(), 6 * 8);

        thr1.join();
        thr2.join();
        return 0;
    }
    return -1;
}

/**
 * Method for private computation of the ReLU function.
 * @param proxy
 * @param x - vector of variables for which to compute ReLU
 * @param size - size of vector x
 * @return vector of resulting values for each position in x, must be deleted if not needed anymore.
 */
uint64_t* ReLU(Party *const proxy, const uint64_t *const x, uint64_t size){
    uint64_t K = (RING_SIZE >> 1); // N is the ring size - 1 = 2^64 -1
    if (proxy->GetPRole() == proxy1 || proxy->GetPRole() == proxy2) {
        uint64_t commonValues[2*size];
        for (uint64_t & commonValue : commonValues) {
            commonValue = proxy->GenerateCommonRandom() | 0x1; // common values must be odd
        }
        auto* e = new uint64_t [2*size]; // 2 shares per value to compute Relu for
        auto* f = new int[size];
        auto* g = new int[size];
        auto* h = new int[size];
        auto* t = new uint64_t [size];
        for (uint64_t i = 0; i < size; i++){
            // create even random shares: store first all 'e_0', then all 'e_1'
            e[i] = proxy->GenerateRandom() & EVEN_MASK;
            e[i+size] = proxy->GenerateRandom() & EVEN_MASK;

            // init f, g, and h per value in x
            f[i] = proxy->GenerateCommonRandom() & 0x1;
            g[i] = proxy->GenerateCommonRandom() & 0x1;
            h[i] = proxy->GenerateCommonRandom() & 0x1;

            // make the shares more random by adding i to e_i for S_i:
            // is 0*size or 1*size --> add pRole to either according e_0 or e_1
            e[i + f[i]*size] += proxy->GetPRole();

            t[i] = x[i] & K; // get first L-1 bit of the share
        }
        uint64_t* d = ModularConversion(proxy, t, size);
        delete [] t;

        auto* z = new uint64_t [size];
        uint64_t values[5*size];
        for(uint64_t j = 0; j<size; j++){
            z[j] = x[j] - d[j];

            // compute parts of a, b and c:
            values[j*5] = proxy->GetPRole() * f[j] * (K + 1) - z[j];                  // a_0
            values[j*5 + 1] = proxy->GetPRole() * (1 - f[j]) * (K + 1) - z[j];        // a_1
            values[j*5 + 2] = x[j] + e[j + g[j]*size];                              // b
            values[j*5 + 3] = (e[j + h[j]*size]) * commonValues[j*2];               // c_0
            values[j*5 + 4] = (e[j + (1-h[j])*size]) * commonValues[j*2+1];         // c_1
        }
        delete[] e;
        delete[] d;

        // proxy sends a,b and c to helper:
        unsigned char *ptr_out = proxy->GetBuffer1();
        AddValueToCharArray(values, &ptr_out, 5 * size);
        Send(proxy->GetSocketHelper(), proxy->GetBuffer1(), 5 * size * 8);

        // receive fresh share from helper
        Receive(proxy->GetSocketHelper(), proxy->GetBuffer1(),
                6 * size * 8); // order is (ab[0], ab[1], then all (ac[0], ac[1], ac[2], ac[3])

        ptr_out = proxy->GetBuffer1();
        uint64_t ab[2*size];
        uint64_t ac[4*size];
        ConvertToArray(&ptr_out, &ab[0], 2 * size);
        ConvertToArray(&ptr_out, &ac[0], 4 * size);

        for(uint64_t i = 0; i<size; i++) {
            uint64_t em;
            if (g[i] == h[i]) {
                uint64_t r0_inverse = GetModularInverse(commonValues[i * 2]);
                em = ac[i*4 + 2 * f[i]] * r0_inverse;
            } else {
                uint64_t r1_inverse = GetModularInverse(commonValues[i * 2 + 1]);
                em = ac[i*4 + 2 * f[i] + 1] * r1_inverse;
            }

            //uint64_t r0_inverse = GetModularInverse(commonValues[i*3]);
            uint64_t xm = ab[i*2 + f[i]] - em;
            z[i] = x[i] - xm;
        }
        delete [] f;
        delete [] g;
        delete [] h;
        return z;
    }
    else if (proxy->GetPRole() == helper) {
        ModularConversion(proxy, nullptr, size);

        Receive(proxy->GetSocketP1(), proxy->GetBuffer1(), 5 * size * 8);
        Receive(proxy->GetSocketP2(), proxy->GetBuffer2(), 5 * size * 8);
        unsigned char *ptr = proxy->GetBuffer1();
        unsigned char *ptr2 = proxy->GetBuffer2();

        uint64_t reconstructedVals[5*size];
        for (uint i = 0; i < 5*size; i++) {
            reconstructedVals[i] = (ConvertToLong(&ptr) + ConvertToLong(&ptr2));
            if ((i % 5) < 2){
                reconstructedVals[i] /= (K+1);
            }
        }
        // BenchmarkReconstruct values a, b and c obtained from proxy 1 and 2
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
        uint64_t* tmp = ConvertToUint64(Random1dData(proxy, 6 * size), 6 * size);
        for (uint i = 0; i < 2*size; i++) {
            share1[i] = tmp[i];
            share2[i] = ab[i] - tmp[i];

            share1[i+2*size] = tmp[i+2*size];
            share2[i+2*size] = ac[i] - tmp[i+2*size];

            share1[i+4*size] = tmp[i+4*size];
            share2[i+4*size] = ac[i + 2*size] - tmp[i+4*size];
        }

        // send shares to proxy 1 and 2
        unsigned char *ptr_back = proxy->GetBuffer1();
        unsigned char *ptr_back2 = proxy->GetBuffer2();
        AddValueToCharArray(share1, &ptr_back, 6 * size);
        AddValueToCharArray(share2, &ptr_back2, 6 * size);

        thread thr1 = thread(Send, proxy->GetSocketP1(), proxy->GetBuffer1(), 6 * size * 8);
        thread thr2 = thread(Send, proxy->GetSocketP2(), proxy->GetBuffer2(), 6 * size * 8);

        thr1.join();
        thr2.join();

        return nullptr;
    }
    return nullptr;
}


/**
 * Method for private computation of the derivative of the ReLU function.
 * @param proxy
 * @param x - variable x for which to compute ReLU'(x), the derivative of the ReLU function.
 * @return
 */
uint64_t DerivativeReLU(Party *const proxy, uint64_t x){
    uint64_t K = (RING_SIZE >> 1); // N is the ring size - 1 = 2^64 -1
    // K is 2^63 - 1
    uint8_t exchangingBit = 2;
    if (proxy->GetPRole() == proxy1 || proxy->GetPRole() == proxy2) {
        // init
        uint64_t f = proxy->GenerateCommonRandom() & 0x1; // generate common random bit

        uint64_t t = x & K; // get first L-1 bit of the share
        K += 1; // increase K by 1 K is 2^63
        uint64_t d = ModularConversion(proxy, t);
        uint64_t z = x - d;

        // compute parts of a:
        uint64_t values[exchangingBit];
        values[0] = proxy->GetPRole() * f * K - z;       // a_0
        values[1] = proxy->GetPRole() * (1 - f) * K - z; // a_1

        // proxy sends a to helper:
        unsigned char *ptr = proxy->GetBuffer1();
        AddValueToCharArray(values, &ptr, exchangingBit);
        Send(proxy->GetSocketHelper(), proxy->GetBuffer1(), exchangingBit * 8);

        // receive fresh share from helper
        Receive(proxy->GetSocketHelper(), proxy->GetBuffer1(), exchangingBit * 8); // order is share of a[0], a[1]
        ptr = proxy->GetBuffer1();

        z = proxy->GetPRole() - ConvertToLong(&ptr);
        if (f) { // if f is 1  we get the next long value in the buffer.
            z = proxy->GetPRole() - ConvertToLong(&ptr);
        }
        return z;
    }
    else if (proxy->GetPRole() == helper) {
        K += 1;
        ModularConversion(proxy, 0);

        Receive(proxy->GetSocketP1(), proxy->GetBuffer1(), exchangingBit * 8);
        Receive(proxy->GetSocketP2(), proxy->GetBuffer2(), exchangingBit * 8);
        unsigned char *ptr1 = proxy->GetBuffer1();
        unsigned char *ptr2 = proxy->GetBuffer2();

        uint64_t reconstructedVals[exchangingBit];
        reconstructedVals[0] = (ConvertToLong(&ptr1) + ConvertToLong(&ptr2)) / K;
        reconstructedVals[1] = (ConvertToLong(&ptr1) + ConvertToLong(&ptr2)) / K;

        //reassign buffer because ptr1 and ptr2 were incremented by convert2Long calls.
        ptr1 = proxy->GetBuffer1();
        ptr2 = proxy->GetBuffer2();

        uint64_t tmp = proxy->GenerateRandom();
        AddValueToCharArray(tmp, &ptr1);
        AddValueToCharArray(reconstructedVals[0] - tmp, &ptr2);
        tmp = proxy->GenerateRandom();
        AddValueToCharArray(tmp, &ptr1);
        AddValueToCharArray(reconstructedVals[1] - tmp, &ptr2);

        thread thr1 = thread(Send, proxy->GetSocketP1(), proxy->GetBuffer1(), exchangingBit * 8);
        thread thr2 = thread(Send, proxy->GetSocketP2(), proxy->GetBuffer2(), exchangingBit * 8);
        thr1.join();
        thr2.join();

        return 0;
    }
    return -1;
}

/**
 * Method for private computation of the derivative of the ReLU function.
 * @param proxy
 * @param x - vector of variables for which to compute ReLU' (the derivative of the ReLU function).
 * @param size - size of vector x
 * @return vector containing the resulting DerivativeReLU(x), must be deleted if not needed anymore.
 */
uint64_t* DerivativeReLU(Party *const proxy, const uint64_t *const x, uint32_t size){
    uint64_t K = (RING_SIZE >> 1); // N is the ring size - 1 = 2^64 -1
    // K is 2^63 - 1
    uint32_t exchangingBit = 2 * size;
    if (proxy->GetPRole() == proxy1 || proxy->GetPRole() == proxy2) {
        // init
        auto* f = new uint64_t [size];
        auto* t = new uint64_t [size];
        for (uint32_t i = 0; i < size; i++){
            f[i] = proxy->GenerateCommonRandom() & 0x1;
            t[i] = x[i] & K; // get first L-1 bit of the share
        }
        K += 1;
        uint64_t* d = ModularConversion(proxy, t, size);
        delete [] t;
        uint64_t* z = Subtract(x, d, size);
        delete[] d;

        // compute parts of a:
        uint64_t values[exchangingBit];
        for (uint32_t i = 0; i < size; i++){
            values[i * 2] = proxy->GetPRole() * f[i] * K - z[i];         // a_0 of the i-th variable
            values[i * 2 + 1] = proxy->GetPRole() * (1 - f[i]) * K - z[i];   // a_1
        }

        // proxy sends a to helper:
        unsigned char *ptr = proxy->GetBuffer1();
        AddValueToCharArray(values, &ptr, exchangingBit);
        Send(proxy->GetSocketHelper(), proxy->GetBuffer1(), exchangingBit * 8);

        // receive fresh share from helper
        Receive(proxy->GetSocketHelper(), proxy->GetBuffer1(), exchangingBit * 8); // order is share of a[0], a[1]
        ptr = proxy->GetBuffer1();

        for (uint64_t c = 0; c < size; c++){
            z[c] = proxy->GetPRole() - ConvertToLong(&ptr);
            if (f[c]) {  // if f is 1  we get a1.
                z[c] = proxy->GetPRole() - ConvertToLong(&ptr);
            }
            else{// increase ptr so that in next iteration not a1 is read:
                (*ptr)+=7;
            }
        }
        delete[] f;
        return z;
    }
    else if (proxy->GetPRole() == helper) {
        K += 1;
        ModularConversion(proxy, nullptr, size);

        Receive(proxy->GetSocketP1(), proxy->GetBuffer1(), exchangingBit * 8);
        Receive(proxy->GetSocketP2(), proxy->GetBuffer2(), exchangingBit * 8);
        unsigned char *ptr1 = proxy->GetBuffer1();
        unsigned char *ptr2 = proxy->GetBuffer2();

        uint64_t reconstructedVals[exchangingBit];
        for (uint64_t v = 0; v < exchangingBit; v++) {
            reconstructedVals[v] = (ConvertToLong(&ptr1) + ConvertToLong(&ptr2)) / K;
        }
        //reassign buffer because ptr1 and ptr2 were incremented
        ptr1 = proxy->GetBuffer1();
        ptr2 = proxy->GetBuffer2();

        uint64_t* tmp = ConvertToUint64(Random1dData(proxy, exchangingBit), exchangingBit);  // values for proxy1
        AddValueToCharArray(tmp, &ptr1, exchangingBit);
        uint64_t * share = Subtract(reconstructedVals, tmp, exchangingBit);
        AddValueToCharArray(share, &ptr2, exchangingBit);

        thread thr1 = thread(Send, proxy->GetSocketP1(), proxy->GetBuffer1(), exchangingBit * 8);
        thread thr2 = thread(Send, proxy->GetSocketP2(), proxy->GetBuffer2(), exchangingBit * 8);
        thr1.join();
        thr2.join();
        delete [] share;

        return nullptr;
    }
    return nullptr;
}

// LAYER helper FUNCTIONS:

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
uint64_t ***Increase(
    const uint64_t *const *const *const input,
    uint32_t channel,
    uint32_t height,
    uint32_t width,
    uint32_t k_dim,
    uint32_t stride
) {
    uint32_t k_size = k_dim * k_dim;
    uint32_t last_row_start = height - k_dim + 1;
    uint32_t last_col_start = width - k_dim + 1;
    uint32_t conv_height = last_row_start/stride;
    uint32_t conv_width = last_col_start/stride;
    // stretch the input for vectorized MatrixVectorMultiply
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
uint64_t **Pad(const uint64_t *const *const input, uint32_t rows, uint64_t cols, uint64_t padding_value, uint32_t padding_size){
    uint32_t padded_row_length = 2*padding_size+cols;
    uint32_t padded_col_length = 2*padding_size+rows;
    auto** padded_input = new uint64_t *[padded_col_length];

    for (uint32_t i = 0; i<padded_col_length; i++){
        padded_input[i] = new uint64_t[padded_row_length];
        // init whole matrix with padding value
        //memset(padded_input[i], padding_value, padded_row_length*8);
        for (int j = 0; j < padding_size; ++j) {
            padded_input[i][j] = padding_value;                     //left padding
            padded_input[i][rows+padding_size + j] = padding_value; //right padding
        }
        if(i >= padding_size && i < (padding_size + rows)){
            //memcpy(padded_input + i*padded_row_length + padding_size, input[i - padding_size], cols * 8);        // copy values
            for (int j = padding_size; j < padded_row_length-padding_size; ++j) {
                padded_input[i][j] = input[i-padding_size][j-padding_size];
            }
        }
        else{
            // pad the top and bottom rows
            for (int j = padding_size; j < padded_row_length-padding_size; ++j) {
                padded_input[i][j] = padding_value;
            }
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
uint64_t * Flatten(const uint64_t *const *const *const images, uint32_t i_height, uint32_t i_width, uint32_t i_number){
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
 * Implements the function of a convolutional layer (ConvolutionalLayer) using ReLU as activation function
 * and then Maxpool with a 2x2 filter if according parameter is set.
 *
 * @param proxy
 * @param input data on which convolution is performed using the provided kernels.
 * The input is supposed to be in 3D shape with WxHxC.
 * @param i_channel number of channels of the input; must be > 0
 * @param i_height size of input in one dimension per i_channel; must be > 0
 * @param i_width size of input in the other dimension; must be > 0
 * @param kernel vector containing all kernel to be used for convolution.
 * The parameter kernel has shape i_channel x output_channel x k_dim * k_dim.
 * So each channel of a kernel is represented by a vector of length k_dim * k_dim.
 * @param k_dim dimension of the symmetric kernel; must be > 0
 * @param output_channel number of kernels with length k_dim * k_dim; must be > 0
 * For each kernel there will be one output channel in the result, where the value at the according location is
 * the sum of the single multiplications between input and kernel of same channel number.
 * @param stride step size of each kernel to shift per iteration; must be > 0
 * @param max_win_height defines window heigth of maxpool, which will be realized after Relu activation.
 * A max_win_height <= 0 indicates that no maxpool shall be performed after Relu activation.
 * @param max_win_width defines window width of maxpool, which will be realized after Relu activation.
 * A max_win_width <= 0 indicates that no maxpool shall be performed after Relu activation.
 * In order for maxpool to be performed it must hold: (max_win_width > 0) and (max_win_height > 0)
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
uint64_t*** ConvolutionalLayer(
    Party *const proxy,
    const uint64_t *const *const *const input,
    uint32_t i_channel,
    uint32_t i_height,
    uint32_t i_width,
    const uint64_t *const *const *const kernel,
    uint32_t k_dim,
    uint32_t output_channel,
    uint32_t stride,
    uint32_t max_win_height,
    uint32_t max_win_width,
    const uint64_t *const bias,
    bool last_conv = false
){
    uint32_t k_size = k_dim * k_dim;
    auto conv_width = static_cast<uint32_t>(floor((i_width - k_dim + 1) / stride));
    auto conv_height = static_cast<uint32_t>(floor((i_height - k_dim + 1) / stride));
    uint32_t conv_len = conv_width * conv_height;
    bool doMaxpooling = ((max_win_width > 0) && (max_win_height > 0) and (max_win_width + max_win_height) > 2);

    if (proxy->GetPRole() == proxy1 || proxy->GetPRole() == proxy2) {
        uint64_t *** stretched_input = Increase(input, i_channel, i_height, i_width, k_dim, stride); // stretched_input is the same for each kernel
        uint64_t *** conv_input = Transpose(stretched_input, i_channel, conv_len, k_size);
        uint64_t ***conv_result = MatrixMatrixMultiply(proxy, kernel, conv_input, i_channel, output_channel, k_size,
                                                       conv_len);

        uint64_t** summed_channel_conv;
        if(i_channel == 1){
            summed_channel_conv = conv_result[0];
        }
        else{
            summed_channel_conv = Add(proxy, conv_result, i_channel, output_channel, conv_len);
        }
        auto *conv_reshaped = new uint64_t [conv_len*output_channel];
        for (int k = 0; k < output_channel; ++k) {
            //memcpy(conv_reshaped + k*conv_len, summed_channel_conv + k, conv_len*8); //number of bytes == length * 8 because 64 bits == 8 bytes
            for (int i = 0; i < conv_len; ++i) {
                conv_reshaped[k*conv_len + i] = summed_channel_conv[k][i] + bias[k];
            }
        }
        for (int i = 0; i < i_channel; ++i) {
            for (int o = 0; o < output_channel; ++o) {
                delete[] conv_result[i][o];
                if(i == 1){ //if i_channel == 1 this will not be reached and summed_channel_conv has been deleted via con_result
                    delete[] summed_channel_conv[o];
                }
            }
            for (int row = 0; row < conv_len; ++row) {
                delete[] stretched_input[i][row];
            }
            for (int k = 0; k < k_size; ++k) {
                delete[] conv_input[i][k];
            }
            delete[] stretched_input[i];
            delete[] conv_input[i];
            delete[] conv_result[i];
            if(i == 1){ //if i_channel == 1 this will not be reached and summed_channel_conv has been deleted via con_result
                delete[] summed_channel_conv;
            }
        }
        delete[] stretched_input;
        delete[] conv_input;
        delete[] conv_result;

        // ACTIVATION:
        uint64_t* conv_activated = ReLU(proxy, conv_reshaped, conv_len * output_channel);
        uint32_t out_width = conv_width;
        uint32_t out_height = conv_height;
        uint32_t out_len = out_height*out_width;
        if (doMaxpooling) {
            out_height /= max_win_height;
            out_width /= max_win_width;
            out_len = out_height*out_width;
            // Maxpool:
            //conv_activated is like all conv_results staked on top of each other, rows not only conv_height but output_channel*conv_height
            conv_activated = Max(proxy, conv_activated, conv_height * output_channel, conv_width, max_win_height,
                                 max_win_width);
        }

        auto ***conv_layer = new uint64_t **[output_channel];
        for (uint32_t k = 0; k < output_channel; k++) {
            if(last_conv){
                if(k == 0){
                    conv_layer = new uint64_t **[1];
                    conv_layer[0] = new uint64_t*[1];
                    conv_layer[0][0] = new uint64_t[output_channel*out_len];
                }
                //BenchmarkConvolutionalLayer[0][0] = conv_activated;
                for (uint32_t row = 0; row < out_height; row++) {
                    for (uint32_t col = 0; col < out_width; col++) {
                        conv_layer[0][0][k*out_len + row*out_width + col] = conv_activated[k*out_len + row * out_width + col];
                    }
                }
            }
            else{
                // bring result in matrix shape
                conv_layer[k] = new uint64_t *[out_height];
                for (uint32_t row = 0; row < out_height; row++) {
                    conv_layer[k][row] = new uint64_t[out_width];
                    for (uint32_t col = 0; col < out_width; col++) {
                        conv_layer[k][row][col] = conv_activated[k*out_len + row * out_width + col];
                    }
                }
            }
        }
        delete[] conv_reshaped;
        delete[] conv_activated;
        return conv_layer;
    }
    else if (proxy->GetPRole() == helper){
        // convolution:
        MatrixMatrixMultiply(proxy, nullptr, nullptr, 0, i_channel * output_channel * k_size * conv_len, 0, 0);
        // ACTIVATION:
        ReLU(proxy, nullptr, conv_len * output_channel);
        if (doMaxpooling){
            // Maxpool:
            Max(proxy, nullptr, output_channel * conv_height, conv_width, max_win_height, max_win_width);
        }
        return nullptr;
    }
    return nullptr;
}

/**
 * Implements functionality of a fully connected layer (FullyConnectedLayer) without activation function.
 * @param proxy
 * @param input the input vector of length in_size to be fully connected to the output nodes.
 * @param in_size length of the input vector (when input is in more dimensional shape, use the Flattening method Flatten before)
 * @param weights to be used must be of shape node_number x in_size
 * @param node_number number of output nodes of this layer
 * @param bias vector of length node_number. For each output node there is one bias value which is added.
 * @return the output layer that have been computed by using the dot product between input values and weights
 * and the according bias added. It will be of length node_number.
 * return vector must be deleted if not needed anymore.
 */
uint64_t* FullyConnectedLayer(
    Party *const proxy,
    const uint64_t *const input,
    int in_size,
    const uint64_t *const *const weights,
    int node_number,
    uint64_t* bias
){
    if (proxy->GetPRole() == proxy1 || proxy->GetPRole() == proxy2){
        uint64_t *output = MatrixVectorMultiply(proxy, weights, input, node_number, in_size);
        uint64_t *added_bias = Add(proxy, output, bias, node_number);
        delete[] output;
        return added_bias;
    }
    else if (proxy->GetPRole() == helper){
        MatrixVectorMultiply(proxy, nullptr, nullptr, node_number * in_size, 0);
        return nullptr;
    }
    else{
        cout << "Error: unknown proxy Role for fully connected layer." << endl;
        return nullptr;
    }
}

#endif //PPAUC_CNN_H

