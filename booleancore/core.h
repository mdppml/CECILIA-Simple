//
// Created by Seyma Selcan on 30.03.23.
//

#ifndef BOOL_CORE_H
#define BOOL_CORE_H

#include "../core/Party.h"
#include <thread>
#include <mutex>
#include <bitset>


/**
 * This function is for testing boolean subtraction
 * Sorting protocol will no use it directly
 * */

uint64_t *RECB(Party* proxy, uint64_t *a, uint32_t sz) {

    uint64_t *b = new uint64_t[sz];
    if ( proxy->getPRole() == P1 ) {
        unsigned char *ptr = proxy->getBuffer1();
        for (int i = 0; i < sz; i++) {
            addVal2CharArray(a[i], &ptr);
        }
        thread thr1 = thread(Send,proxy->getSocketP2(), proxy->getBuffer1(), sz*8);
        thread thr2 = thread(Receive,proxy->getSocketP2(), proxy->getBuffer2(), sz*8);
        thr1.join();
        thr2.join();

        ptr = proxy->getBuffer2();
        for (int i = 0; i < sz; i++) {
            b[i] = convert2Long(&ptr);
        }

    }
    else if ( proxy->getPRole() == P2) {
        unsigned char *ptr = proxy->getBuffer1();
        for (int i = 0; i < sz; i++) {
            addVal2CharArray(a[i], &ptr );

        }
        thread thr1 = thread(Send,proxy->getSocketP1(), proxy->getBuffer1(), sz*8);
        thread thr2 = thread(Receive,proxy->getSocketP1(), proxy->getBuffer2(), sz*8);
        thr1.join();
        thr2.join();
        ptr = proxy->getBuffer2();
        for (int i = 0; i < sz; i++) {
            b[i] = convert2Long(&ptr);

        }
    }
    for (int i = 0; i < sz; i++) {
        b[i] = (a[i] ^ b[i]) ;
    }
    return b;
}

uint8_t *RECB(Party* proxy, uint8_t *a, uint32_t sz) {
    uint8_t *b = new uint8_t[sz];

    if ( proxy->getPRole() == P1 ) {
        unsigned char *ptr = proxy->getBuffer1();
        (*ptr) = 0;
        uint8_t bit_index = 7;
        for (int i = 0; i < sz; i++) {
            addBit2CharArray(a[i], &ptr, &bit_index);
        }
        thread thr1 = thread(Send,proxy->getSocketP2(), proxy->getBuffer1(), sz/8+1);
        thread thr2 = thread(Receive,proxy->getSocketP2(), proxy->getBuffer2(), sz/8+1);
        thr1.join();
        thr2.join();

        ptr = proxy->getBuffer2();
        bit_index =7;
        for (int i = 0; i < sz; i++) {
            b[i] = convert2Byte(&ptr, &bit_index);
        }

    }
    else if ( proxy->getPRole() == P2) {
        unsigned char *ptr = proxy->getBuffer1();
        (*ptr) = 0;
        uint8_t bit_index = 7;
        for (int i = 0; i < sz; i++) {
            addBit2CharArray(a[i], &ptr, &bit_index);
        }
        thread thr1 = thread(Send,proxy->getSocketP1(), proxy->getBuffer1(), sz/8+1);
        thread thr2 = thread(Receive,proxy->getSocketP1(), proxy->getBuffer2(), sz/8+1);
        thr1.join();
        thr2.join();
        ptr = proxy->getBuffer2();
        bit_index =7;
        for (int i = 0; i < sz; i++) {
            b[i] = convert2Byte(&ptr, &bit_index);
        }

    }
    for (int i = 0; i < sz; i++) {
        b[i] = (a[i] ^ b[i]) ;
    }
    return b;
}

/**
 * @param mt1 3-by-size array whose rows will be a_i, b_i and c_i, respectively
 * @param mt2 3-by-size array whose rows will be a_i, b_i and c_i, respectively
 * @param size the number of multiplication triples that will be generated
 */
void GenerateBoolMultiplicationTriple(Party* proxy, uint8_t *c1, uint32_t sz) {

    for (int i = 0; i < sz; i++) {
        uint8_t a0 = proxy->generateCommonRandomByte();
        uint8_t a1 = proxy->generateCommonRandomByte2();
        uint8_t b0 = proxy->generateCommonRandomByte();
        uint8_t b1 = proxy->generateCommonRandomByte2();
        uint8_t c0=  proxy->generateCommonRandomByte();
        c1[i] = (((a0^a1)&(b0^b1)) ^ c0); //(a0^a1)*(b0+b1) - c0

    }
}
/** Vectorized AND operation for XOR shared numbers
 * @param a first operand in and operation
 * @param b second operand in and operation
 * @param size number of elements in a and b arrays
 * Each bit is represented with 1 byte so values of a[i] (or b[i]) is either 1 or 0
 * Multiplication triples formulation: a^b = c and mt[0]=a,  mt[1]=b, mt[2]:c
 *
 * */
uint8_t *AND(Party* proxy, uint8_t *a, uint8_t *b, uint32_t size) {
    uint32_t sz =(size/8 +1);
    if (proxy->getPRole() == HELPER) {
        uint8_t *c1 = new uint8_t[sz];
        GenerateBoolMultiplicationTriple(proxy, c1, sz);

        unsigned char *ptr_out2 = proxy->getBuffer2();
        (*ptr_out2) = 0;
        for (int j = 0; j < sz; j++) {
            addVal2CharArray(c1[j], &ptr_out2);
        }

        Send( proxy->getSocketP2(), proxy->getBuffer2(), sz);
        delete[] c1;
        return NULL;
    } else if (proxy->getPRole() == P1 || proxy->getPRole() == P2) {
        uint8_t *mt[3];
        mt[0] = new uint8_t[sz*8]; //a
        mt[1] = new uint8_t[sz*8]; //b
        mt[2] = new uint8_t[sz*8]; //c
        uint8_t *concat_e_f = new uint8_t[size * 2];
        if (proxy->getPRole() == P2) {
            Receive(proxy->getSocketHelper(), proxy->getBuffer1(), sz);
            unsigned char *ptr = proxy->getBuffer1();
            for (int i = 0; i < sz; ++i) {
                auto a_trip = proxy->generateCommonRandomByte2();
                auto b_trip = proxy->generateCommonRandomByte2();
                auto c_trip = ptr[i];
                for (int j = 0; j < 8; ++j) {
                    mt[0][i*8+j] = (a_trip>>j)&0x1;
                    mt[1][i*8+j] = (b_trip>>j)&0x1;
                    mt[2][i*8+j] = (c_trip>>j)&0x1;
                }
            }

            for (int i = 0; i < size; ++i) {
                concat_e_f[i] = a[i] ^ mt[0][i];
                concat_e_f[i + size] = b[i] ^ mt[1][i];
            }

        }
        else { //P1
            for (int i = 0; i < sz; ++i) {
                auto a_trip = proxy->generateCommonRandomByte2();
                auto b_trip = proxy->generateCommonRandomByte2();
                auto c_trip = proxy->generateCommonRandomByte2();
                for (int j = 0; j < 8; ++j) {
                    mt[0][i*8+j] = (a_trip>>j)&0x1;
                    mt[1][i*8+j] = (b_trip>>j)&0x1;
                    mt[2][i*8+j] = (c_trip>>j)&0x1;
                }
            }
            for (int i = 0; i < size; ++i) {
                concat_e_f[i] = a[i] ^ mt[0][i];
                concat_e_f[i + size] = b[i] ^ mt[1][i];
            }
        }
        uint8_t *e_f = RECB(proxy, concat_e_f, size * 2);
        uint8_t *e = e_f;
        uint8_t *f = &e_f[size];
        uint8_t *z = new uint8_t[size];
        for (int i = 0; i < size; i++) {
            z[i] = (proxy->getPRole() & e[i] & f[i]) ^ (f[i] & mt[0][i]) ^ (e[i] & mt[1][i]) ^ mt[2][i];
        }
        delete [] concat_e_f;
        for (auto &i : mt) {
            delete[] i;
        }
        return z;
    }
    return NULL;
}


uint64_t *BooleanSubstract(Party* proxy, uint64_t *a, uint64_t *b, uint32_t sz) {
    if (proxy->getPRole() == P1 ||  proxy->getPRole() == P2) {
        uint64_t *result = new uint64_t[sz];
        for (int i = 0; i < sz; i++) {
            result[i] = 0;
        }
        uint8_t *c = new uint8_t[sz];
        uint8_t *a_bit = new uint8_t[sz];
        uint8_t *b_bit = new uint8_t[sz];
        uint8_t *ainvxc = new uint8_t[sz];
        uint8_t *bxc = new uint8_t[sz];
        for (int i = 0; i < sz; ++i) {
            c[i] = 0;
        }
        for (int j = 0; j<64; j++){
            for (int i = 0; i < sz; i++) {
                //c[i] = 0;
                a_bit[i] = (a[i]>>j)&0x1;
                b_bit[i] = (b[i]>>j)&0x1;
                uint8_t a_inverse_bit = 1 - a_bit[i];
                if(proxy->getPRole() == P1) a_inverse_bit = a_bit[i];
                ainvxc[i] = a_inverse_bit^c[i];
                bxc[i] = b_bit[i]^c[i];
            }
            uint8_t * ainvNbxc = AND(proxy, ainvxc, bxc, sz);

            for (int i = 0; i < sz; i++) {
                result[i] = ((b_bit[i]^c[i]^a_bit[i])<<j)^result[i];
                c[i] = ainvNbxc[i]^c[i];
            }
        }

        delete [] c;
        delete [] a_bit;
        delete [] b_bit;
        delete [] ainvxc;
        delete [] bxc;
        return result;
    }
    else if ( proxy->getPRole() == HELPER){
        for (int j = 0; j<64; j++){
            AND(proxy, 0,0, sz);
        }
    }
}

uint64_t *Convert2XOR(Party* proxy, uint64_t *a, uint32_t sz) {

    if ( proxy->getPRole() == HELPER ) {
        auto *ar = new uint64_t[sz];
        unsigned char *ptr1 = proxy->getBuffer1();
        unsigned char *ptr2 = proxy->getBuffer2();

        thread thr1 = thread(Receive,proxy->getSocketP1(), proxy->getBuffer1(), sz * 8);
        thread thr2 = thread(Receive,proxy->getSocketP2(), proxy->getBuffer2(), sz * 8);
        thr1.join();
        thr2.join();
        for (int i = 0; i < sz; i++) {      //Receive (a+r) shares and add them
            ar[i] = convert2Long(&ptr1);
            ar[i] += convert2Long(&ptr2);
        }


        //we need to create shares to send
        ptr1 = proxy->getBuffer1();
        ptr2 = proxy->getBuffer2();
        uint64_t tempShare;
        for (int i = 0; i < sz; i++) {
            tempShare = proxy->generateRandom();
            addVal2CharArray(tempShare, &ptr1);     // XOR share of (a+r) for P0
            addVal2CharArray(ar[i]^tempShare, &ptr2);   //P1 share
        }
        thr1 = thread(Send, proxy->getSocketP1(), proxy->getBuffer1(), sz * 8);
        thr2 = thread( Send, proxy->getSocketP2(), proxy->getBuffer2(), sz * 8);
        thr1.join();
        thr2.join();

        return nullptr;
    }
    else { //P0 or P1
        unsigned char *ptr = proxy->getBuffer1();
        uint64_t *r = new uint64_t[sz];
        uint64_t *ar = new uint64_t[sz];  //a+r_i
        uint64_t *r_i = new uint64_t[sz];
        for (int i = 0; i < sz; ++i) {
            r[i] = proxy->generateCommonRandom();
            r_i[i] = proxy->createShare(r[i]);
            ar[i] = a[i] +r_i[i];
        }
        ptr = proxy->getBuffer1();
        for (int i = 0; i < sz; i++) {
            addVal2CharArray(ar[i], &ptr);
        }

        Send(proxy->getSocketHelper(), proxy->getBuffer1(), sz * 8);  //sent ar to helper

        Receive(proxy->getSocketHelper(),proxy->getBuffer1(),sz*8);   // receive XOR share of (a+r)

        ptr = proxy->getBuffer1();
        for (int i = 0; i < sz; i++) {
            ar[i] = convert2Long(&ptr);
        }

        if (proxy->getPRole() == P1)
            ar = BooleanSubstract(proxy, ar, r, sz);  //eliminate r
        return ar;

    }
}

#endif
