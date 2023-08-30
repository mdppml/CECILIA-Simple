//
// Created by Seyma Selcan on 30.03.23.
//

#ifndef BOOL_CORE_H
#define BOOL_CORE_H

#include "../core/Party.h"
#include <thread>
#include <mutex>
#include <bitset>
#include <climits>

/**
 * This function is for testing boolean subtraction
 * Sorting protocol will no use it directly
 * */

uint64_t *RECB(Party* proxy, uint64_t *a, uint32_t sz) {

    uint64_t *b = new uint64_t[sz];
    if ( proxy->GetPRole() == proxy1 ) {
        unsigned char *ptr = proxy->GetBuffer1();
        for (int i = 0; i < sz; i++) {
            AddValueToCharArray(a[i], &ptr);
        }
        thread thr1 = thread(Send,proxy->GetSocketP2(), proxy->GetBuffer1(), sz*8);
        thread thr2 = thread(Receive,proxy->GetSocketP2(), proxy->GetBuffer2(), sz*8);
        thr1.join();
        thr2.join();

        ptr = proxy->GetBuffer2();
        for (int i = 0; i < sz; i++) {
            b[i] = ConvertToLong(&ptr);
        }

    }
    else if ( proxy->GetPRole() == proxy2) {
        unsigned char *ptr = proxy->GetBuffer1();
        for (int i = 0; i < sz; i++) {
            AddValueToCharArray(a[i], &ptr );

        }
        thread thr1 = thread(Send,proxy->GetSocketP1(), proxy->GetBuffer1(), sz*8);
        thread thr2 = thread(Receive,proxy->GetSocketP1(), proxy->GetBuffer2(), sz*8);
        thr1.join();
        thr2.join();
        ptr = proxy->GetBuffer2();
        for (int i = 0; i < sz; i++) {
            b[i] = ConvertToLong(&ptr);

        }
    }
    for (int i = 0; i < sz; i++) {
        b[i] = (a[i] ^ b[i]) ;
    }
    return b;
}

uint8_t *RECB(Party* proxy, uint8_t *a, uint32_t sz) {
    uint8_t *b = new uint8_t[sz];

    if ( proxy->GetPRole() == proxy1 ) {
        unsigned char *ptr = proxy->GetBuffer1();
        (*ptr) = 0;
        uint8_t bit_index = 7;
        for (int i = 0; i < sz; i++) {
            AddBitToCharArray(a[i], &ptr, &bit_index);
        }
        thread thr1 = thread(Send,proxy->GetSocketP2(), proxy->GetBuffer1(), sz/8+1);
        thread thr2 = thread(Receive,proxy->GetSocketP2(), proxy->GetBuffer2(), sz/8+1);
        thr1.join();
        thr2.join();

        ptr = proxy->GetBuffer2();
        bit_index =7;
        for (int i = 0; i < sz; i++) {
            b[i] = ConvertToByte(&ptr, &bit_index);
        }

    }
    else if ( proxy->GetPRole() == proxy2) {
        unsigned char *ptr = proxy->GetBuffer1();
        (*ptr) = 0;
        uint8_t bit_index = 7;
        for (int i = 0; i < sz; i++) {
            AddBitToCharArray(a[i], &ptr, &bit_index);
        }
        thread thr1 = thread(Send,proxy->GetSocketP1(), proxy->GetBuffer1(), sz/8+1);
        thread thr2 = thread(Receive,proxy->GetSocketP1(), proxy->GetBuffer2(), sz/8+1);
        thr1.join();
        thr2.join();
        ptr = proxy->GetBuffer2();
        bit_index =7;
        for (int i = 0; i < sz; i++) {
            b[i] = ConvertToByte(&ptr, &bit_index);
        }

    }
    for (int i = 0; i < sz; i++) {
        b[i] = (a[i] ^ b[i]) ;
    }
    return b;
}


uint8_t *RECB2(Party* proxy, uint8_t *a, uint32_t sz) {
    uint8_t *b = new uint8_t[sz];

    if ( proxy->GetPRole() == proxy1 ) {
        unsigned char *ptr = proxy->GetBuffer1();
        for (int i = 0; i < sz; i++) {
            AddValueToCharArray(a[i], &ptr);
        }
        thread thr1 = thread(Send,proxy->GetSocketP2(), proxy->GetBuffer1(), sz);
        thread thr2 = thread(Receive,proxy->GetSocketP2(), proxy->GetBuffer2(), sz);
        thr1.join();
        thr2.join();

        ptr = proxy->GetBuffer2();
        for (int i = 0; i < sz; i++) {
            b[i] = ConvertToUint8(&ptr);
        }
    }
    else if ( proxy->GetPRole() == proxy2) {
        unsigned char *ptr = proxy->GetBuffer1();
        for (int i = 0; i < sz; i++) {
            AddValueToCharArray(a[i], &ptr);
        }
        thread thr1 = thread(Send,proxy->GetSocketP1(), proxy->GetBuffer1(), sz);
        thread thr2 = thread(Receive,proxy->GetSocketP1(), proxy->GetBuffer2(), sz);
        thr1.join();
        thr2.join();

        ptr = proxy->GetBuffer2();
        for (int i = 0; i < sz; i++) {
            b[i] = ConvertToUint8(&ptr);
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
        uint8_t a0 = proxy->GenerateCommonRandomByte();
        uint8_t a1 = proxy->GenerateCommonRandomByte2();
        uint8_t b0 = proxy->GenerateCommonRandomByte();
        uint8_t b1 = proxy->GenerateCommonRandomByte2();
        uint8_t c0= proxy->GenerateCommonRandomByte();
        c1[i] = (((a0^a1)&(b0^b1)) ^ c0); //(a0^a1)*(b0+b1) - c0

    }
}


/** Vectorized And operation for XOR shared numbers
 * @param a first operand in and operation
 * @param b second operand in and operation
 * @param size number of elements in a and b arrays
 * Each bit is represented with 1 byte so values of a[i] (or b[i]) is either 1 or 0
 * Multiplication triples formulation: a^b = c and mt[0]=a,  mt[1]=b, mt[2]:c
 *
 * */
uint8_t *And(Party* proxy, uint8_t *a, uint8_t *b, uint32_t size) {
    uint32_t sz =(size/8 +1);
    if (proxy->GetPRole() == helper) {
        uint8_t *c1 = new uint8_t[sz];
        GenerateBoolMultiplicationTriple(proxy, c1, sz);

        unsigned char *ptr_out2 = proxy->GetBuffer2();
        (*ptr_out2) = 0;
        for (int j = 0; j < sz; j++) {
            AddValueToCharArray(c1[j], &ptr_out2);
        }

        Send( proxy->GetSocketP2(), proxy->GetBuffer2(), sz);
        delete[] c1;
        return nullptr;
    } else if (proxy->GetPRole() == proxy1 || proxy->GetPRole() == proxy2) {
        uint8_t *mt[3];
        mt[0] = new uint8_t[sz*8]; //a
        mt[1] = new uint8_t[sz*8]; //b
        mt[2] = new uint8_t[sz*8]; //c
        uint8_t *concat_e_f = new uint8_t[size * 2];
        if (proxy->GetPRole() == proxy2) {
            Receive(proxy->GetSocketHelper(), proxy->GetBuffer1(), sz);
            unsigned char *ptr = proxy->GetBuffer1();
            for (int i = 0; i < sz; ++i) {
                auto a_trip = proxy->GenerateCommonRandomByte2();
                auto b_trip = proxy->GenerateCommonRandomByte2();
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
                auto a_trip = proxy->GenerateCommonRandomByte2();
                auto b_trip = proxy->GenerateCommonRandomByte2();
                auto c_trip = proxy->GenerateCommonRandomByte2();
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
            z[i] = (proxy->GetPRole() & e[i] & f[i]) ^ (f[i] & mt[0][i]) ^ (e[i] & mt[1][i]) ^ mt[2][i];
        }
        delete [] concat_e_f;
        for (auto &i : mt) {
            delete[] i;
        }
        return z;
    }
    return nullptr;
}

uint8_t *AND2(Party* proxy, uint8_t *a, uint8_t *b, uint32_t size) {
    if (proxy->GetPRole() == helper) {
        uint8_t *c1 = new uint8_t[size];
        GenerateBoolMultiplicationTriple(proxy, c1, size);
        unsigned char *ptr_out2 = proxy->GetBuffer2();
        (*ptr_out2) = 0;
        for (int j = 0; j < size; j++) {
            AddValueToCharArray(c1[j], &ptr_out2);
        }
        Send( proxy->GetSocketP2(), proxy->GetBuffer2(), size);
        delete[] c1;
        return nullptr;
    } else if (proxy->GetPRole() == proxy1 || proxy->GetPRole() == proxy2) {
        uint8_t *mt[3];
        mt[0] = new uint8_t[size]; //a
        mt[1] = new uint8_t[size]; //b
        mt[2] = new uint8_t[size]; //c
        uint8_t *concat_e_f = new uint8_t[size * 2];
        if (proxy->GetPRole() == proxy2) {
            Receive(proxy->GetSocketHelper(), proxy->GetBuffer1(), size);
            unsigned char *ptr = proxy->GetBuffer1();
            for (int i = 0; i < size; ++i) {
                mt[0][i] = proxy->GenerateCommonRandomByte2();
                mt[1][i] = proxy->GenerateCommonRandomByte2();
                mt[2][i] = ptr[i];
            }
            for (int i = 0; i < size; ++i) {
                concat_e_f[i] = a[i] ^ mt[0][i];
                concat_e_f[i + size] = b[i] ^ mt[1][i];
            }
        }
        else { //P1
            for (int i = 0; i < size; ++i) {
                mt[0][i] = proxy->GenerateCommonRandomByte2();
                mt[1][i] = proxy->GenerateCommonRandomByte2();
                mt[2][i] = proxy->GenerateCommonRandomByte2();
            }
            for (int i = 0; i < size; ++i) {
                concat_e_f[i] = a[i] ^ mt[0][i];
                concat_e_f[i + size] = b[i] ^ mt[1][i];
            }
        }
        uint8_t *e_f = RECB2(proxy, concat_e_f, size * 2);
        uint8_t *e = e_f;
        uint8_t *f = &e_f[size];
        uint8_t *z = new uint8_t[size];
        uint8_t role = (proxy->GetPRole()<<8)-proxy->GetPRole();
        for (int i = 0; i < size; i++) {
            z[i] = ((role & e[i] & f[i]) ^ (f[i] & mt[0][i]) ^ (e[i] & mt[1][i]) ^ mt[2][i]);
        }
        delete [] concat_e_f;
        for (auto &i : mt) {
            delete[] i;
        }
        return z;
    }
    return nullptr;
}


uint64_t *BooleanSubstract(Party* proxy, uint64_t *a, uint64_t *b, uint32_t sz) {
    if (proxy->GetPRole() == proxy1 ||  proxy->GetPRole() == proxy2) {
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
                if(proxy->GetPRole() == proxy1) a_inverse_bit = a_bit[i];
                ainvxc[i] = a_inverse_bit^c[i];
                bxc[i] = b_bit[i]^c[i];
            }
            uint8_t * ainvNbxc = And(proxy, ainvxc, bxc, sz);

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
    else if ( proxy->GetPRole() == helper){
        for (int j = 0; j<64; j++){
            And(proxy, nullptr, nullptr, sz);
        }
    }
    return nullptr;
}

/**
 * A variation of BooleanSubstract: Data is stored in a more compact way.
 * Instead of representing bits as bytes it directly uses bits
 *
 * */
uint64_t *BooleanSubstract2(Party* proxy, uint64_t *a, uint64_t *b, uint32_t sz) {
    if (proxy->GetPRole() == proxy1 ||  proxy->GetPRole() == proxy2) {
        uint64_t *result = new uint64_t[sz];
        uint8_t *c = new uint8_t[sz];
        uint8_t *a_bit = new uint8_t[sz];
        uint8_t *b_bit = new uint8_t[sz];
        uint32_t sz2 =  sz/8+1;
        uint8_t *ainvxc = new uint8_t[sz2];
        uint8_t *bxc = new uint8_t[sz2];
        for (int i = 0; i < sz; ++i) {
            c[i] = 0;
            result[i] = 0;
        }
        for (int j = 0; j<64; j++){
            uint8_t bit_index1 =7;
            uint8_t bit_index2 =7;
            ainvxc[0] = 0;
            bxc[0] = 0;
            unsigned char *ptr1 = &ainvxc[0];
            unsigned char *ptr2 = &bxc[0];
            for (int i = 0; i < sz; i++) {
                a_bit[i] = (a[i]>>j)&0x1;
                b_bit[i] = (b[i]>>j)&0x1;
                uint8_t a_inverse_bit = 1 - a_bit[i];
                if(proxy->GetPRole() == proxy1) a_inverse_bit = a_bit[i];
                AddBitToCharArray(a_inverse_bit ^ c[i], &ptr1, &bit_index1);
                AddBitToCharArray(b_bit[i] ^ c[i], &ptr2, &bit_index2);
            }
            uint8_t * ainvNbxc = AND2(proxy, ainvxc, bxc, sz2);

            ptr1 = &ainvNbxc[0];
            bit_index1 =7;
            for (int i = 0; i < sz; i++) {
                result[i] = ((b_bit[i]^c[i]^a_bit[i])<<j)^result[i];
                c[i] = ConvertToByte(&ptr1, &bit_index1) ^ c[i];
            }
        }

        delete [] c;
        delete [] a_bit;
        delete [] b_bit;
        delete [] ainvxc;
        delete [] bxc;
        return result;
    }
    else if ( proxy->GetPRole() == helper){
        for (int j = 0; j<64; j++){
            AND2(proxy, 0,0, sz/8+1);
        }
    }
    return nullptr;
}

/**Protocol to convert Arithmetic shares to XOR shares
 * @param a Arithmetic share
 * @param sz number of elements in the share
 * */
uint64_t *Arithmetic2XOR(Party *const proxy, const uint64_t *const a, uint32_t sz) {

    if ( proxy->GetPRole() == helper ) {
        auto *ar = new uint64_t[sz];
        unsigned char *ptr1 = proxy->GetBuffer1();
        unsigned char *ptr2 = proxy->GetBuffer2();

        thread thr1 = thread(Receive,proxy->GetSocketP1(), proxy->GetBuffer1(), sz * 8);
        thread thr2 = thread(Receive,proxy->GetSocketP2(), proxy->GetBuffer2(), sz * 8);
        thr1.join();
        thr2.join();
        for (int i = 0; i < sz; i++) {      //Receive (a+r) shares and add them
            ar[i] = ConvertToLong(&ptr1);
            ar[i] += ConvertToLong(&ptr2);
        }

        //we need to create shares to send
        ptr1 = proxy->GetBuffer1();
        ptr2 = proxy->GetBuffer2();
        uint64_t tempShare;
        for (int i = 0; i < sz; i++) {
            tempShare = proxy->GenerateRandom();
            AddValueToCharArray(tempShare, &ptr1);     // XOR share of (a+r) for P0
            AddValueToCharArray(ar[i]^tempShare, &ptr2);   //P1 share
        }
        thr1 = thread(Send, proxy->GetSocketP1(), proxy->GetBuffer1(), sz * 8);
        thr2 = thread( Send, proxy->GetSocketP2(), proxy->GetBuffer2(), sz * 8);
        thr1.join();
        thr2.join();
        BooleanSubstract2(proxy, 0, 0, sz);
        return nullptr;
    }
    else { //P0 or proxy1
        unsigned char *ptr = proxy->GetBuffer1();
        uint64_t *r = new uint64_t[sz];
        uint64_t *ar = new uint64_t[sz];  //a+r_i
        uint64_t *r_i = new uint64_t[sz];
        uint64_t *r_i_xor = new uint64_t[sz];
        for (int i = 0; i < sz; ++i) {
            r[i] = proxy->GenerateCommonRandom();
            r_i[i] = proxy->CreateShare(r[i]);
            r_i_xor[i] = proxy->createXORShare(r[i]);
            ar[i] = a[i] + r_i[i];
        }
        ptr = proxy->GetBuffer1();
        for (int i = 0; i < sz; i++) {
            AddValueToCharArray(ar[i], &ptr);
        }
        Send(proxy->GetSocketHelper(), proxy->GetBuffer1(), sz * 8);  //sent ar to helper

        Receive(proxy->GetSocketHelper(),proxy->GetBuffer1(),sz*8);   // receive XOR share of (a+r)

        ptr = proxy->GetBuffer1();
        for (int i = 0; i < sz; i++) {
            ar[i] = ConvertToLong(&ptr);
        }

        auto result = BooleanSubstract2(proxy, ar, r_i_xor, sz);

        delete [] r;
        delete [] r_i;
        return result;

    }
}


/**Protocol for converting XOR shares of 64bit values to Arithmetic shares
 * Simply the inverse of Arithmetic2XOR function
 * @param a XOR share
 * @param sz number of elements in the share
 * */
uint64_t *XOR2Arithmetic(Party* proxy, uint64_t *a, uint32_t sz) {

    if ( proxy->GetPRole() == helper ) {
        auto *ar = new uint64_t[sz];
        unsigned char *ptr1 = proxy->GetBuffer1();
        unsigned char *ptr2 = proxy->GetBuffer2();

        BooleanSubstract2(proxy, 0, 0, sz);

        thread thr1 = thread(Receive,proxy->GetSocketP1(), proxy->GetBuffer1(), sz * 8);
        thread thr2 = thread(Receive,proxy->GetSocketP2(), proxy->GetBuffer2(), sz * 8);
        thr1.join();
        thr2.join();
        for (int i = 0; i < sz; i++) {      //Receive (a-r) XOR shares and recreate it
            ar[i] = ConvertToLong(&ptr1);
            ar[i] ^= ConvertToLong(&ptr2);
        }

        //we need to create shares to send
        ptr1 = proxy->GetBuffer1();
        ptr2 = proxy->GetBuffer2();
        uint64_t tempShare;
        for (int i = 0; i < sz; i++) {
            tempShare = proxy->GenerateRandom();
            AddValueToCharArray(tempShare, &ptr1);     // arithmetic share of (a+r) for P0
            AddValueToCharArray(ar[i]-tempShare, &ptr2);   //P1 share
        }
        thr1 = thread(Send, proxy->GetSocketP1(), proxy->GetBuffer1(), sz * 8);
        thr2 = thread( Send, proxy->GetSocketP2(), proxy->GetBuffer2(), sz * 8);
        thr1.join();
        thr2.join();
        return nullptr;
    }
    else { //P0 or proxy1
        unsigned char *ptr = proxy->GetBuffer1();
        uint64_t *r = new uint64_t[sz];
        uint64_t *r_i = new uint64_t[sz];
        uint64_t *r_i_xor = new uint64_t[sz];
        for (int i = 0; i < sz; ++i) {
            r[i] = 0;//proxy->GenerateCommonRandom();
            if(proxy->GetPRole()== proxy1){
                r_i[i] = 0;//proxy->CreateShare(r[i]);
                r_i_xor[i] = 4;//proxy->createXORShare(r[i]);
            }
            else{
                r_i[i] = 0;//proxy->CreateShare(r[i]);
                r_i_xor[i] = 4;//proxy->createXORShare(r[i]);
            }
        }
        auto ar = BooleanSubstract2(proxy, a, r_i_xor, sz);

        ptr = proxy->GetBuffer1();
        for (int i = 0; i < sz; i++) {
            AddValueToCharArray(ar[i], &ptr);
        }
        Send(proxy->GetSocketHelper(), proxy->GetBuffer1(), sz * 8);  //sent ar to helper

        Receive(proxy->GetSocketHelper(),proxy->GetBuffer1(),sz*8);   // receive XOR share of (a+r)

        ptr = proxy->GetBuffer1();
        for (int i = 0; i < sz; i++) {
            ar[i] = ConvertToLong(&ptr)+r_i[i];
        }

        delete [] r;
        delete [] r_i;
        return ar;

    }
}

/**Protocol for converting XOR shares of bits to Arithmetic shares
 * @param a XOR share
 * @param sz number of elements in the share
 * */
uint64_t *XOR2Arithmetic2(Party* proxy, uint8_t *a, uint32_t sz) {
    uint32_t bsz = sz/8 +1;
    if ( proxy->GetPRole() == helper ) {
        auto *a1 = new uint8_t[sz];
        auto *a2 = new uint8_t[sz];

        thread thr1 = thread(Receive,proxy->GetSocketP1(), proxy->GetBuffer1(), bsz*2);//it will receive 2 things from P0
        thread thr2 = thread(Receive,proxy->GetSocketP2(), proxy->GetBuffer2(), bsz);
        thr1.join();
        thr2.join();

        unsigned char *ptr1 = proxy->GetBuffer1();
        unsigned char *ptr2 = proxy->GetBuffer2();

        uint8_t bit_index1 =7;
        uint8_t bit_index2 =7;
        for (int i = 0; i < sz; i++) {
            uint8_t tmp = ConvertToByte(&ptr2, &bit_index2);
            a1[i] = ConvertToByte(&ptr1, &bit_index1) ^ tmp;   //Recreate and store first possibility in a1
            a2[i] = ConvertToByte(&ptr1, &bit_index1) ^ tmp;
        }
        //we need to create shares to send
        ptr1 = proxy->GetBuffer1();
        ptr2 = proxy->GetBuffer2();
        uint64_t tempShare;
        for (int i = 0; i < sz; i++) {
            tempShare = proxy->GenerateCommonRandom();           //P0 share for 1st possibility
            AddValueToCharArray(a1[i]-tempShare, &ptr2);   //P1 share for 1st possibility
            tempShare = proxy->GenerateCommonRandom2();          //P1 share for 2nd possibility
            AddValueToCharArray(a2[i]-tempShare, &ptr1);   //P0 share for 2nd possibility
        }

        Send( proxy->GetSocketP1(), proxy->GetBuffer1(), sz * 8);
        Send( proxy->GetSocketP2(), proxy->GetBuffer2(), sz * 8);
        delete [] a1;
        delete [] a2;

        return nullptr;
    }
    else { //P0 or proxy1
        unsigned char *ptr = proxy->GetBuffer1();
        uint8_t *r = new uint8_t[bsz];
        uint64_t *result = new uint64_t[sz];

        if (proxy->GetPRole() == proxy1) {

            ptr = proxy->GetBuffer1();
            (*ptr) = 0;
            for (int i = 0; i < bsz; ++i) {
                r[i] = proxy->GenerateCommonRandomByte();
                uint8_t bit_index = 7;
                for (int j = 7; j >= 0; j--) {
                        uint8_t bit = ((a[i]>>j)&0x1)^((r[i]>>j)&0x1);
                    AddBitToCharArray(bit, &ptr, &bit_index);
                    AddBitToCharArray(bit ^ 1, &ptr, &bit_index);
                }
            }

            Send(proxy->GetSocketHelper(), proxy->GetBuffer1(), bsz*2);  //sent ar to helper

            auto r1 = new uint64_t[sz];
            for (int i = 0; i < sz; ++i) {
                r1[i] = proxy->GenerateCommonRandom2();  //this will be a share of first possibility
            }

            Receive(proxy->GetSocketHelper(), proxy->GetBuffer1(), sz*8);   // receive Arithmetic share for second possibility

            ptr = proxy->GetBuffer1();
            for (int i = 0; i < sz; i++) {
                auto select = ((r[i/8]>>(7-(i&7))) & 0x1);
                result[i] = (1- select) * r1[i] + select * ConvertToLong(&ptr);      // if select is 0 take the first possibility else second
            }
        }
        else {  //P2
            ptr = proxy->GetBuffer1();
            for (int i = 0; i < bsz; i++) {
                r[i] = proxy->GenerateCommonRandomByte();
            }

            for (int i = 0; i < bsz; i++) {
                AddValueToCharArray(a[i], &ptr);
            }

            Send(proxy->GetSocketHelper(), proxy->GetBuffer1(), bsz);  //sent ar to helper

            auto r1 = new uint64_t[sz];
            for (int i = 0; i < sz; ++i) {
                r1[i] = proxy->GenerateCommonRandom2();  //this will be a share of second possibility
            }

            Receive(proxy->GetSocketHelper(), proxy->GetBuffer1(), sz * 8);   // receive the share of first possibility

            ptr = proxy->GetBuffer1();
            for (int i = 0; i < sz; i++) {
                auto select = ((r[i/8]>>(7-(i&7))) & 0x1);
                result[i] = (1- select) * ConvertToLong(&ptr) + select * r1[i];      // if select is 0 take the first possibility else second
            }

        }
        delete [] r;
        return result;

    }
}

/**Protocol for converting XOR shares of bits to Arithmetic shares (1->23)
 * @param a XOR share
 * @param sz number of elements in the share
 * */
uint64_t *XOR2Arithmetic3(Party* proxy, uint8_t *a, uint32_t sz) {
    uint32_t mask = 0xfffff;
    uint32_t bsz = sz/8 +1;
    if ( proxy->GetPRole() == helper ) {
        auto *a1 = new uint8_t[sz];
        auto *a2 = new uint8_t[sz];

        thread thr1 = thread(Receive,proxy->GetSocketP1(), proxy->GetBuffer1(), bsz*2);//it will receive 2 things from P0
        thread thr2 = thread(Receive,proxy->GetSocketP2(), proxy->GetBuffer2(), bsz);
        thr1.join();
        thr2.join();

        unsigned char *ptr1 = proxy->GetBuffer1();
        unsigned char *ptr2 = proxy->GetBuffer2();

        uint8_t bit_index1 =7;
        uint8_t bit_index2 =7;
        for (int i = 0; i < sz; i++) {
            uint8_t tmp = ConvertToByte(&ptr2, &bit_index2);
            a1[i] = ConvertToByte(&ptr1, &bit_index1) ^ tmp;   //Recreate and store first possibility in a1
            a2[i] = ConvertToByte(&ptr1, &bit_index1) ^ tmp;
        }
        //we need to create shares to send
        ptr1 = proxy->GetBuffer1();
        ptr2 = proxy->GetBuffer2();
        uint64_t tempShare;
        for (int i = 0; i < sz; i++) {
            tempShare = proxy->GenerateCommonRandom()&mask;           //P0 share for 1st possibility
            AddValueToCharArray((a1[i]-tempShare)&mask, &ptr2,3);   //P1 share for 1st possibility
            tempShare = proxy->GenerateCommonRandom2()&mask;          //P1 share for 2nd possibility
            AddValueToCharArray((a2[i]-tempShare)&mask, &ptr1,3);   //P0 share for 2nd possibility

        }

        Send( proxy->GetSocketP1(), proxy->GetBuffer1(), sz * 3);
        Send( proxy->GetSocketP2(), proxy->GetBuffer2(), sz * 3);
        delete [] a1;
        delete [] a2;

        return nullptr;
    }
    else { //P0 or proxy1
        unsigned char *ptr = proxy->GetBuffer1();
        uint8_t *r = new uint8_t[bsz];
        uint64_t *result = new uint64_t[sz];

        if (proxy->GetPRole() == proxy1) {

            ptr = proxy->GetBuffer1();
            (*ptr) = 0;
            uint8_t bit_index = 7;
            for (int i = 0; i < bsz; ++i) {
                r[i] = proxy->GenerateCommonRandomByte();
                for (int j = 7; j >= 0; j--) {
                    uint8_t bit = ((a[i]>>j)&0x1)^((r[i]>>j)&0x1);
                    AddBitToCharArray(bit, &ptr, &bit_index);
                    AddBitToCharArray(bit ^ 1, &ptr, &bit_index);
                }
            }

            Send(proxy->GetSocketHelper(), proxy->GetBuffer1(), bsz*2);  //sent ar to helper

            auto r1 = new uint64_t[sz];
            for (int i = 0; i < sz; ++i) {
                r1[i] = proxy->GenerateCommonRandom2()&mask;  //this will be a share of first possibility
            }

            Receive(proxy->GetSocketHelper(), proxy->GetBuffer1(), sz*3);   // receive Arithmetic share for second possibility

            ptr = proxy->GetBuffer1();
            for (int i = 0; i < sz; i++) {
                auto select = ((r[i/8]>>(7-(i&7))) & 0x1);
                auto tmp = ConvertToLong(&ptr,3);
                result[i] = (1- select) * r1[i] + select * tmp;      // if select is 0 take the first possibility else second
            }
        }
        else {  //P2
            ptr = proxy->GetBuffer1();
            for (int i = 0; i < bsz; i++) {
                r[i] = proxy->GenerateCommonRandomByte();
            }

            for (int i = 0; i < bsz; i++) {
                AddValueToCharArray(a[i], &ptr);
            }

            Send(proxy->GetSocketHelper(), proxy->GetBuffer1(), bsz);  //sent ar to helper

            auto r1 = new uint64_t[sz];
            for (int i = 0; i < sz; ++i) {
                r1[i] = proxy->GenerateCommonRandom2()&mask;  //this will be a share of second possibility
            }

            Receive(proxy->GetSocketHelper(), proxy->GetBuffer1(), sz * 3);   // receive the share of first possibility

            ptr = proxy->GetBuffer1();
            for (int i = 0; i < sz; i++) {
                auto select = ((r[i/8]>>(7-(i&7))) & 0x1);
                auto tmp = ConvertToLong(&ptr,3);
                result[i] = (1- select) * tmp + select * r1[i];      // if select is 0 take the first possibility else second
            }

        }
        delete [] r;
        return result;

    }
}
#endif
