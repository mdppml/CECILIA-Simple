//
// Created by Seyma Selcan on 02.12.22.
//

#ifndef CECILIA_SORT_H
#define CECILIA_SORT_H

#endif //CECILIA_SORT_H

#include "core.h"
#include "../booleancore/core.h"
#include "../utils/flib.h"
#include "bitset"
#include "string.h"


uint64_t* generateF(Party* proxy, const uint64_t* a, uint32_t size) {
    if(proxy->getPRole() == P1) {
        auto *result = new uint64_t[2 * size];
        for(int i = 0; i < size; ++i) {
            result[i] = (1 - a[i]);
            result[i + size] = a[i];
        }
        return result;
    }
    else if(proxy->getPRole() == P2) {
        auto *result = new uint64_t[2 * size];
        for(int i = 0; i < size; ++i) {
            result[i] = (0 - a[i]);
            result[i + size] = a[i];
        }
        return result;
    }
    else return nullptr;
}

uint64_t* generateS(const uint64_t* a, uint32_t size) {
    auto *result = new uint64_t[size];
    result[0] = a[0];
    for(int i = 1; i < size; ++i) {
        result[i] = (result[i-1] + a[i]);
    }
    return result;
}

uint64_t* generateP(const uint64_t* t, uint32_t size) {
    auto *result = new uint64_t[size];
    for(int i = 0; i < size; ++i) {
        result[i] = (t[i] + t[i + size]) ;  //removed %L correct?
    }
    return result;

}
uint64_t* generateF(Party* proxy, const uint64_t* a, uint32_t size, uint64_t mask) {
    if(proxy->getPRole() == P1) {
        auto *result = new uint64_t[2 * size];
        for(int i = 0; i < size; ++i) {
            result[i] = (1 - a[i])&mask;
            result[i + size] = a[i];
        }
        return result;
    }
    else if(proxy->getPRole() == P2) {
        auto *result = new uint64_t[2 * size];
        for(int i = 0; i < size; ++i) {
            result[i] = (0 - a[i])&mask;
            result[i + size] = a[i];
        }
        return result;
    }
    else return nullptr;
}

uint64_t* generateS(const uint64_t* a, uint32_t size, uint64_t mask) {
    auto *result = new uint64_t[size];
    result[0] = a[0];
    for(int i = 1; i < size; ++i) {
        result[i] = (result[i-1] + a[i])&mask;
    }
    return result;
}

uint64_t* generateP(const uint64_t* t, uint32_t size, uint64_t mask) {
    auto *result = new uint64_t[size];
    for(int i = 0; i < size; ++i) {
        result[i] = (t[i] + t[i + size])&mask ;  //removed %L correct?
    }
    return result;

}

uint64_t *generatePermutation(Party *proxy, uint64_t *x, uint32_t size) {
    // paperdaki ilk algoritma, getbitarray in sonucundaki her eleman x olarak giriyor

    if(proxy->getPRole() == P1 || proxy->getPRole() == P2) {
        uint64_t *f = generateF(proxy, x, size);
        uint64_t *s = generateS(f, size * 2);
        uint64_t *t = MUL(proxy, f, s, size * 2);
        uint64_t *p = generateP(t, size);
        delete[] f;
        delete[] s;
        delete[] t;
        return p;
    }
    else{
        MUL(proxy,0,0,size * 2);
        return 0;
    }
}

uint64_t *generatePermutation(Party *proxy, uint64_t *x, uint32_t size, uint32_t ringbits) {
    // paperdaki ilk algoritma, getbitarray in sonucundaki her eleman x olarak giriyor
    uint64_t mask =(1<< ringbits) -1;
    if(proxy->getPRole() == P1 || proxy->getPRole() == P2) {
        uint64_t *f = generateF(proxy, x, size, mask);
        uint64_t *s = generateS(f, size * 2, mask);
        uint64_t *t = MUL(proxy, f, s, size * 2, ringbits);
        uint64_t *p = generateP(t, size, mask);
        delete[] f;
        delete[] s;
        delete[] t;
        return p;
    }
    else{   //HELPER
        MUL(proxy,0,0,size * 2, ringbits);
        return 0;
    }
}

uint64_t *getRandomPermutation(uint64_t randoms[], uint32_t size){
    uint64_t to_permute[size];
    for(int i = 0; i < size; i++)  to_permute[i] = i+1;
    uint64_t* result = new uint64_t[size];
    int last = size-1;
    for(int i = 0; i<size;i++){
        uint64_t index = randoms[i]%(last+1);
        result[i] = to_permute[index];
        to_permute[index] = to_permute[last];
        last--;
    }
    return result;
}

/**ApplyPermutation for Sorting: Given sorting permutation for index i and common permutation,
 * this function computes the permuted/sorted shares of data secretly
 * construct permuted p and v shares
 * @param p sorting permutation
 * @param v data shares (arithmetic)
 * @param pi common permutation
 * @param size number of elements in data
 * */
uint64_t *applyPermutation(Party *proxy, uint64_t *p, uint64_t *v, uint64_t *pi, uint32_t size) {
 //   double tt1, tt2, tt3, tt4, tt5=0;

    // paperdaki ikinci algoritma, sortlanmis v nin sharelarini verecek
    auto* r = new uint64_t[size];
    // construct permuted p and v shares
    auto* pip = new uint64_t[size];
    auto* piv = new uint64_t[size];
    auto* pir = new uint64_t[size];
    auto* pv_inv = new uint64_t[size];
    auto* pr_inv = new uint64_t[size];
    long long x;   //pointer lazım mı?
    uint64_t  n =  0x1fffffffffffffff;

//    using namespace std::chrono;
//    auto microseconds_since_epoch = duration_cast<microseconds>(system_clock::now().time_since_epoch()).count();
//    cout << microseconds_since_epoch << endl;

    if (proxy->getPRole() == P1) {
        x = proxy->generateCommonRandom() & n;
        for (int i = 0; i < size; i++) {
            r[i] = proxy->generateRandom() & n;   //need to be smaller than n
            v[i] += r[i];
            r[i] = (uint64_t) multMod((long long) r[i], x, n);
        }
        for (int i = 0; i < size; i++) {
            pip[i] = p[pi[i]-1];
            piv[i] = v[pi[i]-1];
            pir[i] = r[pi[i]-1];
        }
        unsigned char *ptr = proxy->getBuffer1();
        for (int i = 0; i < size; i++) {
            addVal2CharArray(pip[i], &ptr);
            addVal2CharArray(piv[i], &ptr);
            addVal2CharArray(pir[i], &ptr);
        }
        thread thr1 = thread(Send, proxy->getSocketHelper(), proxy->getBuffer1(), size * 24);
        thread thr2 = thread( Receive, proxy->getSocketHelper(),proxy->getBuffer2(),size * 8); //receives a share from pv_inv
        thr1.join();
        thr2.join();
        ptr = proxy->getBuffer2();
        for (int i = 0; i < size; i++) {
            pv_inv[i] = convert2Long(&ptr);   //they got the pvinv shares but P1 needs to eliminate the effect of r
        }
    } else if (proxy->getPRole() == P2) {
        x = proxy->generateCommonRandom() & n;
        for (int i = 0; i < size; i++) {
            pip[i] = p[pi[i] - 1];
            piv[i] = v[pi[i] - 1];
        }
        unsigned char *ptr = proxy->getBuffer1();
        for (int i = 0; i < size; i++) {
            addVal2CharArray(pip[i], &ptr);
            addVal2CharArray(piv[i], &ptr);
        }
        thread thr1 = thread(Send, proxy->getSocketHelper(), proxy->getBuffer1(), size * 16);
        thread thr2 = thread(Receive, proxy->getSocketHelper(), proxy->getBuffer2(),size * 16); //receives a share from pv_inv
        thr1.join();
        thr2.join();
        ptr = proxy->getBuffer2();
        for (int i = 0; i < size; i++) {
            pv_inv[i] = convert2Long(&ptr);   //they got the pvinv shares but P1 needs to eliminate the effect of r
            pr_inv[i] = convert2Long(&ptr);
        }
        long long xinv = getModularInverse_n(x, n);

        for (int i = 0; i < size; i++) {
            pr_inv[i] = (uint64_t) multMod((long long) pr_inv[i], xinv, n);
            pv_inv[i] = pv_inv[i] - pr_inv[i];

        }
    }
    else { // HELPER
        unsigned char *ptr1 = proxy->getBuffer1();
        unsigned char *ptr2 = proxy->getBuffer2();

        thread thr1 = thread(Receive,proxy->getSocketP1(), proxy->getBuffer1(), size * 24);
        thread thr2 = thread(Receive,proxy->getSocketP2(), proxy->getBuffer2(), size * 16);
        thr1.join();
        thr2.join();
        for (int i = 0; i < size; i++) {
            pip[i] = convert2Long(&ptr1);
            piv[i] = convert2Long(&ptr1);
            pir[i] = convert2Long(&ptr1);
            pip[i] += convert2Long(&ptr2);
            piv[i] += convert2Long(&ptr2);
        }
        for(int i = 0; i < size; i++){  //applying inverse permutation, check
            pv_inv[pip[i]-1] = piv[i];   //-1 because permutation starts from 1
            pr_inv[pip[i]-1] = pir[i];
        }

        //we need to create shares to send
        ptr1 = proxy->getBuffer1();
        ptr2 = proxy->getBuffer2();
        uint64_t tempShare;
        for (int i = 0; i < size; i++) {
            tempShare = proxy->generateRandom();
            addVal2CharArray(tempShare, &ptr1);   //aslında sharelar olmalı

            addVal2CharArray(pv_inv[i] - tempShare, &ptr2);
            addVal2CharArray(pr_inv[i], &ptr2);

        }

        thr1 = thread(Send, proxy->getSocketP1(), proxy->getBuffer1(), size * 8);
        thr2 = thread( Send, proxy->getSocketP2(), proxy->getBuffer2(), size * 16);
        thr1.join();
        thr2.join();

    }
    delete[] r;
    delete[] pip;
    delete[] piv;
    delete[] pir;
    delete[] pr_inv;
    return pv_inv;
}


/**ApplyPermutation for Sorting: Given sorting permutation for index i and common permutation,
 * this function computes the permuted/sorted shares of data secretly
 * construct permuted p and v shares
 * @param p sorting permutation
 * @param v data shares (arithmetic)
 * @param pi common permutation
 * @param size number of elements in data
 * */
uint64_t *applyPermutationN(Party *proxy, uint64_t *p, uint64_t *v, uint64_t *pi, uint32_t size, uint64_t ringbits) {
    //   double tt1, tt2, tt3, tt4, tt5=0;

    // paperdaki ikinci algoritma, sortlanmis v nin sharelarini verecek
    auto bsz = ceil((ringbits+1)/8.0);
    auto mask = (1<< (ringbits+1))-1; //ringbits+1 because number to be permuted 24 bit not 23
    auto mask2 = (1<< (ringbits))-1;
    auto* r = new uint64_t[size];
    // construct permuted p and v shares
    auto* pip = new uint64_t[size];
    auto* piv = new uint64_t[size];
    auto* pir = new uint64_t[size];
    auto* pv_inv = new uint64_t[size];
    auto* pr_inv = new uint64_t[size];
    long long x;   //pointer lazım mı?
    uint64_t  n =  0x7FFFF; //(2^19)-1

//    using namespace std::chrono;
//    auto microseconds_since_epoch = duration_cast<microseconds>(system_clock::now().time_since_epoch()).count();
//    cout << microseconds_since_epoch << endl;

    if (proxy->getPRole() == P1) {
        x = proxy->generateCommonRandom() & n;
        for (int i = 0; i < size; i++) {
            r[i] = proxy->generateRandom() & n;   //need to be smaller than n
            v[i] = (v[i]+r[i])&mask;
            r[i] = (uint64_t) multMod((long long) r[i], x, n);
        }
        for (int i = 0; i < size; i++) {
            pip[i] = p[pi[i]-1];
            piv[i] = v[pi[i]-1];
            pir[i] = r[pi[i]-1];
        }
        unsigned char *ptr = proxy->getBuffer1();
        for (int i = 0; i < size; i++) {
            addVal2CharArray(pip[i], &ptr, bsz);
            addVal2CharArray(piv[i], &ptr, bsz);
            addVal2CharArray(pir[i], &ptr, bsz);
        }
        thread thr1 = thread(Send, proxy->getSocketHelper(), proxy->getBuffer1(), size *bsz*3);
        thread thr2 = thread( Receive, proxy->getSocketHelper(),proxy->getBuffer2(),size * bsz); //receives a share from pv_inv
        thr1.join();
        thr2.join();
        ptr = proxy->getBuffer2();
        for (int i = 0; i < size; i++) {
            pv_inv[i] = convert2Long(&ptr, bsz);   //they got the pvinv shares but P1 needs to eliminate the effect of r
        }
    } else if (proxy->getPRole() == P2) {
        x = proxy->generateCommonRandom() & n;
        for (int i = 0; i < size; i++) {
            pip[i] = p[pi[i] - 1];
            piv[i] = v[pi[i] - 1];
        }
        unsigned char *ptr = proxy->getBuffer1();
        for (int i = 0; i < size; i++) {
            addVal2CharArray(pip[i], &ptr, bsz);
            addVal2CharArray(piv[i], &ptr, bsz);
        }
        thread thr1 = thread(Send, proxy->getSocketHelper(), proxy->getBuffer1(), size * bsz*2);
        thread thr2 = thread(Receive, proxy->getSocketHelper(), proxy->getBuffer2(),size * 2* bsz); //receives a share from pv_inv
        thr1.join();
        thr2.join();
        ptr = proxy->getBuffer2();
        for (int i = 0; i < size; i++) {
            pv_inv[i] = convert2Long(&ptr, bsz);   //they got the pvinv shares but P1 needs to eliminate the effect of r
            pr_inv[i] = convert2Long(&ptr, bsz);
        }
        long long xinv = getModularInverse_n(x, n);

        for (int i = 0; i < size; i++) {
            pr_inv[i] = (uint64_t) multMod((long long) pr_inv[i], xinv, n);
            pv_inv[i] = (pv_inv[i] - pr_inv[i])&mask;

        }
    }
    else { // HELPER
        unsigned char *ptr1 = proxy->getBuffer1();
        unsigned char *ptr2 = proxy->getBuffer2();

        thread thr1 = thread(Receive,proxy->getSocketP1(), proxy->getBuffer1(), size * bsz*3);
        thread thr2 = thread(Receive,proxy->getSocketP2(), proxy->getBuffer2(), size * bsz*2);
        thr1.join();
        thr2.join();
        for (int i = 0; i < size; i++) {
            pip[i] = convert2Long(&ptr1, bsz);
            piv[i] = convert2Long(&ptr1,bsz);
            pir[i] = convert2Long(&ptr1,bsz);
            pip[i] = (pip[i]+convert2Long(&ptr2,bsz))&mask2;
            piv[i] = (piv[i]+convert2Long(&ptr2,bsz))&mask;
        }
        for(int i = 0; i < size; i++){  //applying inverse permutation, check
            pv_inv[pip[i]-1] = piv[i];   //-1 because permutation starts from 1
            pr_inv[pip[i]-1] = pir[i];
        }

        //we need to create shares to send
        ptr1 = proxy->getBuffer1();
        ptr2 = proxy->getBuffer2();
        uint64_t tempShare;
        for (int i = 0; i < size; i++) {
            tempShare = proxy->generateRandom()&mask;
            addVal2CharArray(tempShare, &ptr1, bsz);   //aslında sharelar olmalı
            addVal2CharArray((pv_inv[i] - tempShare)&mask, &ptr2,bsz);
            addVal2CharArray(pr_inv[i], &ptr2, bsz);

        }

        thr1 = thread(Send, proxy->getSocketP1(), proxy->getBuffer1(), size * bsz);
        thr2 = thread( Send, proxy->getSocketP2(), proxy->getBuffer2(), size * bsz*2);
        thr1.join();
        thr2.join();

    }
    delete[] r;
    delete[] pip;
    delete[] piv;
    delete[] pir;
    delete[] pr_inv;
    return pv_inv;
}

/**ApplyPermutation for Sorting (for XOR shared data)
 * @param p sorting permutation
 * @param v data shares (XOR shared)
 * @param pi common permutation
 * @param size number of elements in data
 * */
uint64_t *applyPermutationB(Party *proxy, uint64_t *p, uint64_t *v, uint64_t *pi, uint32_t size) {
    //   double tt1, tt2, tt3, tt4, tt5=0;

    auto* r = new uint64_t[size];
    auto* pip = new uint64_t[size];
    auto* piv = new uint64_t[size];
    auto* pir = new uint64_t[size];
    auto* pv_inv = new uint64_t[size];
    auto* pr_inv = new uint64_t[size];
    long long x;   //pointer lazım mı?
    uint64_t  n =  0x1fffffffffffffff;

//    using namespace std::chrono;
//    auto microseconds_since_epoch = duration_cast<microseconds>(system_clock::now().time_since_epoch()).count();
//    cout << microseconds_since_epoch << endl;
    if (proxy->getPRole() == P1) {
        for (int i = 0; i < size; i++) {
            r[i] = proxy->generateRandom() & n;   //need to be smaller than n
            v[i] ^= r[i];
        }
    }

    if (proxy->getPRole() == P1 || proxy->getPRole() == P2) {
        x = MersenneMod(proxy->generateCommonRandom(), n, 61); // 61 is log n
        for (int i = 0; i < size; i++) { //permute p and v using pi
            pip[i] = p[pi[i]-1];
            piv[i] = v[pi[i]-1];
        }

        unsigned char *ptr = proxy->getBuffer1();

        for (int i = 0; i < size; i++) {
            addVal2CharArray(pip[i], &ptr);
            addVal2CharArray(piv[i], &ptr);
        }

        thread thr1 = thread(Send, proxy->getSocketHelper(), proxy->getBuffer1(), size * 16);
        thread thr2 = thread( Receive, proxy->getSocketHelper(),proxy->getBuffer2(),size * 8); //receives a share from pv_inv
        thr1.join();
        thr2.join();
        if(proxy->getPRole() == P1) {
            for (int i = 0; i < size; i++) {
                r[i] = (uint64_t) multMod((long long) r[i], x, n);      // r' = r*y modn
            }
            for (int i = 0; i < size; i++) {
                pir[i] = r[pi[i]-1];        //permute r'
            }
            ptr = proxy->getBuffer1();
            for (int i = 0; i < size; i++) {
                addVal2CharArray(pir[i], &ptr);
            }
            Send(proxy->getSocketHelper(), proxy->getBuffer1(), size * 8);  //sent pir to helper
        }

        ptr = proxy->getBuffer2();
        for (int i = 0; i < size; i++) {
            pv_inv[i] = convert2Long(&ptr);   //they got the pvinv shares but P2 needs to eliminate the effect of r
        }

        if(proxy->getPRole() == P2){
            Receive(proxy->getSocketHelper(),proxy->getBuffer1(),size*8);
            ptr = proxy->getBuffer1();
            for (int i = 0; i < size; i++) {
                pr_inv[i] = convert2Long(&ptr);
            }

            long long xinv = getModularInverse_n(x, n);
            for (int i = 0; i < size; i++) {
                pr_inv[i] = (uint64_t) multMod((long long) pr_inv[i], xinv, n);
                pv_inv[i] = pv_inv[i]^pr_inv[i]; //P2 eliminates r effect
            }
        }
    }
    else { // HELPER

        unsigned char *ptr1 = proxy->getBuffer1();
        unsigned char *ptr2 = proxy->getBuffer2();

        thread thr1 = thread(Receive,proxy->getSocketP1(), proxy->getBuffer1(), size * 16);
        thread thr2 = thread(Receive,proxy->getSocketP2(), proxy->getBuffer2(), size * 16);
        thr1.join();
        thr2.join();
        for (int i = 0; i < size; i++) {
            pip[i] = convert2Long(&ptr1);
            piv[i] = convert2Long(&ptr1);
            pip[i] += convert2Long(&ptr2);
            piv[i] += convert2Long(&ptr2);
        }
        for (int i = 0; i < size; ++i) {
            cout << pip[i] << endl;
        }

        for(int i = 0; i < size; i++){  //applying inverse permutation, check
            pv_inv[pip[i]-1] = piv[i];   //-1 because permutation starts from 1
        }

        //we need to create shares to send
        ptr1 = proxy->getBuffer1();
        ptr2 = proxy->getBuffer2();
        uint64_t tempShare;
        for (int i = 0; i < size; i++) {
            tempShare = proxy->generateRandom();
            addVal2CharArray(tempShare, &ptr1);   //aslında sharelar olmalı
            addVal2CharArray(pv_inv[i] ^ tempShare, &ptr2);

        }

        thr1 = thread(Send, proxy->getSocketP1(), proxy->getBuffer1(), size * 8);
        thr2 = thread( Send, proxy->getSocketP2(), proxy->getBuffer2(), size * 8);
        thr1.join();
        thr2.join();

        Receive(proxy->getSocketP1(), proxy->getBuffer1(), size * 8);

        ptr1 = proxy->getBuffer1();
        for (int i = 0; i < size; i++) {
            pir[i] = convert2Long(&ptr1);
        }

        for(int i = 0; i < size; i++){  //applying inverse permutation, check (pip)-1 pir
            pr_inv[pip[i]-1] = pir[i];
        }

        ptr1 = proxy->getBuffer1();
        for (int i = 0; i < size; i++) {
            addVal2CharArray(pr_inv[i], &ptr1);
        }

        Send(proxy->getSocketP2(), proxy->getBuffer1(), size * 8);


    }
    delete[] r;
    delete[] pip;
    delete[] piv;
    delete[] pir;
    delete[] pr_inv;
    return pv_inv;
}
uint64_t *SORT(Party *proxy, uint64_t *a, uint32_t size) {  //size = size of array
    int LT= 64;
    if (proxy->getPRole() == HELPER) {
        for(int i = 0; i < LT; ++i) {
            MSB(proxy, 0, size);
            generatePermutation(proxy, 0, size);
            applyPermutation(proxy,0,0, 0, size);
        }
        return 0;
    }
    else {  //P1 or P2
        double tt0 = 0, tt1 = 0, tt2 = 0, tt3=0, tt4=0, tt5=0;
        auto *randoms = new uint64_t[size];
        auto* res = new uint64_t[size];
        auto* to_shift = new uint64_t[size];
        for (int j = 0; j < size ; ++j) {
            res[j] = a[j];
            to_shift[j] = res[j];
        }

        for(int i = 0; i < LT; ++i) {
            auto start = chrono::high_resolution_clock::now();
            for (int k = 0; k < size; ++k) {
                to_shift[k] <<=(LT-1-i);
            }
            auto end = chrono::high_resolution_clock::now();
            tt0 = tt0 + chrono::duration_cast<chrono::nanoseconds>(end - start).count();

            start = chrono::high_resolution_clock::now();
            uint64_t* msb_array = MSB(proxy, to_shift, size);   //obtain msb from shifted array
            end = chrono::high_resolution_clock::now();
            tt1 = tt1 + chrono::duration_cast<chrono::nanoseconds>(end - start).count();

            start = chrono::high_resolution_clock::now();
            uint64_t* perm = generatePermutation(proxy, msb_array, size);  //obtain permutation for msb
            end = chrono::high_resolution_clock::now();
            tt2 = tt2 + chrono::duration_cast<chrono::nanoseconds>(end - start).count();

            start = chrono::high_resolution_clock::now();
            for (int j = 0; j < size ; j++)
                randoms[j] = proxy->generateCommonRandom();
            end = chrono::high_resolution_clock::now();
            tt3 =  tt3 + chrono::duration_cast<chrono::nanoseconds>(end - start).count();

            start = chrono::high_resolution_clock::now();
            uint64_t* pi = getRandomPermutation(randoms, size);
            end = chrono::high_resolution_clock::now();
            tt4 = tt4+ chrono::duration_cast<chrono::nanoseconds>(end - start).count();

            start = chrono::high_resolution_clock::now();
            uint64_t* temp = applyPermutation(proxy,perm,res, pi, size);  //apply the permutation to the current array
            end = chrono::high_resolution_clock::now();
            tt5 = tt5 + chrono::duration_cast<chrono::nanoseconds>(end - start).count();
            for (int j = 0; j < size ; ++j) {
                res[j] = temp[j];  //res always keeps the permuted version
            }

            for (int j = 0; j < size ; ++j) {    //we need to backup res but also have a copy of it to shift and obtain msbs
                to_shift[j] = res[j];
            }
            delete[] temp;
            delete[] pi;
        }
        cout << "to shift\t" << tt0*1e-6 << " msec\n";
        cout << "MSBs\t" << tt1*1e-6 << " msec\n";
        cout << "GenSortPermt\t" << tt2*1e-6 << " msec\n";
        cout << "GenComRan\t" << tt3*1e-6 << " msec\n";
        cout << "GenComPerm\t" << tt4*1e-6 << " msec\n";
        cout << "ApplyPerm\t"<< tt5*1e-6 << " msec\n";
        delete[] randoms;
        delete[] to_shift;
        return res;
    }
    return 0;
}

/** Sorting algorithm that uses XOR share for bit decomposition
 **
 ***/
uint64_t *SORT(Party *proxy, uint64_t *a, uint32_t size, uint32_t ringbits) {  //size = size of array
    int LT= 64;
    int bsz = ceil(size/8.0);
    cout << "ringbits" << ringbits <<  endl;
    if (proxy->getPRole() == HELPER) {

        Arithmetic2XOR(proxy, 0, size);
        cout << "Helper 1" << endl;
        for(int i = 0; i < LT; ++i) {
            XOR2Arithmetic3(proxy, 0, size);
            cout << "Helper 2" << endl;
            generatePermutation(proxy, 0, size, ringbits);
            cout << "Helper 3" << endl;
            applyPermutationN(proxy,0,0, 0, size, ringbits);
            cout << "Helper 4" << endl;
        }
        return 0;
    }
    else {  //P1 or P2
        auto sz2 = size / 8 + 1;
        double tt0 = 0, tt1 = 0, tt2 = 0, tt3 = 0, tt4 = 0, tt5 = 0;
        auto *randoms = new uint64_t[size];
        auto *res = new uint64_t[size];
        auto *permG = new uint64_t[size];   //permutation general: stores the combination of permutations
        auto *dc = new uint8_t[bsz];       //keeps the bits at the current(i) index
        auto *dn = new uint8_t[bsz];       //keeps the bits at the next(i+1) index
        auto tmp = new uint64_t[size];

        auto start = chrono::high_resolution_clock::now();
        auto a_xor = Arithmetic2XOR(proxy, a, size);
        auto end = chrono::high_resolution_clock::now();
        tt0 = chrono::duration_cast<chrono::nanoseconds>(end - start).count();

        for (int j = 0; j < size; ++j) {
            res[j] = a[j];
        }
        for (int i = 0; i < size; ++i) {
            if (proxy->getPRole() ==P1) permG[i] =  i+1;
            else permG[i] =  0;
        }

        uint8_t bit_index = 7;
        uint8_t * ptr = &dc[0];
        for (int j = 0; j < size; ++j) {
            addBit2CharArray(((a[j])&0x1), &ptr, &bit_index);
        }
        for(int i = 0; i < LT; ++i) {
            //get the ith index of all elements and store them in dc
            cout << "In for" << endl;

            cout << "1" << endl;
            auto dca = XOR2Arithmetic3(proxy, dc, size);

            cout << "2" << endl;
            auto permC = generatePermutation(proxy, dca, size, ringbits);
            auto permc_rec = RECN(proxy, permC, size, ringbits);
            for (int j = 0; j < size; ++j) {
                cout << permc_rec[j] << endl;
            }
            cout << "3" << endl;
            // concat dn and permG bits  this is wrong
            for (int j = 0; j < size; ++j) {
                tmp[i] = (permG[j]<<1) + (a[j]>>(i+1))&1;
            }
            cout << "4" << endl;

            for (int j = 0; j < size ; j++)
                randoms[j] = proxy->generateCommonRandom();
            cout << "5" << endl;

            uint64_t* pi = getRandomPermutation(randoms, size);

            auto tmp2 = applyPermutationN(proxy,permC,tmp, pi, size, ringbits);


            cout << "6" << endl;

            bit_index = 7;
            ptr = &dc[0];
            for (int j = 0; j < size; ++j) {
                addBit2CharArray((tmp2[j])&0x1, &ptr, &bit_index);
                permG[j] = tmp2[j]>>1;
            }
            cout << "7" << endl;

        }

        delete[] randoms;
        delete[] a_xor;
        delete[] dc;
        delete[] dn;
        return permG;
    }
    return 0;
}


///** Sorting algorithm that uses XOR share for bit decomposition
// **
// ***/
//uint64_t *SORTv2(Party *proxy, uint64_t *a, uint32_t size) {  //size = size of array
//    int LT= 64;
//    if (proxy->getPRole() == HELPER) {
//
//        Arithmetic2XOR(proxy, 0, size);
//        for(int i = 0; i < LT; ++i) {
//            XOR2Arithmetic(proxy, 0, size);
//            //XOR2Arithmetic2(proxy, 0, size);
//            generatePermutation(proxy, 0, size);
//            applyPermutation(proxy,0,0, 0, size);
//            applyPermutationB(proxy,0,0, 0, size);
//            cout << "Here 4\n";
//        }
//        return 0;
//    }
//    else {  //P1 or P2
//        double tt0 = 0, tt1 = 0, tt2 = 0, tt3=0, tt4=0, tt5=0;
//        auto *randoms = new uint64_t[size];
//        auto* res = new uint64_t[size];
//        auto* to_shift = new uint64_t[size];
//        auto start = chrono::high_resolution_clock::now();
//        auto a_xor = Arithmetic2XOR(proxy, a, size);
//        auto end = chrono::high_resolution_clock::now();
//        tt0 = chrono::duration_cast<chrono::nanoseconds>(end - start).count();
//
//        for (int j = 0; j < size ; ++j) {
//            res[j] = a[j];
//        }
//        for(int i = 0; i < LT; ++i) {
//            for (int j = 0; j < size; ++j) {
//                to_shift[j] = (a_xor[j] >> i) & ((uint64_t)1);
//            }
//            cout << "==========="<< i << "===========" << endl;
//
//            start = chrono::high_resolution_clock::now();
//            //auto a_i = XOR2Arithmetic2(proxy, to_shift, size);  //obtain permutation for msb
//            auto a_i = XOR2Arithmetic(proxy, to_shift, size);  //obtain permutation for msb
//
//            end = chrono::high_resolution_clock::now();
//            tt1 = chrono::duration_cast<chrono::nanoseconds>(end - start).count();
//
//            start = chrono::high_resolution_clock::now();
//            uint64_t* perm = generatePermutation(proxy, a_i, size);  //obtain permutation for msb
//            end = chrono::high_resolution_clock::now();
//            tt2 = tt2 + chrono::duration_cast<chrono::nanoseconds>(end - start).count();
//
//            auto p_rec = REC(proxy, perm, size);
//            for (int j = 0; j < size; ++j) {
//                cout <<  "p_rec[j]" <<  p_rec[j] << endl;
//            }
//            cout  << endl;
//
//            start = chrono::high_resolution_clock::now();
//            for (int j = 0; j < size ; j++)
//                randoms[j] = proxy->generateCommonRandom();
//            end = chrono::high_resolution_clock::now();
//            tt3 =  tt3 + chrono::duration_cast<chrono::nanoseconds>(end - start).count();
//
//            start = chrono::high_resolution_clock::now();
//            uint64_t* pi = getRandomPermutation(randoms, size);
//            end = chrono::high_resolution_clock::now();
//            tt4 = tt4+ chrono::duration_cast<chrono::nanoseconds>(end - start).count();
//            auto a_xor_rec = REC(proxy, res, size);
//            for (int j = 0; j < size; ++j) {
//                cout << "a_xor[j] " << a_xor_rec[j] << endl;
//            }
//            start = chrono::high_resolution_clock::now();
//            uint64_t* temp = applyPermutation(proxy,perm,res, pi, size);  //apply the permutation to the current array
//            auto temp2 = applyPermutationB(proxy,perm,a_xor, pi, size);  //apply the permutation to the current array
//            end = chrono::high_resolution_clock::now();
//            tt5 = tt5 + chrono::duration_cast<chrono::nanoseconds>(end - start).count();
//            a_xor_rec = REC(proxy, temp, size);
//            for (int j = 0; j < size; ++j) {
//                cout << "a_xor_rec " << a_xor_rec[j] << endl;
//            }
//            for (int j = 0; j < size ; ++j) {
//                res[j] = temp[j];  //res always keeps the permuted version
//                a_xor[j] = temp2[j];
//            }
//
//            delete[] temp;
//            delete[] temp2;
//            delete[] pi;
//        }
//        cout << "Arithmetic to XOR\t" << tt0*1e-6 << " msec\n";
//        cout << "MSBs\t" << tt1*1e-6 << " msec\n";
//        cout << "GenSortPermt\t" << tt2*1e-6 << " msec\n";
//        cout << "GenComRan\t" << tt3*1e-6 << " msec\n";
//        cout << "GenComPerm\t" << tt4*1e-6 << " msec\n";
//        cout << "ApplyPerm\t"<< tt5*1e-6 << " msec\n";
//        delete[] randoms;
//        delete[] to_shift;
//        return res;
//    }
//    return 0;
//}


