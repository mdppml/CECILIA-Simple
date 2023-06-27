//
// Created by Seyma Selcan on 02.12.22.
//

#ifndef CECILIA_SORT_H
#define CECILIA_SORT_H

#endif //CECILIA_SORT_H

#include "core.h"
#include "../utils/flib.h"
#include "bitset"
#include "string.h"


uint64_t* generateF(Party *const proxy, const uint64_t *const a, uint32_t size) {
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

uint64_t* generateS(const uint64_t *const a, uint32_t size) {
    auto *result = new uint64_t[size];
    result[0] = a[0];
    for(int i = 1; i < size; ++i) {
        result[i] = (result[i-1] + a[i]);
    }
    return result;
}

uint64_t* generateP(const uint64_t *const t, uint32_t size) {
    auto *result = new uint64_t[size];
    for(int i = 0; i < size; ++i) {
        result[i] = (t[i] + t[i + size]) ;  //removed %L correct?
    }
    return result;
}

uint64_t *generatePermutation(Party *const proxy, const uint64_t *const x, uint32_t size) {
    if(proxy->getPRole() == P1 || proxy->getPRole() == P2) {
        uint64_t *f = generateF(proxy, x, size);
        uint64_t *s = generateS(f, size * 2);
        uint64_t *t = Multiply(proxy, f, s, size * 2);
        uint64_t *p = generateP(t, size);
        delete[] f;
        delete[] s;
        delete[] t;
        return p;
    }
    else{
        Multiply(proxy,0,0,size * 2);
        return 0;
    }
}

uint64_t *getRandomPermutation(const uint64_t *const randoms, uint32_t size){

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

uint64_t *applyPermutation(Party *proxy, const uint64_t *const p, uint64_t *const v, const uint64_t *const pi, uint32_t size) {

    // paperdaki ikinci algoritma, sortlanmis v nin sharelarini verecek
    auto* r = new uint64_t[size];
    // construct permuted p and v shares
    auto* pip = new uint64_t[size];
    auto* piv = new uint64_t[size];
    auto* pir = new uint64_t[size];
    auto* pv_inv = new uint64_t[size];
    auto* pr_inv = new uint64_t[size];
    long long x;   //pointer lazım mı?
    long long n = (long long) (((long long)1 << 61) - 1);
    x = proxy->generateCommonRandom() % n;

    if (proxy->getPRole() == P1) {
        for (int i = 0; i < size; i++) {
            r[i] = proxy->generateRandom() % n;   //need to be smaller than n
            v[i] += r[i];
        }
    }

    if (proxy->getPRole() == P1 || proxy->getPRole() == P2) {

        for (int i = 0; i < size; i++) {
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
                r[i] = (uint64_t) multMod((long long) r[i], x, n);
            }

            for (int i = 0; i < size; i++) {
                pir[i] = r[pi[i]-1];
            }
            ptr = proxy->getBuffer1();
            for (int i = 0; i < size; i++) {
                addVal2CharArray(pir[i], &ptr);
            }
            Send(proxy->getSocketHelper(), proxy->getBuffer1(), size * 8);  //sent pir to helper
        }

        ptr = proxy->getBuffer2();
        for (int i = 0; i < size; i++) {
            pv_inv[i] = convert2Long(&ptr);   //they got the pvinv shares but P1 needs to eliminate the effect of r
        }

        if(proxy->getPRole() == P2){   //P2
            Receive(proxy->getSocketHelper(),proxy->getBuffer1(),size*8);
            ptr = proxy->getBuffer1();
            for (int i = 0; i < size; i++) {
                pr_inv[i] = convert2Long(&ptr);
            }

            long long xinv = getModularInverse_n(x, n);
            for (int i = 0; i < size; i++) {
                pr_inv[i] = (uint64_t) multMod((long long) pr_inv[i], xinv, n);
                pv_inv[i] = pv_inv[i] - pr_inv[i];
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
            addVal2CharArray(pv_inv[i] - tempShare, &ptr2);

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

uint64_t *Sort(Party *const proxy, const uint64_t *const a, uint32_t size) {  //size = size of array
    int LT= 64;
    if (proxy->getPRole() == HELPER) {
        for(int i = 0; i < LT; ++i) {
            MostSignificantBit(proxy, 0, size);
            generatePermutation(proxy, 0, size);
            applyPermutation(proxy,0,0, 0, size);
        }
        return 0;
    }
    else {  //P1 or P2
        double tt1, tt2, tt3, tt4, tt5=0;
        auto *randoms = new uint64_t[size];
        auto* res = new uint64_t[size];
        auto* to_shift = new uint64_t[size];
        for (int j = 0; j < size ; ++j) {
            res[j] = a[j];
            to_shift[j] = res[j];
        }

        for(int i = 0; i < LT; ++i) {
            for (int k = 0; k < size; ++k) {
                to_shift[k] <<=(LT-1-i);
            }
            auto start = chrono::high_resolution_clock::now();
            uint64_t* msb_array = MostSignificantBit(proxy, to_shift, size);   //obtain msb from shifted array
            auto end = chrono::high_resolution_clock::now();
            tt1 = tt1 + chrono::duration_cast<chrono::nanoseconds>(end - start).count();

            start = chrono::high_resolution_clock::now();
            uint64_t* perm = generatePermutation(proxy, msb_array, size);  //obtain permutation for msb
            end = chrono::high_resolution_clock::now();
            tt2 = tt2 +chrono::duration_cast<chrono::nanoseconds>(end - start).count();

            start = chrono::high_resolution_clock::now();
            for (int j = 0; j < size ; j++)
                randoms[j] = proxy->generateCommonRandom();
            end = chrono::high_resolution_clock::now();
            tt3 = tt3 + chrono::duration_cast<chrono::nanoseconds>(end - start).count();

            start = chrono::high_resolution_clock::now();
            uint64_t* pi = getRandomPermutation(randoms, size);
            end = chrono::high_resolution_clock::now();
            tt4 = tt4 + chrono::duration_cast<chrono::nanoseconds>(end - start).count();

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
        cout << "MSBs\t" << tt1*1e-9 << " sec\n";
        cout << "GenSortPermt\t" << tt2*1e-9 << " sec\n";
        cout << "GenComRan\t" << tt3*1e-9 << " sec\n";
        cout << "GenComPerm\t" << tt4*1e-9 << " sec\n";
        cout << "ApplyPerm\t"<< tt5*1e-9 << " sec\n";
        delete[] randoms;
        return res;
    }
    return 0;
}