//
// Created by Seyma Selcan on 02.12.22.
//

#ifndef CECILIA_SORT_H
#define CECILIA_SORT_H

#include "core.h"
#include "../utils/flib.h"
#include "bitset"
#include "string.h"


uint64_t* GenerateF(Party *const proxy, const uint64_t *const a, uint32_t size) {
    if(proxy->GetPRole() == proxy1) {
        auto *result = new uint64_t[2 * size];
        for(int i = 0; i < size; ++i) {
            result[i] = (1 - a[i]);
            result[i + size] = a[i];
        }
        return result;
    }
    else if(proxy->GetPRole() == proxy2) {
        auto *result = new uint64_t[2 * size];
        for(int i = 0; i < size; ++i) {
            result[i] = (0 - a[i]);
            result[i + size] = a[i];
        }
        return result;
    }
    else return nullptr;
}

uint64_t* GenerateS(const uint64_t *const a, uint32_t size) {
    auto *result = new uint64_t[size];
    result[0] = a[0];
    for(int i = 1; i < size; ++i) {
        result[i] = (result[i-1] + a[i]);
    }
    return result;
}

uint64_t* GenerateP(const uint64_t *const t, uint32_t size) {
    auto *result = new uint64_t[size];
    for(int i = 0; i < size; ++i) {
        result[i] = (t[i] + t[i + size]) ;  //removed %L correct?
    }
    return result;
}

uint64_t *GeneratePermutation(Party *const proxy, const uint64_t *const x, uint32_t size) {
    if(proxy->GetPRole() == proxy1 || proxy->GetPRole() == proxy2) {
        uint64_t *f = GenerateF(proxy, x, size);
        uint64_t *s = GenerateS(f, size * 2);
        uint64_t *t = Multiply(proxy, f, s, size * 2);
        uint64_t *p = GenerateP(t, size);
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

uint64_t *GetRandomPermutation(const uint64_t *const randoms, uint32_t size){

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

uint64_t *ApplyPermutation(Party *proxy, const uint64_t *const p, uint64_t *const v, const uint64_t *const pi, uint32_t size) {

    // the second algorithm in the paper will give the shares of sorted v
    auto* r = new uint64_t[size];
    // construct permuted p and v shares
    auto* pip = new uint64_t[size];
    auto* piv = new uint64_t[size];
    auto* pir = new uint64_t[size];
    auto* pv_inv = new uint64_t[size];
    auto* pr_inv = new uint64_t[size];
    long long x;   //pointer lazım mı?
    long long n = (long long) (((long long)1 << 61) - 1);
    x = proxy->GenerateCommonRandom() % n;

    if (proxy->GetPRole() == proxy1) {
        for (int i = 0; i < size; i++) {
            r[i] = proxy->GenerateRandom() % n;   //need to be smaller than n
            v[i] += r[i];
        }
    }

    if (proxy->GetPRole() == proxy1 || proxy->GetPRole() == proxy2) {

        for (int i = 0; i < size; i++) {
            pip[i] = p[pi[i]-1];
            piv[i] = v[pi[i]-1];
        }
        unsigned char *ptr = proxy->GetBuffer1();

        for (int i = 0; i < size; i++) {
            AddValueToCharArray(pip[i], &ptr);
            AddValueToCharArray(piv[i], &ptr);
        }
        thread thr1 = thread(Send, proxy->GetSocketHelper(), proxy->GetBuffer1(), size * 16);
        thread thr2 = thread(Receive, proxy->GetSocketHelper(), proxy->GetBuffer2(), size * 8); //receives a share from pv_inv
        thr1.join();
        thr2.join();
        if(proxy->GetPRole() == proxy1) {

            for (int i = 0; i < size; i++) {
                r[i] = (uint64_t) MultMod((long long) r[i], x, n);
            }

            for (int i = 0; i < size; i++) {
                pir[i] = r[pi[i]-1];
            }
            ptr = proxy->GetBuffer1();
            for (int i = 0; i < size; i++) {
                AddValueToCharArray(pir[i], &ptr);
            }
            Send(proxy->GetSocketHelper(), proxy->GetBuffer1(), size * 8);  //sent pir to helper
        }

        ptr = proxy->GetBuffer2();
        for (int i = 0; i < size; i++) {
            pv_inv[i] = ConvertToLong(&ptr);   //they got the pvinv shares but proxy1 needs to eliminate the effect of r
        }

        if(proxy->GetPRole() == proxy2){   //proxy2
            Receive(proxy->GetSocketHelper(), proxy->GetBuffer1(), size * 8);
            ptr = proxy->GetBuffer1();
            for (int i = 0; i < size; i++) {
                pr_inv[i] = ConvertToLong(&ptr);
            }

            long long xinv = GetModularInverseN(x, n);
            for (int i = 0; i < size; i++) {
                pr_inv[i] = (uint64_t) MultMod((long long) pr_inv[i], xinv, n);
                pv_inv[i] = pv_inv[i] - pr_inv[i];
            }

        }
    }
    else { // helper

        unsigned char *ptr1 = proxy->GetBuffer1();
        unsigned char *ptr2 = proxy->GetBuffer2();

        thread thr1 = thread(Receive, proxy->GetSocketP1(), proxy->GetBuffer1(), size * 16);
        thread thr2 = thread(Receive, proxy->GetSocketP2(), proxy->GetBuffer2(), size * 16);
        thr1.join();
        thr2.join();

        for (int i = 0; i < size; i++) {
            pip[i] = ConvertToLong(&ptr1);
            piv[i] = ConvertToLong(&ptr1);
            pip[i] += ConvertToLong(&ptr2);
            piv[i] += ConvertToLong(&ptr2);
        }

        for(int i = 0; i < size; i++){  //applying inverse permutation, check
            pv_inv[pip[i]-1] = piv[i];   //-1 because permutation starts from 1
        }

        //we need to create shares to send
        ptr1 = proxy->GetBuffer1();
        ptr2 = proxy->GetBuffer2();
        uint64_t tempShare;
        for (int i = 0; i < size; i++) {
            tempShare = proxy->GenerateRandom();
            AddValueToCharArray(tempShare, &ptr1);   //aslında sharelar olmalı
            AddValueToCharArray(pv_inv[i] - tempShare, &ptr2);

        }

        thr1 = thread(Send, proxy->GetSocketP1(), proxy->GetBuffer1(), size * 8);
        thr2 = thread(Send, proxy->GetSocketP2(), proxy->GetBuffer2(), size * 8);
        thr1.join();
        thr2.join();

        Receive(proxy->GetSocketP1(), proxy->GetBuffer1(), size * 8);

        ptr1 = proxy->GetBuffer1();
        for (int i = 0; i < size; i++) {
            pir[i] = ConvertToLong(&ptr1);
        }


        for(int i = 0; i < size; i++){  //applying inverse permutation, check (pip)-1 pir
            pr_inv[pip[i]-1] = pir[i];
        }

        ptr1 = proxy->GetBuffer1();
        for (int i = 0; i < size; i++) {
            AddValueToCharArray(pr_inv[i], &ptr1);
        }

        Send(proxy->GetSocketP2(), proxy->GetBuffer1(), size * 8);


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
    if (proxy->GetPRole() == helper) {
        for(int i = 0; i < LT; ++i) {
            MostSignificantBit(proxy, 0, size);
            GeneratePermutation(proxy, 0, size);
            ApplyPermutation(proxy, 0, 0, 0, size);
        }
        return 0;
    }
    else {  //proxy1 or proxy2
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
            uint64_t* perm = GeneratePermutation(proxy, msb_array, size);  //obtain permutation for msb
            end = chrono::high_resolution_clock::now();
            tt2 = tt2 +chrono::duration_cast<chrono::nanoseconds>(end - start).count();

            start = chrono::high_resolution_clock::now();
            for (int j = 0; j < size ; j++)
                randoms[j] = proxy->GenerateCommonRandom();
            end = chrono::high_resolution_clock::now();
            tt3 = tt3 + chrono::duration_cast<chrono::nanoseconds>(end - start).count();

            start = chrono::high_resolution_clock::now();
            uint64_t* pi = GetRandomPermutation(randoms, size);
            end = chrono::high_resolution_clock::now();
            tt4 = tt4 + chrono::duration_cast<chrono::nanoseconds>(end - start).count();

            start = chrono::high_resolution_clock::now();
            uint64_t* temp = ApplyPermutation(proxy, perm, res, pi, size);  //apply the permutation to the current array
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
#endif //CECILIA_SORT_H