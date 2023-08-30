//
// Created by Seyma Selcan on 02.12.22.
//

#ifndef CECILIA_SORT_H
#define CECILIA_SORT_H

#include "core.h"
#include "../booleancore/core.h"
#include "../utils/flib.h"
#include "bitset"
#include "string.h"

double mul_time = 0.0;
double comp_time = 0.0;
double x2a_time = 0.0;
double gp_time = 0.0;
double grp_time = 0.0;
double app_time = 0.0;
double t_time = 0.0;
double a2x_time = 0.0;

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
uint64_t* GenerateF(Party *const proxy, const uint64_t *const a, uint32_t size, uint64_t mask) {
    if(proxy->GetPRole() ==proxy1) {
        auto *result = new uint64_t[2 * size];
        for(int i = 0; i < size; ++i) {
            result[i] = (1 - a[i])&mask;
            result[i + size] = a[i];
        }
        return result;
    }
    else if(proxy->GetPRole() == proxy2) {
        auto *result = new uint64_t[2 * size];
        for(int i = 0; i < size; ++i) {
            result[i] = (0 - a[i])&mask;
            result[i + size] = a[i];
        }
        return result;
    }
    else return nullptr;
}

uint64_t* GenerateS(const uint64_t *const a, uint32_t size, uint64_t mask) {
    auto *result = new uint64_t[size];
    result[0] = a[0];
    for(int i = 1; i < size; ++i) {
        result[i] = (result[i-1] + a[i])&mask;
    }
    return result;
}

uint64_t* GenerateP(const uint64_t *const t, uint32_t size, uint64_t mask) {
    auto *result = new uint64_t[size];
    for(int i = 0; i < size; ++i) {
        result[i] = (t[i] + t[i + size])&mask ;  //removed %L correct?
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
        Multiply(proxy, nullptr, nullptr, size * 2);
        return nullptr;
    }
}

uint64_t *GeneratePermutation(Party *const proxy, const uint64_t *const x, uint32_t size, uint32_t ringbits) {
    uint64_t mask =(1<< ringbits) -1;
    if(proxy->GetPRole() ==proxy1 || proxy->GetPRole() == proxy2) {
        uint64_t *f = GenerateF(proxy, x, size, mask);
        uint64_t *s = GenerateS(f, size * 2, mask);
        auto start = chrono::high_resolution_clock::now();
        uint64_t *t = MultiplyNarrow(proxy, f, s, size * 2, ringbits);
        auto end = chrono::high_resolution_clock::now();
        mul_time +=
                chrono::duration_cast<chrono::nanoseconds>(end - start).count()*1e-9;
        uint64_t *p = GenerateP(t, size, mask);
        delete[] f;
        delete[] s;
        delete[] t;
        return p;
    }
    else{   //helper
        MultiplyNarrow(proxy, nullptr, nullptr, size * 2, ringbits);
        return nullptr;
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

/**ApplyPermutation for Sorting: Given sorting permutation for index i and common permutation,
 * this function computes the permuted/sorted shares of data secretly
 * construct permuted p and v shares
 * @param p sorting permutation
 * @param v data shares (arithmetic)
 * @param pi common permutation
 * @param size number of elements in data
 * */
uint64_t *ApplyPermutation(Party *proxy, const uint64_t *const p, uint64_t *const v, const uint64_t *const pi, uint32_t size) {
    auto* r = new uint64_t[size];
    // construct permuted p and v shares
    auto* pip = new uint64_t[size];
    auto* piv = new uint64_t[size];
    auto* pir = new uint64_t[size];
    auto* pv_inv = new uint64_t[size];
    auto* pr_inv = new uint64_t[size];
    long long x;   //pointer lazım mı?
    uint64_t  n =  0x1fffffffffffffff;

    if (proxy->GetPRole() ==proxy1) {
        x = proxy->GenerateCommonRandom() & n;
        for (int i = 0; i < size; i++) {
            r[i] = proxy->GenerateRandom() & n;   //need to be smaller than n
            v[i] += r[i];
            r[i] = (uint64_t) MultMod((long long) r[i], x, n);
        }
        for (int i = 0; i < size; i++) {
            pip[i] = p[pi[i]-1];
            piv[i] = v[pi[i]-1];
            pir[i] = r[pi[i]-1];
        }
        unsigned char *ptr = proxy->GetBuffer1();
        for (int i = 0; i < size; i++) {
            AddValueToCharArray(pip[i], &ptr);
            AddValueToCharArray(piv[i], &ptr);
            AddValueToCharArray(pir[i], &ptr);
        }
        thread thr1 = thread(Send, proxy->GetSocketHelper(), proxy->GetBuffer1(), size * 24);
        thread thr2 = thread( Receive, proxy->GetSocketHelper(),proxy->GetBuffer2(),size * 8); //receives a share from pv_inv
        thr1.join();
        thr2.join();
        ptr = proxy->GetBuffer2();
        for (int i = 0; i < size; i++) {
            pv_inv[i] = ConvertToLong(&ptr);   //they got the pvinv shares but proxy1 needs to eliminate the effect of r
        }
    } else if (proxy->GetPRole() == proxy2) {
        x = proxy->GenerateCommonRandom() & n;
        for (int i = 0; i < size; i++) {
            pip[i] = p[pi[i] - 1];
            piv[i] = v[pi[i] - 1];
        }
        unsigned char *ptr = proxy->GetBuffer1();
        for (int i = 0; i < size; i++) {
            AddValueToCharArray(pip[i], &ptr);
            AddValueToCharArray(piv[i], &ptr);
        }
        thread thr1 = thread(Send, proxy->GetSocketHelper(), proxy->GetBuffer1(), size * 16);
        thread thr2 = thread(Receive, proxy->GetSocketHelper(), proxy->GetBuffer2(),size * 16); //receives a share from pv_inv
        thr1.join();
        thr2.join();
        ptr = proxy->GetBuffer2();
        for (int i = 0; i < size; i++) {
            pv_inv[i] = ConvertToLong(&ptr);   //they got the pvinv shares but proxy1 needs to eliminate the effect of r
            pr_inv[i] = ConvertToLong(&ptr);
        }
        long long xinv = GetModularInverseN(x, n);

        for (int i = 0; i < size; i++) {
            pr_inv[i] = (uint64_t) MultMod((long long) pr_inv[i], xinv, n);
            pv_inv[i] = pv_inv[i] - pr_inv[i];

        }
    }
    else { // helper
        unsigned char *ptr1 = proxy->GetBuffer1();
        unsigned char *ptr2 = proxy->GetBuffer2();

        thread thr1 = thread(Receive,proxy->GetSocketP1(), proxy->GetBuffer1(), size * 24);
        thread thr2 = thread(Receive,proxy->GetSocketP2(), proxy->GetBuffer2(), size * 16);
        thr1.join();
        thr2.join();
        for (int i = 0; i < size; i++) {
            pip[i] = ConvertToLong(&ptr1);
            piv[i] = ConvertToLong(&ptr1);
            pir[i] = ConvertToLong(&ptr1);
            pip[i] += ConvertToLong(&ptr2);
            piv[i] += ConvertToLong(&ptr2);
        }
        for(int i = 0; i < size; i++){  //applying inverse permutation, check
            pv_inv[pip[i]-1] = piv[i];   //-1 because permutation starts from 1
            pr_inv[pip[i]-1] = pir[i];
        }

        //we need to create shares to send
        ptr1 = proxy->GetBuffer1();
        ptr2 = proxy->GetBuffer2();
        uint64_t tempShare;
        for (int i = 0; i < size; i++) {
            tempShare = proxy->GenerateRandom();
            AddValueToCharArray(tempShare, &ptr1);   //aslında sharelar olmalı

            AddValueToCharArray(pv_inv[i] - tempShare, &ptr2);
            AddValueToCharArray(pr_inv[i], &ptr2);

        }

        thr1 = thread(Send, proxy->GetSocketP1(), proxy->GetBuffer1(), size * 8);
        thr2 = thread( Send, proxy->GetSocketP2(), proxy->GetBuffer2(), size * 16);
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
 [[maybe_unused]] uint64_t *ApplyPermutationNarrow(
     Party *const proxy,
     const uint64_t *const p,
     uint64_t *const v,
     const uint64_t *const pi,
     uint32_t size,
     uint64_t ringbits
 ) {
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
    uint64_t  n =  0x7FFFF; //(2^19)-1 mersenne prime

    if (proxy->GetPRole() ==proxy1) {
        x = proxy->GenerateCommonRandom() & n;
        for (int i = 0; i < size; i++) {
            r[i] = proxy->GenerateRandom() & n;   //need to be smaller than nu
            auto vpart = ((v[i]>>1)+(r[i]>>1))&mask2;
            auto bpart = (v[i]^r[i])&0x1;
            v[i] = (vpart<<1)+bpart;
            r[i] = (uint64_t) MultMod((long long) r[i], x, n);
        }
        for (int i = 0; i < size; i++) {
            pip[i] = p[pi[i]-1];
            piv[i] = v[pi[i]-1];
            pir[i] = r[pi[i]-1];
        }
        unsigned char *ptr = proxy->GetBuffer1();
        for (int i = 0; i < size; i++) {
            AddValueToCharArray(pip[i], &ptr, bsz);
            AddValueToCharArray(piv[i], &ptr, bsz);
            AddValueToCharArray(pir[i], &ptr, bsz);
        }
        thread thr1 = thread(Send, proxy->GetSocketHelper(), proxy->GetBuffer1(), size *bsz*3);
        thread thr2 = thread( Receive, proxy->GetSocketHelper(),proxy->GetBuffer2(),size * bsz); //receives a share from pv_inv
        thr1.join();
        thr2.join();
        ptr = proxy->GetBuffer2();
        for (int i = 0; i < size; i++) {
            pv_inv[i] = ConvertToLong(&ptr, bsz);   //they got the pvinv shares but proxy1 needs to eliminate the effect of r
        }
    } else if (proxy->GetPRole() == proxy2) {
        x = proxy->GenerateCommonRandom() & n;
        for (int i = 0; i < size; i++) {
            pip[i] = p[pi[i] - 1];
            piv[i] = v[pi[i] - 1];
        }

        unsigned char *ptr = proxy->GetBuffer1();
        for (int i = 0; i < size; i++) {
            AddValueToCharArray(pip[i], &ptr, bsz);
            AddValueToCharArray(piv[i], &ptr, bsz);
        }
        thread thr1 = thread(Send, proxy->GetSocketHelper(), proxy->GetBuffer1(), size * bsz*2);
        thread thr2 = thread(Receive, proxy->GetSocketHelper(), proxy->GetBuffer2(),size * 2* bsz); //receives a share from pv_inv
        thr1.join();
        thr2.join();
        ptr = proxy->GetBuffer2();
        for (int i = 0; i < size; i++) {
            pv_inv[i] = ConvertToLong(&ptr, bsz);   //they got the pvinv shares but proxy1 needs to eliminate the effect of r
            pr_inv[i] = ConvertToLong(&ptr, bsz);
        }
        long long xinv = GetModularInverseN(x, n);

        for (int i = 0; i < size; i++) {
            pr_inv[i] = (uint64_t) MultMod((long long) pr_inv[i], xinv, n);
            pv_inv[i] = ((((pv_inv[i]>>1)-(pr_inv[i]>>1))&mask2)<<1)+((pv_inv[i]^pr_inv[i])&0x1);
        }
    }
    else { // helper
        unsigned char *ptr1 = proxy->GetBuffer1();
        unsigned char *ptr2 = proxy->GetBuffer2();

        thread thr1 = thread(Receive,proxy->GetSocketP1(), proxy->GetBuffer1(), size * bsz*3);
        thread thr2 = thread(Receive,proxy->GetSocketP2(), proxy->GetBuffer2(), size * bsz*2);
        thr1.join();
        thr2.join();
        for (int i = 0; i < size; i++) {
            pip[i] = ConvertToLong(&ptr1, bsz);
            piv[i] = ConvertToLong(&ptr1,bsz);
            pir[i] = ConvertToLong(&ptr1,bsz);
            pip[i] = (pip[i]+ConvertToLong(&ptr2,bsz))&mask2;
            auto tmp = ConvertToLong(&ptr2,bsz);
            piv[i] = ((((piv[i]>>1)+(tmp>>1))&mask2)<<1)+((piv[i]^tmp)&0x1);
        }
        for(int i = 0; i < size; i++){  //applying inverse permutation, check
            pv_inv[pip[i]-1] = piv[i];   //-1 because permutation starts from 1
            pr_inv[pip[i]-1] = pir[i];
        }

        //we need to create shares to send
        ptr1 = proxy->GetBuffer1();
        ptr2 = proxy->GetBuffer2();
        uint64_t tempShare,tmpshare2;
        for (int i = 0; i < size; i++) {
            tempShare = proxy->GenerateRandom()&mask;
            AddValueToCharArray(tempShare, &ptr1, bsz);   //aslında sharelar olmalı
            tmpshare2 =  ((((pv_inv[i]>>1)-(tempShare>>1))&mask2)<<1)+((pv_inv[i]^tempShare)&0x1);
            AddValueToCharArray(tmpshare2, &ptr2,bsz);
            AddValueToCharArray(pr_inv[i], &ptr2, bsz);

        }

        thr1 = thread(Send, proxy->GetSocketP1(), proxy->GetBuffer1(), size * bsz);
        thr2 = thread( Send, proxy->GetSocketP2(), proxy->GetBuffer2(), size * bsz*2);
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

uint8_t *ApplyPermutationNarrow2(
    Party *const proxy,
    const uint64_t *const p,
    uint8_t *const v,
    const uint64_t *const pi,
    uint32_t size,
    uint64_t ringbits
) {
    auto bsz = ceil((ringbits+1)/8.0);
    auto mask = (1<< (ringbits+1))-1; //ringbits+1 because number to be permuted 24 bit not 23
    auto mask2 = (1<< (ringbits))-1;
    // construct permuted p and v shares
    auto* pip = new uint64_t[size];
    auto* pir = new uint8_t[size];
    auto* pr_inv = new uint8_t[size];


    if (proxy->GetPRole() ==proxy1) {
        uint8_t map[128];
        auto* r = new uint8_t[size];
        for (int i=0;i<128;i++){
            map[i] = proxy->GenerateCommonRandomByte() & 0x1;
        }
        for (int i = 0; i < size; i++) {
            r[i] = proxy->GenerateRandomByte()&0x1;
            v[i] = (v[i] + r[i])&0x1;
            int index = proxy->GenerateRandomByte()&0x7f;   //need to be smaller than nu
            while (map[index] != r[i]) {
                index = proxy->GenerateRandomByte() & 0x7f;
            }
            r[i] = (index<<1)+v[i];
        }
        for (int i = 0; i < size; i++) {
            pip[i] = p[pi[i]-1];
            pir[i] = r[pi[i]-1];
        }
        unsigned char *ptr = proxy->GetBuffer1();
        for (int i = 0; i < size; i++) {
            AddValueToCharArray(pip[i], &ptr, bsz);
            AddValueToCharArray(pir[i], &ptr);
        }
        Send(proxy->GetSocketHelper(), proxy->GetBuffer1(), size*(bsz+1));
        Receive(proxy->GetSocketHelper(), proxy->GetBuffer2(),size); //receives a share from pv_inv
        ptr = proxy->GetBuffer2();
        for (int i = 0; i < size; i++) {
            pr_inv[i] = ConvertToUint8(&ptr);   //they got the pvinv shares but proxy1 needs to eliminate the effect of r
        }
        delete []r;
    } else if (proxy->GetPRole() == proxy2) {
        uint8_t map[128];
        for (int i=0;i<128;i++){
            map[i] = proxy->GenerateCommonRandomByte() & 0x1;
        }
        for (int i = 0; i < size; i++) {
            pip[i] = p[pi[i]-1];
            pir[i] = v[pi[i]-1];
        }
        unsigned char *ptr = proxy->GetBuffer1();
        for (int i = 0; i < size; i++) {
            AddValueToCharArray(pip[i], &ptr, bsz);
            AddValueToCharArray(pir[i], &ptr);
        }
        Send(proxy->GetSocketHelper(), proxy->GetBuffer1(), size*(bsz+1));
        Receive(proxy->GetSocketHelper(), proxy->GetBuffer2(),size); //receives a share from pv_inv
        ptr = proxy->GetBuffer2();
        for (int i = 0; i < size; i++) {
            uint8_t tmp = ConvertToUint8(&ptr);
            pr_inv[i] = ((tmp&0x1)-map[tmp>>1])&0x1;
        }
    }
    else { // helper
        thread thr1 = thread(Receive,proxy->GetSocketP1(), proxy->GetBuffer1(), size*(bsz+1));
        thread thr2 = thread(Receive,proxy->GetSocketP2(), proxy->GetBuffer2(), size*(bsz+1));
        thr1.join();
        thr2.join();
        unsigned char *ptr1 = proxy->GetBuffer1();
        unsigned char *ptr2 = proxy->GetBuffer2();
        for (int i = 0; i < size; i++) {
            pip[i] = ConvertToLong(&ptr1, bsz);
            pir[i] = ConvertToUint8(&ptr1);
            pip[i] = (pip[i]+ConvertToLong(&ptr2,bsz))&mask2;
            auto tmp = ConvertToUint8(&ptr2);
            pir[i] = ((((pir[i]>>1)+(tmp>>1))&0x7f)<<1)+((pir[i]^tmp)&0x1);
        }
        for(int i = 0; i < size; i++){  //applying inverse permutation, check
            pr_inv[pip[i]-1] = pir[i];   //-1 because permutation starts from 1
        }
        //we need to create shares to send
        ptr1 = proxy->GetBuffer1();
        ptr2 = proxy->GetBuffer2();
        uint8_t tmpShare,tmpShare2;
        for (int i = 0; i < size; i++) {
            tmpShare = proxy->GenerateRandomByte()&0x1;
            AddValueToCharArray(tmpShare, &ptr1);   //aslında sharelar olmalı
            tmpShare2 =  (pr_inv[i]^tmpShare);
            AddValueToCharArray(tmpShare2, &ptr2);
        }
        thr1 = thread(Send, proxy->GetSocketP1(), proxy->GetBuffer1(), size);
        thr2 = thread( Send, proxy->GetSocketP2(), proxy->GetBuffer2(), size);
        thr1.join();
        thr2.join();
    }
    delete[] pip;
    delete[] pir;
    return pr_inv;
}

uint64_t *ComposePermutations(
    Party *const proxy,
    const uint64_t *const p,
    const uint64_t *const r,
    uint32_t size,
    uint64_t ringbits
) {
    auto bsz = ceil((ringbits)/8.0);
    auto mask = (1<< (ringbits))-1;
    if (proxy->GetPRole() ==proxy1 || proxy->GetPRole() == proxy2){
        uint64_t *randoms = new uint64_t[size];
        uint64_t *pip = new uint64_t[size];
        uint64_t *pr = new uint64_t[size];
        for (int j = 0; j < size; ++j) {
            randoms[j] = proxy->GenerateCommonRandom();
        }
        uint64_t* pi = GetRandomPermutation(randoms, size);

        unsigned char *ptr = proxy->GetBuffer1();
        for (int i = 0; i < size; i++) { //permute p using pi
            pip[i] = p[pi[i]-1];
            AddValueToCharArray(pip[i], &ptr,bsz);
            AddValueToCharArray(r[i], &ptr,bsz);
        }
        thread thr1 = thread(Send, proxy->GetSocketHelper(), proxy->GetBuffer1(), size * bsz * 2);
        thread thr2 = thread(Receive, proxy->GetSocketHelper(), proxy->GetBuffer2(),size * bsz); //receives a share from pv_inv
        thr1.join();
        thr2.join();

        ptr = proxy->GetBuffer2();
        for (int i = 0; i < size; i++) {
            pr[pi[i]-1] = ConvertToLong(&ptr, bsz); // = pipr[i]
        }
        delete [] pip;
        delete [] randoms;
        return pr;
    }else{
        uint64_t *pip = new uint64_t[size];
        uint64_t *r_prime = new uint64_t[size];
        unsigned char *ptr1 = proxy->GetBuffer1();
        unsigned char *ptr2 = proxy->GetBuffer2();
        thread thr1 = thread(Receive,proxy->GetSocketP1(), proxy->GetBuffer1(), size * bsz*2);
        thread thr2 = thread(Receive,proxy->GetSocketP2(), proxy->GetBuffer2(), size * bsz*2);
        thr1.join();
        thr2.join();
        for (int i = 0; i < size; i++) {
            pip[i] = ConvertToLong(&ptr1, bsz);
            r_prime[i] = ConvertToLong(&ptr1,bsz);
            pip[i] = (pip[i]+ConvertToLong(&ptr2,bsz))&mask;
            r_prime[i] = (r_prime[i]+ConvertToLong(&ptr2,bsz))&mask;
        }
        ptr1 = proxy->GetBuffer1();
        ptr2 = proxy->GetBuffer2();
        for (int i = 0; i < size; i++) {
            auto pipr_prime = r_prime[pip[i]-1];
            auto tempShare = proxy->GenerateRandom()&mask;
            AddValueToCharArray(tempShare, &ptr1, bsz);
            AddValueToCharArray((pipr_prime-tempShare)&mask, &ptr2,bsz);
        }
        thr1 = thread(Send, proxy->GetSocketP1(), proxy->GetBuffer1(), size * bsz);
        thr2 = thread( Send, proxy->GetSocketP2(), proxy->GetBuffer2(), size * bsz);
        thr1.join();
        thr2.join();
        delete [] pip;
        delete [] r_prime;
        return nullptr;
    }
    return nullptr;
}


/**ApplyPermutation for Sorting (for XOR shared data)
 * @param p sorting permutation
 * @param v data shares (XOR shared)
 * @param pi common permutation
 * @param size number of elements in data
 * */
[[maybe_unused]] uint64_t *ApplyPermutationBoolean(Party *const proxy, const uint64_t *const p, uint64_t *const v, const uint64_t *const pi, uint32_t size) {
    auto* r = new uint64_t[size];
    auto* pip = new uint64_t[size];
    auto* piv = new uint64_t[size];
    auto* pir = new uint64_t[size];
    auto* pv_inv = new uint64_t[size];
    auto* pr_inv = new uint64_t[size];
    long long x;   //pointer lazım mı?
    uint64_t  n =  0x1fffffffffffffff;
    if (proxy->GetPRole() ==proxy1) {
        for (int i = 0; i < size; i++) {
            r[i] = proxy->GenerateRandom() & n;   //need to be smaller than n
            v[i] ^= r[i];
        }
    }

    if (proxy->GetPRole() ==proxy1 || proxy->GetPRole() == proxy2) {
        x = MersenneMod(proxy->GenerateCommonRandom(), n, 61); // 61 is log n
        for (int i = 0; i < size; i++) { //permute p and v using pi
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
                pir[i] = r[pi[i]-1];        //permute r'
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

uint64_t **ApplyPermutationVectorised(
    Party *const proxy,
    const uint64_t *const p,
    uint64_t *const *const v,
    const uint64_t *const pi,
    uint32_t size,
    uint32_t categories
) {
    auto** r = new uint64_t*[categories];
    auto* pip = new uint64_t[size];
    auto** piv = new uint64_t*[categories];
    auto** pir = new uint64_t*[categories];
    auto** pv_inv = new uint64_t*[categories];
    auto** pr_inv = new uint64_t*[categories];
    long long x;
    long long n = (long long) (((long long)1 << 61) - 1);

    if (proxy->GetPRole() ==proxy1) {
        x = proxy->GenerateCommonRandom() % n;
        for (int i = 0; i < categories; ++i) {
            r[i] = new uint64_t [size];
            for (int j = 0; j < size; j++) {
                r[i][j] = proxy->GenerateRandom() % n;   //need to be smaller than n
                v[i][j] += r[i][j];
                r[i][j] = (uint64_t) MultMod((long long) r[i][j], x, n);
            }
        }

        for (int i = 0; i < size; ++i) {
            pip[i] = p[pi[i]-1]; //permute p with pi
        }

        for (int i = 0; i < categories; ++i) {
            piv[i] = new uint64_t [size];
            pir[i] = new uint64_t [size];
            for (int j = 0; j < size; j++) {
                piv[i][j] = v[i][pi[j]-1];         //permute data
                pir[i][j] = r[i][pi[j]-1];

            }
        }

        unsigned char *ptr = proxy->GetBuffer1();

        for (int i = 0; i < size; i++) {
            AddValueToCharArray(pip[i], &ptr);
        }

        for (int i = 0; i < categories; ++i) {
            for (int j = 0; j < size; j++) {
                AddValueToCharArray(piv[i][j], &ptr);
                AddValueToCharArray(pir[i][j], &ptr);
            }
        }

        thread thr1 = thread(Send, proxy->GetSocketHelper(), proxy->GetBuffer1(), size * (2*categories+1)* 8); //send pip and pivs
        thread thr2 = thread( Receive, proxy->GetSocketHelper(),proxy->GetBuffer2(),size * categories * 8); //receives a share from pv_invs
        thr1.join();
        thr2.join();

        ptr = proxy->GetBuffer2();
        for (int i = 0; i < categories; ++i) {
            pv_inv[i] = new uint64_t [size];
            for (int j = 0; j < size; j++) {
                pv_inv[i][j] = ConvertToLong(&ptr);   //they got the pvinv shares butproxy1 needs to eliminate the effect of r
            }
        }
    }

    else if (proxy->GetPRole() == proxy2) {
        x = proxy->GenerateCommonRandom() % n;
        for (int i = 0; i < size; ++i) {
            pip[i] = p[pi[i]-1]; //permute p with pi
        }
        for (int i = 0; i < categories; ++i) {
            piv[i] = new uint64_t [size];
            for (int j = 0; j < size; j++) {
                piv[i][j] = v[i][pi[j]-1];         //permute data
            }
        }

        unsigned char *ptr = proxy->GetBuffer1();

        for (int i = 0; i < size; i++) {
            AddValueToCharArray(pip[i], &ptr);
        }

        for (int i = 0; i < categories; ++i) {
            for (int j = 0; j < size; j++) {
                AddValueToCharArray(piv[i][j], &ptr);
            }
        }

        thread thr1 = thread(Send, proxy->GetSocketHelper(), proxy->GetBuffer1(), size * (categories+1)* 8); //send pip and pivs
        thread thr2 = thread( Receive, proxy->GetSocketHelper(),proxy->GetBuffer2(),size * categories * 16); //receives a share from pv_invs
        thr1.join();
        thr2.join();

        ptr = proxy->GetBuffer2();
        for (int i = 0; i < categories; ++i) {
            pv_inv[i] = new uint64_t [size];
            pr_inv[i]=new uint64_t [size];
            for (int j = 0; j < size; j++) {
                pv_inv[i][j] = ConvertToLong(&ptr);   //they got the pvinv shares butproxy1 needs to eliminate the effect of r
                pr_inv[i][j] = ConvertToLong(&ptr);

            }
        }
        long long xinv = GetModularInverseN(x, n);

        for (int i = 0; i < categories; ++i) {
            for (int j = 0; j < size; j++) {
                pr_inv[i][j] = (uint64_t) MultMod((long long) pr_inv[i][j], xinv, n);
                pv_inv[i][j] = pv_inv[i][j] - pr_inv[i][j];                }
        }

    }
    else { // helper

        unsigned char *ptr1 = proxy->GetBuffer1();
        unsigned char *ptr2 = proxy->GetBuffer2();
        thread thr1 = thread(Receive,proxy->GetSocketP1(), proxy->GetBuffer1(), size * (2*categories+1)* 8);
        thread thr2 = thread(Receive,proxy->GetSocketP2(), proxy->GetBuffer2(), size *  (categories+1)* 8);
        thr1.join();
        thr2.join();
        for (int i = 0; i < size; i++) {
            pip[i] = ConvertToLong(&ptr1);
            pip[i] += ConvertToLong(&ptr2);
        }

        for (int i = 0; i < categories; ++i) {
            piv[i] = new uint64_t [size];
            pir[i] = new uint64_t [size];
            for (int j = 0; j < size; j++) {
                piv[i][j] = ConvertToLong(&ptr1);
                pir[i][j] = ConvertToLong(&ptr1);
                piv[i][j] += ConvertToLong(&ptr2);
            }
        }

        for (int i = 0; i < categories; ++i) {
            pv_inv[i]= new uint64_t [size];
            pr_inv[i]= new uint64_t [size];
            for (int j = 0; j < size; j++) {
                pv_inv[i][pip[j]-1] = piv[i][j];   //-1 because permutation starts from 1
                pr_inv[i][pip[j]-1] = pir[i][j];
            }
        }

        //we need to create shares to send
        ptr1 = proxy->GetBuffer1();
        ptr2 = proxy->GetBuffer2();
        uint64_t tempShare;
        for (int i = 0; i < categories; ++i) {
            for (int j = 0; j < size; j++) {
                tempShare = proxy->GenerateRandom();
                AddValueToCharArray(tempShare, &ptr1);   //aslında sharelar olmalı
                AddValueToCharArray(pv_inv[i][j] - tempShare, &ptr2);
                AddValueToCharArray(pr_inv[i][j], &ptr2);
            }
        }

        thr1 = thread(Send, proxy->GetSocketP1(), proxy->GetBuffer1(), size * categories* 8);
        thr2 = thread( Send, proxy->GetSocketP2(), proxy->GetBuffer2(), size * categories* 16);

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

// Sort algo that works with MostSignificantBit on 64-bit integers
uint64_t *Sort(Party *const proxy, const uint64_t *const a, uint32_t size) {  //size = size of array
    int LT= 64;
    if (proxy->GetPRole() == helper) {
        for(int i = 0; i < LT; ++i) {
            MostSignificantBit(proxy, nullptr, size);
            GeneratePermutation(proxy, nullptr, size);
            ApplyPermutation(proxy, nullptr, nullptr, nullptr, size);
        }
        return nullptr;
    }
    else {  //P1 or proxy2
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
            uint64_t* msb_array = MostSignificantBit(proxy, to_shift, size);   //obtain msb from shifted array
            end = chrono::high_resolution_clock::now();
            tt1 = tt1 + chrono::duration_cast<chrono::nanoseconds>(end - start).count();

            start = chrono::high_resolution_clock::now();
            uint64_t* perm = GeneratePermutation(proxy, msb_array, size);  //obtain permutation for msb
            end = chrono::high_resolution_clock::now();
            tt2 = tt2 + chrono::duration_cast<chrono::nanoseconds>(end - start).count();

            start = chrono::high_resolution_clock::now();
            for (int j = 0; j < size ; j++)
                randoms[j] = proxy->GenerateCommonRandom();
            end = chrono::high_resolution_clock::now();
            tt3 =  tt3 + chrono::duration_cast<chrono::nanoseconds>(end - start).count();

            start = chrono::high_resolution_clock::now();
            uint64_t* pi = GetRandomPermutation(randoms, size);
            end = chrono::high_resolution_clock::now();
            tt4 = tt4+ chrono::duration_cast<chrono::nanoseconds>(end - start).count();

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
    return nullptr;
}

/** Sorting algorithm that uses XOR share for bit decomposition
 ** This algo works with ApplyPermutationNarrow2 and has some errors currently in SortNarrow
 ***/
uint64_t *SortNarrow(Party *const proxy, const uint64_t *const a, uint32_t size, uint32_t ringbits) {  //size = size of array
    int LT= 64;
    int bsz = ceil(size/8.0);
    if (proxy->GetPRole() == helper) {
        ArithmeticToXor(proxy, nullptr, size);
        XorToArithmetic3(proxy, nullptr, size);
        GeneratePermutation(proxy, nullptr, size, ringbits);
        for(int i = 1; i < LT ; ++i) {
            ApplyPermutationNarrow2(proxy, nullptr, nullptr, nullptr, size, ringbits);
            XorToArithmetic3(proxy, nullptr, size);
            GeneratePermutation(proxy, nullptr, size, ringbits);
            ComposePermutations(proxy, nullptr, nullptr, size, ringbits);
        }
        return nullptr;
    }
    else {  //P1 or proxy2
        auto start1 = chrono::high_resolution_clock::now();
        auto *randoms = new uint64_t[size];
        auto *dc = new uint8_t[bsz];       //keeps the bits at the current(i) index
        auto *dn = new uint8_t[bsz];       //keeps the bits at the next(i+1) index
        auto tmp = new uint8_t[size];

        auto start = chrono::high_resolution_clock::now();
        auto a_xor = ArithmeticToXor(proxy, a, size);
        auto end = chrono::high_resolution_clock::now();
        a2x_time +=
                chrono::duration_cast<chrono::nanoseconds>(end - start).count()*1e-9;

        uint8_t bit_index = 7;
        uint8_t *ptr = &dc[0];
        dc[0] = 0;
        for (int j = 0; j < size; j++) {
            AddBitToCharArray(((a_xor[j]) & 0x1), &ptr, &bit_index);
        }
        auto dca = XorToArithmetic3(proxy, dc, size);
        auto permG = GeneratePermutation(proxy, dca, size, ringbits);
        for(int i = 1; i < LT; ++i) {
            for (int j = 0; j < size; ++j) {
                tmp[j] = (a_xor[j]>>i)&0x1;
                randoms[j] = proxy->GenerateCommonRandom();
            }
            auto start = chrono::high_resolution_clock::now();
            uint64_t* pi = GetRandomPermutation(randoms, size);
            auto end = chrono::high_resolution_clock::now();
            grp_time +=
                    chrono::duration_cast<chrono::nanoseconds>(end - start).count()*1e-9;



            start = chrono::high_resolution_clock::now();
            auto tmp2 = ApplyPermutationNarrow2(proxy, permG, tmp, pi, size, ringbits);
            end = chrono::high_resolution_clock::now();
            app_time +=
                    chrono::duration_cast<chrono::nanoseconds>(end - start).count()*1e-9;


            bit_index = 7;
            ptr = &dc[0];
            dc[0] = 0;
            for (int j = 0; j < size; ++j) {
                AddBitToCharArray((tmp2[j]) & 0x1, &ptr, &bit_index);
            }
            start = chrono::high_resolution_clock::now();
            dca = XorToArithmetic3(proxy, dc, size);
            end = chrono::high_resolution_clock::now();
            x2a_time +=
                    chrono::duration_cast<chrono::nanoseconds>(end - start).count()*1e-9;


            start = chrono::high_resolution_clock::now();
            auto permC = GeneratePermutation(proxy, dca, size, ringbits);
            end = chrono::high_resolution_clock::now();
            gp_time +=
                    chrono::duration_cast<chrono::nanoseconds>(end - start).count()*1e-9;


            start = chrono::high_resolution_clock::now();
            permG = ComposePermutations(proxy, permG, permC, size, ringbits);
            end = chrono::high_resolution_clock::now();
            comp_time +=
                    chrono::duration_cast<chrono::nanoseconds>(end - start).count()*1e-9;
        }

        delete[] randoms;
        delete[] tmp;
        delete[] dc;
        delete[] dn;
        auto end1 = chrono::high_resolution_clock::now();
        t_time +=
                chrono::duration_cast<chrono::nanoseconds>(end1 - start1).count()*1e-9;
        return permG;

    }
    return nullptr;
}

/**  Vectorised Sorting Algorithm
 *   It gets a vector of vectors as input. Calculates the sorting permutation based on only one of the vectors.
 *   Then applies the sorting permutation on all vectors.
 *   Just like sorting a table based on column
 * @param a - secret share of the matrix that needs to be sorted
 * @param size - size of the elements.
 * @param columns - number of columns in the table.
 * @param index - index of the pivot column.
 * @return The sorted table
 * */
uint64_t **Sort(
    Party *const proxy,
    const uint64_t *const *const a,
    uint32_t size,
    uint32_t categories,
    uint32_t index
) {  //size = size of array
    int LT= 64;
    if (proxy->GetPRole() == helper) {
        for(int i = 0; i < LT; ++i) {
            MostSignificantBit(proxy, nullptr, size);
            GeneratePermutation(proxy, nullptr, size);
            ApplyPermutationVectorised(proxy, nullptr, nullptr, nullptr, size, categories);
        }
        return nullptr;
    }
    else {  //P1 or proxy2
        double tt1, tt2, tt3, tt4, tt5=0;
        auto *randoms = new uint64_t[size];
        auto** res = new uint64_t*[categories];
        auto* to_shift = new uint64_t[size];

        for (int i = 0; i < categories; ++i) {
            res[i] = new uint64_t[size];
            for (int j = 0; j < size; ++j) {
                res[i][j] = a[i][j];
            }
        }

        for (int j = 0; j < size; ++j) {
            to_shift[j] = res[index][j];
        }
        for(int i = 0; i < LT; ++i) {

            for (int k = 0; k < size; ++k) {
                to_shift[k] <<=(LT-1-i);
            }
            uint64_t* msb_array = MostSignificantBit(proxy, to_shift, size);   //obtain msb from shifted array

            uint64_t* perm = GeneratePermutation(proxy, msb_array, size);  //obtain permutation for msb

            for (int j = 0; j < size ; j++)
                randoms[j] = proxy->GenerateCommonRandom();

            uint64_t* pi = GetRandomPermutation(randoms, size);

            res = ApplyPermutationVectorised(proxy, perm, res, pi, size, categories);  //apply the permutation to the current array
            for (int j = 0; j < size; ++j) {
                to_shift[j] = res[index][j];
            }
            delete[] pi;
        }
        delete[] randoms;
        return res;
    }
    return nullptr;
}

#endif //CECILIA_SORT_H