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
        result[i] = (t[i] + t[i + size])&mask ;
    }
    return result;

}

uint64_t *GeneratePermutation(Party *const proxy, const uint64_t *const x, uint32_t size) {
    if(proxy->GetPRole() == proxy1 || proxy->GetPRole() == proxy2) {
        uint64_t *f = GenerateF(proxy, x, size);
        uint64_t *s = GenerateS(f, size * 2);
        uint64_t *t = Multiply(proxy, f, s, size * 2,0);
        uint64_t *p = GenerateP(t, size);
        delete[] f;
        delete[] s;
        delete[] t;
        return p;
    }
    else{
        Multiply(proxy, nullptr, nullptr, size * 2, 0);
        return nullptr;
    }
}

uint64_t *GeneratePermutation(Party *const proxy, const uint64_t *const x, uint32_t size, uint32_t ringbits) {
    uint64_t mask =(1<< ringbits) -1;
    if(proxy->GetPRole() ==proxy1 || proxy->GetPRole() == proxy2) {
        uint64_t *f = GenerateF(proxy, x, size, mask);
        uint64_t *s = GenerateS(f, size * 2, mask);
        auto start = chrono::high_resolution_clock::now();
        uint64_t *t = MultiplyNarrow(proxy, f, s, size * 2, ringbits, 0);
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
        MultiplyNarrow(proxy, nullptr, nullptr, size * 2, ringbits,0 );
        return nullptr;
    }
}
void GetRandomPermutationInplace(const uint64_t *const randoms, uint32_t size, uint64_t*result){
    uint64_t * to_permute=new std::uint64_t[size];
    for(int i = 0; i < size; i++)  to_permute[i] = i+1;
    int last = size-1;
    for(int i = 0; i<size;i++){
        uint64_t index = randoms[i]%(last+1);
        result[i] = to_permute[index];
        to_permute[index] = to_permute[last];
        last--;
    }
    delete [] to_permute;
}

uint64_t *GetRandomPermutation(const uint64_t *const randoms, uint32_t size){
    uint64_t * to_permute=new std::uint64_t[size];
    for(int i = 0; i < size; i++)  to_permute[i] = i+1;
    uint64_t* result = new uint64_t[size];
    int last = size-1;
    for(int i = 0; i<size;i++){
        uint64_t index = randoms[i]%(last+1);
        result[i] = to_permute[index];
        to_permute[index] = to_permute[last];
        last--;
    }
    delete [] to_permute;
    return result;
}
uint64_t *Arithmetic2Permutation(Party *const proxy, const uint64_t *const p, uint32_t size, uint32_t ringbits) {
    uint64_t mask =(1<< ringbits) -1;
    auto bsz = ceil((ringbits)/8.0);
    uint64_t*  pip= new std::uint64_t[size];
    uint64_t* randoms = new std::uint64_t[size];
    uint64_t*  perm= new std::uint64_t[size];

    if(proxy->GetPRole() == proxy1 || proxy->GetPRole() == proxy2) {
        for (int j = 0; j < size; ++j) {
            randoms[j] = proxy->GenerateCommonRandom();
        }
        uint64_t* pi = GetRandomPermutation(randoms, size);
        delete [] randoms;

        unsigned char *ptr = proxy->GetBuffer1();
        for (int i = 0; i < size; ++i) {
            AddValueToCharArray(p[pi[i]-1], &ptr, bsz);     //permute p with pi and add to buffer to send to helper
        }
        Send( proxy->GetSocketHelper(), proxy->GetBuffer1(), size * bsz); //send pip to helper
        Receive( proxy->GetSocketHelper(),proxy->GetBuffer2(),size * bsz); //receives a share from pv_inv
        unsigned char *ptr2 = proxy->GetBuffer2();

        if (proxy->GetPRole() == proxy1){
            for (int i = 0; i < size; i++) {
                pip[i] = ConvertToLong(&ptr2, bsz);   //receive pip permutation component
            }
            for (int i = 0; i < size; ++i) {
                perm[pi[i]-1] = pip[i];     //permute pip with pi^-1 and add to buffer to get rid of pi
            }
        }
        else{
            for (int i = 0; i < size; i++) {
                perm[i] = ConvertToLong(&ptr2, bsz);   //receive pip permutation component
            }
        }
        return perm;
    }
    else{   //helper

        unsigned char *ptr1 = proxy->GetBuffer1();
        unsigned char *ptr2 = proxy->GetBuffer2();

        thread thr1 = thread(Receive,proxy->GetSocketP1(), proxy->GetBuffer1(), size * bsz);
        thread thr2 = thread(Receive,proxy->GetSocketP2(), proxy->GetBuffer2(), size * bsz);
        thr1.join();
        thr2.join();
        for (int i = 0; i < size; i++) {
            pip[i] = ConvertToLong(&ptr1, bsz);
            pip[i] = (pip[i]+ConvertToLong(&ptr2, bsz))&mask;      //reconstruct pip
        }
        for (int j = 0; j < size; ++j) {
            randoms[j] = (proxy->GenerateRandom());
        }
        uint64_t* e = GetRandomPermutation(randoms, size);      //permutation components for p0
        delete [] randoms;
        for (int i = 0; i < size; ++i) {
            perm[e[i]-1] = pip[i];                  //permutation component for p1
        }

        ptr1 = proxy->GetBuffer1();
        ptr2 = proxy->GetBuffer2();
        uint64_t tempShare;
        for (int i = 0; i < size; i++) {
            AddValueToCharArray(e[i]&mask, &ptr1, bsz);
            AddValueToCharArray(perm[i]&mask, &ptr2, bsz);
        }

        thr1 = thread(Send, proxy->GetSocketP1(), proxy->GetBuffer1(), size * bsz);
        thr2 = thread( Send, proxy->GetSocketP2(), proxy->GetBuffer2(), size * bsz);
        thr1.join();
        thr2.join();
        return nullptr;
    }
    delete [] pip;

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
            pv_inv[i] = ConvertToLong(&ptr);   //they got the pvinv shares but P1 needs to eliminate the effect of r
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
            pv_inv[i] = ConvertToLong(&ptr);   //they got the pvinv shares but P1 needs to eliminate the effect of r
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
    //   double tt1, tt2, tt3, tt4, tt5=0;

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
    long long x;
    uint64_t  n =  0x7FFFF; //(2^19)-1 mersenne prime

//    using namespace std::chrono;
//    auto microseconds_since_epoch = duration_cast<microseconds>(system_clock::now().time_since_epoch()).count();
//    cout << microseconds_since_epoch << endl;

    if (proxy->GetPRole() == proxy1) {
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
            pv_inv[i] = ConvertToLong(&ptr, bsz);   //they got the pvinv shares but P1 needs to eliminate the effect of r
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
            pv_inv[i] = ConvertToLong(&ptr, bsz);   //they got the pvinv shares but P1 needs to eliminate the effect of r
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


    if (proxy->GetPRole() == proxy1) {
        uint8_t map[128];
        auto* r = new uint8_t[size];
        for (int i=0;i<128;i++){
            map[i] = proxy->GenerateCommonRandomByte()&0x1;
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
            pr_inv[i] = ConvertToUint8(&ptr);   //they got the pvinv shares but P1 needs to eliminate the effect of r
        }
        delete []r;
    } else if (proxy->GetPRole() == proxy2) {
        uint8_t map[128];
        for (int i=0;i<128;i++){
            map[i] = proxy->GenerateCommonRandomByte()&0x1;
        }
        for (int i = 0; i < size; i++) {
            pir[i] = v[pi[i]-1];
        }
        unsigned char *ptr = proxy->GetBuffer1();
        for (int i = 0; i < size; i++) {
            AddValueToCharArray(p[i], &ptr, bsz);
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
        uint64_t *pip1 = new uint64_t[size];
        uint64_t *pip2 = new uint64_t[size];
        unsigned char *ptr1 = proxy->GetBuffer1();
        unsigned char *ptr2 = proxy->GetBuffer2();
        for (int i = 0; i < size; i++) {
            pip1[i] = ConvertToLong(&ptr1, bsz);
            pir[i] = ConvertToUint8(&ptr1);
            pip2[i] = ConvertToLong(&ptr2,bsz)&mask2;
            auto tmp = ConvertToUint8(&ptr2);
            pir[i] = ((((pir[i]>>1)+(tmp>>1))&0x7f)<<1)+((pir[i]^tmp)&0x1);
        }

        for (int i = 0; i < size; ++i) {
            pip[i] = pip2[pip1[i]-1];
        }
        delete [] pip1;
        delete [] pip2;
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

/**Compose two permutations a◯b in 1 round
 * Less Communication overhead than ComposePermutations2Rounds
 *@param a is the first permutation's (arithmetic)share
 *@param b is the second permutation's (arithmetic)share
 *@param size is the number of elements in the vector
 * */
uint64_t *ComposePermutations1Round(
        Party *const proxy,
        const uint64_t *const a,
        const uint64_t *const b,
        uint32_t size,
        uint64_t ringbits
) {
    auto bsz = ceil((ringbits)/8.0);
    auto mask = (1<< (ringbits))-1;

    if (proxy->GetPRole() == proxy1 || proxy->GetPRole() == proxy2){
        uint64_t *randoms_t = new uint64_t[size];
        uint64_t *randoms_k = new uint64_t[size];
        auto t = new uint64_t[size];
        uint64_t* k = new uint64_t[size];
        for (int j = 0; j < size; ++j) {
            randoms_t[j] = proxy->GenerateCommonRandom();
            randoms_k[j] = proxy->GenerateCommonRandom();
        }

        thread thr0 = thread(GetRandomPermutationInplace, randoms_t, size, t);
        thread thr1 = thread(GetRandomPermutationInplace, randoms_k, size, k);
        thr1.join();
        thr0.join();
        delete [] randoms_k;
        delete [] randoms_t;

        uint64_t *pia = new uint64_t[size];
        uint64_t *a_perm = new uint64_t[size];
        uint64_t *b_perm = new uint64_t[size];

        unsigned char *ptr2 = proxy->GetBuffer2();
        unsigned char *ptr = proxy->GetBuffer1();

        if (proxy->GetPRole() == proxy1){
            //pia = t◯a_perm and pib = k◯b
            for(int i = 0; i < size; i++){
                pia[i] = a[t[i]-1];
                b_perm[i] = b[k[i]-1];      //this keeps the kb0
            }

            ptr = proxy->GetBuffer1();
            //(t◯a_perm)◯t = pia◯t
            for(int i = 0; i < size; i++){
                a_perm[i] = t[pia[i]-1];
                AddValueToCharArray(a_perm[i], &ptr,bsz);
                AddValueToCharArray(b_perm[i], &ptr,bsz);
            }
            Send( proxy->GetSocketHelper(), proxy->GetBuffer1(), size * bsz * 2);
            Receive( proxy->GetSocketHelper(), proxy->GetBuffer2(),size * bsz); //
            ptr2 = proxy->GetBuffer2();
            for (int i = 0; i < size; i++) {
                pia[t[i]-1] = ConvertToLong(&ptr2, bsz); //receive and get rid of t
            }
        }
        else{
            uint64_t *k_inv = new uint64_t[size];
            //pia = (t^-1)◯a_perm and b_perm = k◯b1
            for(int i = 0; i < size; i++){
                pia[t[i]-1] = a[i];
                b_perm[i] = b[k[i]-1];
                k_inv[k[i]-1] = i+1;
            }
            ptr = proxy->GetBuffer1();
            //(t^-1)◯a_perm◯(k^-1) = pia◯(k^-1)
            for(int i = 0; i < size; i++){
                a_perm[i] = k_inv[pia[i]-1];
                AddValueToCharArray(a_perm[i], &ptr,bsz);
                AddValueToCharArray(b_perm[i], &ptr,bsz);
            }
            Send( proxy->GetSocketHelper(), proxy->GetBuffer1(), size * bsz * 2);
            delete [] k_inv;
            Receive( proxy->GetSocketHelper(), proxy->GetBuffer2(),size * bsz); //
            ptr2 = proxy->GetBuffer2();
            for (int i = 0; i < size; i++) {
                pia[i] = ConvertToLong(&ptr2, bsz); //receive and get rid of t
            }
        }

        delete [] a_perm;
        delete [] b_perm;
        delete [] t;
        delete [] k;
        return pia;


    }else{ //helper
        uint64_t *pia = new uint64_t[size];
        uint64_t *randoms_e = new uint64_t[size];
        uint64_t *e = new uint64_t[size];
        uint64_t *f = new uint64_t[size];
        uint64_t *h = new uint64_t[size];
        unsigned char *ptr1 = proxy->GetBuffer1();
        unsigned char *ptr2 = proxy->GetBuffer2();

        thread thr1 = thread(Receive,proxy->GetSocketP1(), proxy->GetBuffer1(), size * bsz*2);
        thread thr2 = thread(Receive,proxy->GetSocketP2(), proxy->GetBuffer2(), size * bsz*2);
        thr1.join();
        thr2.join();
        ptr1 = proxy->GetBuffer1();
        ptr2 = proxy->GetBuffer2();
        for (int i = 0; i < size; i++) {
            e[i] = ConvertToLong(&ptr1, bsz);   //ta0t
            h[i] = ConvertToLong(&ptr1,bsz);     //kb0
            f[i] = ConvertToLong(&ptr2,bsz);     //(t^-1)a1(k^-1)
            h[i] = (h[i]+ConvertToLong(&ptr2,bsz))&mask;     //kb = kb0+kb1
        }

        for (int i = 0; i < size; i++) {
            pia[i] = f[e[i]-1];   //t(a0)t(t^-1)(a1)(k^-1)
        }

        for (int i = 0; i < size; i++) {
            e[i] = h[pia[i]-1];   //tab
        }

        for (int j = 0; j < size; ++j) {
            randoms_e[j] = proxy->GenerateRandom();
        }

        //Create permutation elements(shares)
        //pia = e◯f and pib= h◯j
        uint64_t* e_ = GetRandomPermutation(randoms_e, size);
        delete [] randoms_e;
        for (int i = 0; i < size; ++i) {
            f[e_[i]-1]= e[i];
        }

        ptr1 = proxy->GetBuffer1();
        ptr2 = proxy->GetBuffer2();
        for (int i = 0; i < size; i++) {
            AddValueToCharArray(e_[i], &ptr1, bsz);
            AddValueToCharArray(f[i], &ptr2,bsz);
        }

        thr1 = thread(Send, proxy->GetSocketP1(), proxy->GetBuffer1(), size * bsz);
        thr2 = thread(Send, proxy->GetSocketP2(), proxy->GetBuffer2(), size * bsz);
        thr1.join();
        thr2.join();
        delete [] pia;
        delete [] e;
        delete [] f;
        delete [] h;
        return nullptr;
    }
    return nullptr;
}

/**Compose two permutations a◯b in 2 rounds
 *@param a is the first permutation's (arithmetic)share
 *@param b is the second permutation's (arithmetic)share
 *@param size is the number of elements in the vector
 * */
uint64_t *ComposePermutations2Rounds(
        Party *const proxy,
        const uint64_t *const a,
        const uint64_t *const b,
        uint32_t size,
        uint64_t ringbits
) {
    auto bsz = ceil((ringbits)/8.0);
    auto mask = (1<< (ringbits))-1;

    if (proxy->GetPRole() == proxy1 || proxy->GetPRole() == proxy2){
        uint64_t *randoms = new uint64_t[size];
        uint64_t *randoms_t = new uint64_t[size];
        uint64_t *randoms_k = new uint64_t[size];
        uint64_t *randoms_m = new uint64_t[size];
        auto t = new uint64_t[size];
        uint64_t* k = new uint64_t[size];
        uint64_t* m = new uint64_t[size];
        for (int j = 0; j < size; ++j) {
            randoms[j] = proxy->GenerateCommonRandom();
            randoms_t[j] = proxy->GenerateCommonRandom();
            randoms_k[j] = proxy->GenerateCommonRandom();
            randoms_m[j] = proxy->GenerateCommonRandom();
        }
        uint64_t* pi = new uint64_t[size];
        thread thr = thread(GetRandomPermutationInplace, randoms, size, pi);
        thread thr0 = thread(GetRandomPermutationInplace, randoms_t, size, t);
        thread thr1 = thread(GetRandomPermutationInplace, randoms_k, size, k);
        thread thr2 = thread(GetRandomPermutationInplace, randoms_m, size, m);
        thr.join();
        thr1.join();
        thr2.join();
        thr0.join();
        delete [] randoms;
        delete [] randoms_k;
        delete [] randoms_t;
        delete [] randoms_m;
        uint64_t *pia = new uint64_t[size];
        uint64_t *pib = new uint64_t[size];
        uint64_t *a_perm = new uint64_t[size];
        uint64_t *b_perm = new uint64_t[size];

        unsigned char *ptr2 = proxy->GetBuffer2();
        unsigned char *ptr = proxy->GetBuffer1();
        for (int i = 0; i < size; i++) { //permute a and b using pi
            pia[i] = a[pi[i]-1];
            pib[i] = b[pi[i]-1];
            AddValueToCharArray(pia[i], &ptr,bsz);
            AddValueToCharArray(pib[i], &ptr,bsz);
        }
        thr1 = thread(Send, proxy->GetSocketHelper(), proxy->GetBuffer1(), size * bsz * 2);
        thr2 = thread(Receive, proxy->GetSocketHelper(), proxy->GetBuffer2(),size * bsz); //receives a share from pv_inv
        thr1.join();
        thr2.join();


        if (proxy->GetPRole() == proxy1){
            uint64_t *randoms_e = new uint64_t[size];
            for (int j = 0; j < size; ++j) {
                randoms_e[j] = proxy->GenerateCommonRandom2();
            }

            //Create permutation elements(shares)
            //pia = e◯f and pib= h◯j
            pia = GetRandomPermutation(randoms_e, size);    //generate pia perm share e
            delete [] randoms_e;
            ptr2 = proxy->GetBuffer2();
            for (int i = 0; i < size; i++) {
                b_perm[pi[i]-1] = ConvertToLong(&ptr2, bsz); // receive pib permutation shares and eliminate pi
                a_perm[pi[i]-1] = pia[i];   //eliminate pi

            }

            //pia = t◯a_perm and pib = k◯b_perm
            for(int i = 0; i < size; i++){
                pia[i] = a_perm[t[i]-1];
                pib[i] = b_perm[k[i]-1];
            }

            ptr = proxy->GetBuffer1();
            //(t◯a_perm)◯t = pia◯t and (k◯b_perm)◯m = pib◯m
            for(int i = 0; i < size; i++){
                a_perm[i] = t[pia[i]-1];
                b_perm[i] = m[pib[i]-1];
                AddValueToCharArray(a_perm[i], &ptr,bsz);
                AddValueToCharArray(b_perm[i], &ptr,bsz);
            }
            Send( proxy->GetSocketHelper(), proxy->GetBuffer1(), size * bsz * 2);

        }
        else{
            uint64_t *randoms_j = new uint64_t[size];
            for (int j = 0; j < size; ++j) {
                randoms_j[j] = proxy->GenerateCommonRandom2();
            }
            //Create permutation elements(shares)
            //pia = e◯f and pib= h◯j
            pib = GetRandomPermutation(randoms_j, size);    //generate pib perm share j
            delete [] randoms_j;

            ptr2 = proxy->GetBuffer2();
            for (int i = 0; i < size; i++) {
                a_perm[i] = ConvertToLong(&ptr2, bsz); // receive pia permutation shares f
            }

            uint64_t *k_inv = new uint64_t[size];
            //pia = (t^-1)◯a_perm and b_perm = (k^-1)◯pib
            for(int i = 0; i < size; i++){
                pia[t[i]-1] = a_perm[i];
                b_perm[m[i]-1] = pib[i];
                k_inv[k[i]-1] = i+1;
            }
            ptr = proxy->GetBuffer1();
            //(t^-1)◯a_perm◯(k^-1) = pia◯(k^-1)
            for(int i = 0; i < size; i++){
                a_perm[i] = k_inv[pia[i]-1];
                AddValueToCharArray(a_perm[i], &ptr,bsz);
                AddValueToCharArray(b_perm[i], &ptr,bsz);
            }
            Send( proxy->GetSocketHelper(), proxy->GetBuffer1(), size * bsz * 2);
            delete [] k_inv;
        }
        Receive( proxy->GetSocketHelper(), proxy->GetBuffer2(),size * bsz); //
        ptr2 = proxy->GetBuffer2();
        for (int i = 0; i < size; i++) {
            pib[t[i]-1] = ConvertToLong(&ptr2, bsz); //receive and get rid of t
        }
        delete [] pia;
        delete [] a_perm;
        delete [] b_perm;
        delete [] pi;
        delete [] t;
        delete [] k;
        delete [] m;
        return pib;


    }else{ //helper
        uint64_t *pia = new uint64_t[size];
        uint64_t *pib = new uint64_t[size];
        uint64_t *randoms_e = new uint64_t[size];
        uint64_t *randoms_j = new uint64_t[size];
        uint64_t *f = new uint64_t[size];
        uint64_t *h = new uint64_t[size];
        unsigned char *ptr1 = proxy->GetBuffer1();
        unsigned char *ptr2 = proxy->GetBuffer2();
        thread thr1 = thread(Receive,proxy->GetSocketP1(), proxy->GetBuffer1(), size * bsz*2);
        thread thr2 = thread(Receive,proxy->GetSocketP2(), proxy->GetBuffer2(), size * bsz*2);
        thr1.join();
        thr2.join();
        for (int i = 0; i < size; i++) {
            pia[i] = ConvertToLong(&ptr1, bsz);
            pib[i] = ConvertToLong(&ptr1,bsz);
            pia[i] = (pia[i]+ConvertToLong(&ptr2,bsz))&mask;     //pia reconstructed
            pib[i] = (pib[i]+ConvertToLong(&ptr2,bsz))&mask;     //pib reconstructed
        }
        for (int j = 0; j < size; ++j) {
            randoms_e[j] = proxy->GenerateCommonRandom();
            randoms_j[j] = proxy->GenerateCommonRandom2();
        }

        //Create permutation elements(shares)
        //pia = e◯f and pib= h◯j
        uint64_t* e = GetRandomPermutation(randoms_e, size);
        uint64_t* j = GetRandomPermutation(randoms_j, size);
        for (int i = 0; i < size; ++i) {
            randoms_j[j[i]-1]= i+1;
        }
        for (int i = 0; i < size; ++i) {
            f[e[i]-1]= pia[i];
            h[i]= randoms_j[pib[i]-1];
        }

        ptr1 = proxy->GetBuffer1();
        ptr2 = proxy->GetBuffer2();
        for (int i = 0; i < size; i++) {
            AddValueToCharArray(h[i], &ptr1, bsz);
            AddValueToCharArray(f[i], &ptr2,bsz);
        }

        thr1 = thread(Send, proxy->GetSocketP1(), proxy->GetBuffer1(), size * bsz);
        thr2 = thread( Send, proxy->GetSocketP2(), proxy->GetBuffer2(), size * bsz);
        thr1.join();
        thr2.join();

        thr1 = thread(Receive,proxy->GetSocketP1(), proxy->GetBuffer1(), size * bsz*2);
        thr2 = thread(Receive,proxy->GetSocketP2(), proxy->GetBuffer2(), size * bsz*2);
        thr1.join();
        thr2.join();

        ptr1 = proxy->GetBuffer1();
        ptr2 = proxy->GetBuffer2();
        for (int i = 0; i < size; i++) {
            e[i] = ConvertToLong(&ptr1, bsz);   //ta0t
            h[i] = ConvertToLong(&ptr1,bsz);     //kb0k
            f[i] = ConvertToLong(&ptr2,bsz);     //(t^-1)a1(k^-1)
            j[i] = ConvertToLong(&ptr2,bsz);     //(k^-1)b1
        }

        for (int i = 0; i < size; i++) {
            pia[i] = f[e[i]-1];   //t(a0)t(t^-1)(a1)(k^-1)
            pib[i] = j[h[i]-1];     //kb0k(k^-1)b1
        }

        for (int i = 0; i < size; i++) {
            e[i] = pib[pia[i]-1];   //tab
        }

        ptr1 = proxy->GetBuffer1();
        ptr2 = proxy->GetBuffer2();
        for (int i = 0; i < size; i++) {
            auto tempShare = proxy->GenerateRandom()&mask;
            AddValueToCharArray(tempShare, &ptr1, bsz);
            AddValueToCharArray((e[i]-tempShare)&mask, &ptr2,bsz);
        }

        Send( proxy->GetSocketP1(), proxy->GetBuffer1(), size * bsz);
        Send( proxy->GetSocketP2(), proxy->GetBuffer2(), size * bsz);

        delete [] pia;
        delete [] e;
        delete [] f;
        delete [] h;
        delete [] j;
        delete [] pib;
        return nullptr;
    }
    return nullptr;
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

    if (proxy->GetPRole() == proxy1 || proxy->GetPRole() == proxy2){
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
        return NULL;
    }
    return NULL;
}


/**ApplyPermutation for Sorting (for XOR shared data)
 * @param p sorting permutation
 * @param v data shares (XOR shared)
 * @param pi common permutation
 * @param size number of elements in data
 * */
[[maybe_unused]] uint64_t *ApplyPermutationBoolean(Party *const proxy, const uint64_t *const p, uint64_t *const v, const uint64_t *const pi, uint32_t size) {
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
    if (proxy->GetPRole() == proxy1) {
        for (int i = 0; i < size; i++) {
            r[i] = proxy->GenerateRandom() & n;   //need to be smaller than n
            v[i] ^= r[i];
        }
    }

    if (proxy->GetPRole() == proxy1 || proxy->GetPRole() == proxy2) {
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
        thread thr2 = thread( Receive, proxy->GetSocketHelper(),proxy->GetBuffer2(),size * 8); //receives a share from pv_inv
        thr1.join();
        thr2.join();
        if(proxy->GetPRole() == proxy1) {
            for (int i = 0; i < size; i++) {
                r[i] = (uint64_t) MultMod((long long) r[i], x, n);      // r' = r*y modn
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
            pv_inv[i] = ConvertToLong(&ptr);   //they got the pvinv shares but P2 needs to eliminate the effect of r
        }

        if(proxy->GetPRole() == proxy2){
            Receive(proxy->GetSocketHelper(),proxy->GetBuffer1(),size*8);
            ptr = proxy->GetBuffer1();
            for (int i = 0; i < size; i++) {
                pr_inv[i] = ConvertToLong(&ptr);
            }

            long long xinv = GetModularInverseN(x, n);
            for (int i = 0; i < size; i++) {
                pr_inv[i] = (uint64_t) MultMod((long long) pr_inv[i], xinv, n);
                pv_inv[i] = pv_inv[i]^pr_inv[i]; //P2 eliminates r effect
            }
        }
    }
    else { // helper

        unsigned char *ptr1 = proxy->GetBuffer1();
        unsigned char *ptr2 = proxy->GetBuffer2();

        thread thr1 = thread(Receive,proxy->GetSocketP1(), proxy->GetBuffer1(), size * 16);
        thread thr2 = thread(Receive,proxy->GetSocketP2(), proxy->GetBuffer2(), size * 16);
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
            AddValueToCharArray(pv_inv[i] ^ tempShare, &ptr2);

        }

        thr1 = thread(Send, proxy->GetSocketP1(), proxy->GetBuffer1(), size * 8);
        thr2 = thread( Send, proxy->GetSocketP2(), proxy->GetBuffer2(), size * 8);
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

uint64_t **ApplyPermutationVectorized(
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

    if (proxy->GetPRole() == proxy1) {
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
                pv_inv[i][j] = ConvertToLong(&ptr);   //they got the pvinv shares but P1 needs to eliminate the effect of r
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
                pv_inv[i][j] = ConvertToLong(&ptr);   //they got the pvinv shares but P1 needs to eliminate the effect of r
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

// Sort algo that works with MostSignificantBit on 64 bit integers
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
    return NULL;
}

/** Sorting algorithm that uses XOR share for bit decomposition
 ** This algo works with ApplyPermutationNarrow2 and has some errors currently in SortNarrow
 ***/
uint64_t *SortNarrow(Party *const proxy, const uint64_t *const a, uint32_t size, uint32_t ringbits) {  //size = size of array
    int LT= 64;
    int bsz = ceil(size/8.0);
    if (proxy->GetPRole() == helper) {
        ArithmeticToXor(proxy, 0, size);
        XorToArithmetic3(proxy, 0, size);
        GeneratePermutation(proxy, 0, size, ringbits);
        Arithmetic2Permutation(proxy, 0, size, ringbits);
        for(int i = 1; i < LT ; ++i) {
            ApplyPermutationNarrow2(proxy, 0, 0, 0, size, ringbits);
            XorToArithmetic3(proxy, 0,  size);
            GeneratePermutation(proxy, 0, size, ringbits);
            ComposePermutations1Round(proxy, 0, 0, size, ringbits);
        }
        return NULL;
    }
    else {  //P1 or P2
        auto start1 = chrono::high_resolution_clock::now();
        auto *randoms = new uint64_t[size];
        auto *dc = new uint8_t[bsz];       //keeps the bits at the current(i) index
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
            AddBitToCharArray(((a_xor[j])&0x1), &ptr, &bit_index);
        }
        auto dca = XorToArithmetic3(proxy, dc, size);
        auto permG = GeneratePermutation(proxy, dca, size, ringbits);
        permG = Arithmetic2Permutation(proxy, permG, size, ringbits);   //Convert arithmetic share to permutation components
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
                AddBitToCharArray((tmp2[j])&0x1, &ptr, &bit_index);
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
            permG = ComposePermutations1Round(proxy, permG, permC, size, ringbits);
            end = chrono::high_resolution_clock::now();
            comp_time +=
                    chrono::duration_cast<chrono::nanoseconds>(end - start).count()*1e-9;
        }

        delete[] randoms;
        delete[] tmp;
        delete[] dc;
        //delete[] dn;
        auto end1 = chrono::high_resolution_clock::now();
        t_time +=
                chrono::duration_cast<chrono::nanoseconds>(end1 - start1).count()*1e-9;
        return permG;

    }
    return NULL;
}

/**  Vectorized Sorting Algorithm
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
            MostSignificantBit(proxy, 0, size);
            GeneratePermutation(proxy, 0, size);
            ApplyPermutationVectorized(proxy, 0, 0, 0, size, categories);
        }
        return 0;
    }
    else {  //P1 or P2
        double tt1, tt2, tt3, tt4, tt5=0;
        auto *randoms = new uint64_t[size];
        auto** res = new uint64_t*[categories];
        auto* to_shift = new uint64_t[size];
        //Get the

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

            res = ApplyPermutationVectorized(proxy, perm, res, pi, size, categories);  //apply the permutation to the current array
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