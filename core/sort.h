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
double perm_f = 0.0;
double perm_s = 0.0;
double perm_p = 0.0;

uint32_t* GenerateF(Party *const proxy, uint32_t *a, uint32_t size, uint32_t mask) {
    if(proxy->GetPRole() ==proxy1) {
        auto *result = new uint32_t[2 * size];
        memcpy(result+size,a,size*4);
        for(int i = 0; i < size; i++) {
            result[i]= (1 - a[i])&mask;
        }
        return result;
    }
    else if(proxy->GetPRole() == proxy2) {
        auto *result = new uint32_t[2 * size];
        memcpy(result+size,a,size*4);
        for(int i = 0; i < size; i++) {
            result[i]= (0 - a[i])&mask;
        }
        return result;
    }
    else return nullptr;
}

uint32_t* GenerateS(const uint32_t *const a, uint32_t size, uint32_t mask) {
    auto *result = new uint32_t[size];
    result[0] = a[0];
    for(int i = 1; i < size; i++) {
        result[i] = (result[i-1] + a[i])&mask;
    }
    return result;
}

uint32_t *GeneratePermutation(Party *const proxy, uint32_t *x, uint32_t size, uint32_t ringbits) {
    uint64_t mask =(1<< ringbits) -1;
    if(proxy->GetPRole() ==proxy1 || proxy->GetPRole() == proxy2) {
        uint32_t* x_ptr = x;
        auto start = chrono::high_resolution_clock::now();
        uint32_t *f = GenerateF(proxy, x, size, mask);
        auto end = chrono::high_resolution_clock::now();
        mul_time +=
                chrono::duration_cast<chrono::nanoseconds>(end - start).count()*1e-9;
        uint32_t *s = GenerateS(f, size * 2, mask);
        uint32_t *p = MultiplexNarrow(proxy,s,s+size,f+size,size,ringbits);

        delete[] f;
        delete[] s;
        return p;
    }
    else{   //helper
        MultiplexNarrow(proxy,nullptr,nullptr,nullptr,size,ringbits);
	return nullptr;
    }
}

void GetRandomPermutation(uint32_t * to_permute, uint32_t ind1, uint32_t ind2, uint32_t size, uint32_t ringbits= 20){
    for(int i = 0; i < size; i++)  to_permute[i] = i+1;
    auto mask = (1<<(ringbits-1))-1;
    ind1 = ind1 & mask;
    ind2 = ind2 & mask;
    for(int i = 0; i < size; i++){
        auto value = to_permute[ind1];
        to_permute[ind1] = to_permute[i];
        to_permute[i] = value;
        ind1 = (ind1 + ind2) & mask;
    }
}

uint32_t *Arithmetic2Permutation(Party *const proxy, const uint32_t *const p, uint32_t size, uint32_t ringbits) {
    uint64_t mask =(1<< ringbits) -1;
    int bsz = (ringbits)/8 + 1;
    uint32_t*  pip= new std::uint32_t[size];
    uint32_t*  perm= new std::uint32_t[size];

    if(proxy->GetPRole() == proxy1 || proxy->GetPRole() == proxy2) {
        uint64_t rnd = proxy->GenerateCommonRandom();
        uint64_t rnd1 = rnd >> 32;
        uint64_t rnd2 = rnd &  0xffffffff;
        uint32_t* pi = new uint32_t[size];
        GetRandomPermutation(pi, rnd1, rnd2, size, ringbits);

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
        delete[] pi;
        return perm;
    }
    else{   //helper

        unsigned char *ptr1 = proxy->GetBuffer1();
        unsigned char *ptr2 = proxy->GetBuffer2();
        uint64_t rnd = proxy->GenerateRandom();
        uint64_t rnd1 = rnd >> 32;
        uint64_t rnd2 = rnd &  0xffffffff;
        thread thr1 = thread(Receive,proxy->GetSocketP1(), proxy->GetBuffer1(), size * bsz);
        thread thr2 = thread(Receive,proxy->GetSocketP2(), proxy->GetBuffer2(), size * bsz);
        uint32_t *e = new uint32_t[size];
        thread thr3 = thread(GetRandomPermutation, e, rnd1, rnd2, size, ringbits);
        thr1.join();
        thr2.join();
        thr3.join();
        for (int i = 0; i < size; i++) {
            pip[i] = ConvertToLong(&ptr1, bsz);
            pip[i] = (pip[i]+ConvertToLong(&ptr2, bsz))&mask;      //reconstruct pip
        }
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


uint8_t *ApplyPermutation(
        Party *const proxy,
        const uint32_t *const p,
        uint8_t *const v,
        const uint32_t *const pi,
        uint32_t size,
        uint64_t ringbits
) {
    auto bsz = ceil((ringbits)/8.0);
    // construct permuted p and v shares
    auto* pip = new uint32_t[size];
    auto* pir = new uint8_t[size];
    auto* pr_inv = new uint8_t[size];

    if (proxy->GetPRole() == proxy1) {
        uint8_t map[128];
        auto* r = new uint8_t[size];
        for (int i=0;i<128;i++){
            map[i] = proxy->GenerateCommonRandomByte()&0x1;
        }
        for (int i = 0; i < size; i++) {
            uint8_t tmp = proxy->GenerateRandomByte();
            r[i] = tmp&0x1;
            v[i] = (v[i] + r[i])&0x1;
            int index = tmp>>1;   //need to be smaller than nu
            while (map[index] != r[i]) {
                index = (index+1) & 0x7f;
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
        uint32_t *pip1 = new uint32_t[size];
        uint32_t *pip2 = new uint32_t[size];
        unsigned char *ptr1 = proxy->GetBuffer1();
        unsigned char *ptr2 = proxy->GetBuffer2();
        for (int i = 0; i < size; i++) {
            pip1[i] = ConvertToLong(&ptr1, bsz);
            pir[i] = ConvertToUint8(&ptr1);
            pip2[i] = ConvertToLong(&ptr2,bsz);
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

uint32_t *ComposePermutations(
        Party *const proxy,
        const uint32_t *const a,
        const uint32_t *const b,
        const uint32_t *t,
        uint32_t size,
        uint32_t ringbits
) {
    auto bsz = ceil((ringbits)/8.0);
    auto mask = (1<< (ringbits))-1;

    if (proxy->GetPRole() == proxy1 || proxy->GetPRole() == proxy2){
        uint32_t *randoms_t = new uint32_t[size];

        uint32_t *pia = new uint32_t[size];
        uint32_t *a_perm = new uint32_t[size];
        uint32_t *b_perm = new uint32_t[size];

        unsigned char *ptr2 = proxy->GetBuffer2();
        unsigned char *ptr = proxy->GetBuffer1();

        if (proxy->GetPRole() == proxy1){
            std::thread thr1([pia,a,t,size](){
                for(int i = 0; i < size; i++){
                    pia[i] = a[t[i]-1];
                }
            });
            std::thread thr2([b_perm,b,t,size](){
                for(int i = 0; i < size; i++){
                    b_perm[i] = b[t[i]-1];
                }
            });
            thr1.join();
            thr2.join();

            ptr = proxy->GetBuffer1();
            //(t◯a_perm)◯t = pia◯t
            for(int i = 0; i < size; i++){
                a_perm[i] = t[pia[i]-1];
                AddValueToCharArray(a_perm[i], &ptr,bsz);
                AddValueToCharArray(b_perm[i], &ptr,bsz);
            }
            Send(proxy->GetSocketHelper(), proxy->GetBuffer1(), size * bsz * 2);
            Receive(proxy->GetSocketHelper(), proxy->GetBuffer2(),size * bsz); //
            ptr2 = proxy->GetBuffer2();
            for (int i = 0; i < size; i++) {
                pia[t[i]-1] = ConvertToLong(&ptr2, bsz); //receive and get rid of t
            }
        }
        else{
            uint64_t *k_inv = new uint64_t[size];
            std::thread thr1([pia,a,t,size](){
                for(int i = 0; i < size; i++){
                    pia[t[i]-1] = a[i];
                }
            });
            std::thread thr2([b_perm,b,t,size](){
                for(int i = 0; i < size; i++){
                    b_perm[i] = b[t[i]-1];
                }
            });
            std::thread thr3([k_inv,t,size](){
                for(int i = 0; i < size; i++){
                    k_inv[t[i]-1] = i+1;
                }
            });
            thr1.join();
            thr2.join();
            thr3.join();

            ptr = proxy->GetBuffer1();
            for(int i = 0; i < size; i++){
                a_perm[i] = k_inv[pia[i]-1];
                AddValueToCharArray(a_perm[i], &ptr,bsz);
                AddValueToCharArray(b_perm[i], &ptr,bsz);
            }
            Send(proxy->GetSocketHelper(), proxy->GetBuffer1(), size * bsz * 2);
            Receive(proxy->GetSocketHelper(), proxy->GetBuffer2(),size * bsz); //
            ptr2 = proxy->GetBuffer2();
            for (int i = 0; i < size; i++) {
                pia[i] = ConvertToLong(&ptr2, bsz); //receive and get rid of t
            }
            delete [] k_inv;
        }

        delete [] a_perm;
        delete [] b_perm;
        return pia;


    }else{ //helper
        uint32_t *pia = new uint32_t[size];
        uint32_t *e = new uint32_t[size];
        uint32_t *f = new uint32_t[size];
        uint32_t *h = new uint32_t[size];
        unsigned char *ptr1 = proxy->GetBuffer1();
        unsigned char *ptr2 = proxy->GetBuffer2();

        //Create permutation elements(shares)
        //pia = e◯f and pib= h◯j
        uint32_t *e_ = new uint32_t[size];
        uint64_t rnd = proxy->GenerateRandom();
        uint32_t rnd2 = rnd&0xffffffff;
        uint32_t rnd1 = rnd>>32;


        thread thr1 = thread(Receive,proxy->GetSocketP1(), proxy->GetBuffer1(), size * bsz*2);
        thread thr2 = thread(Receive,proxy->GetSocketP2(), proxy->GetBuffer2(), size * bsz*2);
        thread thr3 = thread(GetRandomPermutation, e_, rnd1, rnd2, size, ringbits);
        thr1.join();
        thr2.join();
        thr3.join();
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
            f[e_[i]-1]=  h[pia[i]-1];
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
}

/** Sorting algorithm that uses XOR share for bit decomposition
 ** This algo works with ApplyPermutation and has some errors currently in Sort
 ***/
uint32_t *Sort(Party *const proxy, const uint64_t *const a, uint32_t size, uint32_t ringbits) {  //size = size of array
    int LT= 64;
    int bsz = (size/8.0)+1;
    if (proxy->GetPRole() == helper) {
        ArithmeticToXor(proxy, 0, size);
        XorToArithmetic3(proxy, 0, size, ringbits);
        GeneratePermutation(proxy, 0, size, ringbits);
        Arithmetic2Permutation(proxy, 0, size, ringbits);
        for(int i = 1; i < LT ; ++i) {
            ApplyPermutation(proxy, 0, 0, 0, size, ringbits);
            XorToArithmetic3(proxy, 0,  size, ringbits);
            GeneratePermutation(proxy, 0, size, ringbits);
            ComposePermutations(proxy, 0, 0, 0, size, ringbits);
        }
        return NULL;
    }
    else {  //P1 or P2
        auto start1 = chrono::high_resolution_clock::now();
        auto *dc = new uint8_t[bsz];       //keeps the bits at the current(i) index
        auto *tmp = new uint8_t[size];

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
        auto dca = XorToArithmetic3(proxy, dc, size, ringbits);
        auto permG = GeneratePermutation(proxy, dca, size, ringbits);

        permG = Arithmetic2Permutation(proxy, permG, size, ringbits);//Convert arithmetic share to permutation components


        for(int i = 1; i < LT; ++i) {
            uint32_t* pi = new uint32_t[size];
            uint32_t* pi1 = new uint32_t[size];
            thread thr1 = thread(GetRandomPermutation, pi, proxy->GenerateCommonRandom(), proxy->GenerateCommonRandom(),
                                 size, ringbits);
            thread thr2 = thread(GetRandomPermutation, pi1, proxy->GenerateCommonRandom(),
                                 proxy->GenerateCommonRandom(), size, ringbits);
            for (int j = 0; j < size; ++j) {
                tmp[j] = (a_xor[j]>>i)&0x1;
            }
            thr1.join();
            thr2.join();

            start = chrono::high_resolution_clock::now();
            auto tmp2 = ApplyPermutation(proxy, permG, tmp, pi, size, ringbits);
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
            dca = XorToArithmetic3(proxy, dc, size, ringbits);
            end = chrono::high_resolution_clock::now();
            x2a_time +=
                    chrono::duration_cast<chrono::nanoseconds>(end - start).count()*1e-9;


            start = chrono::high_resolution_clock::now();
            auto permC = GeneratePermutation(proxy, dca, size, ringbits);
            end = chrono::high_resolution_clock::now();
            gp_time +=
                    chrono::duration_cast<chrono::nanoseconds>(end - start).count()*1e-9;


            start = chrono::high_resolution_clock::now();
            permG = ComposePermutations(proxy, permG, permC, pi1, size, ringbits);
            end = chrono::high_resolution_clock::now();
            comp_time +=
                    chrono::duration_cast<chrono::nanoseconds>(end - start).count()*1e-9;
            delete[] pi;
            delete[] pi1;
            }

        delete[] tmp;
        delete[] dc;
        auto end1 = chrono::high_resolution_clock::now();
        t_time +=
                chrono::duration_cast<chrono::nanoseconds>(end1 - start1).count()*1e-9;
        return permG;

    }
    return NULL;
}

uint32_t *ReconstructNarrowPermutation(Party *const proxy, const uint32_t *const a, uint32_t sz, uint64_t ringbits) {
    auto mask = (1<< ringbits)-1;
    auto bsz = (uint32_t)ceil(ringbits/8.0);
    uint32_t *b = new uint32_t[sz];
    uint32_t *ab = new uint32_t[sz];
    if ( proxy->GetPRole() == proxy1 ) {
        unsigned char *ptr = proxy->GetBuffer1();
        for (int i = 0; i < sz; i++) {
            AddValueToCharArray(a[i], &ptr,bsz);
        }
        thread thr1 = thread(Send,proxy->GetSocketP2(), proxy->GetBuffer1(), sz*bsz);
        thread thr2 = thread(Receive,proxy->GetSocketP2(), proxy->GetBuffer2(), sz*bsz);
        thr1.join();
        thr2.join();

        ptr = proxy->GetBuffer2();
        for (int i = 0; i < sz; i++) {
            b[i]=ConvertToLong(&ptr,bsz);
        }
        for (int i = 0; i < sz; i++) {
            ab[i] = (b[a[i]-1]);
        }

    } else if ( proxy->GetPRole() == proxy2) {
        unsigned char *ptr = proxy->GetBuffer1();
        for (int i = 0; i < sz; i++) {
            AddValueToCharArray(a[i], &ptr,bsz);
        }
        thread thr1 = thread(Send,proxy->GetSocketP1(), proxy->GetBuffer1(), sz*bsz);
        thread thr2 = thread(Receive,proxy->GetSocketP1(), proxy->GetBuffer2(), sz*bsz);
        thr2.join();
        thr1.join();
        ptr = proxy->GetBuffer2();
        for (int i = 0; i < sz; i++) {
            b[i]=ConvertToLong(&ptr,bsz);
        }
        for (int i = 0; i < sz; i++) {
            ab[i] = (a[b[i]-1]) ;
        }
    }
    delete [] b;
    return ab;
}
