//
// Created by Mete Akgun on 03.07.20.
//

#ifndef PML_PARTY_H
#define PML_PARTY_H


#include <stdio.h>
#include <string.h> //strlen
#include <errno.h>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <cmath>
using namespace std;

#include "../utils/constant.h"
#include "../utils/flib.h"
#include "../utils/connection.h"

class Party {
public:

    Party(role r, uint16_t hport=7777, const string hip="127.0.0.1", uint16_t cport=8888,const string cip="127.0.0.1") {
        p_role = r;

        if (p_role == P1) {
            socket_helper = connect2helper(cip, hip, cport, hport, r);
            socket_p2 = open_P1(cip, cport);
        } else if (p_role == P2) {
            socket_helper = connect2helper(cip, hip, cport, hport + 1, r);
            socket_p1 = connect2P1(cip, cport);
        } else if (p_role == HELPER) {
            int client_socket[2];
            open_helper(hip, hport, hport + 1, client_socket);

            socket_p1 = client_socket[0];
            socket_p2 = client_socket[1];
        }
        SetCommonSeed();

        // compute the number of bits for exponential
        // for positive power
        bool flag = false;
        int used_bits = 0;
        // for negative power
        bool neg_flag = false;
        int zero_bits = 0;
        for(int i = 7; i >= 0; i--) {
            // positive power
            int n_bits_by_i = ceil(log2(exp(pow(2, i))));
            if(L_BIT - 1 - 2 * FRAC - n_bits_by_i - used_bits >= 0) {
                used_bits += n_bits_by_i;
                max_power += pow(2, i);
            }
            if(!flag && L_BIT - 2 - 2 * FRAC - n_bits_by_i >= 0) {
                flag = true;
                n_bits += i + 1;
            }

            // negative power
            int n_zero_bits_by_i = abs(floor(log2(exp(pow(2, i))))); // how many bits are zero after the comma
            if(zero_bits + n_zero_bits_by_i < FRAC - 1) {
                zero_bits += n_zero_bits_by_i;
                min_power -= pow(2, i);

                if(!neg_flag) {
                    neg_flag = true;
                    neg_n_bits += i + 1;
                }
            }
        }
    }

    ~Party() {
        if (p_role == P1) {
            close(socket_helper);
            close(socket_p2);
        } else if (p_role == P2) {
            close(socket_helper);
            close(socket_p1);
        } else if (p_role == HELPER) {
            close(socket_p1);
            close(socket_p2);
        }
    }

    void SetCommonSeed() {
        if (p_role == P1) {
            srand(time(NULL) + 10000);
            common_seed = rand();
            unsigned char *ptr = &buffer[0];
            addVal2CharArray((uint64_t) common_seed, &ptr);
            Send(socket_p2, buffer, 8);
        } else if (p_role == P2) {
            srand(time(NULL) + 20000);
            Receive(socket_p1, buffer, 8);
            unsigned char *ptr = &buffer[0];
            common_seed = convert2Long(&ptr);
        } else if (p_role == HELPER)
            srand(time(NULL) + 30000);
        seed = rand();
    }

    uint64_t generateRandom() {
        srand(seed);
        uint64_t val = 0;
        for (int i = 3; i >= 0; i -= 1) {
            uint64_t a = rand() & 0xffff;
            val = val ^ (a << (i * 16));
        }
        seed += 4;
        return val;
    }

    uint64_t generateCommonRandom() {
        srand(common_seed);
        uint64_t val = 0;
        for (int i = 3; i >= 0; i -= 1) {
            uint64_t a = rand() & 0xffff;
            val = val ^ (a << (i * 16));
        }
        common_seed += 4;
        return val;
    }

    uint64_t createShare(double val){
        uint64_t v = convert2uint64(val);
        uint64_t share;
        if (p_role ==P1) {
            share = generateCommonRandom();
        }
        else{
            share = v - generateCommonRandom();
        }
        return share;
    }

    uint64_t* createShare(double *val, uint32_t sz){
        uint64_t *v = convert2uint64(val,sz);
        uint64_t *share = new uint64_t[sz];
        for (uint32_t i=0;i<sz;i++){
            if (p_role ==P1) {
                share[i] = generateCommonRandom();
            }
            else{
                share[i] = v[i] - generateCommonRandom();
            }
        }

        return share;
    }

    uint64_t** createShare(double **val, uint32_t n_row, uint32_t n_col){
        uint64_t **v = convert2uint64(val, n_row, n_col);
        uint64_t **share = new uint64_t*[n_row];
        for (uint32_t i = 0; i < n_row; i++){
            share[i] = new uint64_t[n_col];
            for(uint32_t j = 0; j < n_col; j++) {
                if (p_role ==P1) {
                    share[i][j] = generateCommonRandom();
                }
                else{
                    share[i][j] = v[i][j] - generateCommonRandom();
                }
            }
        }
        return share;
    }

    int ReadByte() {
        if (p_role == HELPER) {
            Receive(socket_p1, buffer, 1);
            return (int) buffer[0];
        } else
            return -1;
    }

    uint32_t ReadInt() {
        if (p_role == HELPER) {
            Receive(socket_p1, buffer, 4);
            unsigned char *ptr = &buffer[0];
            return convert2Int(&ptr);
        } else
            return 0;
    }

//    void SendBytes(op o, uint32_t sz = 0,int L1 = 0) {
//        if (p_role == P1) {
//            unsigned char *ptr = &buffer[0];
//            size_t s = 1;
//            addVal2CharArray((uint8_t) o, &ptr);
//            if (sz != 0) {
//                addVal2CharArray((uint32_t) sz, &ptr);
//                s += 4;
//            }
//            if (L1 != 0) {
//                addVal2CharArray((uint8_t) L1, &ptr);
//                s += 1;
//            }
//            Send(socket_helper, buffer, s);
//        }
//    }

    void SendBytes(op o, uint32_t sz = 0, uint32_t sz2 = 0) {
        if (p_role == P1) {
            unsigned char *ptr = &buffer[0];
            size_t s = 1;
            addVal2CharArray((uint8_t) o, &ptr);
            if (sz != 0) {
                addVal2CharArray((uint32_t) sz, &ptr);
                s += 4;
            }
            if (sz2 != 0) {
                addVal2CharArray((uint32_t) sz2, &ptr);
                s += 4;
            }
            Send(socket_helper, buffer, s);
        }
    }

    void PrintBytes() {
        PBytes();
    }

    role getPRole() const {
        return p_role;
    }

    int getSocketP1() const {
        return socket_p1;
    }

    int getSocketP2() const {
        return socket_p2;
    }

    int getSocketHelper() const {
        return socket_helper;
    }

    uint8_t *getBuffer1(){
        return buffer;
    }

    uint8_t *getBuffer2(){
        return buffer2;
    }

    int getNBits() {
        return n_bits;
    }

    void setNBits(int nBits) {
        n_bits = nBits;
    }

    double getMaxPower() const {
        return max_power;
    }

    void setMaxPower(double maxPower) {
        max_power = maxPower;
    }

    int getNegNBits() const {
        return neg_n_bits;
    }

    void setNegNBits(int negNBits) {
        neg_n_bits = negNBits;
    }

    double getMinPower() const {
        return min_power;
    }

    void setMinPower(double minPower) {
        min_power = minPower;
    }
private:
    role p_role;
    int socket_p1,socket_p2,socket_helper;
    uint8_t buffer[BUFFER_SIZE];
    uint8_t buffer2[BUFFER_SIZE];
    uint32_t common_seed;
    uint32_t seed;
    int n_bits = FRAC; // number of bits of a value to consider in the exponential computation
    int neg_n_bits = FRAC; // number of bits of a negative value to consider in the exponential computation
    double max_power = 0;
    double min_power = 0;
};



#endif //PML_PARTY_H
