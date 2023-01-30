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
#include "../utils/AES_CTR_RBG.h"
#include <cryptopp/secblock.h>
using CryptoPP::FixedSizeSecBlock;
#include <cryptopp/osrng.h>
using CryptoPP::OS_GenerateRandomBlock;

class Party {
private:
    double kk = 0.0;
public:

    Party(role r, uint16_t hport=7777, const string hip="127.0.0.1", uint16_t cport=8888,const string cip="127.0.0.1") {
        p_role = r;
        if (p_role == P1) {
            connect2helper(hip, hport,socket_helper);
            open_P0(cip, cport,socket_p1);
        } else if (p_role == P2) {
            connect2helper(hip, hport + SCKNUM,socket_helper);
            connect2P0(cip, cport,socket_p0);
        } else if (p_role == HELPER) {
            open_helper(hip, hport, socket_p0,socket_p1);
        }
        initialiseRBG();

        // pre-compute the truncation mask for negative values based on FRAC
        neg_truncation_mask = (((uint64_t) 1 << FRAC) - 1) << (L_BIT - FRAC);

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
            closesocket(socket_helper);
            closesocket(socket_p1);
        } else if (p_role == P2) {
            closesocket(socket_helper);
            closesocket(socket_p0);
        } else if (p_role == HELPER) {
            closesocket(socket_p0);
            closesocket(socket_p1);
        }
    }

    void initialiseRBG() {
        if (p_role == P1) {
            OS_GenerateRandomBlock(false, buffer, 32);
            Send(&socket_p1[0], buffer, 32);
            common_rbg = new AES_CTR_RBG(buffer, 32);
            common_rbg->GenerateBlock(common_random_buffer, BUFFER_SIZE);
        } else if (p_role == P2) {
            Receive(&socket_p0[0], buffer, 32);
            unsigned char *ptr = &buffer[0];
            common_rbg = new AES_CTR_RBG(ptr, 32);
            common_rbg->GenerateBlock(common_random_buffer, BUFFER_SIZE);
        }
        rbg = new AES_CTR_RBG();
        rbg->GenerateBlock(random_buffer, BUFFER_SIZE);
    }


    uint64_t generateRandom() {
        if (used_random_bytes + 7 >= BUFFER_SIZE) {
            rbg->GenerateBlock(random_buffer, BUFFER_SIZE);
            used_random_bytes = 0;
        }
        uint64_t val = *(uint64_t *)(random_buffer + used_random_bytes);
        used_random_bytes += 8;
        return val;
    }
    uint8_t generateRandomByte() {
        if (used_random_bytes >= BUFFER_SIZE) {
            rbg->GenerateBlock(random_buffer, BUFFER_SIZE);
            used_random_bytes = 0;
        }
        uint8_t val = random_buffer[used_random_bytes];
        used_random_bytes++;
        return val;
    }

    uint64_t generateCommonRandom() {
        if (used_common_random_bytes + 7 >= BUFFER_SIZE) {
            common_rbg->GenerateBlock(common_random_buffer, BUFFER_SIZE);
            used_common_random_bytes = 0;
        }
        uint64_t val = *(uint64_t *)(common_random_buffer + used_common_random_bytes);
        used_common_random_bytes += 8;
        return val;
    }

    uint8_t generateCommonRandomByte() {
        if (used_common_random_bytes >= BUFFER_SIZE) {
            common_rbg->GenerateBlock(common_random_buffer, BUFFER_SIZE);
            used_common_random_bytes = 0;
        }
        uint8_t val = common_random_buffer[used_common_random_bytes];
        used_common_random_bytes++;
        return val;
    }

    uint64_t createShare(uint64_t val){
        uint64_t share;
        if (p_role ==P1) {
            share = generateCommonRandom();
        }
        else{
            share = val - generateCommonRandom();
        }
        return share;
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
        delete[] v;
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
        delete[] v;
        return share;
    }

    int ReadByte() {
        if (p_role == HELPER) {
            Rcv(socket_p0[0], buffer, 1);
            return (int) buffer[0];
        } else
            return -1;
    }

    uint32_t ReadInt() {
        if (p_role == HELPER) {
            Rcv(socket_p0[0], buffer, 4);
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

    /**
     * Sends at least the operational op to the helper.
     * Additional parameters for calling the according operation might be send.
     * @param o the operation to be performed, the operation is one of the building blocks defined as constants in constant.h
     * @param params additional parameters to be send.
     * By default Null, but if additional parameters are desired to send, use a vector holding all those parameters.
     * @param size the size of params. Default is 0 but if params is not NULL, size matches the number of parameters stored in params.
     */
    void SendBytes(op o, uint32_t *params = NULL, uint32_t size = 0) {
        // this will be called by both parties P1 and P2 but must only be executed for one of those,
        // otherwise, the helper will receive the same information twice (or in the worst case not read from P2 and mess up future results)
        if (p_role == P1) {
            unsigned char *ptr = &buffer[0];
            size_t s = 1;
            addVal2CharArray((uint8_t) o, &ptr);
            if (params != NULL && size > 0) {
                for(uint32_t i = 0; i<size; i++){
                    addVal2CharArray((uint32_t) params[i], &ptr);
                    s += 4; // one 32 bit value requires 2^4 bits; if we were to store 64 bit values --> += 8
                }
            }
            Snd(socket_helper[0], buffer, s);
        }
    }


    void PrintBytes() {
        PBytes();
    }

    void PrintPaperFriendly(double time_taken) {
        cout << "Paper\t" << (bytesSend / 1e6) << "\t" << (bytesReceived / 1e6) << "\t" << fixed << time_taken << setprecision(9) << endl;
    }

    role getPRole() const {
        return p_role;
    }

    int *getSocketP1(){
        return socket_p0;
    }

    int *getSocketP2(){
        return socket_p1;
    }

    int *getSocketHelper(){
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

    uint64_t getNegTruncationMask() const {
        return neg_truncation_mask;
    }
private:
    role p_role;
    int socket_p0[SCKNUM],socket_p1[SCKNUM],socket_helper[SCKNUM];
    uint8_t buffer[BUFFER_SIZE];
    uint8_t buffer2[BUFFER_SIZE];
    CryptoPP::byte random_buffer[BUFFER_SIZE];
    CryptoPP::byte common_random_buffer[BUFFER_SIZE];
    size_t used_random_bytes = 0;
    size_t used_common_random_bytes = 0;
    AES_CTR_RBG* common_rbg;
    AES_CTR_RBG* rbg;
    int n_bits = FRAC; // number of bits of a value to consider in the exponential computation
    int neg_n_bits = FRAC; // number of bits of a negative value to consider in the exponential computation
    double max_power = 0;
    double min_power = 0;
    uint64_t neg_truncation_mask = 0;
};



#endif //PML_PARTY_H
