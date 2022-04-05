//
// Created by Mete Akgun on 28.12.21.
//

#ifndef CORE_H
#define CORE_H

#include "Party.h"
#include "../utils/test_functions.h"
#include <thread>
#include <mutex>

uint64_t REC(Party* proxy, uint64_t a, uint64_t mask=RING_N) {

    uint64_t b;
    if ( proxy->getPRole() == P1) {
        unsigned char *ptr = proxy->getBuffer1();
        addVal2CharArray(a, &ptr);
        Send(proxy->getSocketP2(), proxy->getBuffer1(), 8);
        Receive(proxy->getSocketP2(), proxy->getBuffer1(), 8);
        ptr = proxy->getBuffer1();
        b = convert2Long(&ptr);

    } else if ( proxy->getPRole() == P2) {
        unsigned char *ptr = proxy->getBuffer1();
        addVal2CharArray(a, &ptr);
        Send(proxy->getSocketP1(), proxy->getBuffer1(), 8);
        Receive(proxy->getSocketP1(), proxy->getBuffer1(), 8);
        ptr = proxy->getBuffer1();
        b = convert2Long(&ptr);
    }
    return (a + b) & mask;
}

uint64_t *REC(Party* proxy, uint64_t *a, uint32_t sz, uint64_t mask=RING_N) {

    uint64_t *b = new uint64_t[sz];
    if ( proxy->getPRole() == P1 ) {
        unsigned char *ptr = proxy->getBuffer1();
        for (int i = 0; i < sz; i++) {
            addVal2CharArray(a[i], &ptr);
        }
        Send(proxy->getSocketP2(), proxy->getBuffer1(), sz * 8);
        Receive(proxy->getSocketP2(), proxy->getBuffer1(), sz * 8);
        ptr = proxy->getBuffer1();
        for (int i = 0; i < sz; i++) {
            b[i] = convert2Long(&ptr);
        }

    } else if ( proxy->getPRole() == P2) {
        unsigned char *ptr = proxy->getBuffer1();
        for (int i = 0; i < sz; i++) {
            addVal2CharArray(a[i], &ptr);
        }
        Send(proxy->getSocketP1(), proxy->getBuffer1(), sz * 8);
        Receive(proxy->getSocketP1(), proxy->getBuffer1(), sz * 8);
        ptr = proxy->getBuffer1();
        for (int i = 0; i < sz; i++) {
            b[i] = convert2Long(&ptr);
        }
    }
    for (int i = 0; i < sz; i++) {
        b[i] = (a[i] + b[i]) & mask;
    }
    return b;
}

uint64_t** REC(Party *proxy, uint64_t **a, uint32_t n_row, uint32_t n_col) {
    // Reconstruct a secret shared 2D array
    uint64_t **b = new uint64_t*[n_row];
    if (proxy->getPRole() == P1) {
        unsigned char *ptr = proxy->getBuffer1();
        for (int i = 0; i < n_row; i++) {
            b[i] = new uint64_t[n_col];
            for( int j = 0; j < n_col; j++) {
                addVal2CharArray(a[i][j], &ptr);
            }
        }
        Send(proxy->getSocketP2(), proxy->getBuffer1(), n_row * n_col * 8);
        Receive(proxy->getSocketP2(), proxy->getBuffer1(), n_row * n_col * 8);
        ptr = proxy->getBuffer1();
        for (int i = 0; i < n_row; i++) {
            for(int j = 0; j < n_col; j++) {
                b[i][j] = convert2Long(&ptr);
            }
        }

    } else if (proxy->getPRole() == P2) {
        unsigned char *ptr = proxy->getBuffer1();
        for (int i = 0; i < n_row; i++) {
            for( int j = 0; j < n_col; j++) {
                addVal2CharArray(a[i][j], &ptr);
            }
        }
        Send(proxy->getSocketP1(), proxy->getBuffer1(), n_row * n_col * 8);
        Receive(proxy->getSocketP1(), proxy->getBuffer1(), n_row * n_col * 8);
        ptr = proxy->getBuffer1();
        for (int i = 0; i < n_row; i++) {
            b[i] = new uint64_t[n_col];
            for( int j = 0; j < n_col; j++) {
                b[i][j] = convert2Long(&ptr);
            }
        }
    }
    for (int i = 0; i < n_row; i++) {
        for( int j = 0; j < n_col; j++) {
            b[i][j] = (a[i][j] + b[i][j]);
        }
    }
    return b;
}

uint64_t ADD(Party* proxy, uint64_t a, uint64_t b) {
    return a + b;
}

void GenerateMultiplicationTriple(Party* proxy, uint64_t **mt1, uint64_t **mt2, uint32_t size) {
    // Input(s)
    // mt1 and mt2: 3-by-size array whose rows will be a_i, b_i and c_i, respectively
    // size: the number of multiplication triples that will be generated
    // Return(s)
    // None

    srand(time(NULL));
    for (int i = 0; i < size; i++) {
        uint64_t tmp_a = proxy->generateRandom();
        uint64_t tmp_b = proxy->generateRandom();
        uint64_t tmp_c = tmp_a * tmp_b; // mod operation here?

        // a
        mt1[0][i] = proxy->generateRandom();
        mt2[0][i] = tmp_a - mt1[0][i];

        // b
        mt1[1][i] = proxy->generateRandom();
        mt2[1][i] = tmp_b - mt1[1][i];

        // c
        mt1[2][i] = proxy->generateRandom();
        mt2[2][i] = tmp_c - mt1[2][i];

        // cout << mt1[0][i] << " " << mt2[0][i] << " " <<  mt1[1][i] << " " << mt2[1][i] << " " << mt1[2][i] << " " << mt2[2][i] << endl;
    }
}

uint64_t MUX(Party* proxy, uint64_t x, uint64_t y, uint64_t b) {

    if ( proxy->getPRole() == P1) {
        uint64_t r1 = proxy->generateCommonRandom(), r2 = proxy->generateCommonRandom(), r3 = proxy->generateCommonRandom(), r4 = proxy->generateCommonRandom();
        //uint64_t s = proxy->generateRandom();
        //x = x + s ;
        //y = y + s;
        uint64_t m1 = (b * (x - y)) - (r3*(x-y)) - (r3*r4) - (r2*b);
        //m1 = m1 >> FRAC;
        uint64_t m2 = b + r1;
        uint64_t m3 = x - y + r4;
        unsigned char *ptr = proxy->getBuffer1();
        addVal2CharArray(m2, &ptr);
        addVal2CharArray(m3, &ptr);
        Send(proxy->getSocketHelper(), proxy->getBuffer1(), 16);
        Receive(proxy->getSocketHelper(), proxy->getBuffer1(), 8);
        ptr = proxy->getBuffer1();
        m1 = m1 + convert2Long(&ptr);
        m1 = m1 >> FRAC;
        x = x - m1;


        return x;
    } else if ( proxy->getPRole() == P2) {
        uint64_t r1 = proxy->generateCommonRandom(), r2 = proxy->generateCommonRandom(), r3 = proxy->generateCommonRandom(), r4 = proxy->generateCommonRandom();
        //uint64_t s = proxy->generateRandom();
        //x = x + s;
        //y = y + s;
        uint64_t m1 = (b * (x - y)) - (r1*(x-y)) - (r1*r2) - (r4*b);
        //m1 = -1 * ((-1 * m1) >> FRAC);
        uint64_t m2 = x - y + r2;
        uint64_t m3 = b + r3;
        unsigned char *ptr = proxy->getBuffer1();
        //addVal2CharArray(m1, &ptr);
        addVal2CharArray(m2, &ptr);
        addVal2CharArray(m3, &ptr);
        Send(proxy->getSocketHelper(), proxy->getBuffer1(), 16);
        Receive(proxy->getSocketHelper(), proxy->getBuffer1(), 8);
        ptr = proxy->getBuffer1();
        m1 = m1 + convert2Long(&ptr);
        m1 = -1 * ((-1 * m1) >> FRAC);
        x = x - m1 ;

        return x;
    } else if ( proxy->getPRole() == HELPER) {
        Receive(proxy->getSocketP1(), proxy->getBuffer1(), 16);
        Receive(proxy->getSocketP2(), proxy->getBuffer2(), 16);
        unsigned char *ptr = proxy->getBuffer1();
        unsigned char *ptr2 = proxy->getBuffer2();
        //uint64_t m1 = convert2Long(&ptr);
        uint64_t m2 = convert2Long(&ptr);
        uint64_t m3 = convert2Long(&ptr);
        //uint64_t m4 = convert2Long(&ptr2);
        uint64_t m5 = convert2Long(&ptr2);
        uint64_t m6 = convert2Long(&ptr2);

        uint64_t m = /*m1 + m4 -*/ (m2 * m5) + (m3 * m6);

        uint64_t m1 = proxy->generateRandom();
        m2 = m - m1;

        ptr = proxy->getBuffer1();
        addVal2CharArray(m1, &ptr);
        Send(proxy->getSocketP1(), proxy->getBuffer1(), 8);
        ptr2 = proxy->getBuffer2();
        addVal2CharArray(m2, &ptr2);
        Send(proxy->getSocketP2(), proxy->getBuffer2(), 8);
        return 0;
    }
    return -1;
}

uint64_t* MUX(Party* proxy, uint64_t *x, uint64_t *y, uint64_t *b, uint32_t sz) {
    if ( proxy->getPRole() == P1){
        unsigned char *ptr = proxy->getBuffer1();
        uint64_t *res = new uint64_t[sz];
        uint64_t *m1 = new uint64_t[sz];
        for (uint32_t i = 0; i < sz; i++) {
            uint64_t r1=proxy->generateCommonRandom(), r2=proxy->generateCommonRandom(), r3=proxy->generateCommonRandom(), r4=proxy->generateCommonRandom();

            m1[i] = (b[i] * (x[i] - y[i])) - (r2*b[i]) - (r3*(x[i] - y[i])) - (r3*r4);
            uint64_t m2 = b[i] + r1;
            uint64_t m3 = x[i] - y[i] + r4;

            addVal2CharArray(m2,&ptr);
            addVal2CharArray(m3,&ptr);
        }
        Send(proxy->getSocketHelper(), proxy->getBuffer1(), sz*16);
        Receive(proxy->getSocketHelper(), proxy->getBuffer1(), sz*8);
        ptr = proxy->getBuffer1();
        for (uint32_t i = 0; i < sz; i++) {
            res[i] = m1[i] + convert2Long(&ptr);
            res[i] = res[i] >> FRAC;
            res[i] = x[i] - res[i];
        }
        delete [] m1;
        return res;

    }else if ( proxy->getPRole() == P2){
        unsigned char *ptr = proxy->getBuffer1();
        uint64_t *res = new uint64_t[sz];
        uint64_t *m1 = new uint64_t[sz];
        for (uint32_t i = 0; i < sz; i++) {
            uint64_t r1=proxy->generateCommonRandom(), r2=proxy->generateCommonRandom(), r3=proxy->generateCommonRandom(), r4=proxy->generateCommonRandom();

            m1[i] = (b[i] * (x[i] - y[i])) - (r1*(x[i] - y[i])) - (r1*r2) - (r4*b[i]);
            uint64_t m2 = x[i] - y[i] + r2;
            uint64_t m3 = b[i] + r3;

            addVal2CharArray(m2,&ptr);
            addVal2CharArray(m3,&ptr);

        }
        Send(proxy->getSocketHelper(), proxy->getBuffer1(), sz*16);
        Receive(proxy->getSocketHelper(), proxy->getBuffer1(), sz*8);
        ptr = proxy->getBuffer1();
        for (uint32_t i = 0; i < sz; i++) {
            res[i] = m1[i] + convert2Long(&ptr);
            res[i] = -1 * ((-1 * res[i]) >> FRAC);
            res[i] = x[i] - res[i];
        }
        delete [] m1;
        return res;

    }else if ( proxy->getPRole() == HELPER){
        thread thr1 = thread(Receive,proxy->getSocketP1(), proxy->getBuffer1(), sz*16);
        thread thr2 = thread(Receive,proxy->getSocketP2(), proxy->getBuffer2(), sz*16);
        thr1.join();
        thr2.join();

        //Receive(proxy->getSocketP1(), proxy->getBuffer1(), sz*24);
        unsigned char *ptr = proxy->getBuffer1();
        unsigned char *ptr_out = proxy->getBuffer1();

        //Receive(proxy->getSocketP2(), proxy->getBuffer2(), sz*24);
        unsigned char *ptr2 = proxy->getBuffer2();
        unsigned char *ptr_out2 = proxy->getBuffer2();
        for (uint32_t i = 0; i < sz; i++) {
            uint64_t m2 = convert2Long(&ptr);
            uint64_t m3 = convert2Long(&ptr);

            uint64_t m5 = convert2Long(&ptr2);
            uint64_t m6 = convert2Long(&ptr2);

            uint64_t m = (m2 * m5) + (m3 * m6);
            m2 = proxy->generateRandom();
            m3 = m-m2;
            addVal2CharArray(m2,&ptr_out);
            addVal2CharArray(m3,&ptr_out2);
        }
        thr1 = thread(Send,proxy->getSocketP1(), proxy->getBuffer1(), sz*8);
        thr2 = thread(Send,proxy->getSocketP2(), proxy->getBuffer2(), sz*8);
        thr1.join();
        thr2.join();
        return NULL;
    }
    return NULL;
}

uint8_t PCB(Party* proxy, uint64_t a, uint8_t *b, int L1) {
    // b is a boolean share
    // a is reconstructed value
    // check b>a

    if ( proxy->getPRole() == P1 ||  proxy->getPRole() == P2) {
        uint8_t w_sum = 0;
        for (int i = L1 - 1; i >= 0; i--) {
            uint8_t a_bit = bit(a, i);
            int k = i;
            uint8_t w = mod((b[k] +  proxy->getPRole() * a_bit - 2 * a_bit * b[k]) % LP, LP);
            proxy->getBuffer1()[k] =
                    (mod(( proxy->getPRole() * a_bit - b[k] +  proxy->getPRole() + w_sum), LP) * (proxy->generateCommonRandom() % (LP - 1) + 1)) %
                    LP;
            w_sum = (w_sum + w) % LP;
        }
        for (int i = 0; i < L1; i++) {
            int ind1 = (proxy->generateCommonRandom() % L1);
            int ind2 = (proxy->generateCommonRandom() % L1);
            uint8_t tmp = proxy->getBuffer1()[ind1];
            proxy->getBuffer1()[ind1] = proxy->getBuffer1()[ind2];
            proxy->getBuffer1()[ind2] = tmp;
        }
        Send(proxy->getSocketHelper(), proxy->getBuffer1(), L1);
        Receive(proxy->getSocketHelper(), proxy->getBuffer1(), 1);
        uint8_t r = proxy->getBuffer1()[0];
        return r;

    } else if ( proxy->getPRole() == HELPER) {
        Receive(proxy->getSocketP1(), proxy->getBuffer1(), L1);
        Receive(proxy->getSocketP2(), proxy->getBuffer2(), L1);
        unsigned char *ptr_out = proxy->getBuffer1();
        unsigned char *ptr_out2 = proxy->getBuffer2();
        uint8_t res = 0;
        for (int i = 0; i < L1; i++) {
            proxy->getBuffer1()[i] = (proxy->getBuffer1()[i] + proxy->getBuffer2()[i]) % LP;
            if (((int) proxy->getBuffer1()[i]) == 0) {
                res = 1;
                break;
            }
        }
        uint8_t res1 = proxy->generateRandom() % 2;
        uint8_t res2 = res ^res1;

        addVal2CharArray(res1, &ptr_out);
        addVal2CharArray(res2, &ptr_out2);
        Send(proxy->getSocketP1(), proxy->getBuffer1(), 1);
        Send(proxy->getSocketP2(), proxy->getBuffer2(), 1);
        return 0;
    }
    return -1;
}

uint8_t *PCB(Party* proxy, uint64_t *a, uint8_t *b, uint32_t sz, int L1) {
    // b is a boolean share
    // a is reconstructed value
    // check b>a

    if ( proxy->getPRole() == P1 ||  proxy->getPRole() == P2) {
        for (int j = 0; j < sz; j++) {
            int jk = j * L1;
            uint8_t w_sum = 0;
            for (int i = L1 - 1; i >= 0; i--) {
                uint8_t a_bit = bit(a[j], i);
                int k = jk + i;
                uint8_t w = mod((b[k] +  proxy->getPRole() * a_bit - 2 * a_bit * b[k]) % LP, LP);
                proxy->getBuffer1()[k] =
                        (mod(( proxy->getPRole() * a_bit - b[k] +  proxy->getPRole() + w_sum), LP) * (proxy->generateCommonRandom() % (LP - 1) + 1)) %
                        LP;
                w_sum = (w_sum + w) % LP;
            }
            for (int i = 0; i < L1; i++) {
                int ind1 = (proxy->generateCommonRandom() % L1) + jk;
                int ind2 = (proxy->generateCommonRandom() % L1) + jk;
                uint8_t tmp = proxy->getBuffer1()[ind1];
                proxy->getBuffer1()[ind1] = proxy->getBuffer1()[ind2];
                proxy->getBuffer1()[ind2] = tmp;
            }
        }
        Send(proxy->getSocketHelper(), proxy->getBuffer1(), sz * L1);
        Receive(proxy->getSocketHelper(), proxy->getBuffer1(), sz);
        uint8_t *r = new uint8_t[sz];
        for (int i = 0; i < sz; i++)
            r[i] = proxy->getBuffer1()[i];
        return r;

    } else if ( proxy->getPRole() == HELPER) {
        Receive(proxy->getSocketP1(), proxy->getBuffer1(), sz * L1);
        Receive(proxy->getSocketP2(), proxy->getBuffer2(), sz * L1);
        unsigned char *ptr_out = proxy->getBuffer1();
        unsigned char *ptr_out2 = proxy->getBuffer2();
        for (int j = 0; j < sz; j++) {
            int jk = j * L1;
            uint8_t res = 0;
            for (int i = 0; i < L1; i++) {
                proxy->getBuffer1()[jk + i] = (proxy->getBuffer1()[jk + i] + proxy->getBuffer2()[jk + i]) % LP;
                if (((int) proxy->getBuffer1()[jk + i]) == 0) {
                    res = 1;
                    break;
                }
            }
            uint8_t res1 = proxy->generateRandom() % 2;
            uint8_t res2 = res ^res1;

            addVal2CharArray(res1, &ptr_out);
            addVal2CharArray(res2, &ptr_out2);
        }
        Send(proxy->getSocketP1(), proxy->getBuffer1(), sz);
        Send(proxy->getSocketP2(), proxy->getBuffer2(), sz);
        return NULL;
    }
    return NULL;
}

uint64_t MOC(Party* proxy, uint64_t x) {
    // x is a value in the ring 2^63.
    if ( proxy->getPRole() == P1 ||  proxy->getPRole() == P2) {
        uint64_t z_1;
        uint64_t ya;
        uint8_t yb[L_BIT - 1];
        uint8_t w;
        Receive(proxy->getSocketHelper(), proxy->getBuffer1(), (8 + L_BIT));
        unsigned char *ptr = proxy->getBuffer1();
        // an arithmetic share of y from the helper  : ya
        ya = convert2Long(&ptr);
        // boolean shares of each bit of y: yb
        convert2Array(&ptr, &yb[0], L_BIT - 1);
        // a boolean share of whether addition of arithmetic shares (ya_1+ya_2) of y wraps around in the ring 2^63 :w
        w = (*ptr);
        // z1 is x+ya
        z_1 = (x + ya) & N1_MASK;
        // z is the reconstruction of z_1. From z, P1 or P2 can not learn x because x is masked with ya.
        uint64_t z = REC(proxy, z_1,N1_MASK);
        // computes y>z where yb are boolean shares of each bit of y
        uint8_t wc = PCB(proxy, z, yb, L_BIT - 1);

        // w=1 means there is an overflow
        // wc=1 means there is an overflow
        // w = w ^ wc determines whether there is an overflow
        w = w ^ wc;
        if ( proxy->getPRole() == P1 && z_1 > z)
            // if the addition of arithmetic shares of z wraps around P1 adds 2^63 the arithmetic share of z.
            z_1 = z_1 + N1;
        // removing y and the overflow (if there is one) from z
        z_1 = z_1 - (ya + w * N1);
        return z_1;
    }
    else if ( proxy->getPRole() == HELPER) {
        //cout << "start helper MOC" << endl;
        unsigned char *ptr_out = proxy->getBuffer1();
        unsigned char *ptr_out2 = proxy->getBuffer2();
        // helper picks a random number in the ring 2^63
        uint64_t y = proxy->generateRandom() & N1_MASK;
        // helper creates two shares for y in the ring 2^63: ya_1 and ya_2
        uint64_t ya_1 = proxy->generateRandom() & N1_MASK;
        uint64_t ya_2 = (y - ya_1) & N1_MASK;

        // adding ya_1 and ya_2 to proxy->getBuffer1() and proxy->getBuffer2() respectively.
        addVal2CharArray(ya_1, &ptr_out);
        addVal2CharArray(ya_2, &ptr_out2);

        // helper creates two boolean shares for each bit of y : yb_1 and yb_2
        // writing yb_1 and yb_2 to proxy->getBuffer1() and proxy->getBuffer2() respectively.
        for (int j = 0; j < L_BIT - 1; j++) {
            uint8_t k = (y >> j) & 0x1;
            uint8_t yb_1 = proxy->generateRandom() % LP;
            uint8_t yb_2 = mod(k - yb_1, LP);
            addVal2CharArray(yb_1, &ptr_out);
            addVal2CharArray(yb_2, &ptr_out2);
        }

        // if ya_1+ya_2 wraps around 2^63 w=1.
        uint8_t w = 0;
        if (ya_1 > y)
            w = 1;

        // creating two boolean shares of w : w_1 and w_2
        uint8_t w_1 = proxy->generateRandom() % 2;
        uint8_t w_2 = w ^w_1;
        // writing w_1 and w_2 to proxy->getBuffer1() and proxy->getBuffer2() respectively.
        addVal2CharArray(w_1, &ptr_out);
        addVal2CharArray(w_2, &ptr_out2);

        // sending values to P1 and P2
        thread thr1 = thread(Send,proxy->getSocketP1(), proxy->getBuffer1(), (8 + L_BIT));
        thread thr2 = thread(Send,proxy->getSocketP2(), proxy->getBuffer2(), (8 + L_BIT));
        thr1.join();
        thr2.join();

        // P1 and P2 will call PrivateCompareBool
        uint8_t tmp = PCB(proxy, 0, 0, L_BIT - 1);
        return 0;
    }
    return -1;
}

uint64_t *MOC(Party* proxy, uint64_t *x, uint32_t sz) {
    // x is a value in the ring 2^63.
    // sz is the length of x.
    if ( proxy->getPRole() == P1 ||  proxy->getPRole() == P2) {
        uint64_t *z_1 = new uint64_t[sz];
        uint64_t *ya = new uint64_t[sz];
        uint8_t *yb = new uint8_t[sz * (L_BIT - 1)];
        uint8_t *w = new uint8_t[sz];
        Receive(proxy->getSocketHelper(), proxy->getBuffer1(), sz * (8 + L_BIT));
        unsigned char *ptr = proxy->getBuffer1();
        for (int i = 0; i < sz; i++) {
            ya[i] = convert2Long(&ptr);
            convert2Array(&ptr, &yb[i * (L_BIT - 1)], L_BIT - 1);
            w[i] = (*ptr);
            ptr++;
            z_1[i] = (x[i] + ya[i]) & N1_MASK;
        }
        uint64_t *z = REC(proxy, z_1, sz, N1_MASK);
        uint8_t *wc = PCB(proxy, z, yb, sz, L_BIT - 1);

        for (int i = 0; i < sz; i++) {
            w[i] = w[i] ^ wc[i];
            if ( proxy->getPRole() == P1 && z_1[i] > z[i])
                z_1[i] = z_1[i] + N1;
            z_1[i] = (z_1[i] - (ya[i] + w[i] * N1));
        }
        delete[] ya;
        delete[] yb;
        delete[] w;
        return z_1;
    }
    else if ( proxy->getPRole() == HELPER) {
        unsigned char *ptr_out = proxy->getBuffer1();
        unsigned char *ptr_out2 = proxy->getBuffer2();
        for (int i = 0; i < sz; i++) {
            uint64_t y = proxy->generateRandom() & N1_MASK;
            uint64_t ya_1 = proxy->generateRandom() & N1_MASK;
            uint64_t ya_2 = (y - ya_1) & N1_MASK;
            addVal2CharArray(ya_1, &ptr_out);
            addVal2CharArray(ya_2, &ptr_out2);
            for (int j = 0; j < L_BIT - 1; j++) {
                uint8_t k = (y >> j) & 0x1;
                uint8_t yb_1 = proxy->generateRandom() % LP;
                uint8_t yb_2 = mod(k - yb_1, LP);
                addVal2CharArray(yb_1, &ptr_out);
                addVal2CharArray(yb_2, &ptr_out2);
            }
            uint8_t w = 0;
            if (ya_1 > y)
                w = 1;
            uint8_t w_1 = proxy->generateRandom() % 2;
            uint8_t w_2 = w ^w_1;
            addVal2CharArray(w_1, &ptr_out);
            addVal2CharArray(w_2, &ptr_out2);
        }

        thread thr1 = thread(Send,proxy->getSocketP1(), proxy->getBuffer1(), sz * (8 + L_BIT));
        thread thr2 = thread(Send,proxy->getSocketP2(), proxy->getBuffer2(), sz * (8 + L_BIT));
        thr1.join();
        thr2.join();
        // P1 and P2 will call MPrivateCompareBool
        uint8_t *tmp = PCB(proxy, 0, 0, sz, L_BIT - 1);
        return NULL;
    }
    return NULL;
}

uint64_t MSB(Party *proxy, uint64_t x) {
    if ( proxy->getPRole() == P1 ||  proxy->getPRole() == P2) {
        uint64_t d_k = x & N1_MASK;
        uint64_t d = MOC(proxy, d_k);
        uint64_t m = x - d;
        unsigned char *ptr = proxy->getBuffer1();

        ptr = proxy->getBuffer1();
        uint64_t added_random;
        added_random = proxy->generateCommonRandom()%2;
        if ( proxy->getPRole() == P1) {
            if (added_random == 0) {
                addVal2CharArray(m, &ptr);
                addVal2CharArray(m + N1, &ptr);
            }else{
                addVal2CharArray(m + N1, &ptr);
                addVal2CharArray(m, &ptr);
            }
        }
        else {
            addVal2CharArray(m, &ptr);
            addVal2CharArray(m, &ptr);
        }

        Send(proxy->getSocketHelper(), proxy->getBuffer1(), 8 * 2);
        Receive(proxy->getSocketHelper(), proxy->getBuffer1(),  8 * 2);
        ptr = proxy->getBuffer1();

        if (added_random == 0){
            m = convert2Long(&ptr);
            ptr += 8;
        }else{
            ptr += 8;
            m = convert2Long(&ptr);
        }

        return m;

    } else if ( proxy->getPRole() == HELPER) {

        unsigned char *ptr_out = proxy->getBuffer1();
        unsigned char *ptr_out2 = proxy->getBuffer2();
        uint64_t d = MOC(proxy, 0);


        thread thr1 = thread(Receive,proxy->getSocketP1(), proxy->getBuffer1(), 8 * 2);
        thread thr2 = thread(Receive,proxy->getSocketP2(), proxy->getBuffer2(), 8 * 2);
        thr1.join();
        thr2.join();

        unsigned char *ptr = proxy->getBuffer1();
        ptr_out = proxy->getBuffer1();

        unsigned char *ptr2 = proxy->getBuffer2();
        ptr_out2 = proxy->getBuffer2();

        uint64_t v1 = convert2Long(&ptr2) + convert2Long(&ptr);
        uint64_t v2 = convert2Long(&ptr2) + convert2Long(&ptr);
        v1 = convert2uint64((double)(v1/N1));
        v2 = convert2uint64((double)(v2/N1));
        uint64_t s1 = proxy->generateRandom();
        uint64_t s2 = v1 - s1;
        addVal2CharArray(s1, &ptr_out);
        addVal2CharArray(s2, &ptr_out2);
        s1 = proxy->generateRandom();
        s2 = v2 - s1;
        addVal2CharArray(s1, &ptr_out);
        addVal2CharArray(s2, &ptr_out2);

        thr1 = thread(Send,proxy->getSocketP1(), proxy->getBuffer1(),  8 * 2);
        thr2 = thread(Send,proxy->getSocketP2(), proxy->getBuffer2(), 8 * 2);
        thr1.join();
        thr2.join();
        return 0;
    }
    return -1;
}

uint64_t *MSB(Party* proxy, uint64_t *x, uint32_t sz) {
    if ( proxy->getPRole() == P1 ||  proxy->getPRole() == P2) {
        uint64_t* d_k = new uint64_t[sz];
        uint64_t *m = new uint64_t[sz];
        for (int i = 0; i < sz; i++) {
            d_k[i] = x[i] & N1_MASK;
        }
        uint64_t* d = MOC(proxy, d_k,sz);
        delete [] d_k;
        for (int i = 0; i < sz; i++) {
            m[i] = x[i] - d[i];
        }
        unsigned char *ptr = proxy->getBuffer1();

        ptr = proxy->getBuffer1();
        uint64_t *added_random = new uint64_t[sz];
        for (int i = 0; i < sz; i++) {

            added_random[i] = (proxy->generateCommonRandom()%2);
            if ( proxy->getPRole() == P1) {
                if (added_random[i] == 0) {
                    addVal2CharArray(m[i], &ptr);
                    addVal2CharArray(m[i] + N1, &ptr);
                }else{
                    addVal2CharArray(m[i] + N1, &ptr);
                    addVal2CharArray(m[i], &ptr);
                }
            }
            else {
                addVal2CharArray(m[i], &ptr);
                addVal2CharArray(m[i], &ptr);
            }
        }
        Send(proxy->getSocketHelper(), proxy->getBuffer1(), sz * 8 * 2);
        Receive(proxy->getSocketHelper(), proxy->getBuffer1(), sz * 8 * 2);
        ptr = proxy->getBuffer1();

        for (int i = 0; i < sz; i++) {
            if (added_random[i] == 0){
                m[i] = convert2Long(&ptr);
                ptr += 8;
            }else{
                ptr += 8;
                m[i] = convert2Long(&ptr);
            }
        }

        delete[] added_random;
        return m;

    } else if ( proxy->getPRole() == HELPER) {

        unsigned char *ptr_out = proxy->getBuffer1();
        unsigned char *ptr_out2 = proxy->getBuffer2();
        uint64_t* d = MOC(proxy,0,sz);

        thread thr1 = thread(Receive,proxy->getSocketP1(), proxy->getBuffer1(), sz * 8 * 2);
        thread thr2 = thread(Receive,proxy->getSocketP2(), proxy->getBuffer2(), sz * 8 * 2);
        thr1.join();
        thr2.join();

        unsigned char *ptr = proxy->getBuffer1();
        ptr_out = proxy->getBuffer1();

        unsigned char *ptr2 = proxy->getBuffer2();
        ptr_out2 = proxy->getBuffer2();

        for (uint32_t i = 0; i < sz; i++) {
            uint64_t v1 = convert2Long(&ptr2) + convert2Long(&ptr);
            uint64_t v2 = convert2Long(&ptr2) + convert2Long(&ptr);
            v1 = convert2uint64((double)(v1/N1));
            v2 = convert2uint64((double)(v2/N1));
            uint64_t s1 = proxy->generateRandom();
            uint64_t s2 = v1 - s1;
            addVal2CharArray(s1, &ptr_out);
            addVal2CharArray(s2, &ptr_out2);
            s1 = proxy->generateRandom();
            s2 = v2 - s1;
            addVal2CharArray(s1, &ptr_out);
            addVal2CharArray(s2, &ptr_out2);
        }
        thr1 = thread(Send,proxy->getSocketP1(), proxy->getBuffer1(), sz * 8 * 2);
        thr2 = thread(Send,proxy->getSocketP2(), proxy->getBuffer2(), sz * 8 * 2);
        thr1.join();
        thr2.join();
        return NULL;
    }
    return NULL;
}

uint64_t *CMP(Party* proxy, uint64_t *x, uint64_t *y,uint32_t sz) {
    if ( proxy->getPRole() == P1 ||  proxy->getPRole() == P2) {
        uint64_t* diff = new uint64_t[sz];
        for (int i = 0; i < sz; i++) {
            diff[i] = x[i] - y[i];
        }
        uint64_t* m = MSB(proxy, diff,sz);
        for (int i = 0; i < sz; i++) {
            m[i] =  (proxy->getPRole()<<FRAC) - m[i];
        }
        return m;
    }else if ( proxy->getPRole() == HELPER) {
        uint64_t* m = MSB(proxy, 0,sz);
        return NULL;
    }
    return NULL;
}

uint64_t CMP(Party* proxy, uint64_t x, uint64_t y) {
    if ( proxy->getPRole() == P1 ||  proxy->getPRole() == P2) {
        uint64_t diff = x - y;
        return  (proxy->getPRole()<<FRAC) - MSB(proxy, diff);
    }else if ( proxy->getPRole() == HELPER) {
        uint64_t m = MSB(proxy,0);
        return 0;
    }
    return -1;
}

uint64_t MUL(Party* proxy, uint64_t a, uint64_t b) {
    /*
     * Input(s)
     * a: a share of the first multiplicand - uint64_t
     * b: a share of the second multiplicand - uint64_t
     *
     * Output(s)
     * Returns the share of the multiplication of a and b - uint64_t
     */
    if(DEBUG_FLAG >= 1)
        cout << "************************************************************\nNF_MUL is called" << endl;
    if (proxy->getPRole() == HELPER) {
        uint64_t *mt1[3];
        uint64_t *mt2[3];
        for (int i = 0; i < 3; i++) {
            mt1[i] = new uint64_t[1];
            mt2[i] = new uint64_t[1];
        }
        GenerateMultiplicationTriple(proxy, mt1, mt2, 1);

        // send the multiplication triples to P1
        unsigned char *ptr_out = proxy->getBuffer1();
        for (auto &i : mt1) {
            addVal2CharArray(i[0], &ptr_out);
        }

        // addVal2CharArray(mt1, &ptr_out, 3, size); // a special method is needed here!
        Send(proxy->getSocketP1(), proxy->getBuffer1(), 3 * 8);

        // send the multiplication triples to P2
        unsigned char *ptr_out2 = proxy->getBuffer2();
        for (auto &i : mt2) {
            addVal2CharArray(i[0], &ptr_out2);
        }
        // addVal2CharArray(mt2, &ptr_out2, 3, size);
        Send(proxy->getSocketP2(), proxy->getBuffer2(), 3 * 8);

        for (int i = 0; i < 3; i++) {
            delete[] mt1[i];
            delete[] mt2[i];
        }
        if(DEBUG_FLAG >= 1)
            cout << "Returning from NF_MUL...\n************************************************************" << endl;
        return 0;

    } else if (proxy->getPRole() == P1 || proxy->getPRole() == P2) {
        Receive(proxy->getSocketHelper(), proxy->getBuffer1(), 3 * 8);
        unsigned char *ptr = proxy->getBuffer1();
        uint64_t mt[3];
        for (auto &i : mt) {
            i = convert2Long(&ptr);
        }

        uint64_t e_f[2];
        e_f[0] = a - mt[0];
        e_f[1] = b - mt[1];

//        uint64_t e_shr = a - mt[0];
//        uint64_t f_shr = b - mt[1];

//        uint64_t e = Reconstruct(e_shr);
//        uint64_t f = Reconstruct(f_shr);

        uint64_t* rec_e_f = REC(proxy,e_f, 2);

//        uint64_t z = proxy->getPRole() * e * f + f * mt[0] + e * mt[1] + mt[2];
        uint64_t z = proxy->getPRole() * rec_e_f[0] * rec_e_f[1] + rec_e_f[1] * mt[0] + rec_e_f[0] * mt[1] + mt[2];

        // restore the fractional part - refer to SecureNN for more details
        if (proxy->getPRole() == P1) {
            z = z >> FRAC;
        } else {
            z = -1 * ((-1 * z) >> FRAC);
        }

        delete [] rec_e_f;

        if(DEBUG_FLAG >= 1)
            cout << "Returning from NF_MUL...\n************************************************************" << endl;
        return z;
    } else {
        return -1;
    }
}

uint64_t *PMUL(Party* proxy, uint64_t *a, uint64_t *b, uint32_t size) {
    if(DEBUG_FLAG >= 1)
        cout << "************************************************************\nPMNF_MUL is called" << endl;
    if (proxy->getPRole() == HELPER) {
        uint64_t *mt1[3];
        uint64_t *mt2[3];
        for (int i = 0; i < 3; i++) {
            mt1[i] = new uint64_t[size];
            mt2[i] = new uint64_t[size];
        }

        GenerateMultiplicationTriple(proxy, mt1, mt2, size);

        // send the multiplication triples to P1
        unsigned char *ptr_out = proxy->getBuffer1();
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < size; j++) {
                addVal2CharArray(mt1[i][j], &ptr_out);
            }
        }
        // addVal2CharArray(mt1, &ptr_out, 3, size); // a special method is needed here!
        Send(proxy->getSocketP1(), proxy->getBuffer1(), size * 3 * 8);

        // send the multiplication triples to P2
        unsigned char *ptr_out2 = proxy->getBuffer2();
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < size; j++) {
                addVal2CharArray(mt2[i][j], &ptr_out2);
            }
        }

        // addVal2CharArray(mt2, &ptr_out2, 3, size);
        Send(proxy->getSocketP2(), proxy->getBuffer2(), size * 3 * 8);

        for (int i = 0; i < 3; i++) {
            delete[] mt1[i];
            delete[] mt2[i];
        }

        Receive(proxy->getSocketP1(), proxy->getBuffer1(), 1);
        Receive(proxy->getSocketP2(), proxy->getBuffer2(), 1);
        if(DEBUG_FLAG >= 1)
            cout << "Returning from PMNF_MUL...\n************************************************************" << endl;
        return 0;

    } else if (proxy->getPRole() == P1 || proxy->getPRole() == P2) {
        //total_mul += size;
        Receive(proxy->getSocketHelper(), proxy->getBuffer1(), size * 3 * 8);

        unsigned char *ptr = proxy->getBuffer1();
        // uint64_t **mt = new uint64_t*[3];
        uint64_t *mt[3];
        for (int i = 0; i < 3; i++) {
            mt[i] = new uint64_t[size];
            for (int j = 0; j < size; j++) {
                mt[i][j] = convert2Long(&ptr);
            }
        }

        // concatenated form of e and f shares
        uint64_t* concat_e_f = new uint64_t[size * 2];
        for(int i = 0; i < size; i++) {
            concat_e_f[i] = a[i] - mt[0][i];
            concat_e_f[i + size] = b[i] - mt[1][i];
        }

        uint64_t* e_f = REC(proxy, concat_e_f, size * 2);

        uint64_t *e = e_f;
        uint64_t *f = &e_f[size];

        uint64_t *z = new uint64_t[size];
        for (int i = 0; i < size; i++) {
            z[i] = proxy->getPRole() * e[i] * f[i] + f[i] * mt[0][i] + e[i] * mt[1][i] + mt[2][i];
//            cout << i << ": " << z[i] << endl;
            if (proxy->getPRole() == P1) {
                z[i] = z[i] >> FRAC;
            } else {
                z[i] = -1 * ((-1 * z[i]) >> FRAC);
            }
        }

        delete [] e_f;
        for (auto &i : mt) {
            delete[] i;
        }
        proxy->getBuffer1()[0] = 0;
        Send(proxy->getSocketHelper(), proxy->getBuffer1(), 1);
        if(DEBUG_FLAG >= 1)
            cout << "Returning from PMNF_MUL...\n************************************************************" << endl;
        return z;
    } else {
        return nullptr;
    }
}

uint64_t *MUL(Party* proxy, uint64_t *a, uint64_t *b, uint32_t size) {
    /*
     * Input(s)
     * a: one of the vectors of shares of the multiplicands - uint64_t vector
     * b: the other vector of shares of the multiplicands - uint64_t vector
     * size: the size of the vectors a and b
     *
     * Output(s)
     * Returns an uint64_t vector containing the share of the result of the multiplication
     */
    cout << "MUL:Size: " << size << endl;
    if(DEBUG_FLAG >= 1)
        cout << "************************************************************\nMNF_MUL is called" << endl;
    if (proxy->getPRole() == HELPER) {
        int filled_size = 0;
        size_t partial_size = MAXMUL;
        while (filled_size < size) {
            if ((size - filled_size) < MAXMUL) {
                partial_size = (size - filled_size);
            }
            PMUL(proxy,0, 0, partial_size);
            filled_size += partial_size;
        }
    } else if (proxy->getPRole() == P1 || proxy->getPRole() == P2) {

        /*if(DEBUG_FLAG >= 2)
            printValue("MNF_MUL size", (double)size);*/

        uint64_t *result = new uint64_t[size];
        int filled_size = 0;
        size_t partial_size = MAXMUL;
        while (filled_size < size) {
            if ((size - filled_size) < MAXMUL) {
                partial_size = (size - filled_size);
            }
            uint64_t *partial_result = PMUL(proxy, a, b, partial_size);
            std::copy(partial_result, partial_result + partial_size, result + filled_size);
            delete[] partial_result;
            filled_size += partial_size;
            a += partial_size;
            b += partial_size;
        }
        if(DEBUG_FLAG >= 1)
            cout << "Returning from MNF_MUL...\n************************************************************" << endl;
        return result;
    }
    if(DEBUG_FLAG >= 1)
        cout << "Returning from MNF_MUL...\n************************************************************" << endl;
    return nullptr;
}

uint64_t EXP(Party* proxy, uint64_t a) {
    /*
     * Note that this function considers only the specific number of least significant bits not to cause overflow.
     * This is different for positive and negative powers.
     *
     * Input(s)
     * a: the value that will be used as the power of exp
     * n_bits: the number of bits that one can need to cover the largest possible power
     *
     * Output(s)
     * Returns the arithmetic secret share of exp(a)
     */
    int p_role = proxy->getPRole();
    int n_bits = proxy->getNBits();
    int neg_n_bits = proxy->getNegNBits();
    if (p_role == P1 || p_role == P2) {
        // compute the absolute of the input value
        uint64_t msb_a = MSB(proxy, a);
        uint64_t abs_a = ((uint64_t) 0) - a;

        // selection of the correct contribution from each bit of the input value based on the msb of the input value
        uint64_t *pos_e_contributions = new uint64_t[n_bits + 1];
        uint64_t *neg_e_contributions = new uint64_t[n_bits + 1];
        uint64_t *one_contributions = new uint64_t[n_bits];
        uint64_t *repeated_msb_a = new uint64_t[n_bits + 1];

        pos_e_contributions[0] = a;
        neg_e_contributions[0] = abs_a;
        repeated_msb_a[0] = msb_a;

        for (int i = n_bits - 1; i >= 0; i--) {
            pos_e_contributions[n_bits - i] = p_role * convert2uint64(exp(pow(2, i - FRAC)));
            if (i > neg_n_bits - 1) {
                neg_e_contributions[n_bits - i] = p_role * (((uint64_t) 1) << FRAC);
            } else {
                neg_e_contributions[n_bits - i] = p_role * convert2uint64(1.0 / exp(pow(2, i - FRAC)));
            }
            one_contributions[(n_bits - 1) - i] = p_role * (((uint64_t) 1) << FRAC);
            repeated_msb_a[n_bits - i] = msb_a;
        }

        uint64_t *e_contributions = MUX(proxy, pos_e_contributions, neg_e_contributions, repeated_msb_a, n_bits + 1);
        uint64_t new_a = e_contributions[0];
        e_contributions = &e_contributions[1];
        uint64_t horse = REC(proxy, new_a);

        // arrange all the shifted versions of the input value for MMSB
        uint64_t *partial_a = new uint64_t[n_bits];
        uint64_t tmp = new_a << (L_BIT - n_bits);
        for (int i = 0; i < n_bits; i++) {
            partial_a[i] = tmp << i;
        }

        uint64_t *bit_shares = MSB(proxy, partial_a, n_bits);
        uint64_t *at = REC(proxy, bit_shares, n_bits);

        // selection of the contribution of the bits of the value
        uint64_t *contributions = MUX(proxy, one_contributions, e_contributions, bit_shares, n_bits);

        // binary-tree-based multiplication of the contributions into the exponential
        uint64_t res = 1;
        int current_size = n_bits;
        bool flag = false;
        uint64_t remaining = 0;

        for (int i = 0; i < (int) ceil(log2(n_bits)); i++) {
            uint64_t *tmp1 = contributions;
            uint64_t *tmp2 = &contributions[current_size / 2];

            if (current_size % 2 == 1) {
                if (!flag) {
                    remaining = contributions[current_size - 1];
                    flag = true;
                } else {
                    tmp1 = new uint64_t[(current_size + 1) / 2];
                    tmp2 = new uint64_t[(current_size + 1) / 2];

                    size_t partial_size = current_size / 2;

                    copy(contributions, contributions + partial_size + 1, tmp1);
                    copy(contributions + partial_size + 1, contributions + current_size, tmp2);

                    tmp2[current_size / 2] = remaining;

                    current_size++;

                    flag = false;
                }
            }
            contributions = MUL(proxy, tmp1, tmp2, current_size / 2);
            current_size /= 2;
        }
        return contributions[0];
    }
    else if (p_role == HELPER) {
        MSB(proxy, 0);
        MUX(proxy, 0, 0, 0, n_bits + 1);
        MSB(proxy, 0, n_bits);
        MUX(proxy, 0, 0, 0, n_bits);

        int current_size = n_bits;
        bool flag = false;
        for (int i = 0; i < (int) ceil(log2(n_bits)); i++) {
            if (current_size % 2 == 1) {
                if (!flag) {
                    flag = true;
                } else {
                    current_size++;
                    flag = false;
                }
            }
            MUL(proxy, 0, 0, current_size / 2);
            current_size /= 2;
        }

        return 0;
    }
    return 0;
}

uint64_t* EXP(Party* proxy, uint64_t *a, uint32_t size) {
    /*
     * Note that this function considers only the specific number of least significant bits not to cause overflow.
     * This is different for positive and negative powers.
     *
     * Input(s)
     * a: the vector of values that will be used as the power of exp
     * n_bits: the number of bits that one can need to cover the largest possible power
     *
     * Output(s)
     * Returns the arithmetic secret share of exp(a)
     */
    int p_role = proxy->getPRole();
    int n_bits = proxy->getNBits();
    int neg_n_bits = proxy->getNegNBits();

    if (p_role == P1 || p_role == P2) {
        // compute the absolute of the input value
        uint64_t* msb_a = MSB(proxy, a, size);
        uint64_t* abs_a = new uint64_t[size];
        for(uint32_t i = 0; i < size; i++) {
            abs_a[i] = ((uint64_t) 0) - a[i];
        }

        // compute the possible contribution of positive and negative values
        uint64_t* pec = new uint64_t[n_bits];
        uint64_t* nec = new uint64_t[n_bits];
        if(p_role == P2) {
            for (int i = n_bits - 1; i >= 0; i--) {
                pec[n_bits - i - 1] = convert2uint64(exp(pow(2, i - FRAC)));
                if (i > neg_n_bits - 1) {
                    nec[n_bits - i - 1] = (((uint64_t) 1) << FRAC);
                } else {
                    nec[n_bits - i - 1] = convert2uint64(1.0 / exp(pow(2, i - FRAC)));
                }
            }
        }
        else {
            for (int i = n_bits - 1; i >= 0; i--) {
                pec[n_bits - i - 1] = 0;
                nec[n_bits - i - 1] = 0;
            }
        }

        // selection of the correct contribution from each bit of the input value based on the msb of the input value
        uint64_t *pos_e_contributions = new uint64_t[size * (n_bits + 1)]; // if the value is positive
        uint64_t *neg_e_contributions = new uint64_t[size * (n_bits + 1)]; // if the value is negative
        uint64_t *one_contributions = new uint64_t[size * n_bits]; // if the bit is zero regardless of the sign
        uint64_t *repeated_msb_a = new uint64_t[size * (n_bits + 1)]; // this will be used as a selection bit for the contributions of all bits

        for(uint32_t i = 0; i < size; i++) {
            pos_e_contributions[i * (n_bits + 1)] = a[i];
            neg_e_contributions[i * (n_bits + 1)] = abs_a[i];
            repeated_msb_a[i * (n_bits + 1)] = msb_a[i];
            for (int bi = 0; bi < n_bits; bi++) {
                pos_e_contributions[(i * (n_bits + 1)) + bi + 1] = pec[bi];
                neg_e_contributions[(i * (n_bits + 1)) + bi + 1] = nec[bi];
                one_contributions[(i * n_bits) + bi] = p_role * (((uint64_t) 1) << FRAC);
                repeated_msb_a[(i * (n_bits + 1)) + bi + 1] = msb_a[i];
            }
        }
        uint64_t *e_contributions = MUX(proxy, pos_e_contributions, neg_e_contributions, repeated_msb_a, size * (n_bits + 1));

        uint64_t* new_a = new uint64_t[size];
        uint64_t* selected_e_contributions = new uint64_t[size * n_bits];
        for(uint32_t i = 0; i < size; i++) {
            new_a[i] = e_contributions[i * (n_bits + 1)];
            for(uint32_t j = 0; j < n_bits; j++) {
                selected_e_contributions[(i * n_bits) + j] = e_contributions[(i * (n_bits + 1)) + j + 1];
            }
        }

        // arrange all the shifted versions of the input value for MMSB
        uint64_t *partial_a = new uint64_t[size * n_bits];
        for(uint32_t i = 0; i < size; i++) {
            for (uint32_t j = 0; j < n_bits; j++) {
                partial_a[(i * n_bits) + j] = new_a[i] << (L_BIT - n_bits + j);
            }
        }

        // get secret shared form of the bits of the values that could contribute into the result
        uint64_t *bit_shares = MSB(proxy, partial_a, size * n_bits);

        // selection of the contribution of the bits of the value
        uint64_t *contributions = MUX(proxy, one_contributions, selected_e_contributions, bit_shares, size * n_bits);

        // binary-tree-based multiplication of the contributions into the exponential
        int cs = n_bits;
        bool flag = false;
        uint64_t* remaining = new uint64_t[size];
        uint64_t *tmp1, *tmp2;
        for (int j = 0; j < (int) ceil(log2(n_bits)); j++) {
            tmp1 = contributions;
            tmp2 = &contributions[cs / 2];

            if (cs % 2 == 1) {
                if (!flag) {
                    for(uint32_t i = 0; i < size; i++){
                        remaining[i] = contributions[(i * cs) + cs - 1];
                    }

                    tmp1 = new uint64_t[size * (cs / 2)];
                    tmp2 = new uint64_t[size * (cs / 2)];

                    for( uint32_t i = 0; i < size; i++) {
                        copy(contributions + (i * cs), contributions + (i * cs) + (cs / 2), tmp1 + (i * (cs / 2)));
                        copy(contributions + (i * cs) + (cs / 2), contributions + (i * cs) + 2 * (cs / 2), tmp2 + (i * (cs / 2)));
                    }

                    flag = true;
                } else {
                    tmp1 = new uint64_t[size * ((cs + 1) / 2)];
                    tmp2 = new uint64_t[size * ((cs + 1) / 2)];

                    size_t partial_size = cs / 2;

                    for(uint32_t i = 0; i < size; i++) {
                        copy(contributions + (i * cs), contributions + (i * cs) + ((cs + 1) / 2), tmp1 + (i * ((cs + 1) / 2)));
                        copy(contributions + (i * cs) + ((cs + 1) / 2), contributions + ((i + 1) * cs), tmp2 + (i * ((cs + 1) / 2)));
                        tmp2[(i + 1) * ((cs + 1) / 2) - 1] = remaining[i];
                    }

                    cs++;
                    flag = false;
                }
            }
            else {
                tmp1 = new uint64_t[size * (cs / 2)];
                tmp2 = new uint64_t[size * (cs / 2)];

                for( uint32_t i = 0; i < size; i++) {
                    copy(contributions + (i * cs), contributions + (i * cs) + (cs / 2), tmp1 + (i * (cs / 2)));
                    copy(contributions + (i * cs) + (cs / 2), contributions + ((i + 1) * cs), tmp2 + (i * (cs / 2)));
                }
            }
            delete [] contributions;
            contributions = MUL(proxy, tmp1, tmp2, size * (cs / 2));

            delete [] tmp1;
            delete [] tmp2;

            cs /= 2;
        }

        // deleting dynamically allocated arrays
        delete [] msb_a;
        delete [] abs_a;
        delete [] pec;
        delete [] nec;
        delete [] pos_e_contributions;
        delete [] neg_e_contributions;
        delete [] one_contributions;
        delete [] repeated_msb_a;
        delete [] e_contributions;
        delete [] partial_a;
        delete [] bit_shares;
        delete [] remaining;

        return contributions;
    }
    else if ( p_role == HELPER) {
        MSB(proxy, 0, size);
        MUX(proxy, 0, 0, 0, size * (n_bits + 1));
        MSB(proxy, 0, size * n_bits);
        MUX(proxy, 0, 0, 0, size * n_bits);

        int current_size = n_bits;
        bool flag = false;
        for (int i = 0; i < (int) ceil(log2(n_bits)); i++) {
            if (current_size % 2 == 1) {
                if (!flag) {
                    flag = true;
                } else {
                    current_size++;
                    flag = false;
                }
            }
            MUL(proxy, 0, 0, size * (current_size / 2));
            current_size /= 2;
        }

        return 0;
    }
    else {
        return nullptr;
    }
}

uint64_t* PSUM(Party* proxy, uint64_t *a, uint32_t size, uint32_t d) {
    // PartialSum: sum the elements of each section separately
    // Input(s)
    // x: the vector on which we perform the partial summation
    // size: the length of the vector x
    // d: the size of the part that we will use to partition x
    int p_role = proxy->getPRole();
    if(p_role == P1 || p_role == P2) {
        uint64_t *ps_x = new uint64_t[size / d];
        for(uint32_t base = 0; base < size; base += d) {
            uint64_t tmp = 0;
            for(uint32_t i = 0; i < d; i++) {
                tmp += a[base + i];
            }
            ps_x[base / d] = tmp;
        }
        return ps_x;
    }
    else {
        return NULL;
    }
}

uint64_t DP(Party* proxy, uint64_t *a, uint64_t *b, uint32_t size) {
    // This function computes the dot product of two single arithmetically shared vectors.
    // Input(s)
    // x and y: vectors of the given size
    // size: the lenght of the vectors
    // Return(s)
    // The dot product of x and y as a single value
    int p_role = proxy->getPRole();
    if(p_role == P1 || p_role == P2) {
        // compute elementwise multiplication of vectors
        uint64_t *ew_xy = MUL(proxy, a, b, size);
        // sum the result of the multiplications to obtain the dot product
        uint64_t res = 0;
        for(uint32_t i = 0; i < size; i++) {
            res += ew_xy[i];
        }
        delete [] ew_xy;

        return res;
    }
    else if(p_role == HELPER) {
        MUL(proxy, 0, 0, size);
        return 0;
    }
    else {
        return -1;
    }
}

uint64_t* DP(Party* proxy, uint64_t *a, uint64_t *b, uint32_t size, uint32_t d) {
    // This function computes the dot product of arithmetically shared vectors, which are formed by vectors.
    // Input(s)
    // x and y: vectors formed by vectors of given size
    // size: the lenght of the vectors
    // d: the size of the partial vectors forming the main vectors
    // Return(s)
    // Dot product of the given (size/d) vectors as a vector of (size / d)
    int p_role = proxy->getPRole();
    if(p_role == P1 || p_role == P2) {
        // compute elementwise multiplication of vectors
        uint64_t *ew_xy = MUL(proxy, a, b, size);
        // sum the vectors in the main vector
        uint64_t *dp_shr = PSUM(proxy, ew_xy, size, d);

        delete [] ew_xy;

        return dp_shr;
    }
    else if(p_role == HELPER) {
        MUL(proxy, 0, 0, size);
        return NULL;
    }
    else {
        return NULL;
    }

}

uint64_t** MATMATMUL(Party* proxy, uint64_t **a, uint64_t **b, uint32_t a_row, uint32_t a_col, uint32_t b_col) {
    /*
     * Perform multiplication of matrices a and b. The function assumes that the number of columns of a equals to
     * the number of rows of b.
     *
     * Input(s)
     * a: two dimensional matrix of size a_row-by-a_col
     * b: two dimensional matrix of size a_col-by-b_col
     *
     * Output(s)
     * Returns a matrix of size a_row-by-b_col
     */
    int p_role = proxy->getPRole();
    if (p_role == P1 || p_role == P2) {
        // form a single vector for each matrices such that all required multiplications can be performed in one go
        uint32_t size = a_row * a_col * b_col;
        uint64_t *concat_a = new uint64_t[size];
        uint64_t *concat_b = new uint64_t[size];

        for (uint32_t i = 0; i < size; i++) {
            concat_a[i] = a[i / (a_col * b_col)][i % a_col];
            concat_b[i] = b[i % a_col][(i % (a_col * b_col)) / a_col];
        }

        uint64_t *tmp = MUL(proxy, concat_a, concat_b, size);

        // recover the resulting matrix
        uint64_t **res = new uint64_t *[a_row];
        uint32_t ind = 0;
        uint64_t tmp_sum;
        for (uint32_t i = 0; i < a_row; i++) {
            res[i] = new uint64_t[b_col];
            for (uint32_t j = 0; j < b_col; j++) {
                tmp_sum = 0;
                for (uint32_t k = ind; k < ind + a_col; k++) {
                    tmp_sum += tmp[k];
                }
                ind += a_col;
                res[i][j] = tmp_sum;
            }
        }

        delete[] concat_a;
        delete[] concat_b;
        delete[] tmp;

        return res;
    }
    else if( p_role == HELPER) {
        // note that a_row is the required size of the multiplication that will be performed in MATMATMUL
        MUL(proxy, NULL, NULL, a_row);
        return NULL;
    }
    else {
        return NULL;
    }
}

uint64_t*** MATMATMUL(Party* proxy, uint64_t*** a, uint64_t*** b, uint32_t n_matrices, uint32_t a_row, uint32_t a_col, uint32_t b_col) {
    /*
     * Perform several multiplications of matrices of size a_row-by-a_col and a_col-by-b_col stored in a and b.
     *
     * Input(s)
     * a: three dimensional matrix of size n_matrices-by-a_row-by-a_col
     * b: three dimensional matrix of size n_matrices-by-a_col-by-b_col
     *
     * Output(s)
     * Returns a matrix of size n_matrices-by-a_row-by-b_col
     */
    int p_role = proxy->getPRole();
    if (p_role == P1 || p_role == P2) {
        // form a single vector for each matrices such that all required multiplications can be performed in one go
        uint32_t size = n_matrices * a_row * a_col * b_col;
        uint32_t size2 = a_row * a_col * b_col;
        uint64_t *concat_a = new uint64_t[size];
        uint64_t *concat_b = new uint64_t[size];

        for(uint32_t n = 0; n < n_matrices; n++) {
            for (uint32_t i = 0; i < size2; i++) {
                concat_a[size2 * n + i] = a[n][i / (a_col * b_col)][i % a_col];
                concat_b[size2 * n + i] = b[n][i % a_col][(i % (a_col * b_col)) / a_col];
            }
        }

        uint64_t *tmp = MUL(proxy, concat_a, concat_b, size);

        // recover the resulting matrix
        uint64_t ***res = new uint64_t **[n_matrices];
        uint32_t ind = 0;
        uint64_t tmp_sum;
        for(uint32_t n = 0; n < n_matrices; n++) {
            res[n] = new uint64_t*[a_row];
            for (uint32_t i = 0; i < a_row; i++) {
                res[n][i] = new uint64_t[b_col];
                for (uint32_t j = 0; j < b_col; j++) {
                    tmp_sum = 0;
                    for (uint32_t k = ind; k < ind + a_col; k++) {
                        tmp_sum += tmp[k];
                    }
                    ind += a_col;
                    res[n][i][j] = tmp_sum;
                }
            }
        }

        delete[] concat_a;
        delete[] concat_b;
        delete[] tmp;

        return res;
    }
    else if( p_role == HELPER) {
        if(a_row <= 0) {
            throw invalid_argument("core::MATMATMUL-Helper: The given size is " + to_string(a_row) + ". It has to be positive integer.");
        }
        // note that a_row is the required size of the multiplication that will be performed in MATMATMUL
        MUL(proxy, 0, 0, a_row);
        return NULL;
    }
    else {
        return nullptr;
    }
}

uint64_t* MATVECMUL(Party* proxy, uint64_t **a, uint64_t *b, uint32_t a_row, uint32_t a_col) {
    /*
     * Perform multiplication of matrice a and vector b. The function assumes that the number of columns of a equals to
     * the length of b.
     *
     * Input(s)
     * a: two dimensional matrix of size a_row-by-a_col
     * b: a vector of size a_col
     *
     * Output(s)
     * Returns a vector of length a_row
     */
    int p_role = proxy->getPRole();
    if (p_role == P1 || p_role == P2) {
        // form a single vector for each matrices such that all required multiplications can be performed in one go
        uint32_t size = a_row * a_col;
        uint64_t *concat_a = new uint64_t[size];
        uint64_t *concat_b = new uint64_t[size];

        for (uint32_t i = 0; i < size; i++) {
            concat_a[i] = a[i / a_col][i % a_col];
            concat_b[i] = b[i % a_col];
        }

        uint64_t *tmp = MUL(proxy, concat_a, concat_b, size);

        // recover the resulting vector
        uint64_t *res = new uint64_t [a_row];
        uint32_t ind = 0;
        uint64_t tmp_sum;
        for (uint32_t i = 0; i < a_row; i++) {
            tmp_sum = 0;
            for (uint32_t k = ind; k < ind + a_col; k++) {
                tmp_sum += tmp[k];
            }
            ind += a_col;
            res[i] = tmp_sum;
        }

        delete[] concat_a;
        delete[] concat_b;
        delete[] tmp;

        return res;
    }
    else if(p_role == HELPER) {
        // note that a_row is the required size of the multiplication that will be performed in MATVECMUL
        return MUL(proxy, NULL, NULL, a_row);
    }
    else {
        return nullptr;
    }
}

uint64_t** MATVECMUL(Party* proxy, uint64_t ***a, uint64_t **b, uint32_t n_matrices, uint32_t a_row, uint32_t a_col) {
    /*
     * Perform n_matrices multiplications of matrices of size a_row-by-a_col and vectors of size a_col stored in a and
     * b, respectively.
     *
     * Input(s)
     * a: three dimensional matrix of size n_matrices-by-a_row-by-a_col
     * b: two dimensional matrix of size n_matrices-by-a_col
     *
     * Output(s)
     * Returns a two-dimensional matrix of size n_matrices-by-a_row
     */
    int p_role = proxy->getPRole();
    if (p_role == P1 || p_role == P2) {
        // form a single vector for each matrices such that all required multiplications can be performed in one go
        uint32_t size = n_matrices * a_row * a_col;
        uint32_t size2 = a_row * a_col;
        uint64_t *concat_a = new uint64_t[size];
        uint64_t *concat_b = new uint64_t[size];

        for(uint32_t n = 0; n < n_matrices; n++) {
            for (uint32_t i = 0; i < size2; i++) {
                concat_a[size2 * n + i] = a[n][i / a_col][i % a_col];
                concat_b[size2 * n + i] = b[n][i % a_col];
            }
        }

        uint64_t *tmp = MUL(proxy, concat_a, concat_b, size);

        // recover the resulting vector
        uint64_t **res = new uint64_t*[n_matrices];
        uint32_t ind = 0;
        uint64_t tmp_sum;
        for(uint32_t n = 0; n < n_matrices; n++) {
            res[n] = new uint64_t[a_row];
            for (uint32_t i = 0; i < a_row; i++) {
                tmp_sum = 0;
                for (uint32_t k = ind; k < ind + a_col; k++) {
                    tmp_sum += tmp[k];
                }
                ind += a_col;
                res[n][i] = tmp_sum;
            }
        }

        delete[] concat_a;
        delete[] concat_b;
        delete[] tmp;

        return res;
    }
    else if( p_role == HELPER) {
        // note that a_row is the required size of the multiplication that will be performed in MATVECMUL
        MUL(proxy, NULL, NULL, a_row);
        return NULL;
    }
    else {
        return nullptr;
    }
}

uint64_t DIV_NN(Party* proxy, uint64_t a, uint64_t b) {
    uint32_t l = 64;
    if (proxy->getPRole()  == P1 || proxy->getPRole()  == P2){
        uint64_t zero = 0;
        uint64_t u = proxy->createShare(0);
        uint64_t z;
        for (uint8_t i = 0; i < l; i++){
            zero |= 1UL << i;                 // set the i-th bit
            bool uBit = (u >> i) & 1U;
            uint64_t commonShare_w = proxy->createShare(0);

            bool bitToSet = a - uBit - zero * b + commonShare_w;
            z ^= (-bitToSet ^ z) & (1UL << i);   // set the i-th bit of z
            zero &= -(1UL << i);                 // clear the i-th bit
        }


        return z;

    }
    else if (proxy->getPRole() == HELPER) {
        thread thr1 = thread(Receive,proxy->getSocketP1(), proxy->getBuffer1(), 16);
        thread thr2 = thread(Receive,proxy->getSocketP2(), proxy->getBuffer2(), 16);
        thr1.join();
        thr2.join();
// Receive(socket_p1, buffer, 16);
// Receive(socket_p2, buffer2, 16);
        unsigned char *ptr = &proxy->getBuffer1()[0];
        unsigned char *ptr2 = &proxy->getBuffer2()[0];


        uint64_t a1, a2, b1, b2;
        a1 = convert2Long(&ptr);
        b1 = convert2Long(&ptr);
        a2 = convert2Long(&ptr2);
        b2 = convert2Long(&ptr2);
        uint64_t a = a1+a2;
        uint64_t b = b1+b2;


        uint64_t c = ((a*1.0)/b)*PRECISION;

        uint64_t c1 = proxy->generateRandom();
        uint64_t c2 = c - c1;

        ptr = &proxy->getBuffer1()[0];
        addVal2CharArray(c1, &ptr);
        thr1 = thread(Send,proxy->getSocketP1(), proxy->getBuffer1(), 8);
//Send(socket_p1, buffer, 8);
        ptr2 = &proxy->getBuffer2()[0];
        addVal2CharArray(c2, &ptr2);
        thr2 = thread(Send, proxy->getSocketP2(), proxy->getBuffer2(), 8);
        thr1.join();
        thr2.join();
//Send(socket_p2, buffer2, 8);

        return 0;

    }
}

uint64_t DIV(Party* proxy, uint64_t a, uint64_t b) {
    role p_role = proxy->getPRole();
    if (p_role == HELPER) {
        thread thr1 = thread(Receive, proxy->getSocketP1(), proxy->getBuffer1(), 16);
        thread thr2 = thread(Receive, proxy->getSocketP2(), proxy->getBuffer2(), 16);
        thr1.join();
        thr2.join();
// Receive(socket_p1, buffer, 16);
// Receive(socket_p2, buffer2, 16);
        unsigned char *ptr = &proxy->getBuffer1()[0];
        unsigned char *ptr2 = &proxy->getBuffer2()[0];


        uint64_t a1, a2, b1, b2;
        a1 = convert2Long(&ptr);
        b1 = convert2Long(&ptr);
        a2 = convert2Long(&ptr2);
        b2 = convert2Long(&ptr2);
        uint64_t a = a1 + a2;
        uint64_t b = b1 + b2;


        uint64_t c = ((a * 1.0) / b) * PRECISION;

        uint64_t c1 = proxy->generateRandom();
        uint64_t c2 = c - c1;

        ptr = &proxy->getBuffer1()[0];
        addVal2CharArray(c1, &ptr);
        thr1 = thread(Send, proxy->getSocketP1(), proxy->getBuffer1(), 8);
//Send(socket_p1, buffer, 8);
        ptr2 = &proxy->getBuffer2()[0];
        addVal2CharArray(c2, &ptr2);
        thr2 = thread(Send, proxy->getSocketP2(), proxy->getBuffer2(), 8);
        thr1.join();
        thr2.join();
//Send(socket_p2, buffer2, 8);

        return 0;

    } else if (p_role == P1 || p_role == P2) {
        uint64_t d = proxy->generateCommonRandom() & MAX_SAMPLE_MASK;
        uint64_t c = proxy->generateCommonRandom() % d;

//da+cb/db = a/b + c/d
        uint64_t x = d * a + c * b;
        uint64_t y = d * b;

        unsigned char *ptr = &proxy->getBuffer1()[0];
        uint8_t sendingBytesNumber = 16;
        if (p_role == P1) {
            addVal2CharArray((uint8_t) CORE_DIV, &ptr);
            sendingBytesNumber++; // increase because also the command byte is sent.
        }
        addVal2CharArray(x, &ptr);
        addVal2CharArray(y, &ptr);
        Send(proxy->getSocketHelper(), proxy->getBuffer1(), sendingBytesNumber);
        Receive(proxy->getSocketHelper(), proxy->getBuffer1(), 8);
        ptr = &proxy->getBuffer1()[0];
        uint64_t z = convert2Long(&ptr);
        if (p_role == P1) {
            uint64_t res = (((c * 1.0) / d) * PRECISION);
            z = z - res;
        }
        cout << "z = " << REC(proxy, z) << endl;
        return z;
    } else {
        return -1;
    }
}

#endif //CORE_H


