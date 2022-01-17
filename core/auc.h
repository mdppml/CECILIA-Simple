//
// Created by Mete Akgun on 28.12.21.
//

#ifndef PPAUC_AUC_H
#define PPAUC_AUC_H

#include "core.h"



uint64_t* AUCMSB(Party* proxy, uint64_t *x, uint32_t sz) {
    if (proxy->getPRole() == P1 || proxy->getPRole() == P2) {
        uint64_t *z_1 = new uint64_t[sz];
        uint64_t *ya = new uint64_t[sz];
        uint8_t *yb = new uint8_t[sz * (L - 1)];
        uint8_t *w = new uint8_t[sz];
        Receive(proxy->getSocketHelper(), proxy->getBuffer1(), sz * (8 + L));
        unsigned char *ptr = proxy->getBuffer1();
        for (int i = 0; i < sz; i++) {
            ya[i] = convert2Long(&ptr);
            convert2Array(&ptr, &yb[i * (L - 1)], L - 1);
            w[i] = (*ptr);
            ptr++;
            z_1[i] = ((x[i] & N1_MASK) + ya[i]) & N1_MASK;
        }
        uint64_t *z = REC(proxy,z_1, sz,N1_MASK);
        uint8_t *wc = PCB(proxy,z, yb, sz, L - 1);

        for (int i = 0; i < sz; i++) {
            w[i] = w[i] ^ wc[i];
            if (proxy->getPRole() == P1 && z_1[i] > z[i])
                z_1[i] = z_1[i] + N1;
            z_1[i] = x[i] - (z_1[i] - (ya[i] + w[i] * N1));
        }

        delete[] ya;
        delete[] yb;
        delete[] w;
        return z_1;

    } else if (proxy->getPRole() == HELPER) {

        unsigned char *ptr_out = proxy->getBuffer1();
        unsigned char *ptr_out2 = proxy->getBuffer2();
        for (int i = 0; i < sz; i++) {
            uint64_t y = proxy->generateRandom() & N1_MASK;
            uint64_t ya_1 = proxy->generateRandom() & N1_MASK;
            uint64_t ya_2 = (y - ya_1) & N1_MASK;
            addVal2CharArray(ya_1, &ptr_out);
            addVal2CharArray(ya_2, &ptr_out2);
            for (int j = 0; j < L - 1; j++) {
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

        thread thr1 = thread(Send,proxy->getSocketP1(), proxy->getBuffer1(), sz * (8 + L));
        thread thr2 = thread(Send,proxy->getSocketP2(), proxy->getBuffer2(), sz * (8 + L));
        thr1.join();
        thr2.join();
        //Send(proxy->getSocketP1(), proxy->getBuffer1(), sz * (8 + L));
        //Send(proxy->getSocketP2(), proxy->getBuffer2(), sz * (8 + L));
        // P1 and P2 will call PrivateCompareBool
        uint8_t *tmp = PCB(proxy,0, 0, sz, L - 1);
        return NULL;

    }
}

uint64_t* MRound(Party* proxy, uint64_t *a, uint32_t sz) {
    if (proxy->getPRole() == P1 || proxy->getPRole() == P2) {
        uint64_t *b = new uint64_t[sz];
        unsigned char *ptr = proxy->getBuffer1();
        for (int i = 0; i < sz; i++) {
            addVal2CharArray(a[i], &ptr);
        }
        Send(proxy->getSocketHelper(), proxy->getBuffer1(), sz * 8);
        Receive(proxy->getSocketHelper(), proxy->getBuffer1(), sz * 8);
        ptr = proxy->getBuffer1();
        for (int i = 0; i < sz; i++) {
            b[i] = convert2Long(&ptr);
        }
        return b;
    } else if (proxy->getPRole() == HELPER) {
        Receive(proxy->getSocketP1(), proxy->getBuffer1(), sz * 8);
        unsigned char *ptr = proxy->getBuffer1();
        unsigned char *ptr_out = proxy->getBuffer1();

        Receive(proxy->getSocketP2(), proxy->getBuffer2(), sz * 8);
        unsigned char *ptr2 = proxy->getBuffer2();
        unsigned char *ptr_out2 = proxy->getBuffer2();
        for (uint32_t i = 0; i < sz; i++) {
            uint64_t v = convert2Long(&ptr) + convert2Long(&ptr2);
            if (v != 0)
                v = 1;
            uint64_t v1 = proxy->generateRandom();
            uint64_t v2 = v - v1;
            addVal2CharArray(v1, &ptr_out);
            addVal2CharArray(v2, &ptr_out2);
        }
        Send(proxy->getSocketP1(), proxy->getBuffer1(), sz * 8);
        Send(proxy->getSocketP2(), proxy->getBuffer2(), sz * 8);
        return NULL;
    }
}

uint64_t* MDIVISION(Party* proxy, uint64_t *a, uint64_t *b, uint32_t sz) {
    if (proxy->getPRole() == HELPER) {
        thread thr1 = thread(Receive, proxy->getSocketP1(), proxy->getBuffer1(), sz * 16);
        thread thr2 = thread(Receive, proxy->getSocketP2(), proxy->getBuffer2(), sz * 16);
        thr1.join();
        thr2.join();
        // Receive(proxy->getSocketP1(), proxy->getBuffer1(), 16);
        // Receive(proxy->getSocketP2(), proxy->getBuffer2(), 16);
        unsigned char *ptr = proxy->getBuffer1();
        unsigned char *ptr2 = proxy->getBuffer2();


        unsigned char *ptr_out = proxy->getBuffer1();
        unsigned char *ptr_out2 = proxy->getBuffer2();
        for (int i = 0; i < sz; i++) {
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
            addVal2CharArray(c1, &ptr_out);
            addVal2CharArray(c2, &ptr_out2);
        }

        thr1 = thread(Send, proxy->getSocketP1(), proxy->getBuffer1(), 8 * sz);
        thr2 = thread(Send, proxy->getSocketP2(), proxy->getBuffer2(), 8 * sz);
        thr1.join();
        thr2.join();

        return 0;

    } else if (proxy->getPRole() == P1) {

        uint64_t *d = new uint64_t[sz];
        uint64_t *c = new uint64_t[sz];
        unsigned char *ptr = proxy->getBuffer1();
        for (int i = 0; i < sz; i++) {
            d[i] = proxy->generateCommonRandom() & MAX_SAMPLE_MASK;
            c[i] = proxy->generateCommonRandom() % d[i];

            uint64_t x = d[i] * a[i] + c[i] * b[i];
            uint64_t y = d[i] * b[i];
            addVal2CharArray(x, &ptr);
            addVal2CharArray(y, &ptr);
        }
        Send(proxy->getSocketHelper(), proxy->getBuffer1(), sz * 16);
        Receive(proxy->getSocketHelper(), proxy->getBuffer1(), sz * 8);
        ptr = proxy->getBuffer1();
        uint64_t *z = new uint64_t[sz];
        for (int i = 0; i < sz; i++) {
            z[i] = convert2Long(&ptr);
            uint64_t res = (((c[i] * 1.0) / d[i]) * PRECISION);
            z[i] = z[i] - res;
        }
        delete[]d;
        delete[]c;
        return z;

    } else if (proxy->getPRole() == P2) {
        uint64_t *d = new uint64_t[sz];
        uint64_t *c = new uint64_t[sz];
        unsigned char *ptr = proxy->getBuffer1();
        for (int i = 0; i < sz; i++) {
            d[i] = proxy->generateCommonRandom() & MAX_SAMPLE_MASK;
            c[i] = proxy->generateCommonRandom() % d[i];

            uint64_t x = d[i] * a[i] + c[i] * b[i];
            uint64_t y = d[i] * b[i];
            addVal2CharArray(x, &ptr);
            addVal2CharArray(y, &ptr);
        }
        Send(proxy->getSocketHelper(), proxy->getBuffer1(), sz * 16);
        Receive(proxy->getSocketHelper(), proxy->getBuffer1(), sz * 8);
        ptr = proxy->getBuffer1();
        uint64_t *z = new uint64_t[sz];
        for (int i = 0; i < sz; i++) {
            z[i] = convert2Long(&ptr);
        }
        delete[]d;
        delete[]c;
        return z;
    }
}


uint64_t DIVISION(Party* proxy, uint64_t a, uint64_t b) {

    if (proxy->getPRole() == HELPER) {
        thread thr1 = thread(Receive,proxy->getSocketP1(), proxy->getBuffer1(), 16);
        thread thr2 = thread(Receive,proxy->getSocketP2(), proxy->getBuffer2(), 16);
        thr1.join();
        thr2.join();
        // Receive(proxy->getSocketP1(), proxy->getBuffer1(), 16);
        // Receive(proxy->getSocketP2(), proxy->getBuffer2(), 16);
        unsigned char *ptr = proxy->getBuffer1();
        unsigned char *ptr2 = proxy->getBuffer2();


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

        ptr = proxy->getBuffer1();
        addVal2CharArray(c1, &ptr);
        thr1 = thread(Send,proxy->getSocketP1(), proxy->getBuffer1(), 8);
        //Send(proxy->getSocketP1(), proxy->getBuffer1(), 8);
        ptr2 = proxy->getBuffer2();
        addVal2CharArray(c2, &ptr2);
        thr2 = thread(Send,proxy->getSocketP2(), proxy->getBuffer2(), 8);
        thr1.join();
        thr2.join();
        //Send(proxy->getSocketP2(), proxy->getBuffer2(), 8);

        return 0;

    } else if (proxy->getPRole() == P1){
        uint64_t d = proxy->generateCommonRandom() & MAX_SAMPLE_MASK;
        uint64_t c = proxy->generateCommonRandom() % d;

        //da+cb/db = a/b + c/d

        uint64_t x = d*a+c*b;
        uint64_t y = d*b;

        unsigned char *ptr = proxy->getBuffer1();
        addVal2CharArray(x, &ptr);
        addVal2CharArray(y, &ptr);
        Send(proxy->getSocketHelper(), proxy->getBuffer1(), 16);
        Receive(proxy->getSocketHelper(), proxy->getBuffer1(), 8);
        ptr = proxy->getBuffer1();
        uint64_t z = convert2Long(&ptr);
        if (proxy->getPRole() == P1) {
            uint64_t res = (((c * 1.0) / d) * PRECISION);
            z = z - res;
        }
        return z;

    } else if (proxy->getPRole() == P2) {
        uint64_t d = proxy->generateCommonRandom() & MAX_SAMPLE_MASK;
        uint64_t c = proxy->generateCommonRandom() % d;

        //da+cb/db = a/b + c/d

        uint64_t x = d*a+c*b;
        uint64_t y = d*b;

        unsigned char *ptr = proxy->getBuffer1();
        addVal2CharArray(x, &ptr);
        addVal2CharArray(y, &ptr);
        Send(proxy->getSocketHelper(), proxy->getBuffer1(), 16);
        Receive(proxy->getSocketHelper(), proxy->getBuffer1(), 8);
        ptr = proxy->getBuffer1();
        uint64_t z = convert2Long(&ptr);
        return z;
    }
}

#endif //PPAUC_AUC_H