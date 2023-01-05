//
// Created by Mete Akgun on 28.12.21.
//

#ifndef PPAUC_AUC_H
#define PPAUC_AUC_H

#include "core.h"

struct prediction {
    uint64_t val;
    uint64_t label;
};

typedef std::deque<prediction> client_data;

uint64_t *AUCMSB(Party *proxy, uint64_t *x, uint32_t sz) {
    if (proxy->getPRole() == P1 || proxy->getPRole() == P2) {
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
            z_1[i] = ((x[i] & N1_MASK) + ya[i]) & N1_MASK;
        }
        uint64_t *z = REC(proxy, z_1, sz, N1_MASK);
        uint8_t *wc = PCB(proxy, z, yb, sz, L_BIT - 1);

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
            uint8_t w_2 = w ^ w_1;
            addVal2CharArray(w_1, &ptr_out);
            addVal2CharArray(w_2, &ptr_out2);
        }

        thread thr1 = thread(Send, proxy->getSocketP1(), proxy->getBuffer1(), sz * (8 + L_BIT));
        thread thr2 = thread(Send, proxy->getSocketP2(), proxy->getBuffer2(), sz * (8 + L_BIT));
        thr1.join();
        thr2.join();
        //Send(proxy->getSocketP1(), proxy->getBuffer1(), sz * (8 + L_BIT));
        //Send(proxy->getSocketP2(), proxy->getBuffer2(), sz * (8 + L_BIT));
        // P1 and P2 wilL_BIT calL_BIT PrivateCompareBool
        uint8_t *tmp = PCB(proxy, 0, 0, sz, L_BIT - 1);
        return NULL;

    }
    return nullptr;
}

uint64_t *MRound(Party *proxy, uint64_t *a, uint32_t sz) {
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
    return nullptr;
}

uint64_t *MDIVISION(Party *proxy, uint64_t *a, uint64_t *b, uint32_t sz) {
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
    return nullptr;
}

uint64_t DIVISION(Party *proxy, uint64_t a, uint64_t b) {

    if (proxy->getPRole() == HELPER) {
        thread thr1 = thread(Receive, proxy->getSocketP1(), proxy->getBuffer1(), 16);
        thread thr2 = thread(Receive, proxy->getSocketP2(), proxy->getBuffer2(), 16);
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
        uint64_t a = a1 + a2;
        uint64_t b = b1 + b2;


        uint64_t c = ((a * 1.0) / b) * PRECISION;

        uint64_t c1 = proxy->generateRandom();
        uint64_t c2 = c - c1;

        ptr = proxy->getBuffer1();
        addVal2CharArray(c1, &ptr);
        thr1 = thread(Send, proxy->getSocketP1(), proxy->getBuffer1(), 8);
        //Send(proxy->getSocketP1(), proxy->getBuffer1(), 8);
        ptr2 = proxy->getBuffer2();
        addVal2CharArray(c2, &ptr2);
        thr2 = thread(Send, proxy->getSocketP2(), proxy->getBuffer2(), 8);
        thr1.join();
        thr2.join();
        //Send(proxy->getSocketP2(), proxy->getBuffer2(), 8);

        return 0;

    } else if (proxy->getPRole() == P1) {
        uint64_t d = proxy->generateCommonRandom() & MAX_SAMPLE_MASK;
        uint64_t c = proxy->generateCommonRandom() % d;

        //da+cb/db = a/b + c/d

        uint64_t x = d * a + c * b;
        uint64_t y = d * b;

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

        uint64_t x = d * a + c * b;
        uint64_t y = d * b;

        unsigned char *ptr = proxy->getBuffer1();
        addVal2CharArray(x, &ptr);
        addVal2CharArray(y, &ptr);
        Send(proxy->getSocketHelper(), proxy->getBuffer1(), 16);
        Receive(proxy->getSocketHelper(), proxy->getBuffer1(), 8);
        ptr = proxy->getBuffer1();
        uint64_t z = convert2Long(&ptr);
        return z;
    }
    return -1;
}

uint64_t AUCNOTIE(Party *proxy, client_data *c_data, uint32_t size) {
    if(proxy->getPRole() == P1 || proxy->getPRole() == P2) {
//        uint32_t size = c_data[0].size();
        cout << "Size: " << size << endl;
        uint64_t *labels = new uint64_t[size];
        int i = 0;
        for (prediction n: c_data[0]) {
            labels[i++] = n.label;
        }

        uint64_t TP = 0;
        uint64_t FP = 0;
        uint64_t pre_FP = 0;
        uint64_t numerator = 0;
        uint64_t *mul1 = new uint64_t[size];
        uint64_t *mul2 = new uint64_t[size];
        for (int i = 0; i < size; i++) {
            TP = ADD(proxy, TP, labels[i]);
            FP = proxy->getPRole() * convert2uint64(i)  - TP;
//        cout << i << "\tTP: " << convert2double(REC(proxy, TP)) << "\tFP: " << convert2double(REC(proxy, FP)) << endl;
//        if (proxy->getPRole() == P1)
//            FP = i - TP;
//        else
//            FP = 0 - TP;

            mul1[i] = TP;
            mul2[i] = FP - pre_FP;

            pre_FP = FP;
        }
//    cout << "TP: " << convert2double(REC(proxy, TP)) << "\tFP: " << convert2double(REC(proxy, FP)) << endl;


//    cout << "MUL is being called..." << endl;
        uint64_t *area = MUL(proxy, mul1, mul2, size);
//    cout << "MUL is over" << endl;

        for (int i = 0; i < size; i++) {
            numerator = ADD(proxy, numerator, area[i]);
        }
        delete[] mul1;
        delete[] mul2;
        delete[] area;

        uint64_t FN = TP;
        uint64_t TN = proxy->createShare( (uint64_t) 0);
        TN = proxy->getPRole() * convert2uint64(size) - TP;
//    uint64_t TN = 0;
//    if (proxy->getPRole() == P1)
//        TN = size - TP;
//    else
//        TN = 0 - TP;

//    cout << "MUL is being called... Again!" << endl;
//        cout << convert2double(REC(proxy, TN)) << "\t" << convert2double(REC(proxy, FN)) << endl;
        uint64_t denominator = MUL(proxy, FN, TN);

//    cout << "NORM is being called..." << endl;
        uint64_t num[1] = {numerator};
        uint64_t den[1] = {denominator};
//        cout << convert2double(REC(proxy, numerator)) << "\t" << convert2double(REC(proxy, denominator)) << endl;
        uint64_t *auc = NORM(proxy, num, den, 1);

        delete[] labels;

        return auc[0];
    }
    else if(proxy->getPRole() == HELPER) {
        cout << "Calling AUCNOTIE..." << endl;
        MUL(proxy, 0, 0, size);
        cout << "First MUL is over." << endl;
        MUL(proxy, 0, 0);
        cout << "Second MUL is over." << endl;
        NORM(proxy, 0, 0, 1);
        cout << "Returning from AUCNOTIE..." << endl;
    }
    return -1;
}

//uint64_t *SORT(Party *proxy, client_data *c_data, uint64_t size, int delta, int nstation) {
//    if (proxy->getPRole() == P1 || proxy->getPRole() == P2) {
//        int tmp_delta = delta;
//        while (nstation != 1) {
//            int i = 0;
//            int ns = nstation - (nstation % 2);
//            while (i < ns) {
//                delta = tmp_delta;
//                int fl_index = i;
//                int ll_index = i + 1;
//                if (c_data[i].size() < c_data[i + 1].size()) {
//                    fl_index = i + 1;
//                    ll_index = i;
//                }
//                uint64_t *diff = new uint64_t[2 * c_data[ll_index].size()];
//                uint64_t *mux_val1 = new uint64_t[2 * c_data[ll_index].size()];
//                uint64_t *mux_val2 = new uint64_t[2 * c_data[ll_index].size()];
//                uint64_t *mux_res = new uint64_t[2 * c_data[ll_index].size()];
//                client_data sorted;
//
//
//                // sort
//                while (!c_data[fl_index].empty() && !c_data[ll_index].empty()) {
//                    int diff_size = c_data[fl_index].size();
//                    if (diff_size > c_data[ll_index].size())
//                        diff_size = c_data[ll_index].size();
//
//                    for (int j = 0; j < diff_size; j++) {
//                        diff[j] = c_data[fl_index][j].val - c_data[ll_index][j].val;
//                        mux_val1[j] = c_data[fl_index][j].val;
//                        mux_val2[j] = c_data[ll_index][j].val;
//                        mux_val1[j + diff_size] = c_data[fl_index][j].label;
//                        mux_val2[j + diff_size] = c_data[ll_index][j].label;
//                    }
//
//
////                proxy->SendBytes(MMSB, diff_size);
//                    uint64_t *diff_res = MSB(proxy, diff, diff_size);
//
//                    for (int j = 0; j < diff_size; j++) {
//                        diff[j] = diff_res[j];
//                        diff[j + diff_size] = diff_res[j];
//                    }
//
//
////                    mux_res = proxy->MSelectShare(mux_val1,mux_val2,diff,2*diff_size);
//                    mux_res = MUX(proxy, mux_val1, mux_val2, diff, 2 * diff_size);
//                    for (int j = 0; j < diff_size; j++) {
//                        c_data[fl_index][j].val = mux_res[j];
//                        c_data[fl_index][j].label = mux_res[j + diff_size];
//                    }
//
////                    mux_res = proxy->MSelectShare(mux_val2,mux_val1,diff,2*diff_size);
//                    mux_res = MUX(proxy, mux_val2, mux_val1, diff, 2 * diff_size);
//                    for (int j = 0; j < diff_size; j++) {
//                        c_data[ll_index][j].val = mux_res[j];
//                        c_data[ll_index][j].label = mux_res[j + diff_size];
//                    }
//                    sorted.push_back(c_data[fl_index][0]);
//                    c_data[fl_index].pop_front();
//
//
//                    int min_val = min({delta, (int) c_data[fl_index].size(), (int) c_data[ll_index].size()});
//                    int fl_pop = 1;
//                    int ll_pop = 0;
//
//                    for (int j = 0; j < min_val; j++) {
//                        if (c_data[ll_index].size() == 1) {
//                            delta = 0;
//                            break;
//                        }
//                        if (fl_pop != ll_pop) {
//                            diff[0] = c_data[fl_index][0].val - c_data[ll_index][0].val;
////                            uint64_t *diff_res = proxy->MMSB2(diff, 1);
//                            uint64_t *diff_res = MSB(proxy, diff, 1);
//                            uint64_t cmp = REC(proxy, diff_res[0]);
//                            if (cmp == 0) {
//                                sorted.push_back(c_data[fl_index][0]);
//                                c_data[fl_index].pop_front();
//                                fl_pop++;
//                            } else {
//                                sorted.push_back(c_data[ll_index][0]);
//                                c_data[ll_index].pop_front();
//                                ll_pop++;
//                            }
//                        } else {
//                            sorted.push_back(c_data[fl_index][0]);
//                            c_data[fl_index].pop_front();
//                            fl_pop++;
//                        }
//
//                    }
//
//                }
//                if (c_data[fl_index].size() > 0) {
//                    while (!c_data[fl_index].empty()) {
//                        sorted.push_back(c_data[fl_index][0]);
//                        c_data[fl_index].pop_front();
//                    }
//                } else {
//                    while (!c_data[ll_index].empty()) {
//                        sorted.push_back(c_data[ll_index][0]);
//                        c_data[ll_index].pop_front();
//                    }
//                }
//
//                // merge
//                c_data[i] = sorted;
//                c_data[i + 1].clear();
//
//                // delete
//                delete[] diff;
//                delete[] mux_val1;
//                delete[] mux_val2;
//                delete[] mux_res;
//
//                i += 2;
//            }
//
//            i = 1;
//            c_data[1].clear();
//            while (i < nstation) {
//                if (i % 2 == 0)
//                    c_data[i / 2] = c_data[i];
//                c_data[i].clear();
//                i += 1;
//            }
//            nstation = (nstation / 2) + (nstation % 2);
//
//        }
//    } else if (proxy->getPRole() == HELPER) {
//        int tmp_delta = delta;
//        while (nstation != 1) {
//            int i = 0;
//            int ns = nstation - (nstation % 2);
//            while (i < ns) {
//                delta = tmp_delta;
//                int fl_index = i;
//                int ll_index = i + 1;
//                if (c_data[i].size() < c_data[i + 1].size()) {
//                    fl_index = i + 1;
//                    ll_index = i;
//                }
//                uint64_t *diff = new uint64_t[2 * c_data[ll_index].size()];
//                uint64_t *mux_val1 = new uint64_t[2 * c_data[ll_index].size()];
//                uint64_t *mux_val2 = new uint64_t[2 * c_data[ll_index].size()];
//                uint64_t *mux_res = new uint64_t[2 * c_data[ll_index].size()];
//                client_data sorted;
//
//
//                // sort
//                while (!c_data[fl_index].empty() && !c_data[ll_index].empty()) {
////                    int diff_size = c_data[fl_index].size();
////                    if (diff_size > c_data[ll_index].size())
////                        diff_size = c_data[ll_index].size();
//
//                    // I guess we have to send "diff_size" information in every iteration. There seems to be no way
//                    // for the helper to know it or calculate it without having c_data.
//                    int diff_size = proxy->ReadInt();
//
//                    MSB(proxy, 0, diff_size);
//                    MUX(proxy, 0, 0, 0, 2 * diff_size);
//                    MUX(proxy, 0, 0, 0, 2 * diff_size);
//
//                    int min_val = proxy->ReadInt();
//
////                    int min_val = min({delta, (int) c_data[fl_index].size(), (int) c_data[ll_index].size()});
////                    int fl_pop = 1;
////                    int ll_pop = 0;
//
//                    for (int j = 0; j < min_val; j++) {
//                        if (c_data[ll_index].size() == 1) {
//                            delta = 0;
//                            break;
//                        }
//                        if (fl_pop != ll_pop) {
//                            diff[0] = c_data[fl_index][0].val - c_data[ll_index][0].val;
////                            uint64_t *diff_res = proxy->MMSB2(diff, 1);
//                            uint64_t *diff_res = MSB(proxy, diff, 1);
//                            uint64_t cmp = REC(proxy, diff_res[0]);
//                            if (cmp == 0) {
//                                sorted.push_back(c_data[fl_index][0]);
//                                c_data[fl_index].pop_front();
//                                fl_pop++;
//                            } else {
//                                sorted.push_back(c_data[ll_index][0]);
//                                c_data[ll_index].pop_front();
//                                ll_pop++;
//                            }
//                        } else {
//                            sorted.push_back(c_data[fl_index][0]);
//                            c_data[fl_index].pop_front();
//                            fl_pop++;
//                        }
//
//                    }
//
//                }
//
//                i += 2;
//            }
//
//            i = 1;
//            c_data[1].clear();
//            while (i < nstation) {
//                if (i % 2 == 0)
//                    c_data[i / 2] = c_data[i];
//                c_data[i].clear();
//                i += 1;
//            }
//            nstation = (nstation / 2) + (nstation % 2);
//
//        }
//    }
//    return NULL;
//}

#endif //PPAUC_AUC_H