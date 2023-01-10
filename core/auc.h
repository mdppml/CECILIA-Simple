//
// Created by Mete Akgun on 28.12.21.
//

#ifndef PPAUC_AUC_H
#define PPAUC_AUC_H

#include "core.h"
#include "../utils/auc_utils.h"

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
        cout << "helper MRound" << endl;
        Receive(proxy->getSocketP1(), proxy->getBuffer1(), sz * 8);
        unsigned char *ptr = proxy->getBuffer1();
        unsigned char *ptr_out = proxy->getBuffer1();

        Receive(proxy->getSocketP2(), proxy->getBuffer2(), sz * 8);
        unsigned char *ptr2 = proxy->getBuffer2();
        unsigned char *ptr_out2 = proxy->getBuffer2();
        for (uint32_t i = 0; i < sz; i++) {
            uint64_t v = convert2Long(&ptr) + convert2Long(&ptr2);
//            cout << "v: " << v << endl;
            if (v != 0)
                v = convert2uint64(1);
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

uint64_t ROCNOTIE(Party *proxy, client_data *c_data, uint32_t size) {
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
        cout << "Calling ROCNOTIE..." << endl;
        MUL(proxy, 0, 0, size);
        cout << "First MUL is over." << endl;
        MUL(proxy, 0, 0);
        cout << "Second MUL is over." << endl;
        NORM(proxy, 0, 0, 1);
        cout << "Returning from ROCNOTIE..." << endl;
    }
    return -1;
}

uint64_t ROCWITHTIE(Party *proxy, client_data *c_data, int size) {
    if(proxy->getPRole() == P1 || proxy->getPRole() == P2) {

        // calculate confidence
        uint64_t *preds = new uint64_t[size];
        //uint64_t* dummy_indexes = new uint64_t [size/5];
        //int i=0;
        //int j=0;
//        cout << "==================" << endl;
        for (int k = 0; k < size; k++) {
            if (k != (size - 1))
                preds[k] = c_data[0][k].val - c_data[0][k + 1].val; // get the difference of the current and next prediction
            else
                preds[k] = convert2uint64(2); // the last prediction must be in the list.
//            cout << convert2double(REC(proxy, preds[k])) << endl;

            preds[k] = preds[k] * (proxy->generateCommonRandom() & MAXSCALAR); // multiply with random value
            //i++;
            /*if ((proxy->generateCommonRandom()%10) == 0){ // add dummy predictions
                if ((proxy->generateCommonRandom()%2) == 0) {
                    preds[i] = proxy->generateRandom(); // random dummy
                }else{
                    if (proxy->GetRole()) // zero dummy
                        preds[i] = proxy->generateCommonRandom();
                    else
                        preds[i] = 0 - proxy->generateCommonRandom();

                    dummy_indexes[j++] = i++;
                }
            }*/
        }
//        cout << "==================" << endl;

        //get a permutation of preds
        uint64_t *rpreds = MRound(proxy, preds, size);
        delete[] preds;

        //get the reverse permutation of rpreds

        /*j = 0;
        int l = 0;
        for(int k=0;k<i:k++){ // remove dummy values
            if (k != dummy_indexes[j])
                c_data[0][l++].val = rpreds[k];
            else
                j++;
        }*/
        for (int k = 0; k < size; k++) { // remove dummy values
            c_data[0][k].val = rpreds[k];
        }
        delete[] rpreds;

//        for (int i = 0; i < 1; i++) {
//            cout << "Station : " << i << endl;
//            for (prediction n: c_data[i]) {
//                cout << convert2double(REC(proxy, n.val)) << "\t" << convert2double(REC(proxy, n.label)) << endl;
//            }
//            cout << "[";
//            for (prediction n: c_data[i]) {
//                cout << convert2double(REC(proxy, n.label)) << ", ";
//            }
//            cout << "]" << endl;
//        }


        // calculate AUC
        uint64_t *labels = new uint64_t[size];
        int ind = 0;

        for (prediction n: c_data[0]) {
            labels[ind++] = n.label;
        }

        uint64_t TP = 0;
        uint64_t FP = 0;
        uint64_t pre_FP = 0;
        uint64_t pre_TP = 0;
        uint64_t numerator = 0;
        uint64_t numerator2 = 0;
        uint64_t *areas;
        uint64_t *parts;
        for (int i = 0; i < size; i++) {
            TP = ADD(proxy, TP, labels[i]);
            FP = proxy->getPRole() * convert2uint64(i) - TP;

//            uint64_t area = MUL(proxy, pre_TP, FP - pre_FP);
//
//            area = MUL(proxy, area, c_data[0][i].val);
//            numerator = ADD(proxy, numerator, area);
//
//            area = MUL(proxy, TP - pre_TP, FP - pre_FP);
//
//            area = MUL(proxy, area, c_data[0][i].val);
//            numerator2 = ADD(proxy, numerator2, area);

            areas = MUL(proxy, (uint64_t[]) {pre_TP, TP - pre_TP}, (uint64_t[]) {FP - pre_FP, FP - pre_FP}, 2);
            parts = MUL(proxy, (uint64_t[]) {c_data[0][i].val, c_data[0][i].val}, areas, 2);
            numerator = ADD(proxy, numerator, parts[0]);
            numerator2 = ADD(proxy, numerator2, parts[1]);

            uint64_t *tmp = MUX(proxy, (uint64_t[]) {pre_FP, pre_TP}, (uint64_t[]) {FP, TP},
                                (uint64_t[]) {c_data[0][i].val, c_data[0][i].val}, 2);
            pre_FP = tmp[0];
            pre_TP = tmp[1];
            //pre_TP = proxy->SelectShare(pre_TP,TP,c_data[0][i].val);

            delete [] areas;
            delete [] parts;
        }

        uint64_t FN = TP;
        uint64_t TN = proxy->getPRole() * convert2uint64(size) - TP;

        uint64_t denominator = MUL(proxy, FN, TN);
        denominator = denominator * 2;

        numerator = 2 * numerator + numerator2;

    //    uint64_t auc = proxy->DIVISION(numerator, denominator);
        uint64_t auc = NORM(proxy, (uint64_t[]) {numerator}, (uint64_t[]) {denominator}, 1)[0];

        delete[] labels;

        return auc;
    }
    else if(proxy->getPRole() == HELPER) {
        MRound(proxy, 0, size);

        for (int i = 0; i < size; i++) {
//            MUL(proxy, 0, 0);
//            MUL(proxy, 0, 0);
//            MUL(proxy, 0, 0);
//            MUL(proxy, 0, 0);
            MUL(proxy, 0, 0, 2);
            MUL(proxy, 0, 0, 2);
            MUX(proxy, 0, 0, 0, 2);
        }
//        cout << "check!!" << endl;

        MUL(proxy, 0, 0);
        NORM(proxy, 0, 0, 1);

//        cout << "check!!!!!" << endl;
    }
    return -1;
}

uint64_t PRCURVE(Party *proxy, client_data *c_data, int size) {
    if (proxy->getPRole() == P1 || proxy->getPRole() == P2) {
        // calculate the confidence
        uint64_t *preds = new uint64_t[size];
        //uint64_t* dummy_indexes = new uint64_t [size/5];
        //int i=0;
        //int j=0;
//        cout << "==================" << endl;
        for (int k = 0; k < size; k++) {
            if (k != (size - 1))
                preds[k] = c_data[0][k].val - c_data[0][k + 1].val; // get the difference of the current and next prediction
            else
                preds[k] = convert2uint64(2); // the last prediction must be in the list.
//            cout << convert2double(REC(proxy, preds[k])) << endl;

            preds[k] = preds[k] * (proxy->generateCommonRandom() & MAXSCALAR); // multiply with random value
            //i++;
            /*if ((proxy->generateCommonRandom()%10) == 0){ // add dummy predictions
                if ((proxy->generateCommonRandom()%2) == 0) {
                    preds[i] = proxy->generateRandom(); // random dummy
                }else{
                    if (proxy->GetRole()) // zero dummy
                        preds[i] = proxy->generateCommonRandom();
                    else
                        preds[i] = 0 - proxy->generateCommonRandom();

                    dummy_indexes[j++] = i++;
                }
            }*/
        }
//        cout << "==================" << endl;

        //get a permutation of preds
//        cout << "check 1 -- size: " << size << endl;
        uint64_t *rpreds = MRound(proxy, preds, size);
//        cout << "check 2" << endl;
        delete[] preds;

        //get the reverse permutation of rpreds

        /*j = 0;
        int l = 0;
        for(int k=0;k<i:k++){ // remove dummy values
            if (k != dummy_indexes[j])
                c_data[0][l++].val = rpreds[k];
            else
                j++;
        }*/
        for (int k = 0; k < size; k++) { // remove dummy values
            c_data[0][k].val = rpreds[k];
        }
        delete[] rpreds;

//        for (int i = 0; i < 1; i++) {
//            cout << "Station : " << i << endl;
//            for (prediction n: c_data[i]) {
//                cout << convert2double(REC(proxy, n.val)) << "\t" << convert2double(REC(proxy, n.label)) << endl;
//            }
//            cout << "[";
//            for (prediction n: c_data[i]) {
//                cout << convert2double(REC(proxy, n.label)) << ", ";
//            }
//            cout << "]" << endl;
//        }


        // calculate the auc of the PR curve
        uint64_t *labels = new uint64_t[size];
        int j = 0;
        for (prediction n: c_data[0]) {
            labels[j++] = n.label;
        }

        uint64_t TP = 0;
//        uint64_t pre_prec = convert2uint64(50);
        uint64_t pre_prec = convert2uint64(0.5);
//        cout << "pre_prec: " << convert2double(REC(proxy, pre_prec)) << endl;
        uint64_t tmp = 0;
        uint64_t pre_reca = 0;
        uint64_t numerator = 0;
        uint64_t numerator2 = 0;
        uint64_t *precs = new uint64_t[size];
        uint64_t *recas = new uint64_t[size];
        for (int i = 0; i < size; i++) {
            TP = ADD(proxy, TP, labels[i]);
            recas[i] = TP;
            tmp = tmp + proxy->getPRole() * convert2uint64(1);
            precs[i] = tmp;
        }
        uint64_t prec;
        uint64_t reca;
//        uint64_t *areas;
        uint64_t *parts;
//        proxy->SendBytes(AUC_MDIV, params, 1);
//        precs = MDIVISION(proxy, recas, precs, size);

//        cout << "check 3" << endl;
        precs = NORM(proxy, recas, precs, size);
//        print1DArray("Precs", convert2double(REC(proxy, precs, size), size), size);
//        print1DArray("Recas", convert2double(REC(proxy, recas, size), size), size);
//        cout << "check 4" << endl;
        for (int i = 0; i < size; i++) {
            prec = precs[i];
            reca = recas[i];
//            proxy->SendBytes(CORE_MUL);
//            uint64_t area = MUL(proxy, pre_prec, reca - pre_reca);
//            proxy->SendBytes(CORE_MUL);
//            area = MUL(proxy, area, c_data[0][i].val);
//            numerator = ADD(proxy, numerator, area);
//
//            proxy->SendBytes(CORE_MUL);
//            area = MUL(proxy, reca - pre_reca, prec - pre_prec);
//            proxy->SendBytes(CORE_MUL);
//            area = MUL(proxy, area, c_data[0][i].val);
//            numerator2 = ADD(proxy, numerator2, area);

//            cout << "check 5.1." << i << endl;
//            cout << "reca: " << convert2double(REC(proxy, reca)) << endl;
//            cout << "pre_reca: " << convert2double(REC(proxy, pre_reca)) << endl;
//            cout << "prec: " << convert2double(REC(proxy, prec)) << endl;
//            cout << "pre_prec: " << convert2double(REC(proxy, pre_prec)) << endl;
//            cout << "reca - pre_reca: " << convert2double(REC(proxy, reca - pre_reca)) << endl;
//            cout << "prec - pre_prec: " << convert2double(REC(proxy, prec - pre_prec)) << endl;
//            cout << "c_data.val: " << convert2double(REC(proxy, c_data[0][i].val)) << endl;
            uint64_t *areas = MUL(proxy, (uint64_t[]) {pre_prec, reca - pre_reca}, (uint64_t[]) {reca - pre_reca, prec - pre_prec}, 2);
//            print1DArray("MUL 1: areas", convert2double(REC(proxy, areas, 2), 2), 2);
//            cout << "at " << i << endl;
            parts = MUL(proxy, (uint64_t[]) {c_data[0][i].val, c_data[0][i].val}, areas, 2);
//            print1DArray("MUL 2: parts", convert2double(REC(proxy, parts, 2), 2), 2);
//            cout << "kafa " << i << endl;
            numerator = ADD(proxy, numerator, parts[0]);
            numerator2 = ADD(proxy, numerator2, parts[1]);

//            cout << "check 5.2." << i << endl;
            uint64_t *tmp2 = MUX(proxy, (uint64_t[]) {pre_prec, pre_reca}, (uint64_t[]) {prec, reca},
                                 (uint64_t[]) {c_data[0][i].val, c_data[0][i].val}, 2);
//            print1DArray("MUX result", convert2double(REC(proxy, tmp2, 2), 2), 2);
            pre_prec = tmp2[0];
            pre_reca = tmp2[1];

//            cout << "check 5.3." << i << endl;

//            delete [] areas;
//            delete [] parts;
//            delete [] tmp2;
        }

        numerator = 2 * numerator + numerator2;
        uint64_t denominator = 2 * TP;

//        proxy->SendBytes(AUC_TDIV);
//        uint64_t prc = DIVISION(proxy, numerator, denominator);
//        cout << "check 6" << endl;
        uint64_t prc = NORM(proxy, (uint64_t []) {numerator}, (uint64_t []) {denominator}, 1)[0];
//        cout << "check 7" << endl;

        delete [] precs;
        delete [] recas;
        delete[] labels;

        return prc;
    }
    else if (proxy->getPRole() == HELPER) {
//        cout << "check 1 - size: " << size << endl;
        MRound(proxy, 0, size);
//        cout << "check 2" << endl;
        NORM(proxy, 0, 0, size);
//        cout << "check 3" << endl;
        for (int i = 0; i < size; i++) {
//            cout << "check 4.1." << i << endl;
            MUL(proxy, 0, 0, 2);
            MUL(proxy, 0, 0, 2);
//            cout << "check 4.2." << i << endl;
            MUX(proxy, 0, 0, 0, 2);
//            cout << "check 4.3." << i << endl;
        }
//        cout << "check 5" << endl;
        NORM(proxy, 0, 0, 1);
//        cout << "check 6" << endl;
    }
    return -1;
}

void SORT(Party *proxy, client_data *c_data, int nstation, int delta) {
    if (proxy->getPRole() == P1 || proxy->getPRole() == P2) {
        int tmp_delta = delta;
        while (nstation != 1) {
//        cout << "nstation: " << nstation << endl;
            int i = 0;
            int ns = nstation - (nstation % 2);
            while (i < ns) {
//            cout << "i: " << i << endl;
                delta = tmp_delta;
                int fl_index = i;
                int ll_index = i + 1;
                if (c_data[i].size() < c_data[i + 1].size()) {
                    fl_index = i + 1;
                    ll_index = i;
                }
                uint64_t *diff = new uint64_t[2 * c_data[ll_index].size()];
                uint64_t *mux_val1 = new uint64_t[2 * c_data[ll_index].size()];
                uint64_t *mux_val2 = new uint64_t[2 * c_data[ll_index].size()];
                uint64_t *mux_res = new uint64_t[2 * c_data[ll_index].size()];
                client_data sorted;

                uint64_t tmp_sample_size[100];
                client_data data[2];
                // sort
                int cnt = 0;
                bool flag = false;
                while (!c_data[fl_index].empty() && !c_data[ll_index].empty()) {
//                cout << "fl_index: " << fl_index << "\tll_index: " << ll_index << endl;
//                cout << "***************************************" << endl;
//                tmp_sample_size[0] = c_data[fl_index].size();
//                tmp_sample_size[1] = c_data[ll_index].size();
//                data[0] = c_data[fl_index];
//                data[1] = c_data[ll_index];
//                print_data(proxy, 2, tmp_sample_size, data);
//                cout << "-----------------------------------" << endl;

//                cout << "Count: " << cnt << endl;
                    cnt++;
                    int diff_size = c_data[fl_index].size();
                    if (diff_size > c_data[ll_index].size()) {
                        diff_size = c_data[ll_index].size();
                        flag = true;
                    }
//                cout << "check 1" << endl;

                    for (int j = 0; j < diff_size; j++) {
                        diff[j] = c_data[fl_index][j].val - c_data[ll_index][j].val;
                        mux_val1[j] = c_data[fl_index][j].val;
                        mux_val2[j] = c_data[ll_index][j].val;
                        mux_val1[j + diff_size] = c_data[fl_index][j].label;
                        mux_val2[j + diff_size] = c_data[ll_index][j].label;
                    }
//                cout << "check 2" << endl;

                    uint32_t params[1] = {(uint32_t) diff_size};
                    proxy->SendBytes(CORE_MMSB, params, 1);
                    uint64_t *diff_res = MSB(proxy, diff, diff_size);
//                cout << "check 3" << endl;

                    for (int j = 0; j < diff_size; j++) {
                        diff[j] = diff_res[j];
                        diff[j + diff_size] = diff_res[j];
                    }

//                cout << "check 4" << endl;
//                mux_res = proxy->MSelectShare(mux_val1, mux_val2, diff, 2 * diff_size);
                    params[0] = 2 * diff_size;
                    proxy->SendBytes(CORE_MMUX, params, 1);
                    mux_res = MUX(proxy, mux_val1, mux_val2, diff, 2 * diff_size);
                    for (int j = 0; j < diff_size; j++) {
                        c_data[fl_index][j].val = mux_res[j];
                        c_data[fl_index][j].label = mux_res[j + diff_size];
                    }
//                cout << "check 5" << endl;

                    params[0] = 2 * diff_size;
                    proxy->SendBytes(CORE_MMUX, params, 1);
                    mux_res = MUX(proxy, mux_val2, mux_val1, diff, 2 * diff_size);
//                cout << "check 5.5" << endl;
                    for (int j = 0; j < diff_size; j++) {
                        c_data[ll_index][j].val = mux_res[j];
                        c_data[ll_index][j].label = mux_res[j + diff_size];
                    }

//                cout << "Moving...";
//                if(flag) {
//                    int n_iter = c_data[ll_index].size() - c_data[fl_index].size();
//                    for(int t = 0; t < n_iter; t++) {
//                        c_data[fl_index].push_back(c_data[ll_index][diff_size]);
//                        c_data[ll_index].erase(c_data[ll_index].begin() + diff_size);
//                    }
//
////                    for(int t = 0; t < c_data[ll_index].size() - c_data[fl_index].size(); t++) {
////                        c_data[ll_index].pop_back();
////                    }
//                }
//                cout << " is done!" << endl;

                    sorted.push_back(c_data[fl_index][0]);
                    c_data[fl_index].pop_front();
//                cout << "check 6" << endl;


                    int min_val = min({delta, (int) c_data[fl_index].size(), (int) c_data[ll_index].size()});
                    int fl_pop = 1;
                    int ll_pop = 0;

//                cout << "########### min_val: " << min_val << endl;

                    for (int j = 0; j < min_val; j++) {
//                    cout << j << endl;
                        if (c_data[ll_index].size() == 1) {
                            delta = 0;
                            break;
                        }
                        if (fl_pop != ll_pop) {
                            diff[0] = c_data[fl_index][0].val - c_data[ll_index][0].val;
                            params[0] = 1;
                            proxy->SendBytes(CORE_MMSB, params, 1);
                            uint64_t *diff_res = MSB(proxy, diff, 1);
                            uint64_t cmp = REC(proxy, diff_res[0]);
                            if (cmp == 0) {
                                sorted.push_back(c_data[fl_index][0]);
                                c_data[fl_index].pop_front();
                                fl_pop++;
                            } else {
                                sorted.push_back(c_data[ll_index][0]);
                                c_data[ll_index].pop_front();
                                ll_pop++;
                            }
                        } else {
                            sorted.push_back(c_data[fl_index][0]);
                            c_data[fl_index].pop_front();
                            fl_pop++;
                        }

                    }

//                cout << "-----------------------------------" << endl;
//                cout << "After the loop" << endl;
//                tmp_sample_size[0] = c_data[fl_index].size();
//                tmp_sample_size[1] = c_data[ll_index].size();
//                data[0] = c_data[fl_index];
//                data[1] = c_data[ll_index];
//                print_data(proxy, 2, tmp_sample_size, data);

                    if(c_data[ll_index].size() > c_data[fl_index].size()) {
                        int tmp = ll_index;
                        ll_index = fl_index;
                        fl_index = tmp;
                    }

//                cout << "-----------------------------------" << endl;
//                cout << "After the exchanging" << endl;
//                tmp_sample_size[0] = c_data[fl_index].size();
//                tmp_sample_size[1] = c_data[ll_index].size();
//                data[0] = c_data[fl_index];
//                data[1] = c_data[ll_index];
//                print_data(proxy, 2, tmp_sample_size, data);
//                cout << "***************************************" << endl;

//                cout << "Over" << endl;

                }

//            cout << "Inner while loop is over" << endl;

                if (c_data[fl_index].size() > 0) {
//                cout << "if statement" << endl;
                    while (!c_data[fl_index].empty()) {
                        sorted.push_back(c_data[fl_index][0]);
                        c_data[fl_index].pop_front();
                    }
                } else {
//                cout << "else statement" << endl;
                    while (!c_data[ll_index].empty()) {
                        sorted.push_back(c_data[ll_index][0]);
                        c_data[ll_index].pop_front();
                    }
                }

//            cout << "if statement is over" << endl;

                // merge
                c_data[i] = sorted;
                c_data[i + 1].clear();

//            cout << "merging is over" << endl;

                // delete
                delete[] diff;
                delete[] mux_val1;
                delete[] mux_val2;
                delete[] mux_res;

                i += 2;

            }

//        cout << "inner while loop is over" << endl;

            uint64_t tmp_sample_size[100];
            tmp_sample_size[0] = c_data[0].size();
            i = 1;
            c_data[1].clear();
            while (i < nstation) {
//            cout << "last while loop: " << i << endl;
                if (i % 2 == 0) {
                    tmp_sample_size[i / 2] = c_data[i].size();
                    c_data[i / 2] = c_data[i];
                }
                c_data[i].clear();
                i += 1;
            }
//        cout << "last while loop is over" << endl;
            nstation = (nstation / 2) + (nstation % 2);

//        cout << "=========================================" << endl;
//        print_data(proxy, nstation, tmp_sample_size);
//        cout << "=========================================" << endl;

        }
    }
    else if(proxy->getPRole() == HELPER) {

    }
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