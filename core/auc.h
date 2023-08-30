//
// Created by Mete Akgun on 28.12.21.
//

#ifndef PPAUC_AUC_H
#define PPAUC_AUC_H

#include "core.h"
#include "../utils/auc_utils.h"

uint64_t *AucMostSignificantBit(Party *proxy, uint64_t *x, uint32_t sz) {
    if (proxy->GetPRole() == proxy1 || proxy->GetPRole() == proxy2) {
        uint64_t *z_1 = new uint64_t[sz];
        uint64_t *ya = new uint64_t[sz];
        uint8_t *yb = new uint8_t[sz * (L_BIT - 1)];
        uint8_t *w = new uint8_t[sz];
        Receive(proxy->GetSocketHelper(), proxy->GetBuffer1(), sz * (8 + L_BIT));
        unsigned char *ptr = proxy->GetBuffer1();
        for (int i = 0; i < sz; i++) {
            ya[i] = ConvertToLong(&ptr);
            ConvertToArray(&ptr, &yb[i * (L_BIT - 1)], L_BIT - 1);
            w[i] = (*ptr);
            ptr++;
            z_1[i] = ((x[i] & N1_MASK) + ya[i]) & N1_MASK;
        }
        uint64_t *z = Reconstruct(proxy, z_1, sz, N1_MASK);
        uint8_t *wc = PrivateCompareBool(proxy, z, yb, sz, L_BIT - 1);

        for (int i = 0; i < sz; i++) {
            w[i] = w[i] ^ wc[i];
            if (proxy->GetPRole() == proxy1 && z_1[i] > z[i])
                z_1[i] = z_1[i] + N1;
            z_1[i] = x[i] - (z_1[i] - (ya[i] + w[i] * N1));
        }

        delete[] ya;
        delete[] yb;
        delete[] w;
        return z_1;

    } else if (proxy->GetPRole() == helper) {

        unsigned char *ptr_out = proxy->GetBuffer1();
        unsigned char *ptr_out2 = proxy->GetBuffer2();
        for (int i = 0; i < sz; i++) {
            uint64_t y = proxy->GenerateRandom() & N1_MASK;
            uint64_t ya_1 = proxy->GenerateRandom() & N1_MASK;
            uint64_t ya_2 = (y - ya_1) & N1_MASK;
            AddValueToCharArray(ya_1, &ptr_out);
            AddValueToCharArray(ya_2, &ptr_out2);
            for (int j = 0; j < L_BIT - 1; j++) {
                uint8_t k = (y >> j) & 0x1;
                uint8_t yb_1 = proxy->GenerateRandom() % LP;
                uint8_t yb_2 = Mod(k - yb_1, LP);
                AddValueToCharArray(yb_1, &ptr_out);
                AddValueToCharArray(yb_2, &ptr_out2);
            }
            uint8_t w = 0;
            if (ya_1 > y)
                w = 1;
            uint8_t w_1 = proxy->GenerateRandom() % 2;
            uint8_t w_2 = w ^ w_1;
            AddValueToCharArray(w_1, &ptr_out);
            AddValueToCharArray(w_2, &ptr_out2);
        }

        thread thr1 = thread(Send, proxy->GetSocketP1(), proxy->GetBuffer1(), sz * (8 + L_BIT));
        thread thr2 = thread(Send, proxy->GetSocketP2(), proxy->GetBuffer2(), sz * (8 + L_BIT));
        thr1.join();
        thr2.join();
        uint8_t *tmp = PrivateCompareBool(proxy, 0, 0, sz, L_BIT - 1);
        return NULL;

    }
    return nullptr;
}

uint64_t *Round(Party *proxy, uint64_t *a, uint32_t sz) {
    if (proxy->GetPRole() == proxy1 || proxy->GetPRole() == proxy2) {
        uint64_t *b = new uint64_t[sz];
        unsigned char *ptr = proxy->GetBuffer1();
        for (int i = 0; i < sz; i++) {
            AddValueToCharArray(a[i], &ptr);
        }
        Send(proxy->GetSocketHelper(), proxy->GetBuffer1(), sz * 8);
        Receive(proxy->GetSocketHelper(), proxy->GetBuffer1(), sz * 8);
        ptr = proxy->GetBuffer1();
        for (int i = 0; i < sz; i++) {
            b[i] = ConvertToLong(&ptr);
        }
        return b;
    } else if (proxy->GetPRole() == helper) {
        cout << "helper Round" << endl;
        thread thr1 = thread(Receive, proxy->GetSocketP1(), proxy->GetBuffer1(), sz * 8);
        thread thr2 = thread(Receive, proxy->GetSocketP2(), proxy->GetBuffer2(), sz * 8);
        thr1.join();
        thr2.join();
        unsigned char *ptr = proxy->GetBuffer1();
        unsigned char *ptr_out = proxy->GetBuffer1();
        unsigned char *ptr2 = proxy->GetBuffer2();
        unsigned char *ptr_out2 = proxy->GetBuffer2();
        for (uint32_t i = 0; i < sz; i++) {
            uint64_t v = ConvertToLong(&ptr) + ConvertToLong(&ptr2);
            if (v != 0)
                v = ConvertToUint64(1);
            uint64_t v1 = proxy->GenerateRandom();
            uint64_t v2 = v - v1;
            AddValueToCharArray(v1, &ptr_out);
            AddValueToCharArray(v2, &ptr_out2);
        }
        thr1 = thread(Send, proxy->GetSocketP1(), proxy->GetBuffer1(), sz * 8);
        thr2 = thread(Send, proxy->GetSocketP2(), proxy->GetBuffer2(), sz * 8);
        thr1.join();
        thr2.join();
        return NULL;
    }
    return nullptr;
}

uint64_t *AucDivide(Party *proxy, uint64_t *a, uint64_t *b, uint32_t sz) {
    if (proxy->GetPRole() == helper) {
        thread thr1 = thread(Receive, proxy->GetSocketP1(), proxy->GetBuffer1(), sz * 16);
        thread thr2 = thread(Receive, proxy->GetSocketP2(), proxy->GetBuffer2(), sz * 16);
        thr1.join();
        thr2.join();
        unsigned char *ptr = proxy->GetBuffer1();
        unsigned char *ptr2 = proxy->GetBuffer2();


        unsigned char *ptr_out = proxy->GetBuffer1();
        unsigned char *ptr_out2 = proxy->GetBuffer2();
        for (int i = 0; i < sz; i++) {
            uint64_t a1, a2, b1, b2;
            a1 = ConvertToLong(&ptr);
            b1 = ConvertToLong(&ptr);
            a2 = ConvertToLong(&ptr2);
            b2 = ConvertToLong(&ptr2);
            uint64_t a = a1 + a2;
            uint64_t b = b1 + b2;
            uint64_t c = ((a * 1.0) / b) * PRECISION;

            uint64_t c1 = proxy->GenerateRandom();
            uint64_t c2 = c - c1;
            AddValueToCharArray(c1, &ptr_out);
            AddValueToCharArray(c2, &ptr_out2);
        }

        thr1 = thread(Send, proxy->GetSocketP1(), proxy->GetBuffer1(), 8 * sz);
        thr2 = thread(Send, proxy->GetSocketP2(), proxy->GetBuffer2(), 8 * sz);
        thr1.join();
        thr2.join();

        return 0;

    } else if (proxy->GetPRole() == proxy1) {

        uint64_t *d = new uint64_t[sz];
        uint64_t *c = new uint64_t[sz];
        unsigned char *ptr = proxy->GetBuffer1();
        for (int i = 0; i < sz; i++) {
            d[i] = proxy->GenerateCommonRandom() & MAX_SAMPLE_MASK;
            c[i] = proxy->GenerateCommonRandom() % d[i];

            uint64_t x = d[i] * a[i] + c[i] * b[i];
            uint64_t y = d[i] * b[i];
            AddValueToCharArray(x, &ptr);
            AddValueToCharArray(y, &ptr);
        }
        Send(proxy->GetSocketHelper(), proxy->GetBuffer1(), sz * 16);
        Receive(proxy->GetSocketHelper(), proxy->GetBuffer1(), sz * 8);
        ptr = proxy->GetBuffer1();
        uint64_t *z = new uint64_t[sz];
        for (int i = 0; i < sz; i++) {
            z[i] = ConvertToLong(&ptr);
            uint64_t res = (((c[i] * 1.0) / d[i]) * PRECISION);
            z[i] = z[i] - res;
        }
        delete[]d;
        delete[]c;
        return z;

    } else if (proxy->GetPRole() == proxy2) {
        uint64_t *d = new uint64_t[sz];
        uint64_t *c = new uint64_t[sz];
        unsigned char *ptr = proxy->GetBuffer1();
        for (int i = 0; i < sz; i++) {
            d[i] = proxy->GenerateCommonRandom() & MAX_SAMPLE_MASK;
            c[i] = proxy->GenerateCommonRandom() % d[i];

            uint64_t x = d[i] * a[i] + c[i] * b[i];
            uint64_t y = d[i] * b[i];
            AddValueToCharArray(x, &ptr);
            AddValueToCharArray(y, &ptr);
        }
        Send(proxy->GetSocketHelper(), proxy->GetBuffer1(), sz * 16);
        Receive(proxy->GetSocketHelper(), proxy->GetBuffer1(), sz * 8);
        ptr = proxy->GetBuffer1();
        uint64_t *z = new uint64_t[sz];
        for (int i = 0; i < sz; i++) {
            z[i] = ConvertToLong(&ptr);
        }
        delete[]d;
        delete[]c;
        return z;
    }
    return nullptr;
}

uint64_t AucDivide(Party *proxy, uint64_t a, uint64_t b) {

    if (proxy->GetPRole() == helper) {
        thread thr1 = thread(Receive, proxy->GetSocketP1(), proxy->GetBuffer1(), 16);
        thread thr2 = thread(Receive, proxy->GetSocketP2(), proxy->GetBuffer2(), 16);
        thr1.join();
        thr2.join();
        unsigned char *ptr = proxy->GetBuffer1();
        unsigned char *ptr2 = proxy->GetBuffer2();


        uint64_t a1, a2, b1, b2;
        a1 = ConvertToLong(&ptr);
        b1 = ConvertToLong(&ptr);
        a2 = ConvertToLong(&ptr2);
        b2 = ConvertToLong(&ptr2);
        uint64_t a = a1 + a2;
        uint64_t b = b1 + b2;


        uint64_t c = ((a * 1.0) / b) * PRECISION;

        uint64_t c1 = proxy->GenerateRandom();
        uint64_t c2 = c - c1;

        ptr = proxy->GetBuffer1();
        AddValueToCharArray(c1, &ptr);
        thr1 = thread(Send, proxy->GetSocketP1(), proxy->GetBuffer1(), 8);
        ptr2 = proxy->GetBuffer2();
        AddValueToCharArray(c2, &ptr2);
        thr2 = thread(Send, proxy->GetSocketP2(), proxy->GetBuffer2(), 8);
        thr1.join();
        thr2.join();

        return 0;

    } else if (proxy->GetPRole() == proxy1) {
        uint64_t d = proxy->GenerateCommonRandom() & MAX_SAMPLE_MASK;
        uint64_t c = proxy->GenerateCommonRandom() % d;
        uint64_t x = d * a + c * b;
        uint64_t y = d * b;
        unsigned char *ptr = proxy->GetBuffer1();
        AddValueToCharArray(x, &ptr);
        AddValueToCharArray(y, &ptr);
        Send(proxy->GetSocketHelper(), proxy->GetBuffer1(), 16);
        Receive(proxy->GetSocketHelper(), proxy->GetBuffer1(), 8);
        ptr = proxy->GetBuffer1();
        uint64_t z = ConvertToLong(&ptr);
        if (proxy->GetPRole() == proxy1) {
            uint64_t res = (((c * 1.0) / d) * PRECISION);
            z = z - res;
        }
        return z;

    } else if (proxy->GetPRole() == proxy2) {
        uint64_t d = proxy->GenerateCommonRandom() & MAX_SAMPLE_MASK;
        uint64_t c = proxy->GenerateCommonRandom() % d;
        uint64_t x = d * a + c * b;
        uint64_t y = d * b;

        unsigned char *ptr = proxy->GetBuffer1();
        AddValueToCharArray(x, &ptr);
        AddValueToCharArray(y, &ptr);
        Send(proxy->GetSocketHelper(), proxy->GetBuffer1(), 16);
        Receive(proxy->GetSocketHelper(), proxy->GetBuffer1(), 8);
        ptr = proxy->GetBuffer1();
        uint64_t z = ConvertToLong(&ptr);
        return z;
    }
    return -1;
}

uint64_t RocNoTie(Party *proxy, client_data *c_data, uint32_t size) {
    if(proxy->GetPRole() == proxy1 || proxy->GetPRole() == proxy2) {
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
            TP = Add(proxy, TP, labels[i]);
            FP = proxy->GetPRole() * ConvertToUint64(i) - TP;
            mul1[i] = TP;
            mul2[i] = FP - pre_FP;
            pre_FP = FP;
        }
        uint64_t *area = Multiply(proxy, mul1, mul2, size);
        for (int i = 0; i < size; i++) {
            numerator = Add(proxy, numerator, area[i]);
        }
        delete[] mul1;
        delete[] mul2;
        delete[] area;
        uint64_t FN = TP;
        uint64_t TN = proxy->CreateShare((uint64_t) 0);
        TN = proxy->GetPRole() * ConvertToUint64(size) - TP;
        uint64_t denominator = Multiply(proxy, FN, TN);
        uint64_t num[1] = {numerator};
        uint64_t den[1] = {denominator};
        uint64_t *auc = Normalise(proxy, num, den, 1);
        delete[] labels;
        return auc[0];
    }
    else if(proxy->GetPRole() == helper) {
        Multiply(proxy, 0, 0, size);
        Multiply(proxy, 0, 0);
        Normalise(proxy, 0, 0, 1);
    }
    return -1;
}

uint64_t RocWithTie(Party *proxy, client_data *c_data, int size) {
    if(proxy->GetPRole() == proxy1 || proxy->GetPRole() == proxy2) {

        // calculate confidence
        uint64_t *preds = new uint64_t[size];
        for (int k = 0; k < size; k++) {
            if (k != (size - 1))
                preds[k] = c_data[0][k].val - c_data[0][k + 1].val; // get the difference of the current and next prediction
            else
                preds[k] = ConvertToUint64(2); // the last prediction must be in the list.
            preds[k] = preds[k] * (proxy->GenerateCommonRandom() & MAX_SCALAR); // multiply with random value
        }
        //get a permutation of preds
        uint64_t *rpreds = Round(proxy, preds, size);
        delete[] preds;
        //get the reverse permutation of rpreds
        for (int k = 0; k < size; k++) { // remove dummy values
            c_data[0][k].val = rpreds[k];
        }
        delete[] rpreds;
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
            TP = Add(proxy, TP, labels[i]);
            FP = proxy->GetPRole() * ConvertToUint64(i) - TP;
            uint64_t a[2] = {pre_TP, TP - pre_TP};
            uint64_t b[2] = {FP - pre_FP, FP - pre_FP};
            areas = Multiply(proxy, a, b, 2);
            a[0] = a[1] = c_data[0][i].val;
            parts = Multiply(proxy, a, areas, 2);
            numerator = Add(proxy, numerator, parts[0]);
            numerator2 = Add(proxy, numerator2, parts[1]);
            a[0] = pre_FP;
            a[1] = pre_TP;
            uint64_t y[2] = {FP, TP};
            b[0] = b[1] = c_data[0][i].val;
            uint64_t *tmp = Multiplex(proxy, a, y,
                                      b, 2);
            pre_FP = tmp[0];
            pre_TP = tmp[1];
            delete [] areas;
            delete [] parts;
        }

        uint64_t FN = TP;
        uint64_t TN = proxy->GetPRole() * ConvertToUint64(size) - TP;

        uint64_t denominator = Multiply(proxy, FN, TN);
        denominator = denominator * 2;

        numerator = 2 * numerator + numerator2;
        uint64_t a[1] = {numerator};
        uint64_t b[1] = {denominator};
        uint64_t auc = Normalise(proxy, a, b, 1)[0];

        delete[] labels;

        return auc;
    }
    else if(proxy->GetPRole() == helper) {
        Round(proxy, 0, size);

        for (int i = 0; i < size; i++) {
            Multiply(proxy, 0, 0, 2);
            Multiply(proxy, 0, 0, 2);
            Multiplex(proxy, 0, 0, 0, 2);
        }
        Multiply(proxy, 0, 0);
        Normalise(proxy, 0, 0, 1);
    }
    return -1;
}

uint64_t PrCurve(Party *proxy, client_data *c_data, int size) {
    if (proxy->GetPRole() == proxy1 || proxy->GetPRole() == proxy2) {
        // calculate the confidence
        uint64_t *preds = new uint64_t[size];
        for (int k = 0; k < size; k++) {
            if (k != (size - 1))
                preds[k] = c_data[0][k].val - c_data[0][k + 1].val; // get the difference of the current and next prediction
            else
                preds[k] = ConvertToUint64(2); // the last prediction must be in the list.

            preds[k] = preds[k] * (proxy->GenerateCommonRandom() & MAX_SCALAR); // multiply with random value
        }
        //get a permutation of preds
        uint64_t *rpreds = Round(proxy, preds, size);
        delete[] preds;
        //get the reverse permutation of rpreds
        for (int k = 0; k < size; k++) { // remove dummy values
            c_data[0][k].val = rpreds[k];
        }
        delete[] rpreds;
        // calculate the auc of the PR curve
        uint64_t *labels = new uint64_t[size];
        int j = 0;
        for (prediction n: c_data[0]) {
            labels[j++] = n.label;
        }

        uint64_t TP = 0;
        uint64_t pre_prec = ConvertToUint64(0.5);
        uint64_t tmp = 0;
        uint64_t pre_reca = 0;
        uint64_t numerator = 0;
        uint64_t numerator2 = 0;
        uint64_t *precs = new uint64_t[size];
        uint64_t *recas = new uint64_t[size];
        for (int i = 0; i < size; i++) {
            TP = Add(proxy, TP, labels[i]);
            recas[i] = TP;
            tmp = tmp + proxy->GetPRole() * ConvertToUint64(1);
            precs[i] = tmp;
        }
        uint64_t prec;
        uint64_t reca;
        uint64_t *parts;
        precs = Normalise(proxy, recas, precs, size);
        for (int i = 0; i < size; i++) {
            prec = precs[i];
            reca = recas[i];
            uint64_t a[2] = {pre_prec, reca - pre_reca};
            uint64_t b[2] = {reca - pre_reca, prec - pre_prec};
            uint64_t *areas = Multiply(proxy, a, b, 2);
            a[0] = c_data[0][i].val;
            a[1] = c_data[0][i].val;
            parts = Multiply(proxy, a, areas, 2);
            numerator = Add(proxy, numerator, parts[0]);
            numerator2 = Add(proxy, numerator2, parts[1]);
            uint64_t x[2] = {pre_prec, pre_reca};
            uint64_t y[2] = {prec, reca};
            b[0] = c_data[0][i].val;
            b[1] = c_data[0][i].val;
            uint64_t *tmp2 = Multiplex(proxy, x, y, b, 2);
            pre_prec = tmp2[0];
            pre_reca = tmp2[1];
        }

        numerator = 2 * numerator + numerator2;
        uint64_t denominator = 2 * TP;
        uint64_t a[1] = {numerator};
        uint64_t b[1] = {denominator};
        uint64_t prc = Normalise(proxy, a, b, 1)[0];
        delete [] precs;
        delete [] recas;
        delete[] labels;

        return prc;
    }
    else if (proxy->GetPRole() == helper) {
        Round(proxy, 0, size);
        Normalise(proxy, 0, 0, size);
        for (int i = 0; i < size; i++) {
            Multiply(proxy, 0, 0, 2);
            Multiply(proxy, 0, 0, 2);
            Multiplex(proxy, 0, 0, 0, 2);
        }
        Normalise(proxy, 0, 0, 1);
    }
    return -1;
}

void Sort(Party *proxy, client_data *c_data, int nstation, int delta) {
    if (proxy->GetPRole() == proxy1 || proxy->GetPRole() == proxy2) {
        int cnt = 0;
        int tmp_delta = delta;
        while (nstation != 1) {
            int i = 0;
            int ns = nstation - (nstation % 2);
            while (i < ns) {
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
                bool flag = false;
                while (!c_data[fl_index].empty() && !c_data[ll_index].empty()) {
                    cnt++;
                    int diff_size = c_data[fl_index].size();
                    if (diff_size > c_data[ll_index].size()) {
                        diff_size = c_data[ll_index].size();
                        flag = true;
                    }

                    for (int j = 0; j < diff_size; j++) {
                        diff[j] = c_data[fl_index][j].val - c_data[ll_index][j].val;
                        mux_val1[j] = c_data[fl_index][j].val;
                        mux_val2[j] = c_data[ll_index][j].val;
                        mux_val1[j + diff_size] = c_data[fl_index][j].label;
                        mux_val2[j + diff_size] = c_data[ll_index][j].label;
                    }
                    uint32_t params[1] = {(uint32_t) diff_size};
                    proxy->SendBytes(coreVectorisedMostSignificantBit, params, 1);
                    cnt++;
                    uint64_t *diff_res = MostSignificantBit(proxy, diff, diff_size);
                    for (int j = 0; j < diff_size; j++) {
                        diff[j] = diff_res[j];
                        diff[j + diff_size] = diff_res[j];
                    }
                    params[0] = 2 * diff_size;
                    proxy->SendBytes(coreVectorisedMultiplex, params, 1);
                    cnt++;
                    mux_res = Multiplex(proxy, mux_val1, mux_val2, diff, 2 * diff_size);
                    for (int j = 0; j < diff_size; j++) {
                        c_data[fl_index][j].val = mux_res[j];
                        c_data[fl_index][j].label = mux_res[j + diff_size];
                    }
                    params[0] = 2 * diff_size;
                    proxy->SendBytes(coreVectorisedMultiplex, params, 1);
                    cnt++;
                    mux_res = Multiplex(proxy, mux_val2, mux_val1, diff, 2 * diff_size);
                    for (int j = 0; j < diff_size; j++) {
                        c_data[ll_index][j].val = mux_res[j];
                        c_data[ll_index][j].label = mux_res[j + diff_size];
                    }
                    sorted.push_back(c_data[fl_index][0]);
                    c_data[fl_index].pop_front();
                    int min_val = min({delta, (int) c_data[fl_index].size(), (int) c_data[ll_index].size()});
                    int fl_pop = 1;
                    int ll_pop = 0;
                    for (int j = 0; j < min_val; j++) {
                        if (c_data[ll_index].size() == 1) {
                            delta = 0;
                            break;
                        }
                        if (fl_pop != ll_pop) {
                            diff[0] = c_data[fl_index][0].val - c_data[ll_index][0].val;
                            params[0] = 1;
                            proxy->SendBytes(coreVectorisedMostSignificantBit, params, 1);
                            cnt++;
                            uint64_t *diff_res = MostSignificantBit(proxy, diff, 1);
                            uint64_t cmp = Reconstruct(proxy, diff_res[0]);
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
                    if(c_data[ll_index].size() > c_data[fl_index].size()) {
                        int tmp = ll_index;
                        ll_index = fl_index;
                        fl_index = tmp;
                    }
                }
                if (c_data[fl_index].size() > 0) {
                    while (!c_data[fl_index].empty()) {
                        sorted.push_back(c_data[fl_index][0]);
                        c_data[fl_index].pop_front();
                    }
                } else {
                    while (!c_data[ll_index].empty()) {
                        sorted.push_back(c_data[ll_index][0]);
                        c_data[ll_index].pop_front();
                    }
                }
                // merge
                c_data[i] = sorted;
                c_data[i + 1].clear();
                // delete
                delete[] diff;
                delete[] mux_val1;
                delete[] mux_val2;
                delete[] mux_res;

                i += 2;

            }
            uint64_t tmp_sample_size[100];
            tmp_sample_size[0] = c_data[0].size();
            i = 1;
            c_data[1].clear();
            while (i < nstation) {
                if (i % 2 == 0) {
                    tmp_sample_size[i / 2] = c_data[i].size();
                    c_data[i / 2] = c_data[i];
                }
                c_data[i].clear();
                i += 1;
            }
            nstation = (nstation / 2) + (nstation % 2);
        }
        cout<<"cnt "<<cnt<<endl;
    }
    else if(proxy->GetPRole() == helper) {

    }

}
#endif //PPAUC_AUC_H