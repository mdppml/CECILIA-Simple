#include <cstdlib>
#include <iostream>
#include <sstream>
#include <deque>
#include <chrono>
#include <iomanip>
#include <assert.h>
#include "../../utils/parse_options.h"
#include "../../core/auc.h"
#include <algorithm>

using namespace std;
int nstation;
uint64_t sample_size[100];

client_data *c_data;

void del() {
    delete[] c_data;
}

void calc_auc(Party *proxy) {
    uint32_t size = c_data[0].size();
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
    uint32_t params[1] = {size};
    proxy->SendBytes(coreVectorisedMultiply, params, 1);
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

    cout << ConvertToDouble(Reconstruct(proxy, TN)) << "\t" << ConvertToDouble(Reconstruct(proxy, FN)) << endl;
    proxy->SendBytes(coreMultiply);
    uint64_t denominator = Multiply(proxy, FN, TN);

    uint64_t num[1] = {numerator};
    uint64_t den[1] = {denominator};
    cout << ConvertToDouble(Reconstruct(proxy, numerator)) << "\t" << ConvertToDouble(Reconstruct(proxy, denominator)) << endl;
    params[0] = 1;
    proxy->SendBytes(coreNormalise, params, 1);
    uint64_t *auc = Normalise(proxy, num, den, 1);

    auc = Reconstruct(proxy, auc, 1);
    cout << "AUC :\t" << ConvertToDouble(auc[0]) << endl;
    delete[] labels;
}


void calc_auc_v2(Party *proxy) {
    uint32_t size = c_data[0].size();
    uint32_t params[1] = {size};
    proxy->SendBytes(aucRocNoTie, params, 1);
    uint64_t auc = RocNoTie(proxy, c_data, size);
    cout << "AUC :\t" << ConvertToDouble(Reconstruct(proxy, auc)) << endl;
}

void sort_data(Party *proxy, int nstation, int delta) {
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
            int cnt = 0;
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
                uint64_t *diff_res = MostSignificantBit(proxy, diff, diff_size);
                for (int j = 0; j < diff_size; j++) {
                    diff[j] = diff_res[j];
                    diff[j + diff_size] = diff_res[j];
                }
                params[0] = 2 * diff_size;
                proxy->SendBytes(coreVectorisedMultiplex, params, 1);
                mux_res = Multiplex(proxy, mux_val1, mux_val2, diff, 2 * diff_size);
                for (int j = 0; j < diff_size; j++) {
                    c_data[fl_index][j].val = mux_res[j];
                    c_data[fl_index][j].label = mux_res[j + diff_size];
                }
                params[0] = 2 * diff_size;
                proxy->SendBytes(coreVectorisedMultiplex, params, 1);
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
}

int main(int argc, char *argv[]) {
    uint8_t role = atoi(argv[1]);
    uint16_t cport = atoi(argv[2]);
    string caddress(argv[3]);
    uint16_t hport = atoi(argv[4]);
    string haddress(argv[5]);
    int delta = atoi(argv[6]);
    string ss(argv[7]);

    nstation = 2;

    bool file_flag = false;
    if (!IsPathExist(ss)) {
        if (ss != "") {
            int i = 0;
            cout << ss << endl;
            stringstream sss(ss);
            while (sss.good()) {
                string substr;
                getline(sss, substr, ',');
                if (i == 0) {
                    nstation = stoi(substr);
                    if (nstation <= 2 && nstation > 1000) {
                        cout << "Number of stations must be between 2 and 1000." << endl;
                    }
                } else {
                    sample_size[i - 1] = stoi(substr);
                    if (sample_size[i - 1] <= 0) {
                        cout << "Number of samples must be greater than 0" << endl;
                        exit(0);
                    }
                }
                i++;
                if ((i - 1) == (nstation))
                    break;
            }
            if ((i - 1) != (nstation)) {
                cout << "Missing sample size" << endl;
                exit(0);
            }
        } else {
            nstation = 20;
            for (int i = 0; i < nstation; i++)
                sample_size[i] = 20000;
        }
    } else {
        file_flag = true;
    }


    if (cport != 0) {
        assert(cport < 1 << (sizeof(uint16_t) * 8));
    }

    if (hport != 0) {
        assert(hport < 1 << (sizeof(uint16_t) * 8));
    }

    Party *proxy;
    if (role == 0)
        proxy = new Party(proxy1, hport, haddress, cport, caddress);
    else
        proxy = new Party(proxy2, hport, haddress, cport, caddress);


    // determine which type of data is going to be used -- real or synthetic
    if (file_flag) {
        cout << "real data" << endl;
        FileData(ss, c_data, nstation, sample_size);
    } else {
        cout << "random data" << endl;
        RandomData(proxy, c_data, nstation, sample_size);
    }

    cout << "Number of parties: " << nstation << endl;
    for(int i = 0; i < nstation; i++)
        cout << sample_size[i] << "\t" << endl;
    auto start = chrono::high_resolution_clock::now();
    Sort(proxy, c_data, nstation, delta);
    calc_auc_v2(proxy);

    ios_base::sync_with_stdio(false);
    auto end = chrono::high_resolution_clock::now();
    double time_taken =
            chrono::duration_cast<chrono::nanoseconds>(end - start).count();
    time_taken *= 1e-9;
    cout << "Time taken by program is : " << fixed
         << time_taken << setprecision(9);
    cout << " sec" << endl;

    proxy->PrintPaperFriendly(time_taken);

    proxy->SendBytes(coreEnd);
    del();
    cout << "*****************************" << endl;
    return 0;
}
