#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <deque>
#include <chrono>
#include <iomanip>
#include <assert.h>
//#include "../../core/Party.h"
#include "../../utils/parse_options.h"
#include "../../core/auc.h"
//#include "../../core/core.h"
#include <sys/stat.h>
#include <dirent.h>
#include <algorithm>

using namespace std;
int nstation;
uint64_t sample_size[100];

client_data *c_data;
//
//bool IsPathExist(const std::string &s) {
//    struct stat buf;
//    return (stat(s.c_str(), &buf) == 0);
//}
//
//void ReadDirectory(const std::string &name, vector<string> &v) {
//    DIR *dirp = opendir(name.c_str());
//    struct dirent *dp;
//    while ((dp = readdir(dirp)) != NULL) {
//        v.push_back(dp->d_name);
//    }
//    closedir(dirp);
//}
//
//void RandomData(Party *proxy, int nstation, uint64_t *sample_size) {
//    srand(100);
//    c_data = new client_data[nstation];
//    for (int i = 0; i < nstation; i++) {
//        uint64_t tmp[MAX_SAMPLE_SIZE];
//        for (int j = 0; j < sample_size[i]; j++) {
////            tmp[j] = (rand() % 10000) + 1;
//            tmp[j] = proxy->GenerateCommonRandom() & MAX_SCALAR;
//        }
//
//        SortValues(tmp, sample_size[i]);
//        for (int j = 0; j < sample_size[i]; j++) {
//            uint64_t l = proxy->CreateShare(ConvertToUint64(rand() % 2));
//            c_data[i].push_back({proxy->CreateShare(tmp[j]), l});
//
////            if (Role == 0)
////                c_data[i].push_back({(uint64_t) rand(), (uint64_t) rand()});
////            else
////                c_data[i].push_back({tmp[j] - rand(), l - rand(),});
//
//        }
//    }
//
//}
//
//// Do we need to update FileData function considering the changes that we made in the number format?
//void FileData(string path) {
//    vector<string> f_list;
//    ReadDirectory(path, f_list);
//    c_data = new client_data[f_list.size()];
//    int f_index = 0;
//    for (string file: f_list) {
//        if (file != "." && file != "..") {
//            file = path + "/" + file;
//            ifstream ip(file.c_str());
//            cout << file << endl;
//            if (ip.is_open()) {
//                cout << "ddd" << endl;
//                string tmp;
//                int s_size = 0;
//                while (ip.good()) {
//                    getline(ip, tmp, ',');
//                    if (tmp.empty())
//                        break;
//                    char *end;
//                    uint64_t label = strtoull(tmp.c_str(), &end, 10);
//                    getline(ip, tmp, '\n');
//                    uint64_t pred = strtoull(tmp.c_str(), &end, 10);
//                    s_size++;
//                }
//                sample_size[f_index] = s_size;
//                f_index++;
//            }
//            ip.close();
//        }
//    }
//    nstation = f_index;
//}
//
//void PrintData(int nstation, uint64_t *sample_size) {
//    for (int i = 0; i < nstation; i++) {
//        cout << "Station : " << i << endl;
//        for (prediction n: c_data[i]) {
//            cout << n.val << " " << n.label << endl;
//        }
//    }
//}
//
//void PrintData(Party *proxy, int nstation, uint64_t *sample_size, client_data *data = c_data) {
//    for (int i = 0; i < nstation; i++) {
//        cout << "Station : " << i << endl;
//        for (prediction n: data[i]) {
//            cout << ConvertToDouble(REC(proxy, n.val)) << "\t" << ConvertToDouble(REC(proxy, n.label)) << endl;
//        }
//        cout << "[";
//        for (prediction n: data[i]) {
//            cout << ConvertToDouble(REC(proxy, n.label)) << ", ";
//        }
//        cout << "]" << endl;
//    }
//}

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
        TP = ADD(proxy, TP, labels[i]);
        FP = proxy->GetPRole() * ConvertToUint64(i) - TP;

        mul1[i] = TP;
        mul2[i] = FP - pre_FP;

        pre_FP = FP;
    }
//    cout << "TP: " << ConvertToDouble(REC(proxy, TP)) << "\tFP: " << ConvertToDouble(REC(proxy, FP)) << endl;


//    cout << "MUL is being called..." << endl;
    uint32_t params[1] = {size};
    proxy->SendBytes(coreVectorisedMultiply, params, 1);
    uint64_t *area = MUL(proxy, mul1, mul2, size);
//    cout << "MUL is over" << endl;

    for (int i = 0; i < size; i++) {
        numerator = ADD(proxy, numerator, area[i]);
    }
    delete[] mul1;
    delete[] mul2;
    delete[] area;

    uint64_t FN = TP;
    uint64_t TN = proxy->CreateShare((uint64_t) 0);
    TN = proxy->GetPRole() * ConvertToUint64(size) - TP;

//    cout << "MUL is being called... Again!" << endl;
    cout << ConvertToDouble(REC(proxy, TN)) << "\t" << ConvertToDouble(REC(proxy, FN)) << endl;
    proxy->SendBytes(coreMultiply);
    uint64_t denominator = MUL(proxy, FN, TN);

//    cout << "NORM is being called..." << endl;
    uint64_t num[1] = {numerator};
    uint64_t den[1] = {denominator};
    cout << ConvertToDouble(REC(proxy, numerator)) << "\t" << ConvertToDouble(REC(proxy, denominator)) << endl;
    params[0] = 1;
    proxy->SendBytes(coreNormalise, params, 1);
    uint64_t *auc = NORM(proxy, num, den, 1);

//    cout << "REC is being called..." << endl;
    auc = REC(proxy, auc, 1);
    cout << "AUC :\t" << ConvertToDouble(auc[0]) << endl;
    delete[] labels;
}


void calc_auc_v2(Party *proxy) {
    uint32_t size = c_data[0].size();
    uint32_t params[1] = {size};
    proxy->SendBytes(aucRocNoTie, params, 1);
    uint64_t auc = RocNoTie(proxy, c_data, size);
    cout << "AUC :\t" << ConvertToDouble(REC(proxy, auc)) << endl;
}

void sort_data(Party *proxy, int nstation, int delta) {
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
//                PrintData(proxy, 2, tmp_sample_size, data);
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
                proxy->SendBytes(coreVectorisedMostSignificantBit, params, 1);
                uint64_t *diff_res = MSB(proxy, diff, diff_size);
//                cout << "check 3" << endl;

                for (int j = 0; j < diff_size; j++) {
                    diff[j] = diff_res[j];
                    diff[j + diff_size] = diff_res[j];
                }

//                cout << "check 4" << endl;
//                mux_res = proxy->MSelectShare(mux_val1, mux_val2, diff, 2 * diff_size);
                params[0] = 2 * diff_size;
                proxy->SendBytes(coreVectorisedMultiplex, params, 1);
                mux_res = MUX(proxy, mux_val1, mux_val2, diff, 2 * diff_size);
                for (int j = 0; j < diff_size; j++) {
                    c_data[fl_index][j].val = mux_res[j];
                    c_data[fl_index][j].label = mux_res[j + diff_size];
                }
//                cout << "check 5" << endl;

                params[0] = 2 * diff_size;
                proxy->SendBytes(coreVectorisedMultiplex, params, 1);
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
                        proxy->SendBytes(coreVectorisedMostSignificantBit, params, 1);
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
//                PrintData(proxy, 2, tmp_sample_size, data);

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
//                PrintData(proxy, 2, tmp_sample_size, data);
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
//        PrintData(proxy, nstation, tmp_sample_size);
//        cout << "=========================================" << endl;

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

//    cout << "=========== Data of each station ==================" << endl;
//    PrintData(proxy, nstation, sample_size, true);
//    cout << "================================" << endl;

    auto start = chrono::high_resolution_clock::now();

//    cout << "Sorting..." << endl;
//    sort_data(proxy, nstation, delta);
    Sort(proxy, c_data, nstation, delta);

//    cout << "Data after sorting: " << endl;
//    uint64_t total_n_samples = 0;
//    for(int i = 0; i < nstation; i++) {
//        total_n_samples += sample_size[i];
//    }
//    uint64_t tmp_size[1] = {total_n_samples};
//    PrintData(proxy, 1, tmp_size, true);

//    cout << "... is done!" << endl;
//    calc_auc(proxy);
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

    //PrintData(nstation,sample_size,proxy);

    proxy->SendBytes(coreEnd);
    proxy->PrintBytes();
    del();
    cout << "*****************************" << endl;
    return 0;
}
