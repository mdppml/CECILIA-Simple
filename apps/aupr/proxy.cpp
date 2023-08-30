#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <deque>
#include <chrono>
#include <iomanip>
#include <assert.h>
#include "../../core/auc.h"
#include "../../utils/parse_options.h"
#include <sys/stat.h>
#include <dirent.h>
using namespace std;
int nstation;
uint64_t sample_size[1000];

client_data* c_data;

void del(){
    delete [] c_data;
}

void calc_auc_v2(Party *proxy) {
    uint32_t size = c_data[0].size();
    uint32_t params[1] = {size};
    proxy->SendBytes(aucPrCurve, params, 1);
    uint64_t aupr = PrCurve(proxy, c_data, size);
    cout << "AUPR :\t" << ConvertToDouble(Reconstruct(proxy, aupr)) << endl;
}

int main(int argc, char* argv[]) {
    uint8_t role = atoi(argv[1]);
    uint16_t cport = atoi(argv[2]);
    string caddress(argv[3]);
    uint16_t hport = atoi(argv[4]);
    string haddress(argv[5]);
    int delta = atoi(argv[6]);
    string ss(argv[7]);
    cout<<ss<<endl;

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
    }else{
        file_flag = true;
    }

    if (cport != 0) {
        assert(cport < 1 << (sizeof(uint16_t) * 8));
    }

    if (hport != 0) {
        assert(hport < 1 << (sizeof(uint16_t) * 8));
    }

    Party *proxy;
    if (role==0)
        proxy = new Party(proxy1, hport, haddress, cport, caddress);
    else
        proxy = new Party(proxy2, hport, haddress, cport, caddress);

    if (file_flag) {
        FileData(ss, c_data, nstation, sample_size);
    }
    else{
        RandomData(proxy, c_data, nstation, sample_size);
    }

    cout << "Number of parties: " << nstation << endl;
    for(int i = 0; i < nstation; i++)
        cout << sample_size[i] << "\t" << endl;

    auto start = chrono::high_resolution_clock::now();
    Sort(proxy, c_data, nstation, delta);
    auto mid = chrono::high_resolution_clock::now();

    calc_auc_v2(proxy);

    ios_base::sync_with_stdio(false);
    auto end = chrono::high_resolution_clock::now();
    double time_taken = chrono::duration_cast<chrono::nanoseconds>(end - start).count();
    time_taken *= 1e-9;
    cout << "Time taken by program is : " << fixed << time_taken << setprecision(9) << " sec" << endl;

    proxy->PrintPaperFriendly(time_taken);

    time_taken = chrono::duration_cast<chrono::nanoseconds>(mid - start).count();
    time_taken *= 1e-9;
    cout << "Time taken by Sort is : " << fixed << time_taken << setprecision(9) << " sec" << endl;

    time_taken = chrono::duration_cast<chrono::nanoseconds>(end - mid).count();
    time_taken *= 1e-9;
    cout << "Time taken by AUC computation is : " << fixed << time_taken << setprecision(9) << " sec" << endl;

    proxy->SendBytes(coreEnd);
    PrintBytes();
    del();
    cout<<"*****************************"<<endl;
    return 0;
}
