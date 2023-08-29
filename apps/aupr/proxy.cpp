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

//bool IsPathExist(const std::string &s)
//{
//    struct stat buf;
//    return (stat (s.c_str(), &buf) == 0);
//}
//void read_directory(const std::string& name, vector<string>& v)
//{
//    DIR* dirp = opendir(name.c_str());
//    struct dirent * dp;
//    while ((dp = readdir(dirp)) != NULL) {
//        v.push_back(dp->d_name);
//    }
//    closedir(dirp);
//}
//void random_data(int nstation, uint64_t* sample_size,int role){
//    srand(100);
//    c_data = new client_data[nstation];
//    int ind = 10;
//    for (int i=0;i<nstation;i++){
//        uint64_t tmp[MAXSAMPLESIZE];
//        for (int j=0;j<sample_size[i];j++){
//            tmp[j] = (rand()%10000)+1;
//        }
//
//        sort_values(tmp,sample_size[i]);
//        for (int j=0;j<sample_size[i];j++){
//            uint64_t l = rand()%2;
//
//            if (role == 0)
//                c_data[i].push_back({(uint64_t)rand(),(uint64_t)rand()});
//            else
//                c_data[i].push_back({tmp[j]-rand(),l-rand()});
//
//        }
//    }
//
//}
//void file_data(string path){
//    vector<string> f_list;
//    read_directory(path,f_list);
//    c_data = new client_data[f_list.size()];
//    int f_index = 0;
//    for(string file : f_list) {
//        if (file != "." && file != "..") {
//            file = path + "/" + file;
//            ifstream ip(file.c_str());
//            cout<<file<<endl;
//            if (ip.is_open()){
//                cout<<"ddd"<<endl;
//                string tmp;
//                int s_size = 0;
//                while (ip.good()){
//                    getline(ip,tmp,',');
//                    if (tmp.empty())
//                        break;
//                    char* end;
//                    uint64_t label = strtoull( tmp.c_str(), &end,10 );
//                    getline(ip,tmp,'\n');
//                    uint64_t pred = strtoull( tmp.c_str(), &end,10 );
//                    c_data[f_index].push_back({pred,label});
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
//void print_data(int nstation, uint64_t* sample_size){
//    for (int i=0;i<nstation;i++){
//        cout<<"Station : "<<i<<endl;
//        for(prediction n : c_data[i]) {
//            cout<<n.val<<" "<<n.label<<endl;
//        }
//    }
//}
//void print_data(int nstation, uint64_t* sample_size, Party* proxy){
//    for (int i=0;i<nstation;i++){
//        cout<<"Station : "<<i<<endl;
//        for(prediction n : c_data[i]) {
//            cout<<REC(proxy,n.val)<<"\t"<<REC(proxy,n.label)<<endl;
//        }
//    }
//}

void del(){
    delete [] c_data;
}
//void calc_auc(Party* proxy){
//    uint32_t size = c_data[0].size();
//    uint64_t* labels = new uint64_t [size];
//    int i=0;
//    for(prediction n : c_data[0]) {
//        labels[i++] = n.label;
//    }
//
//    uint64_t TP = 0;
//    uint64_t pre_prec = 50;
//    uint64_t tmp = 0;
//    uint64_t pre_reca = 0;
//    uint64_t numerator = 0;
//    uint64_t numerator2 = 0;
//    uint64_t *precs = new uint64_t[size];
//    uint64_t *recas = new uint64_t[size];
//    for (int i=0;i<size;i++) {
//        TP = ADD(proxy,TP, labels[i]);
//        recas[i] = TP;
//        if (proxy->getPRole() == P1)
//            tmp = tmp + 1;
//        precs[i] = tmp;
//    }
//    uint64_t prec;
//    uint64_t reca;
//    proxy->SendBytes(AUC_MDIV,size);
//    precs = MDIVISION(proxy,recas,precs,size);
//    for (int i=0;i<size;i++){
//        prec = precs[i];
//        reca = recas[i];
//        proxy->SendBytes(CORE_MUL);
//        uint64_t area = MUL(proxy,pre_prec,reca-pre_reca);
//        proxy->SendBytes(CORE_MUL);
//        area = MUL(proxy,area,c_data[0][i].val);
//        numerator = ADD(proxy,numerator,area);
//
//        proxy->SendBytes(CORE_MUL);
//        area = MUL(proxy,reca-pre_reca,prec-pre_prec);
//        proxy->SendBytes(CORE_MUL);
//        area = MUL(proxy,area,c_data[0][i].val);
//        numerator2 = ADD(proxy,numerator2,area);
//
//        proxy->SendBytes(CORE_MMUX,2);
//        uint64_t *tmp2 = MUX(proxy,(uint64_t[]){pre_prec,pre_reca},(uint64_t[]){prec,reca},(uint64_t[]){c_data[0][i].val,c_data[0][i].val},2);
//        pre_prec = tmp2[0];
//        pre_reca = tmp2[1];
//
//
//
//
//    }
//
//    numerator = 2*numerator + numerator2;
//    uint64_t denominator = 2*TP;
//
//    proxy->SendBytes(AUC_TDIV);
//    uint64_t prc = DIVISION(proxy,numerator,denominator);
//
//    prc = REC(proxy,prc);
//
//    cout<<"PRC :\t"<<prc<<endl;
//
//
//    delete [] labels;
//}

void calc_auc_v2(Party *proxy) {
    uint32_t size = c_data[0].size();
    uint32_t params[1] = {size};
    proxy->SendBytes(AUC_PR, params, 1);
    uint64_t aupr = PRCURVE(proxy, c_data, size);
    cout << "AUPR :\t" << convert2double(Reconstruct(proxy, aupr)) << endl;
}

//void calc_confidence(Party* proxy){
//    uint32_t size = c_data[0].size();
//    uint64_t* preds= new uint64_t [size];
//    //uint64_t* dummy_indexes = new uint64_t [size/5];
//    //int i=0;
//    //int j=0;
//    for(int k=0;k<size;k++){
//        if (k!=(size-1))
//            preds[k] = c_data[0][k].val-c_data[0][k+1].val; // get the difference of the current and next prediction
//        else
//            preds[k] = 2; // the last prediction must be in the list.
//
//        preds[k]= preds[k]*proxy->generateCommonRandom(); // multiply with random value
//        //i++;
//        /*if ((proxy->generateCommonRandom()%10) == 0){ // add dummy predictions
//            if ((proxy->generateCommonRandom()%2) == 0) {
//                preds[i] = proxy->generateRandom(); // random dummy
//            }else{
//                if (proxy->GetRole()) // zero dummy
//                    preds[i] = proxy->generateCommonRandom();
//                else
//                    preds[i] = 0 - proxy->generateCommonRandom();
//
//                dummy_indexes[j++] = i++;
//            }
//        }*/
//    }
//
//    //get a permutation of preds
//    proxy->SendBytes(AUC_MROU, size);
//    uint64_t* rpreds = MRound(proxy,preds,size);
//    delete [] preds;
//
//    //get the reverse permutation of rpreds
//
//    /*j = 0;
//    int l = 0;
//    for(int k=0;k<i:k++){ // remove dummy values
//        if (k != dummy_indexes[j])
//            c_data[0][l++].val = rpreds[k];
//        else
//            j++;
//    }*/
//    for(int k=0;k<size;k++){ // remove dummy values
//        c_data[0][k].val = rpreds[k];
//    }
//    delete [] rpreds;
//
//}

//void sort_data(int nstation, Party* proxy, int delta){
//
//    int tmp_delta = delta;
//    while (nstation!=1) {
//        int i=0;
//        int ns = nstation-(nstation%2);
//        while (i<ns){
//            delta = tmp_delta;
//            int fl_index=i;
//            int ll_index=i+1;
//            if (c_data[i].size() < c_data[i+1].size()){
//                fl_index=i+1;
//                ll_index=i;
//            }
//            uint64_t *diff = new uint64_t[2*c_data[ll_index].size()];
//            uint64_t *mux_val1 = new uint64_t[2*c_data[ll_index].size()];
//            uint64_t *mux_val2 = new uint64_t[2*c_data[ll_index].size()];
//            uint64_t *mux_res = new uint64_t[2*c_data[ll_index].size()];
//            client_data sorted;
//
//
//            // sort
//            while (!c_data[fl_index].empty() && !c_data[ll_index].empty()){
//                int diff_size = c_data[fl_index].size();
//                if(diff_size > c_data[ll_index].size())
//                    diff_size = c_data[ll_index].size();
//
//                for (int j=0;j<diff_size;j++){
//                    diff[j] = c_data[fl_index][j].val - c_data[ll_index][j].val;
//                    mux_val1[j] = c_data[fl_index][j].val;
//                    mux_val2[j] = c_data[ll_index][j].val;
//                    mux_val1[j+diff_size] = c_data[fl_index][j].label;
//                    mux_val2[j+diff_size] = c_data[ll_index][j].label;
//                }
//
//
//                proxy->SendBytes(CORE_MMSB, diff_size);
//                uint64_t * diff_res = MSB(proxy,diff,diff_size);
//
//                for (int j=0;j<diff_size;j++){
//                    diff[j] = diff_res[j];
//                    diff[j+diff_size] = diff_res[j];
//                }
//
//
//                proxy->SendBytes(CORE_MMUX, 2*diff_size);
//                mux_res = MUX(proxy,mux_val1,mux_val2,diff,2*diff_size);
//                for (int j=0;j<diff_size;j++){
//                    c_data[fl_index][j].val = mux_res[j];
//                    c_data[fl_index][j].label = mux_res[j+diff_size];
//                }
//
//                proxy->SendBytes(CORE_MMUX, 2*diff_size);
//                mux_res = MUX(proxy,mux_val2,mux_val1,diff,2*diff_size);
//                for (int j=0;j<diff_size;j++){
//                    c_data[ll_index][j].val = mux_res[j];
//                    c_data[ll_index][j].label = mux_res[j+diff_size];
//                }
//                sorted.push_back(c_data[fl_index][0]);
//                c_data[fl_index].pop_front();
//
//
//
//                int min_val = min({delta, (int)c_data[fl_index].size(), (int)c_data[ll_index].size()});
//                int fl_pop = 1;
//                int ll_pop = 0;
//
//                for (int j=0;j<min_val;j++){
//                    if (c_data[ll_index].size() == 1){
//                        delta = 0;
//                        break;
//                    }
//                    if (fl_pop != ll_pop) {
//                        diff[0] = c_data[fl_index][0].val - c_data[ll_index][0].val;
//                        proxy->SendBytes(AUC_MSB, 1);
//                        uint64_t *diff_res = AUCMSB(proxy,diff, 1);
//                        uint64_t cmp = REC(proxy,diff_res[0]);
//                        if (cmp == 0) {
//                            sorted.push_back(c_data[fl_index][0]);
//                            c_data[fl_index].pop_front();
//                            fl_pop++;
//                        } else {
//                            sorted.push_back(c_data[ll_index][0]);
//                            c_data[ll_index].pop_front();
//                            ll_pop++;
//                        }
//                    }
//                    else{
//                        sorted.push_back(c_data[fl_index][0]);
//                        c_data[fl_index].pop_front();
//                        fl_pop++;
//                    }
//
//                }
//
//            }
//            if (c_data[fl_index].size() > 0){
//                while (!c_data[fl_index].empty()){
//                    sorted.push_back(c_data[fl_index][0]);
//                    c_data[fl_index].pop_front();
//                }
//            }else{
//                while (!c_data[ll_index].empty()){
//                    sorted.push_back(c_data[ll_index][0]);
//                    c_data[ll_index].pop_front();
//                }
//            }
//
//            // merge
//            c_data[i] = sorted;
//            c_data[i+1].clear();
//
//            // delete
//            delete[] diff;
//            delete[] mux_val1;
//            delete[] mux_val2;
//            delete[] mux_res;
//
//            i+=2;
//        }
//
//        i=1;
//        c_data[1].clear();
//        while (i<nstation){
//            if (i%2 == 0)
//                c_data[i/2] = c_data[i];
//            c_data[i].clear();
//            i+=1;
//        }
//        nstation = (nstation/2)+(nstation%2);
//
//    }
//}

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
        proxy = new Party(P1,hport, haddress, cport, caddress);
    else
        proxy = new Party(P2,hport, haddress, cport, caddress);

    if (file_flag) {
        file_data(ss, c_data, nstation, sample_size);
    }
    else{
        random_data(proxy, c_data, nstation, sample_size);
    }

    cout << "Number of parties: " << nstation << endl;
    for(int i = 0; i < nstation; i++)
        cout << sample_size[i] << "\t" << endl;

    auto start = chrono::high_resolution_clock::now();

//    sort_data(nstation,proxy, delta);
//    calc_confidence(proxy);
//    calc_auc(proxy);

//    cout << "Data after sorting: " << endl;
    SORT(proxy, c_data, nstation, delta);
    auto mid = chrono::high_resolution_clock::now();
//    cout<<"Confidence mapping was calculated"<<endl;
//    cout << "======================" << endl;
//    uint64_t total_n_samples = 0;
//    for(int i = 0; i < nstation; i++) {
//        total_n_samples += sample_size[i];
//    }
//    uint64_t tmp_size[1] = {total_n_samples};
//    print_data(proxy, 1, tmp_size, c_data, true);
//    cout << "======================" << endl;

    calc_auc_v2(proxy);
//    cout<<"AUC was calculated"<<endl;

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

    proxy->SendBytes(CORE_END);
    proxy->PrintBytes();
    del();
    cout<<"*****************************"<<endl;
    return 0;
}
