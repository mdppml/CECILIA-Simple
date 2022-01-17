#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <deque>
#include <chrono>
#include <iomanip>
#include "../../core/Party.h"
#include "../../utils/parse_options.h"
#include "llib.h"
#include <unistd.h>
#include <sys/stat.h>
#include <dirent.h>
using namespace std;
int nstation;
uint64_t sample_size[1000];

struct prediction{
    uint64_t val;
    uint64_t label;
};

typedef std::deque<prediction> client_data;
client_data* c_data;



bool IsPathExist(const std::string &s)
{
    struct stat buf;
    return (stat (s.c_str(), &buf) == 0);
}
void read_directory(const std::string& name, vector<string>& v)
{
    DIR* dirp = opendir(name.c_str());
    struct dirent * dp;
    while ((dp = readdir(dirp)) != NULL) {
        v.push_back(dp->d_name);
    }
    closedir(dirp);
}

void random_data(int nstation, uint64_t* sample_size,int role){
    srand(100);
    c_data = new client_data[nstation];
    for (int i=0;i<nstation;i++){
        uint64_t tmp[MAXSAMPLESIZE];
        for (int j=0;j<sample_size[i];j++){
            tmp[j] = (rand()%10000)+1;
        }

        sort_values(tmp,sample_size[i]);
        for (int j=0;j<sample_size[i];j++){
            uint64_t l = rand()%2;
            if (role == 0)
                c_data[i].push_back({(uint64_t)rand(),(uint64_t)rand()});
            else
                c_data[i].push_back({tmp[j]-rand(),l-rand()});


        }
    }
}

void file_data(string path){
    vector<string> f_list;
    read_directory(path,f_list);
    c_data = new client_data[f_list.size()];
    int f_index = 0;
    for(string file : f_list) {
        if (file != "." && file != "..") {
            file = path + "/" + file;
            ifstream ip(file.c_str());
            cout<<file<<endl;
            if (ip.is_open()){
                cout<<"ddd"<<endl;
                string tmp;
                int s_size = 0;
                while (ip.good()){
                    getline(ip,tmp,',');
                    if (tmp.empty())
                        break;
                    char* end;
                    uint64_t label = strtoull( tmp.c_str(), &end,10 );
                    getline(ip,tmp,'\n');
                    uint64_t pred = strtoull( tmp.c_str(), &end,10 );
                    c_data[f_index].push_back({pred,label});
                    s_size++;
                }
                sample_size[f_index] = s_size;
                f_index++;
            }
            ip.close();
        }
    }
    nstation = f_index;
}

void print_data(int nstation, uint64_t* sample_size, Party* proxy){
    for (int i=0;i<nstation;i++){
        cout<<"Station : "<<i<<endl;
        for(prediction n : c_data[i]) {
            cout<<proxy->Reconstruct(n.val)<<"\t"<<proxy->Reconstruct(n.label)<<endl;
        }
    }
}

void calc_auc(Party* proxy){
    uint32_t size = c_data[0].size();
    uint64_t* labels = new uint64_t [size];
    int i=0;

    for(prediction n : c_data[0]) {
        labels[i++] = n.label;
    }

    uint64_t TP = 0;
    uint64_t FP = 0;
    uint64_t pre_FP = 0;
    uint64_t pre_TP = 0;
    uint64_t numerator = 0;
    uint64_t numerator2 = 0;
    for (int i=0;i<size;i++) {
        TP = proxy->ADD(TP, labels[i]);
        if (proxy->GetRole()==P1)
            FP = i - TP;
        else
            FP = 0 - TP;

        proxy->SendBytes(MUL);
        uint64_t area = proxy->MUL(pre_TP,FP-pre_FP);
        proxy->SendBytes(MUL);
        area = proxy->MUL(area,c_data[0][i].val);
        numerator = proxy->ADD(numerator,area);

        proxy->SendBytes(MUL);
        area = proxy->MUL(TP-pre_TP,FP-pre_FP);
        proxy->SendBytes(MUL);
        area = proxy->MUL(area,c_data[0][i].val);
        numerator2 = proxy->ADD(numerator2,area);


        uint64_t *tmp = proxy->MSelectShare((uint64_t[]){pre_FP,pre_TP},(uint64_t[]){FP,TP},(uint64_t[]){c_data[0][i].val,c_data[0][i].val},2);
        pre_FP = tmp[0];
        pre_TP = tmp[1];
        //pre_TP = proxy->SelectShare(pre_TP,TP,c_data[0][i].val);

    }


    uint64_t FN = TP;
    uint64_t TN = 0;

    if (proxy->GetRole()==P1)
        TN = size - TP;
    else
        TN = 0 - TP;


    proxy->SendBytes(MUL);
    uint64_t denominator = proxy->MUL(FN,TN);
    denominator = denominator*2;

    numerator = 2*numerator+numerator2;

    //proxy->SendBytes(TDIV);
    uint64_t auc = proxy->DIVISION(numerator,denominator);

    auc = proxy->Reconstruct(auc);
    cout<<"AUC :\t"<<auc<<endl;

    delete [] labels;
}
void calc_confidence(Party* proxy){
    uint32_t size = c_data[0].size();
    uint64_t* preds= new uint64_t [size];
    //uint64_t* dummy_indexes = new uint64_t [size/5];
    //int i=0;
    //int j=0;
    for(int k=0;k<size;k++){
        if (k!=(size-1))
            preds[k] = c_data[0][k].val-c_data[0][k+1].val; // get the difference of the current and next prediction
        else
            preds[k] = 2; // the last prediction must be in the list.

        preds[k]= preds[k]*proxy->generateCommonRandom(); // multiply with random value
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

    //get a permutation of preds
    proxy->SendBytes(MROU, size);
    uint64_t* rpreds = proxy->MRound(preds,size);
    delete [] preds;

    //get the reverse permutation of rpreds

    /*j = 0;
    int l = 0;
    for(int k=0;k<i:k++){ // remove dummy values
        if (k != dummy_indexes[j])
            c_data[0][l++].val = rpreds[k];
        else
            j++;
    }*/
    for(int k=0;k<size;k++){ // remove dummy values
        c_data[0][k].val = rpreds[k];
    }
    delete [] rpreds;

}
void sort_data(int nstation, Party* proxy, int delta) {

    int tmp_delta = delta;
    while (nstation!=1) {
        int i=0;
        int ns = nstation-(nstation%2);
        while (i<ns){
            delta = tmp_delta;
            int fl_index=i;
            int ll_index=i+1;
            if (c_data[i].size() < c_data[i+1].size()){
                fl_index=i+1;
                ll_index=i;
            }
            uint64_t *diff = new uint64_t[2*c_data[ll_index].size()];
            uint64_t *mux_val1 = new uint64_t[2*c_data[ll_index].size()];
            uint64_t *mux_val2 = new uint64_t[2*c_data[ll_index].size()];
            uint64_t *mux_res = new uint64_t[2*c_data[ll_index].size()];
            client_data sorted;


            // sort
            while (!c_data[fl_index].empty() && !c_data[ll_index].empty()){
                int diff_size = c_data[fl_index].size();
                if(diff_size > c_data[ll_index].size())
                    diff_size = c_data[ll_index].size();

                for (int j=0;j<diff_size;j++){
                    diff[j] = c_data[fl_index][j].val - c_data[ll_index][j].val;
                    mux_val1[j] = c_data[fl_index][j].val;
                    mux_val2[j] = c_data[ll_index][j].val;
                    mux_val1[j+diff_size] = c_data[fl_index][j].label;
                    mux_val2[j+diff_size] = c_data[ll_index][j].label;
                }


                proxy->SendBytes(MMSB, diff_size);
                uint64_t * diff_res = proxy->MMSB(diff,diff_size);

                for (int j=0;j<diff_size;j++){
                    diff[j] = diff_res[j];
                    diff[j+diff_size] = diff_res[j];
                }


                mux_res = proxy->MSelectShare(mux_val1,mux_val2,diff,2*diff_size);
                for (int j=0;j<diff_size;j++){
                    c_data[fl_index][j].val = mux_res[j];
                    c_data[fl_index][j].label = mux_res[j+diff_size];
                }

                mux_res = proxy->MSelectShare(mux_val2,mux_val1,diff,2*diff_size);
                for (int j=0;j<diff_size;j++){
                    c_data[ll_index][j].val = mux_res[j];
                    c_data[ll_index][j].label = mux_res[j+diff_size];
                }
                sorted.push_back(c_data[fl_index][0]);
                c_data[fl_index].pop_front();



                int min_val = min({delta, (int)c_data[fl_index].size(), (int)c_data[ll_index].size()});
                int fl_pop = 1;
                int ll_pop = 0;

                for (int j=0;j<min_val;j++){
                    if (c_data[ll_index].size() == 1){
                        delta = 0;
                        break;
                    }
                    if (fl_pop != ll_pop) {
                        diff[0] = c_data[fl_index][0].val - c_data[ll_index][0].val;
                        proxy->SendBytes(MMSB2, 1);
                        uint64_t *diff_res = proxy->MMSB2(diff, 1);
                        uint64_t cmp = proxy->Reconstruct(diff_res[0]);
                        if (cmp == 0) {
                            sorted.push_back(c_data[fl_index][0]);
                            c_data[fl_index].pop_front();
                            fl_pop++;
                        } else {
                            sorted.push_back(c_data[ll_index][0]);
                            c_data[ll_index].pop_front();
                            ll_pop++;
                        }
                    }
                    else{
                        sorted.push_back(c_data[fl_index][0]);
                        c_data[fl_index].pop_front();
                        fl_pop++;
                    }

                }

            }
            if (c_data[fl_index].size() > 0){
                while (!c_data[fl_index].empty()){
                    sorted.push_back(c_data[fl_index][0]);
                    c_data[fl_index].pop_front();
                }
            }else{
                while (!c_data[ll_index].empty()){
                    sorted.push_back(c_data[ll_index][0]);
                    c_data[ll_index].pop_front();
                }
            }

            // merge
            c_data[i] = sorted;
            c_data[i+1].clear();

            // delete
            delete[] diff;
            delete[] mux_val1;
            delete[] mux_val2;
            delete[] mux_res;

            i+=2;
        }

        i=1;
        c_data[1].clear();
        while (i<nstation){
            if (i%2 == 0)
                c_data[i/2] = c_data[i];
            c_data[i].clear();
            i+=1;
        }
        nstation = (nstation/2)+(nstation%2);

    }
}
int main(int argc, char* argv[]) {
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
            stringstream sss(ss);
            while (sss.good()) {
                string substr;
                getline(sss, substr, ',');
                if (i == 0) {
                    nstation = stoi(substr);
                    if (nstation <= 2 || nstation > 1000) {
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
                sample_size[i] = 2000;
        }
    }else{
        file_flag = true;
    }

    Party *proxy;
    if (role == 0)
        proxy = new Party(P1, hport, haddress, cport, caddress);
    else
        proxy = new Party(P2, hport, haddress, cport, caddress);

    if (file_flag) {
        file_data(ss);
    }
    else{
        random_data(nstation,sample_size,role);
    }

    auto start = chrono::high_resolution_clock::now();
    sort_data(nstation,proxy, delta);
    cout<<"Data was sorted"<<endl;
    calc_confidence(proxy);
    cout<<"Confidence mapping was calculated"<<endl;
    calc_auc(proxy);
    cout<<"AUC was calculated"<<endl;
    ios_base::sync_with_stdio(false);
    auto end = chrono::high_resolution_clock::now();
    double time_taken =
            chrono::duration_cast<chrono::nanoseconds>(end - start).count();
    time_taken *= 1e-9;
    cout << "Time taken by program is : " << fixed
         << time_taken << setprecision(9);
    cout << " sec" << endl;


    proxy->SendBytes(END);
    proxy->PrintBytes();
    delete [] c_data;

    delete proxy;
    cout<<"*****************************"<<endl;
    return 0;
}
