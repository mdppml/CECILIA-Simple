#ifndef CECILIA_AUC_UTILS_H
#define CECILIA_AUC_UTILS_H

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <deque>
#include <unistd.h>
#include <sys/stat.h>
#include <dirent.h>
#include "llib.h"

using namespace std;

struct prediction {
    uint64_t val;
    uint64_t label;
};

typedef std::deque<prediction> client_data;

bool IsPathExist(const std::string &s)
{
    struct stat buf;
    return (stat (s.c_str(), &buf) == 0);
}

void ReadDirectory(const std::string& name, vector<string>& v)
{
    DIR* dirp = opendir(name.c_str());
    struct dirent * dp;
    while ((dp = readdir(dirp)) != NULL) {
        v.push_back(dp->d_name);
    }
    closedir(dirp);
}

/**
 * Generate random prediction scores and corresponding random labels. Generate duplicate predictions scores explicitly
 * if flag is true.
 * @param proxy
 * @param c_data
 * @param nstation
 * @param sample_size
 * @param flag: generate duplicate samples if it is true
 */
void RandomData(Party *proxy, client_data *&c_data, int nstation, uint64_t *sample_size, bool flag = false) {
    srand(100);
    c_data = new client_data[nstation];
    for (int i = 0; i < nstation; i++) {
        uint64_t tmp[MAX_SAMPLE_SIZE];
        for (int j = 0; j < sample_size[i]; j++) {
//            tmp[j] = (rand() % 10000) + 1;
            if(j % 2 == 1 && flag)
                tmp[j] = tmp[j - 1];
            else
                tmp[j] = proxy->GenerateCommonRandom() & MAX_SCALAR;
        }

        SortValues(tmp, sample_size[i]);
        for (int j = 0; j < sample_size[i]; j++) {
            uint64_t l = proxy->CreateShare(ConvertToUint64(rand() % 2));
            c_data[i].push_back({proxy->CreateShare(tmp[j]), l});

//            if (Role == 0)
//                c_data[i].push_back({(uint64_t) rand(), (uint64_t) rand()});
//            else
//                c_data[i].push_back({tmp[j] - rand(), l - rand(),});

        }
    }

}

void FileData(string path, client_data *&c_data, int nstation, uint64_t sample_size[]){
    vector<string> f_list;
    ReadDirectory(path, f_list);
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

/**
 * Print the data in the stations. It prints "transfer friendly" version to copy and paste to generate numpy array
 * if the corresponding transfer_friendly flag is true.
 * @param proxy
 * @param nstation
 * @param sample_size
 * @param data
 * @param transfer_friendly: print transfer friendly version to Python numpy array if it is true
 */
void PrintData(Party *proxy, int nstation, uint64_t *sample_size, client_data *data, bool transfer_friendly = false) {
    for (int i = 0; i < nstation; i++) {
        cout << "Station : " << i << endl;
        if(!transfer_friendly) {
            for (prediction n: data[i]) {
                cout << ConvertToDouble(Reconstruct(proxy, n.val)) << "\t" << ConvertToDouble(
                        Reconstruct(proxy, n.label)) << endl;
            }
        }
        else {
            cout << "[";
            for (int n = 0; n < data[i].size(); n++) {
                cout << ConvertToDouble(Reconstruct(proxy, data[i][n].val));
                if(n != data[i].size() - 1) {
                    cout << ", ";
                }
            }
            cout << "]" << endl;
            cout << "[";
            for (int n = 0; n < data[i].size(); n++) {
                cout << ConvertToDouble(Reconstruct(proxy, data[i][n].label));
                if(n != data[i].size() - 1) {
                    cout << ", ";
                }
            }
            cout << "]" << endl;
        }
    }
}

#endif //CECILIA_AUC_UTILS_H
