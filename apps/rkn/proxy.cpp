//
// Created by aburak on 24.03.22.
//

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <deque>
#include <chrono>
#include <iomanip>
#include "../../core/rkn.h"
#include "../../utils/parse_options.h"
//#include "llib.h"
#include <unistd.h>
#include <sys/stat.h>
#include <dirent.h>
#include <cassert>

using namespace std;
int nstation;
uint64_t sample_size[1000];

struct prediction {
    uint64_t val;
    uint64_t label;
    uint64_t flag;
};

typedef std::deque<prediction> client_data;
client_data *c_data;

uint64_t **nf_data;


int main(int argc, char *argv[]) {
    return 0;
}
