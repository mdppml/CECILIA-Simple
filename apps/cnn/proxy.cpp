#include <cstdlib>
#include <iostream>
#include <deque>
#include <chrono>
#include <iomanip>
#include "../../core/cnn.h"
using namespace std;


int main(int argc, char* argv[]) {
    uint8_t role = atoi(argv[1]);

    Party *proxy;
    if (role==0)
        proxy = new Party(P1,7777, "127.0.0.1", 8888, "127.0.0.1");
    else
        proxy = new Party(P2,7777, "127.0.0.1", 8888, "127.0.0.1");


    // CNN inference pipeline, call functions sequantially for inference (Matrix MUL, MAXPOOL, RELU, vs...)


    proxy->SendBytes(CORE_END);
    proxy->PrintBytes();
    return 0;
}
