#include <iostream>
#include "../../core/Party.h"

int main(int argc, char* argv[]) {
    uint16_t port = atoi(argv[2]);
    string address(argv[1]);

    Party *helper = new Party(HELPER,port,address);
    while (1){
        int op = helper->ReadByte();
        if (op == MMSB){
            int sz = helper->ReadInt();
            helper->MMSB(0,sz);
        }else if (op == MMSB2){
            int sz = helper->ReadInt();
            helper->MMSB2(0,sz);
        }else if (op == MUL){
            helper->MUL(0,0);
        }else if (op == MMUL){
            int sz = helper->ReadInt();
            helper->MMUL(0,0,sz);
        }else if (op == MSSA){
            int sz = helper->ReadInt();
            helper->MSelectShare(0,0,0,sz);
        }else if (op == TDIV){
            helper->DIVISION(0,0);
        }
        else if (op == END)
            break;
    }
    helper->PrintBytes();
    cout<<"*****************************"<<endl;
    return 0;
}
