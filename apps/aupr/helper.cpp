#include <iostream>
#include "../../core/auc.h"

int main(int argc, char* argv[]) {
    uint16_t port = atoi(argv[2]);
    string address(argv[1]);

    Party *helper = new Party(HELPER,port,address);
    while (1){
        int op = helper->ReadByte();
        if (op == CORE_MMSB){
            int sz = helper->ReadInt();
            MSB(helper,0,sz);
        }else if (op == AUC_MSB){
            int sz = helper->ReadInt();
            AUCMSB(helper,0,sz);
        }else if (op == CORE_MUL){
            MUL(helper,0,0);
        }else if (op == CORE_MMUL){
            int sz = helper->ReadInt();
            MUL(helper,0,0,sz);
        }else if (op == AUC_MROU){
            int sz = helper->ReadInt();
            MRound(helper,0,sz);
        }else if (op == CORE_MUX){
            MUX(helper,0,0,0);
        }else if (op == CORE_MMUX){
            int sz = helper->ReadInt();
            MUX(helper,0,0,0,sz);
        }
        else if (op == AUC_TDIV){
            DIVISION(helper,0,0);
        }
        else if (op == AUC_MDIV){
            int sz = helper->ReadInt();
            MDIVISION(helper,0,0,sz);
        }
        else if (op == CORE_END)
            break;
    }
    helper->PrintBytes();
    cout<<"*****************************"<<endl;
    return 0;
}
