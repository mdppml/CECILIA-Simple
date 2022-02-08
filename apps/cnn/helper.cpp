#include <iostream>
#include "../../core/core.h"
#include "../../core/cnn.h"

int main(int argc, char* argv[]) {
    uint16_t port = atoi(argv[2]);
    string address(argv[1]);

    Party *helper = new Party(HELPER,port,address);
    while (1){
        int op = helper->ReadByte();
        if (op == CNN_MAX){
            int matrix_size = helper->ReadInt();
            MAX(helper, 0, matrix_size);
        }
        else if (op == CNN_MMAX){
            int mmaxParams = helper->ReadInt();
            uint16_t matrix_size = (mmaxParams >> 16);
            uint16_t window_size = (mmaxParams & 0b0000000011111111);
            MAX(helper, nullptr, matrix_size, window_size);
        }
        else if (op == CNN_RST){
            int sz = helper->ReadInt();
            MSB(helper,0,sz);
        }
        else if (op == CORE_MMSB){
            int sz = helper->ReadInt();
            MSB(helper,0,sz);
        }
        else if (op == CORE_MSB){
            MSB(helper,0);
        }
        else if (op == CORE_MC){
            MOC(helper,0);
        }
        else if (op == CORE_MMC){
            int sz = helper->ReadInt();
            MOC(helper,0,sz);
        }
        else if (op == CORE_CMP){
            CMP(helper,0,0);
        }
        else if (op == CORE_MCMP){
            int sz = helper->ReadInt();
            CMP(helper,0,0,sz);
        }else if (op == CORE_MUX){
            MUX(helper,0, 0, 0);
        }
        else if (op == CORE_MMUX){
            int sz = helper->ReadInt();
            MUX(helper,0, 0, 0, sz);
        }else if (op == CORE_MUL){
            MUL(helper,0, 0);
        }
        else if (op == CORE_MMUL){
            int sz = helper->ReadInt();
            MUL(helper,0, 0, sz);
        }
        else if (op == CORE_END)
            break;
    }
    helper->PrintBytes();
    return 0;
}
