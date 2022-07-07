#include <iostream>
#include "../../core/core.h"
#include "../../core/cnn.h"

int main(int argc, char* argv[]) {
    if (argc < 2){
        cout << "Calling Helper without specifying address (1) and port (2) is not possible." << endl;
        return 1;
    }
    string address(argv[1]);
    uint16_t port = atoi(argv[2]);

    Party *helper = new Party(HELPER,port,address);
    while (1){
        int op = helper->ReadByte();
        if (op == CNN_MAX){
            int matrix_size = helper->ReadInt();
            MAX(helper, 0, matrix_size);
        }
        else if (op == CNN_MMAX){
            uint32_t *mmaxParams = new uint32_t [3];
            mmaxParams[0] = helper->ReadInt();
            mmaxParams[1] = helper->ReadInt();
            mmaxParams[2] = helper->ReadInt();
            MAX(helper, nullptr, mmaxParams[0], mmaxParams[1], mmaxParams[2]);
        }
        else if (op == CNN_RELU){
            RELU(helper,0);
        }
        else if (op == CNN_DRLU){
            DRELU(helper,0);
        }
        else if (op == CNN_CL){
            uint32_t *params = new uint32_t [6];
            params[0] = helper->ReadInt();
            params[1] = helper->ReadInt();
            params[2] = helper->ReadInt();
            params[3] = helper->ReadInt();
            params[4] = helper->ReadInt();
            params[5] = helper->ReadInt();
            CL(helper, nullptr, params[0], params[1], params[2], nullptr, params[3], params[4], params[5]);
        }
        else if (op == CNN_FCL){
            uint32_t *params = new uint32_t [3];
            params[0] = helper->ReadInt();
            params[1] = helper->ReadInt(); //TODO this can be calculated here.
            params[2] = helper->ReadInt();
            FCL(helper, nullptr, params[0], params[1], nullptr, params[2]);
        }
        else if (op == CORE_MSB){
            MSB(helper,0);
        }
        else if (op == CORE_MMSB){
            int sz = helper->ReadInt();
            MSB(helper,0,sz);
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
        else if (op == CORE_DIV){
            DIV(helper,0, 0);
        }
        else if (op == CORE_END) {
            cout << "programm was ended. " << endl;
            break;
        }
    }
    helper->PrintBytes();
    return 0;
}
