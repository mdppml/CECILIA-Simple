#include <iostream>
#include "../../core/core.h"
#include "../../core/cnn.h"

int main(int argc, char* argv[]) {
    if (argc < 3){
        cout << "Calling Helper without specifying address (1) and port (2) is not possible." << endl;
        return -1;
    }
    uint16_t port = atoi(argv[2]);
    string address(argv[1]);

    Party *helper = new Party(HELPER,port,address);
    while (1){
        int op = helper->ReadByte();
        cout << "Operation: " << op << endl;
        if (op == CNN_MAX){
            int matrix_size = helper->ReadInt();
            MAX(helper, 0, matrix_size);
        }
        else if (op == CNN_MMAX){
            int mmaxParams = helper->ReadInt();
            uint16_t mRows = (mmaxParams >> 48);
            uint16_t mCols = (mmaxParams >> 32);
            uint16_t window_size = (mmaxParams & 0b0000000011111111);
            MAX(helper, nullptr, mRows, mCols, window_size);
        }
        else if (op == CNN_RELU){
            RELU(helper,0);
        }
        else if (op == CNN_DRLU){
            DRELU(helper,0);
        }
        else if (op == CNN_CL){
            unsigned char *ptr = helper->getBuffer1();
            Receive(helper->getSocketP2(), helper->getBuffer1(), 4 * 8);
            uint64_t params [4];
            convert2Array(&ptr, &params[0], 4);
            cout << "HELPER CL: i_dim = " << params[0] << ", k_dim= " << params[1] << ", k_number = " << params[2] << ", stride= " << params[3] << endl;
            CL(helper, nullptr, params[0], nullptr, params[1], params[2], params[3]);
            cout << "finished CL1" << endl;
        }
        else if (op == CNN_CL2){
            cout << "operation CL2 received." << endl;
            unsigned char *ptr = helper->getBuffer1();
            Receive(helper->getSocketP2(), helper->getBuffer1(), 4 * 8);
            uint64_t params [4];
            convert2Array(&ptr, &params[0], 4);
            cout << "HELPER CL2: i_dim = " << params[0] << ", i_number= " << params[1] << ", k_dim = " << params[2] << ", stride= " << params[3] << endl;
            //CL(helper, nullptr, params[0], params[1], nullptr, params[2], params[3]);
            cout << "finished CL2" << endl;
        }
        else if (op == CNN_FCL){
            cout << "operation FCL received." << endl;
            unsigned char *ptr = helper->getBuffer1();
            Receive(helper->getSocketP2(), helper->getBuffer1(), 3 * 8);
            uint64_t params [3];
            convert2Array(&ptr, &params[0], 3);
            cout << "HELPER FCL: i_dim = " << params[0] << ", i_number= " << params[1] << ", node_number = " << params[2] << endl;
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
