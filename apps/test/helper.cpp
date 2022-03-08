#include <iostream>
#include "../../core/core.h"
#include "../../core/cnn.h"

int main(int argc, char* argv[]) {
    string address(argv[1]);
    uint16_t port = atoi(argv[2]);

    Party *helper = new Party(HELPER,port,address);
    while (1){
        int op = helper->ReadByte();
        if (op == CORE_MMSB){
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
        else if (op == CNN_MAX){
            cout << "MAX was called..." << endl;
            int matrix_size = helper->ReadInt();
            MAX(helper,nullptr, matrix_size);
            cout << "MAX finished." << endl;
        }
        else if (op == CNN_MMAX){
            cout << "MMAX was called..." << endl;
            uint8_t *buffer = helper->getBuffer1();
            uint8_t *buffer2 = helper->getBuffer2();

            thread thr1 = thread(Receive,helper->getSocketP1(), buffer, 8 * 4);
            thread thr2 = thread(Receive,helper->getSocketP2(), buffer2, 8 * 4);
            thr1.join();
            thr2.join();

            uint64_t mRows = buffer[0];
            uint64_t mCols = buffer[1];
            uint64_t window_size = buffer[2];
            if (mRows == buffer2[0] and mCols == buffer2[1] and window_size == buffer2[2] and buffer[3] == buffer2[3]){
                cout << "mRows = " << mRows << ", mCols = " << mCols << ", wSize = " << window_size << endl;
                MAX(helper, nullptr, mRows, mCols, window_size);
            }
            else{
                cout << "ERROR: received mmax parameters were not matching each other..." << endl;
            }
        }
        else if (op == CNN_RELU){
            RELU(helper, 0);
        }
        else if (op == CNN_DRLU){
            DRELU(helper, 0);
        }
        else if (op == CORE_END)
            break;
    }
    helper->PrintBytes();
    return 0;
}
