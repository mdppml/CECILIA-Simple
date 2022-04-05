//
// Created by aburak on 24.03.22.
//

#include <iostream>
#include "../../core/rkn.h"
#include "../../core/cnn.h"
#include "../../utils/constant.h"
#include "../../utils/parse_options.h"

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
        else if (op == CORE_EXP) {
            EXP(helper, 0);
        }
        else if (op == CORE_MEXP) {
            int sz = helper->ReadInt();
            EXP(helper, 0, sz);
        }
        else if (op == CNN_MAX){
            int matrix_size = helper->ReadInt();
            MAX(helper,nullptr, matrix_size);
        }
        else if (op == CNN_MMAX){
            cout << "MMAX was called..." << endl;

            Receive(helper->getSocketP1(), helper->getBuffer1(), 8 * 4);
            Receive(helper->getSocketP2(), helper->getBuffer2(), 8 * 4);

            unsigned char* ptr = helper->getBuffer1();
            unsigned char* ptr2 = helper->getBuffer2();
            uint64_t mmaxParams [4];
            bool areParamsMatching = true;
            for (uint8_t i = 0; i < 4; i++){
                mmaxParams[i] = convert2Long(&ptr);
                uint64_t p2 = convert2Long(&ptr2);
                if (mmaxParams[i] != p2){
                    cout << "Parameters from P0 and P1 must match... " << convert2double(mmaxParams[i]) << " != " << convert2double(p2) << endl;
                    areParamsMatching = false;
                }
            }
            if (mmaxParams[0] > 0 and mmaxParams[1] > 0 and mmaxParams[2] > 0 and mmaxParams[2] <= mmaxParams[0] and mmaxParams[2] <= mmaxParams[1] and areParamsMatching){
                MAX(helper, nullptr, mmaxParams[0], mmaxParams[1], mmaxParams[2]);
                cout << "finished MMAX" << endl;
            }
            else{
                cout << "ERROR: received mmax parameters were not matching each other or were not in valid range..." << endl;
            }
        }
        else if (op == CNN_RELU){
            RELU(helper, 0);
        }
        else if (op == CNN_DRLU){
            DRELU(helper, 0);
        }
        else if( op == RKN_EIG) {
            int sz = helper->ReadInt();
            EIG(helper, sz);
        }
        else if( op == RKN_MEIG) {
            int n_gms = helper->ReadInt();
            int sz = helper->ReadInt();
            EIG(helper, n_gms, sz);
        }
        else if (op == CORE_END)
            break;
    }
    helper->PrintBytes();
    return 0;
}
