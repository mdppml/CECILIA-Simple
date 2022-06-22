#include <iostream>
#include "../../core/core.h"
#include "../../core/cnn.h"
#include "../../core/rkn.h"

int main(int argc, char* argv[]) {
    string address(argv[1]);
    uint16_t port = atoi(argv[2]);

    Party *helper = new Party(HELPER,port,address);
    while (1){
        int op = helper->ReadByte();
        if (op == CORE_MMSB){
            cout << "CORE_MMSB" << endl;
            int sz = helper->ReadInt();
            MSB(helper,0,sz);
        }
        else if (op == CORE_MSB){
            cout << "CORE_MSB" << endl;
            MSB(helper,0);
        }
        else if (op == CORE_MC){
            cout << "CORE_MC" << endl;
            MOC(helper,0);
        }
        else if (op == CORE_MMC){
            cout << "CORE_MMC" << endl;
            int sz = helper->ReadInt();
            MOC(helper,0,sz);
        }
        else if (op == CORE_CMP){
            cout << "CORE_CMP" << endl;
            CMP(helper,0,0);
        }
        else if (op == CORE_MCMP){
            cout << "CORE_MCMP" << endl;
            int sz = helper->ReadInt();
            CMP(helper,0,0,sz);
        }else if (op == CORE_MUX){
            cout << "CORE_MUX" << endl;
            MUX(helper,0, 0, 0);
        }
        else if (op == CORE_MMUX){
            cout << "CORE_MMUX" << endl;
            int sz = helper->ReadInt();
            MUX(helper,0, 0, 0, sz);
        }else if (op == CORE_MUL){
            cout << "CORE_MUL" << endl;
            MUL(helper,0, 0);
        }
        else if (op == CORE_MMUL){
            cout << "CORE_MMUL" << endl;
            int sz = helper->ReadInt();
            MUL(helper,0, 0, sz);
        }
        else if( op == CORE_DP) {
            cout << "CORE_DP" << endl;
            int sz = helper->ReadInt();
            DP(helper, 0, 0, sz);
        }
        else if( op == CORE_MDP) {
            cout << "CORE_MDP" << endl;
            int sz = helper->ReadInt();
            DP(helper, 0, 0, sz, 0);
        }
        else if (op == CORE_EXP) {
            cout << "CORE_EXP" << endl;
            EXP(helper, 0);
        }
        else if (op == CORE_MEXP) {
            cout << "CORE_MEXP" << endl;
            int sz = helper->ReadInt();
            EXP(helper, 0, sz);
        }
        else if (op == CORE_MATMATMUL) {
            cout << "CORE_MATMATMUL" << endl;
            int sz = helper->ReadInt();
            // note that a_row is the required size of the multiplication that will be performed in MATMATMUL
            MATMATMUL(helper, 0, 0, sz, 0, 0);
        }
        else if (op == CORE_MMATMATMUL) {
            cout << "CORE_MMATMATMUL" << endl;
            int sz = helper->ReadInt();
            // note that a_row is the required size of the multiplication that will be performed in MATMATMUL
            MATMATMUL(helper, 0, 0, 0, sz, 0, 0);
        }
        else if (op == CORE_MATVECMUL) {
            cout << "CORE_MATVECMUL" << endl;
            int sz = helper->ReadInt();
            // note that a_row is the required size of the multiplication that will be performed in MATVECMUL
            MATVECMUL(helper, 0, 0, sz, 0);
        }
        else if (op == CORE_MMATVECMUL) {
            cout << "CORE_MMATVECMUL" << endl;
            int sz = helper->ReadInt();
            // note that a_row is the required size of the multiplication that will be performed in MATVECMUL
            MATVECMUL(helper, 0, 0, 0, sz, 0);
        }
        else if (op == CNN_MAX){
            cout << "CNN_MAX" << endl;
            int matrix_size = helper->ReadInt();
            MAX(helper,nullptr, matrix_size);
        }
        else if (op == CNN_MMAX){
            cout << "CNN_MMAX" << endl;

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
            cout << "CNN_RELU" << endl;
            RELU(helper, 0);
        }
        else if (op == CNN_MRELU){
            cout << "CNN_MRELU" << endl;
            uint64_t size = helper->ReadInt();
            RELU(helper, 0, size);
        }
        else if (op == CNN_DRLU){
            cout << "CNN_DRLU" << endl;
            DRELU(helper, 0);
        }
        else if (op == CNN_MDRLU){
            cout << "CNN_MDRLU" << endl;
            int size = helper->ReadInt();
            DRELU(helper, nullptr, size);
        }
        else if (op == CNN_CL){
            cout << "CNN_CL" << endl;
            unsigned char *ptr = helper->getBuffer1();
            Receive(helper->getSocketP1(), helper->getBuffer1(), 4*8);
            uint64_t params[4];
            convert2Array(&ptr, &params[0], 4);
            CL(helper, nullptr, params[0], nullptr, params[1], params[2], params[3]);
        }
        else if( op == RKN_GM2KM) {
            cout << "RKN_GM2KM" << endl;
            int n_gms = helper->ReadInt();
            int sz = helper->ReadInt();
            GM2KM(helper, 0, 0, n_gms, sz);
        }
        else if( op == RKN_INVSQRT) {
            cout << "RKN_INVSQRT" << endl;
            int sz = helper->ReadInt();
            INVSQRT(helper, 0, sz);
        }
        else if( op == RKN_MINVSQRT) {
            cout << "RKN_MINVSQRT" << endl;
            int n_gms = helper->ReadInt();
            int sz = helper->ReadInt();
            INVSQRT(helper, 0, n_gms, (uint32_t) sz);
        }
        else if( op == RKN_ITER) {
            cout << "RKN_ITER" << endl;
            int size1 = helper->ReadInt();
            int size2 = helper->ReadInt();
            RKN_ITERATION(helper, 0, 0, 0, size1, size2, 0, 0, 0);
        }
        else if (op == CORE_DIV){
            DIV(helper, 0, 0);
        }
        else if (op == CORE_END)
            break;
    }
    helper->PrintBytes();
    return 0;
}
