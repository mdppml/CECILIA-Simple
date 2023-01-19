//
// Created by noah on 21/06/22.
//
#include <iostream>
#include "../core/core.h"
#include "../core/cnn.h"
#include "../core/rkn.h"
#include "../core/auc.h"


int main(int argc, char* argv[]) {
    if (argc != 3) {
        clog << "Error: The program requires exactly two arguments; the IP address of the helper and the port it is listening on." << endl;
        return 1;
    }
    string address(argv[1]);
    uint16_t port = strtol(argv[2], nullptr, 10);

    auto *helper = new Party(HELPER,port,address);
    bool keep_looping = true;
    uint32_t sz, n_gms, size1, size2;
    uint64_t params [9];
    op operation;
    auto start = chrono::high_resolution_clock::now();
    while (keep_looping){
        operation = static_cast<op>(helper->ReadByte());
        switch(operation) {
            case CORE_MMSB:
                sz = helper->ReadInt();
                MSB(helper,nullptr,sz);
                break;
            case CORE_MSB:
                MSB(helper,0);
                break;
            case CORE_MC:
                MOC(helper,0);
                break;
            case CORE_MMC:
                sz = helper->ReadInt();
                MOC(helper,nullptr,sz);
                break;
            case CORE_CMP:
                CMP(helper,0,0);
                break;
            case CORE_MCMP:
                sz = helper->ReadInt();
                CMP(helper,nullptr,nullptr,sz);
                break;
            case CORE_MUX:
                MUX(helper,0, 0, 0);
                break;
            case CORE_MMUX:
                sz = helper->ReadInt();
                MUX(helper,nullptr, nullptr, nullptr, sz);
                break;
            case CORE_MUL:
                MUL(helper,0, 0);
                break;
            case CORE_MMUL:
                sz = helper->ReadInt();
                MUL(helper,nullptr, nullptr, sz);
                break;
            case CORE_DP:
                sz = helper->ReadInt();
                DP(helper, nullptr, nullptr, sz);
                break;
            case CORE_MDP:
                sz = helper->ReadInt();
                DP(helper, nullptr, nullptr, sz, 0);
                break;
            case CORE_EXP:
                EXP(helper, 0);
                break;
            case CORE_MEXP:
                sz = helper->ReadInt();
                EXP(helper, nullptr, sz);
                break;
            case CORE_MATMATMUL:
                sz = helper->ReadInt();
                // note that a_row is the required size of the multiplication that will be performed in MATMATMUL
                MATMATMUL(helper, nullptr, nullptr, sz, 0, 0);
                break;
            case CORE_MMATMATMUL:
                sz = helper->ReadInt();
                // note that a_row is the required size of the multiplication that will be performed in MATMATMUL
                MATMATMUL(helper, nullptr, nullptr, 0, sz, 0, 0);
                break;
            case CORE_MATVECMUL:
                sz = helper->ReadInt();
                // note that a_row is the required size of the multiplication that will be performed in MATVECMUL
                MATVECMUL(helper, nullptr, nullptr, sz, 0);
                break;
            case CORE_MMATVECMUL:
                sz = helper->ReadInt();
                // note that a_row is the required size of the multiplication that will be performed in MATVECMUL
                MATVECMUL(helper, nullptr, nullptr, 0, sz, 0);
                break;
            case CNN_MAX:
                sz = helper->ReadInt();
                MAX(helper,nullptr, sz);
                break;
            case CNN_MMAX:
                params[0] = helper->ReadInt();
                params[1] = helper->ReadInt();
                params[2] = helper->ReadInt();
                params[3] = helper->ReadInt();
                if (
                        params[0] > 0
                        and params[1] > 0
                        and params[2] > 0
                        and params[2] <= params[0]
                        and params[3] <= params[1]
                    ){
                    MAX(helper, nullptr, params[0], params[1], params[2], params[3], true);
                }
                else{
                    cout << "ERROR: received mmax parameters were not in valid range..." << endl;
                }
                break;
            case CNN_ARGMAX:
                sz = helper->ReadInt();
                ARGMAX(helper,nullptr, sz);
                break;
            case CNN_RELU:
                RELU(helper, 0);
                break;
            case CNN_MRELU:
                sz = helper->ReadInt();
                RELU(helper, nullptr, sz);
                break;
            case CNN_DRLU:
                DRELU(helper, 0);
                break;
            case CNN_MDRLU:
                sz = helper->ReadInt();
                DRELU(helper, nullptr, sz);
                break;
            case CNN_CL:
                params[0] = helper->ReadInt();
                params[1] = helper->ReadInt();
                params[2] = helper->ReadInt();
                params[3] = helper->ReadInt();
                params[4] = helper->ReadInt();
                params[5] = helper->ReadInt();
                params[6] = helper->ReadInt();
                params[7] = helper->ReadInt();
                params[8] = helper->ReadInt();
                CL(helper, nullptr, params[0], params[1], params[2], nullptr, params[3], params[4], params[5], params[6],
                   params[7], nullptr, params[8]);
                break;
            case CNN_FCL:
                params[0] = helper->ReadInt();
                params[1] = helper->ReadInt();
                FCL(helper, nullptr, params[0], nullptr, params[1], nullptr);
                break;
            case RKN_GM2KM:
                n_gms = helper->ReadInt();
                sz = helper->ReadInt();
                GM2KM(helper, nullptr, 0, n_gms, sz);
                break;
            case RKN_INVSQRT:
                sz = helper->ReadInt();
                INVSQRT(helper, nullptr, sz);
                break;
            case RKN_MINVSQRT:
                n_gms = helper->ReadInt();
                sz = helper->ReadInt();
                INVSQRT(helper, nullptr, n_gms, sz);
                break;
            case RKN_ITER:
                size1 = helper->ReadInt();
                size2 = helper->ReadInt();
                RKN_ITERATION(helper, nullptr, nullptr, nullptr, size1, size2, 0, 0, 0);
                break;
            case CORE_DIV:
                DIV(helper, 0, 0);
                break;
            case CORE_MDIV:
                size1 = helper->ReadInt();
                DIV(helper, 0, 0, size1);
                break;
            case CORE_MNORM:
                size1 = helper->ReadInt();
                NORM(helper, 0, 0, size1);
                break;
            case AUC_MSB:
                sz = helper->ReadInt();
                AUCMSB(helper,nullptr,sz);
                break;
            case AUC_TDIV:
                DIVISION(helper,0,0);
                break;
            case AUC_MROU:
                sz = helper->ReadInt();
                MRound(helper,nullptr,sz);
                break;
            case AUC_MDIV:
                sz = helper->ReadInt();
                MDIVISION(helper,nullptr,nullptr,sz);
                break;
            case AUC_ROCNOTIE:
                sz = helper->ReadInt();
                cout << "Sz: " << sz << endl;
                ROCNOTIE(helper, 0, sz);
                break;
            case AUC_ROCWITHTIE:
                sz = helper->ReadInt();
                ROCWITHTIE(helper, 0, sz);
                break;
            case AUC_PR:
                sz = helper->ReadInt();
                PRCURVE(helper, 0, sz);
                break;
            case RKN_EIG:
                sz = helper->ReadInt();
                EIG(helper, sz);
                break;
            case RKN_MEIG:
                n_gms = helper->ReadInt();
                sz = helper->ReadInt();
                EIG(helper, n_gms, sz);
                break;
            case CORE_END:
                keep_looping = false;
                break;
        }
    }
    auto end = chrono::high_resolution_clock::now();
    helper->PrintBytes();

    double time_taken = chrono::duration_cast<chrono::nanoseconds>(end - start).count();
    time_taken *= 1e-9;
    helper->PrintPaperFriendly(time_taken);
    delete helper;
    return 0;
}