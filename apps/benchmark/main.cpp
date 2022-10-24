//
// Created by noah on 21/06/22.
//
#include <iostream>
#include "Party_BM.h"

static const int REQUIRED_ARGS_COUNT = 14;

static void print_error() {
    clog
        << "Required arguments not provided. These are the required (and optional) arguments:\n"
        << "1: role (P1, P2 or HELPER)\n"
        << "2 & 3: helper's port and IP\n"
        << "4 & 5: P1's port and IP\n"
        << "6: number of operations (for vectorised operations)\n"
        << "7: vector length (for operations on vectors)\n"
        << "8 & 9: matrix size (x, y) (gram matrix has size x*x)\n"
        << "10: window size (used in MAXPOOL)\n"
        << "11: kernel size (used in CNN)\n"
        << "12: number of kernel (used in CNN)\n"
        << "13: number of cycles (how often to run each function for the benchmark)\n"
        << "(14+: which functions to run (e.g. MUL, MUX), without this, all functions are run)" << endl;
    exit(EXIT_FAILURE);
}

int suppress_stdout() {
    fflush(stdout);

    int stdout_file_descriptor = dup(1);
    int null_file_descriptor = open("/dev/null", O_WRONLY);
    // check null_file_descriptor for error omitted
    dup2(null_file_descriptor, 1);
    close(null_file_descriptor);

    return stdout_file_descriptor;
}

void resume_stdout(int file_descriptor) {
    fflush(stdout);
    dup2(file_descriptor, 1);
    close(file_descriptor);
}


int main(int argc, char* argv[]) {
    if (argc < REQUIRED_ARGS_COUNT) {
        print_error();
    }
    //parse arguments:
    string role_string(argv[1]);
    role proxy_role;
    if (role_string == "P1") {
        proxy_role = P1;
    } else if (role_string == "P2") {
        proxy_role = P2;
    } else if (role_string == "HELPER") {
        proxy_role = HELPER;
    } else {
        print_error();
    }
    uint16_t helper_port = strtol(argv[2], nullptr, 10);
    string helper_ip(argv[3]);
    uint16_t p1_port = strtol(argv[4], nullptr, 10);
    string p1_ip(argv[5]);
    int op_count = stoi(argv[6], nullptr, 10);
    int vector_length = stoi(argv[7]);
    int matrix_x = stoi(argv[8]);
    int matrix_y = stoi(argv[9]);
    int window_size = stoi(argv[10]);
    int kernel_size = stoi(argv[11]);
    int kernel_count = stoi(argv[12]);
    int cycle_count = stoi(argv[13]);
    // create output formatting variables:
    int delay;
    string padding;
    switch(proxy_role) {
        case HELPER:
            delay = 2;
            padding = "";
            break;
        case P1:
            delay = 0;
            padding = "    ";
            break;
        case P2:
            delay = 1;
            padding = "    ";
            break;
    }
    // initialise proxy:
    int file_descriptor = suppress_stdout();
    Party_BM* proxy = (proxy_role == HELPER) ?
            new Party_BM(helper_port, helper_ip, op_count, vector_length, matrix_x, matrix_y, window_size, kernel_size, kernel_count, cycle_count) :
            new Party_BM(proxy_role, helper_port, helper_ip, p1_port, p1_ip, op_count, vector_length, matrix_x, matrix_y, window_size, kernel_size, kernel_count, cycle_count);
    resume_stdout(file_descriptor);
    // parse/obtain function names to run:
    string* functions;
    int function_count = argc-REQUIRED_ARGS_COUNT;
    if (function_count == 0) {
        tie(function_count, functions) = proxy->get_all_function_names();
    } else {
        functions = new string[function_count];
        for (int i = 0; i < function_count; i++) {
            functions[i] = argv[i+REQUIRED_ARGS_COUNT];
        }
    }
    // run benchmark:
    double cpu_time, real_time;
    for (int i = 0; i < function_count; i++) {
        file_descriptor = suppress_stdout();
        tie(cpu_time, real_time) = proxy->benchmark(functions[i]);
        resume_stdout(file_descriptor);
        if (cpu_time != -1) {
            if (proxy_role == P1) {
                cout << "\n" << functions[i] << "\nReal time:       " << real_time << " ms" << endl;
            }
            if (proxy_role == HELPER && (functions[i] == "REC" || functions[i] == "createShare")) {
                // skips helper print for functions that don't use it
                continue;
            }
            sleep(delay);
            cout << role_string << " CPU time: " << padding << cpu_time << " ms" << endl;
            cout << "Bytes: " << (bytesSend + bytesReceived)/cycle_count << endl;
            sleep(3-delay);
            bytesSend = 0;
            bytesReceived = 0;
        }
    }
    delete[] functions;
    delete proxy;
}