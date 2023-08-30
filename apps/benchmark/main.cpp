//
// Created by noah on 21/06/22.
//
#include <iostream>
#include "PartyBm.h"

static const int REQUIRED_ARGS_COUNT = 12;

static void PrintError() {
    clog
        << "Required arguments not provided. These are the required (and optional) arguments:\n"
        << "1: Role (proxy1, proxy2 or helper)\n"
        << "2 & 3: helper's port and IP\n"
        << "4 & 5: proxy1's port and IP\n"
        << "6: vector length (for operations on vectors)\n"
        << "7 & 8: matrix size (x, y) (gram matrix has size x*x)\n"
        << "9: window size (used in MAXPOOL)\n"
        << "10: kernel size (used in CNN)\n"
        << "11: number of kernels (used in CNN)\n"
        << "12: repeats (how often to run each function for the benchmark)\n"
        << "(13+: which functions to run (e.g. Multiply, Multiplex), without this, all functions are run)" << endl;
    exit(EXIT_FAILURE);
}

int SuppressStdOut() {
    fflush(stdout);

    int stdout_file_descriptor = dup(1);
    int null_file_descriptor = open("/dev/null", O_WRONLY);
    // check null_file_descriptor for error omitted
    dup2(null_file_descriptor, 1);
    close(null_file_descriptor);

    return stdout_file_descriptor;
}

void ResumeStdOut(int file_descriptor) {
    fflush(stdout);
    dup2(file_descriptor, 1);
    close(file_descriptor);
}


int main(int argc, char* argv[]) {
    if (argc < REQUIRED_ARGS_COUNT) {
        PrintError();
    }
    //parse arguments:
    string role_string(argv[1]);
    Role proxy_role;
    if (role_string == "proxy1") {
        proxy_role = proxy1;
    } else if (role_string == "proxy2") {
        proxy_role = proxy2;
    } else if (role_string == "helper") {
        proxy_role = helper;
    } else {
        PrintError();
    }
    uint16_t helper_port = strtol(argv[2], nullptr, 10);
    string helper_ip(argv[3]);
    uint16_t p1_port = strtol(argv[4], nullptr, 10);
    string p1_ip(argv[5]);
    int vector_length = stoi(argv[6]);
    int matrix_x = stoi(argv[7]);
    int matrix_y = stoi(argv[8]);
    int window_size = stoi(argv[9]);
    int kernel_size = stoi(argv[10]);
    int kernel_count = stoi(argv[11]);
    int repeats = stoi(argv[12]);
    // create output formatting variables:
    int delay;
    string padding;
    switch(proxy_role) {
        case helper:
            delay = 2;
            padding = "";
            break;
        case proxy1:
            delay = 0;
            padding = "    ";
            break;
        case proxy2:
            delay = 1;
            padding = "    ";
            break;
    }
    // initialise proxy:
    int file_descriptor = SuppressStdOut();
    PartyBm* proxy = (proxy_role == helper) ?
                     new PartyBm(helper_port, helper_ip, vector_length, matrix_x, matrix_y, window_size, kernel_size, kernel_count, repeats) :
                     new PartyBm(proxy_role, helper_port, helper_ip, p1_port, p1_ip, vector_length, matrix_x, matrix_y, window_size, kernel_size, kernel_count, repeats);
    ResumeStdOut(file_descriptor);
    // parse/obtain function names to run:
    string* functions;
    int function_count = argc-REQUIRED_ARGS_COUNT-1;
    if (function_count == 0) {
        tie(function_count, functions) = proxy->GetAllFunctionNames();
    } else {
        functions = new string[function_count];
        for (int i = 0; i < function_count; i++) {
            functions[i] = argv[i+REQUIRED_ARGS_COUNT+1];
        }
    }
    // run benchmark:
    double cpu_time, real_time;
    for (int i = 0; i < function_count; i++) {
        bytes_sent = 0;
        bytes_received = 0;
        file_descriptor = SuppressStdOut();
        tie(cpu_time, real_time) = proxy->Benchmark(functions[i]);
        ResumeStdOut(file_descriptor);
        if (cpu_time != -1) {
            if (proxy_role == proxy1) {
                cout << "\n" << functions[i] << "\nReal time:       " << real_time << " ms" << endl;
            }
            if (proxy_role == helper && (functions[i] == "Reconstruct" || functions[i] == "CreateShare")) {
                // skips helper print for functions that don't use it
                continue;
            }
            sleep(delay);
            cout << role_string << " CPU time: " << padding << cpu_time << " ms" << endl;
            cout << role_string << " bytes:    " << padding << (bytes_sent + bytes_received) / repeats << endl;
            sleep(3-delay);
        }
    }
    delete[] functions;
    delete proxy;
    return 0;
}