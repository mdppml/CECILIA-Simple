#include <cstdlib>
#include <iostream>
#include "Test.h"
#include "TestCore.h"
#include "../../core/Party.h"

static void PrintError() {
    std::clog
        << "Required arguments not provided. These are the required (and optional) arguments:\n"
        << "1: Role (proxy1, proxy2 or helper)\n"
        << "2 & 3: helper's port and IP\n"
        << "4 & 5: proxy1's port and IP\n" << std::endl;
}

static void TestMethods(const Test *const test, Role proxy_role) {
    bool all_working = true;
    bool is_broken;
    for (const auto& method : *test->GetMethods()) {
        is_broken = (test->*method.reference)();
        if (is_broken and (proxy_role == proxy2)) {
            all_working = false;
            std::cout << method.name << " is not working correctly." << std::endl;
        }
    }
    if (all_working and (proxy_role == proxy2)) {
        std::cout << "All tested functions of " << test->name_ << " are working!" << std::endl;
    }
    if (proxy_role != proxy2) {
        std::cout << "Tests for " << test->name_ << " finished execution. Check proxy2 for the results." << std::endl;
    }
}

int main(int argc, char* argv[]) {
    const int REQUIRED_ARGS_COUNT = 5;
    if (argc < REQUIRED_ARGS_COUNT) {
        PrintError();
    }
     //parse arguments:
    std::string role_string(argv[1]);
    Role proxy_role;
    if (role_string == "proxy1") {
        proxy_role = proxy1;
    } else if (role_string == "proxy2") {
        proxy_role = proxy2;
    } else if (role_string == "helper") {
        proxy_role = helper;
    } else {
        PrintError();
        exit(EXIT_FAILURE);
    }
    uint16_t helper_port = strtol(argv[2], nullptr, 10);
    std::string helper_ip(argv[3]);
    uint16_t proxy1_port = strtol(argv[4], nullptr, 10);
    std::string proxy1_ip(argv[5]);

    //Core
    std::shared_ptr<Party> proxy(new Party(proxy_role, helper_port, helper_ip, proxy1_port, proxy1_ip));
    TestCore core(proxy, 5);
    TestMethods(&core, proxy_role);
}
