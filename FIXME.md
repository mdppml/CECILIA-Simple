# PPAURORA

Document the problems to think about at a later state here. If particular file(s) are affected by the problem, mention those at the beginning. The style for problem documentation is:
##### <affected file(s)>
-> problem description

### Computation Specification
##### Programm start
-> define a buffer size depending on machine the program is running on
##### Files processing vectors (multiple elements)
-> what is the maximum number of elements a matrix can contain without getting overflow problems during computations?


### Structure
##### proxy.cpp in apps/test
-> Should we include different folders for test of different applications or include all apps_tests there?
##### utils
-> Creating separate files for addVal2..., convert... and other functions in flib.h ?

### Current Limitations
##### cnn.cpp in core/test
-> MMAX can handle only symmetric matrices for resorting by now (calculation of how many columns to process is not for all cases)