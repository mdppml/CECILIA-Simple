# PPAURORA

PPAURORA is a C++ project for the private computation of AUROC and AUPRC using secure three-party computation.

## Installation

No installation is required.

Make sure to clone the repository using the "--recurse-submodules" or "--recurse" flag to initialise the submodules as well.

If you already have a local version of this repository without submodules, use the command "git submodule update --init --recursive" to initialise the submodules.

## Compiling

### AUROC

#### helper

```bash
c++ -std=gnu++17 -pthread -W -O3 examples/auroc/proxy.cpp core/Party.cpp core/Party.h utils/constant.h utils/parse_options.cpp utils/parse_options.h utils/connection.h utils/flib.h examples/auroc/llib.h -o proxy_auroc
```

#### proxy

```bash
c++ -std=gnu++17 -pthread -W -O3 examples/auroc/helper.cpp core/Party.cpp core/Party.h utils/constant.h utils/parse_options.cpp utils/parse_options.h utils/connection.h utils/flib.h -o helper_auroc
```

### AUROC WITH TIE

#### helper

```bash
c++ -std=gnu++17 -pthread -W -O3 examples/auroctie/proxy.cpp core/Party.cpp core/Party.h utils/constant.h utils/parse_options.cpp utils/parse_options.h utils/connection.h utils/flib.h examples/auroctie/llib.h -o proxy_auroctie
```

#### proxy

```bash
c++ -std=gnu++17 -pthread -W -O3 examples/auroctie/helper.cpp core/Party.cpp core/Party.h utils/constant.h utils/parse_options.cpp utils/parse_options.h utils/connection.h utils/flib.h -o helper_auroctie
```

### AUPRC

#### helper

```bash
c++ -std=gnu++17 -pthread -W -O3 examples/aupr/proxy.cpp core/Party.cpp core/Party.h utils/constant.h utils/parse_options.cpp utils/parse_options.h utils/connection.h utils/flib.h examples/aupr/llib.h -o proxy_aupr
```

#### proxy

```bash
c++ -std=gnu++17 -pthread -W -O3 examples/aupr/helper.cpp core/Party.cpp core/Party.h utils/constant.h utils/parse_options.cpp utils/parse_options.h utils/connection.h utils/flib.h -o helper_aupr
```

## Usage

```bash
./helper_auroc <ip of helper> <port of helper>
./proxy_auroc role <port of proxy 1> <ip of proxy 1> <port of helper> <ip of helper> <delta> <input>
./proxy_auroc role <port of proxy 1> <ip of proxy 1> <port of helper> <ip of helper> <delta> <input>
```

- input = #input parties,#samples of the first input party,#samples of the second input party,...,#samples of the last input party
- delta = delta is a number that specifies how many selections are made after shuffling

```bash
./helper_auroc "172.31.43.235" 7777
./proxy_auroc 0 8888 "127.0.0.1" 7777 "127.0.0.1" 10 "8,1000,1000,1000,1000,1000,1000,1000,1000"
./proxy_auroc 1 8888 "127.0.0.1" 7777 "127.0.0.1" 10 "8,1000,1000,1000,1000,1000,1000,1000,1000"
```


## CNN
### Window size 2; 4x4 matrix
```bash
./helper_cnn "127.0.0.1" 7777
./proxy_cnn 0 8888 "127.0.0.1" 7777 "127.0.0.1" 2 "0.0,-1.2,3.3,2.0,0.0,-1.6,5.5,3.9,1.3,-1.3,0.2,4.3,2.7,3.7,0.8,-1.1" 4 4
./proxy_cnn 1 8888 "127.0.0.1" 7777 "127.0.0.1" 2 "0.0,-1.2,3.3,2.0,0.0,-1.6,5.5,3.9,1.3,-1.3,0.2,4.3,2.7,3.7,0.8,-1.1" 4 4
``` 
### Window size 2; 28x28 matrix
```bash
./helper_cnn "127.0.0.1" 7777
./proxy_cnn 0 8888 "127.0.0.1" 7777 "127.0.0.1" 2 "2, 4, 7, 2, 7, 5, 5, 8, 6, 2, 9, 3, 1, 6, 6, 2, 0, 6, 0, 2, 6, 4, 4, 7, 5, 8, 9, 1, 0, 4, 8, 5, 1, 1, 7, 2, 7, 9, 5, 2, 3, 5, 4, 1, 8, 0, 7, 8, 8, 3, 0, 0, 6, 3, 5, 5, 4, 8, 3, 9, 8, 4, 2, 8, 9, 8, 6, 5, 9, 5, 9, 0, 6, 2, 7, 6, 1, 7, 8, 9, 2, 7, 0, 9, 2, 9, 7, 0, 4, 8, 2, 7, 2, 5, 0, 2, 1, 2, 8, 6, 3, 5, 2, 2, 2, 2, 0, 4, 3, 3, 5, 1, 8, 3, 5, 8, 6, 3, 5, 5, 6, 2, 9, 7, 0, 4, 3, 3, 3, 9, 5, 4, 8, 4, 2, 5, 8, 6, 8, 5, 3, 0, 7, 5, 9, 4, 2, 1, 3, 0, 4, 2, 1, 1, 6, 3, 4, 6, 3, 5, 4, 2, 7, 8, 5, 5, 1, 0, 1, 2, 4, 0, 3, 7, 9, 1, 8, 1, 0, 3, 4, 4, 2, 8, 9, 5, 7, 5, 4, 4, 9, 5, 5, 9, 0, 8, 0, 6, 8, 9, 7, 8, 1, 6, 5, 7, 9, 1, 9, 0, 7, 0, 2, 8, 8, 1, 5, 1, 7, 1, 8, 4, 6, 2, 5, 2, 8, 3, 3, 6, 2, 0, 8, 0, 3, 8, 4, 7, 2, 3, 2, 0, 1, 6, 3, 7, 3, 5, 2, 9, 5, 9, 3, 2, 1, 3, 5, 0, 3, 1, 6, 1, 0, 7, 0, 6, 6, 7, 3, 4, 9, 4, 5, 6, 4, 3, 3, 3, 3, 8, 1, 2, 3, 5, 5, 4, 0, 8, 0, 8, 1, 6, 0, 1, 0, 2, 0, 9, 5, 3, 0, 0, 3, 2, 3, 0, 0, 2, 9, 2, 8, 1, 9, 6, 8, 9, 7, 1, 9, 9, 0, 2, 9, 7, 5, 6, 0, 1, 6, 1, 9, 6, 5, 2, 0, 5, 5, 3, 9, 4, 8, 1, 2, 3, 4, 0, 5, 6, 4, 5, 6, 1, 5, 3, 7, 6, 0, 8, 4, 0, 1, 5, 6, 0, 3, 2, 6, 3, 2, 8, 4, 9, 5, 8, 2, 6, 4, 8, 4, 6, 2, 3, 3, 3, 3, 5, 2, 4, 9, 1, 1, 9, 4, 2, 5, 5, 8, 6, 2, 3, 2, 6, 2, 1, 1, 1, 6, 5, 4, 8, 6, 9, 5, 4, 1, 5, 6, 9, 5, 9, 6, 4, 0, 5, 1, 5, 5, 4, 1, 7, 3, 4, 3, 6, 2, 6, 8, 5, 6, 0, 5, 6, 9, 6, 8, 5, 9, 7, 5, 3, 0, 9, 5, 0, 3, 6, 3, 8, 4, 7, 3, 3, 2, 7, 1, 7, 0, 3, 1, 5, 6, 7, 6, 6, 4, 7, 8, 8, 9, 9, 9, 8, 6, 3, 9, 5, 8, 9, 6, 8, 9, 1, 7, 3, 7, 8, 0, 7, 9, 7, 8, 2, 5, 8, 0, 9, 8, 7, 0, 2, 4, 5, 6, 5, 6, 2, 2, 2, 8, 3, 0, 6, 6, 4, 7, 4, 5, 7, 4, 2, 5, 0, 7, 2, 9, 7, 0, 3, 9, 4, 1, 6, 5, 4, 9, 3, 8, 3, 7, 4, 5, 2, 0, 7, 7, 6, 9, 4, 5, 7, 0, 4, 0, 2, 1, 1, 0, 9, 4, 4, 8, 8, 4, 3, 8, 1, 4, 1, 4, 2, 1, 5, 1, 6, 9, 7, 1, 3, 3, 1, 7, 9, 0, 5, 9, 3, 6, 6, 2, 9, 2, 3, 8, 5, 0, 3, 9, 7, 4, 4, 3, 5, 2, 1, 2, 2, 0, 5, 8, 1, 5, 9, 9, 3, 7, 0, 0, 7, 6, 1, 4, 3, 3, 0, 2, 9, 4, 8, 9, 1, 3, 7, 6, 7, 8, 7, 2, 4, 4, 9, 0, 8, 7, 6, 1, 4, 7, 8, 6, 0, 1, 8, 4, 5, 3, 5, 9, 1, 5, 4, 0, 4, 1, 8, 0, 3, 7, 4, 5, 9, 2, 9, 4, 5, 4, 3, 7, 3, 9, 7, 8, 6, 4, 8, 1, 2, 9, 1, 3, 5, 9, 9, 4, 7, 0, 2, 9, 8, 5, 3, 0, 9, 3, 1, 9, 1, 2, 8, 1, 5, 3, 9, 8, 1, 2, 0, 7, 2, 1, 3, 7, 3, 0, 1, 7, 2, 0, 3, 7, 0, 8, 3, 2, 5, 8, 1, 4, 4, 3, 2, 7, 5, 1, 0, 9, 8, 6, 8, 1, 4, 5, 7, 9, 2, 7, 4, 6, 9, 4, 2, 1, 7, 2, 4, 8, 6, 8, 2, 4, 7, 5, 6, 7, 7" 28 1
./proxy_cnn 1 8888 "127.0.0.1" 7777 "127.0.0.1" 2 "2, 4, 7, 2, 7, 5, 5, 8, 6, 2, 9, 3, 1, 6, 6, 2, 0, 6, 0, 2, 6, 4, 4, 7, 5, 8, 9, 1, 0, 4, 8, 5, 1, 1, 7, 2, 7, 9, 5, 2, 3, 5, 4, 1, 8, 0, 7, 8, 8, 3, 0, 0, 6, 3, 5, 5, 4, 8, 3, 9, 8, 4, 2, 8, 9, 8, 6, 5, 9, 5, 9, 0, 6, 2, 7, 6, 1, 7, 8, 9, 2, 7, 0, 9, 2, 9, 7, 0, 4, 8, 2, 7, 2, 5, 0, 2, 1, 2, 8, 6, 3, 5, 2, 2, 2, 2, 0, 4, 3, 3, 5, 1, 8, 3, 5, 8, 6, 3, 5, 5, 6, 2, 9, 7, 0, 4, 3, 3, 3, 9, 5, 4, 8, 4, 2, 5, 8, 6, 8, 5, 3, 0, 7, 5, 9, 4, 2, 1, 3, 0, 4, 2, 1, 1, 6, 3, 4, 6, 3, 5, 4, 2, 7, 8, 5, 5, 1, 0, 1, 2, 4, 0, 3, 7, 9, 1, 8, 1, 0, 3, 4, 4, 2, 8, 9, 5, 7, 5, 4, 4, 9, 5, 5, 9, 0, 8, 0, 6, 8, 9, 7, 8, 1, 6, 5, 7, 9, 1, 9, 0, 7, 0, 2, 8, 8, 1, 5, 1, 7, 1, 8, 4, 6, 2, 5, 2, 8, 3, 3, 6, 2, 0, 8, 0, 3, 8, 4, 7, 2, 3, 2, 0, 1, 6, 3, 7, 3, 5, 2, 9, 5, 9, 3, 2, 1, 3, 5, 0, 3, 1, 6, 1, 0, 7, 0, 6, 6, 7, 3, 4, 9, 4, 5, 6, 4, 3, 3, 3, 3, 8, 1, 2, 3, 5, 5, 4, 0, 8, 0, 8, 1, 6, 0, 1, 0, 2, 0, 9, 5, 3, 0, 0, 3, 2, 3, 0, 0, 2, 9, 2, 8, 1, 9, 6, 8, 9, 7, 1, 9, 9, 0, 2, 9, 7, 5, 6, 0, 1, 6, 1, 9, 6, 5, 2, 0, 5, 5, 3, 9, 4, 8, 1, 2, 3, 4, 0, 5, 6, 4, 5, 6, 1, 5, 3, 7, 6, 0, 8, 4, 0, 1, 5, 6, 0, 3, 2, 6, 3, 2, 8, 4, 9, 5, 8, 2, 6, 4, 8, 4, 6, 2, 3, 3, 3, 3, 5, 2, 4, 9, 1, 1, 9, 4, 2, 5, 5, 8, 6, 2, 3, 2, 6, 2, 1, 1, 1, 6, 5, 4, 8, 6, 9, 5, 4, 1, 5, 6, 9, 5, 9, 6, 4, 0, 5, 1, 5, 5, 4, 1, 7, 3, 4, 3, 6, 2, 6, 8, 5, 6, 0, 5, 6, 9, 6, 8, 5, 9, 7, 5, 3, 0, 9, 5, 0, 3, 6, 3, 8, 4, 7, 3, 3, 2, 7, 1, 7, 0, 3, 1, 5, 6, 7, 6, 6, 4, 7, 8, 8, 9, 9, 9, 8, 6, 3, 9, 5, 8, 9, 6, 8, 9, 1, 7, 3, 7, 8, 0, 7, 9, 7, 8, 2, 5, 8, 0, 9, 8, 7, 0, 2, 4, 5, 6, 5, 6, 2, 2, 2, 8, 3, 0, 6, 6, 4, 7, 4, 5, 7, 4, 2, 5, 0, 7, 2, 9, 7, 0, 3, 9, 4, 1, 6, 5, 4, 9, 3, 8, 3, 7, 4, 5, 2, 0, 7, 7, 6, 9, 4, 5, 7, 0, 4, 0, 2, 1, 1, 0, 9, 4, 4, 8, 8, 4, 3, 8, 1, 4, 1, 4, 2, 1, 5, 1, 6, 9, 7, 1, 3, 3, 1, 7, 9, 0, 5, 9, 3, 6, 6, 2, 9, 2, 3, 8, 5, 0, 3, 9, 7, 4, 4, 3, 5, 2, 1, 2, 2, 0, 5, 8, 1, 5, 9, 9, 3, 7, 0, 0, 7, 6, 1, 4, 3, 3, 0, 2, 9, 4, 8, 9, 1, 3, 7, 6, 7, 8, 7, 2, 4, 4, 9, 0, 8, 7, 6, 1, 4, 7, 8, 6, 0, 1, 8, 4, 5, 3, 5, 9, 1, 5, 4, 0, 4, 1, 8, 0, 3, 7, 4, 5, 9, 2, 9, 4, 5, 4, 3, 7, 3, 9, 7, 8, 6, 4, 8, 1, 2, 9, 1, 3, 5, 9, 9, 4, 7, 0, 2, 9, 8, 5, 3, 0, 9, 3, 1, 9, 1, 2, 8, 1, 5, 3, 9, 8, 1, 2, 0, 7, 2, 1, 3, 7, 3, 0, 1, 7, 2, 0, 3, 7, 0, 8, 3, 2, 5, 8, 1, 4, 4, 3, 2, 7, 5, 1, 0, 9, 8, 6, 8, 1, 4, 5, 7, 9, 2, 7, 4, 6, 9, 4, 2, 1, 7, 2, 4, 8, 6, 8, 2, 4, 7, 5, 6, 7, 7" 28 1
```

### Window size 3; 3x9 matrix
```bash
./proxy_cnn 0 8888 "127.0.0.1" 7777 "127.0.0.1" 3 "0.0,-1.2,3.3,2.0,0.0,-1.6,5.5,3.9,1.3,-1.3,0.2,4.3,2.7,3.7,0.8,-1.1,0.9,1.6,1.5,5.2,1.4,3.3,0.4,2.5,0.6,4.6,-1.4" 9 3
./proxy_cnn 1 8888 "127.0.0.1" 7777 "127.0.0.1" 3 "0.0,-1.2,3.3,2.0,0.0,-1.6,5.5,3.9,1.3,-1.3,0.2,4.3,2.7,3.7,0.8,-1.1,0.9,1.6,1.5,5.2,1.4,3.3,0.4,2.5,0.6,4.6,-1.4" 9 3
``` 
### Window size 3; 9x9 matrix
```bash
./proxy_cnn 0 8888 "127.0.0.1" 7777 "127.0.0.1" 3 "0.0,-1.2,3.3,2.0,0.0,-1.6,5.5,3.9,1.3,-1.3,0.2,4.3,2.7,3.7,0.8,-1.1,0.9,1.6,1.5,5.2,1.4,3.3,0.4,2.5,0.6,4.6,-1.4,3.2,-1.2,4.5,4.6,0.4,1.9,5.7,2.6,2.9,-1.3,0.4,2.5,2.1,3.6,5.8,-1.9,0.8,4.7,1.9,1.0,1.2,2.5,5.0,1.6,3.2,4.3,1.6,2.2,1.2,3.3,0.5,1.7,1.1,3.1,0.9,2.2,1.3,0.4,0.8,2.1,3.6,4.4,-1.9,-1.0,-1.4,5.7,1.8,2.5,2.3,0.6,0.1,0.9,4.3,-1.8" 9 9
./proxy_cnn 1 8888 "127.0.0.1" 7777 "127.0.0.1" 3 "0.0,-1.2,3.3,2.0,0.0,-1.6,5.5,3.9,1.3,-1.3,0.2,4.3,2.7,3.7,0.8,-1.1,0.9,1.6,1.5,5.2,1.4,3.3,0.4,2.5,0.6,4.6,-1.4,3.2,-1.2,4.5,4.6,0.4,1.9,5.7,2.6,2.9,-1.3,0.4,2.5,2.1,3.6,5.8,-1.9,0.8,4.7,1.9,1.0,1.2,2.5,5.0,1.6,3.2,4.3,1.6,2.2,1.2,3.3,0.5,1.7,1.1,3.1,0.9,2.2,1.3,0.4,0.8,2.1,3.6,4.4,-1.9,-1.0,-1.4,5.7,1.8,2.5,2.3,0.6,0.1,0.9,4.3,-1.8" 9 9
```

### MAXPOOL TESTS:
```bash
./helper_auroc <ip of helper> <port of helper>
./proxy_auroc role <port of proxy 1> <ip of proxy 1> <port of helper> <ip of helper> <randomness> [<m_cols> <m_rows> <w_cols> <w_rows>] <test>
./proxy_auroc role <port of proxy 2> <ip of proxy 2> <port of helper> <ip of helper> <delta> <input>

randomness - defines if random size matrix is generated: 0=no (4 more params are processed), 1=yes 
test - must be one out of 0,1,2 or 3: 
       0 : subtract elementwise
       1 : resort matrix 
       2 : MAX_VAL
       3 : MMAX 
```

```bash
./test_maxpool_helper "127.0.0.1" 7777
./test_maxpool_proxy 1 8888 "127.0.0.1" 7777 "127.0.0.1" 0 1024 1024 2 2
``` 

## License
[MIT](https://choosealicense.com/licenses/mit/)
