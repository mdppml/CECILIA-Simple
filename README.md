# PPAURORA

PPAURORA is a C++ project for the private computation of AUROC and AUPRC using secure three-party computation.

## Installation

No installation is required.

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
       2 : MAX
       3 : MMAX 
```

```bash
./test_maxpool_helper "127.0.0.1" 7777
./test_maxpool_proxy 1 8888 "127.0.0.1" 7777 "127.0.0.1" 0 1024 1024 2 2
``` 

## License
[MIT](https://choosealicense.com/licenses/mit/)
