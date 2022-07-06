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
### Last parameter for proxies defines the model to be used. Can be 0 (Chameleon) or 1 (MiniONN). MiniONNs model parameter are not consistent.
```bash
./helper_cnn "127.0.0.1" 7777
./proxy_cnn 0 8888 "127.0.0.1" 7777 "127.0.0.1" 0
./proxy_cnn 1 8888 "127.0.0.1" 7777 "127.0.0.1" 0
```

```bash
./test_maxpool_helper "127.0.0.1" 7777
./test_maxpool_proxy 1 8888 "127.0.0.1" 7777 "127.0.0.1" 0 1024 1024 2 2
``` 

## License
[MIT](https://choosealicense.com/licenses/mit/)
