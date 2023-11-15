# CECILIA

CECILIA is a three-party computational framework that offers a variety of building blocks to facilitate more complex algorithms in a privacy preserving manner. It is implemented in C++.

## Table of Contents

- [CECILIA](#cecilia)
    - [Table of Contents](#table-of-contents)
    - [Links to the Related Papers](#links-to-the-related-papers)
    - [Installation](#installation)
    - [Building](#building)
        - [Building with CMake](#building-with-cmake)
        - [Building with C++](#building-with-c)
            - [AUROC](#auroc)
            - [AUROC WITH TIE](#auroc-with-tie)
            - [AUPRC](#auprc)
    - [Usage](#usage)
    - [Building Blocks](#building-blocks)
        - [Core](#core)
        - [Boolean Core](#boolean-core)
    - [Applications](#applications)
        - [CNN](#cnn)
        - [RKN](#rkn)
        - [Sorting](#sorting)
        - [Record Linkage](#record-linkage)
    - [Applications to be added soon](#applications-to-be-added-soon)
    - [License](#license)

## Links to the Related Papers

1. [CECILIA: Comprehensive Secure Machine Learning Framework](https://arxiv.org/abs/2202.03023)

2. [ppAURORA: Privacy Preserving Area Under Receiver Operating Characteristic and Precision-Recall Curves](https://arxiv.org/abs/2102.08788)

## Installation

No installation is required.

Make sure to clone the repository using the "--recurse-submodules" or "--recurse" flag to initialise the submodules as well.

If you already have a local version of this repository without submodules, use the command "git submodule update --init --recursive" to initialise the submodules.

The benchmark bash script requires [toxiproxy](https://github.com/Shopify/toxiproxy/releases/latest). Drop the toxiproxy-server and toxiproxy-cli executables into the same directory as the script.

## Building

### Building with CMake

After cloning the repo into directory `CECILIA`, you can build the library `CECILIA` by executing the following commands.

```bash
mkdir build
cd build
```

```bash
cmake -S ../ -DCMAKE_BUILD_TYPE=Release
```

```bash
make
```

After the build completes, the output binaries can be found in `CECILIA/build/` directory.

### Building with C++

If you want to build it with C++, the required comments are given below for `AUROC`, `AUROC with TIE` and `AUPRC` applications.

`Helper` and `Proxy` sides should be compiled seperately.

Please note that these commands may not be up-to-date.

#### AUROC

**_helper_**

```bash
c++ -std=gnu++17 -pthread -W -O3 apps/proxy.cpp core/Party.cpp core/Party.h utils/constant.h utils/parse_options.cpp utils/parse_options.h utils/connection.h utils/flib.h examples/auroc/llib.h -o proxy_auroc
```

**_proxy_**

```bash
c++ -std=gnu++17 -pthread -W -O3 apps/helper.cpp core/Party.cpp core/Party.h utils/constant.h utils/parse_options.cpp utils/parse_options.h utils/connection.h utils/flib.h -o helper_auroc
```

#### AUROC WITH TIE

**_helper_**

```bash
c++ -std=gnu++17 -pthread -W -O3 examples/auroctie/proxy.cpp core/Party.cpp core/Party.h utils/constant.h utils/parse_options.cpp utils/parse_options.h utils/connection.h utils/flib.h examples/auroctie/llib.h -o proxy_auroctie
```

**_proxy_**

```bash
c++ -std=gnu++17 -pthread -W -O3 apps/helper.cpp core/Party.cpp core/Party.h utils/constant.h utils/parse_options.cpp utils/parse_options.h utils/connection.h utils/flib.h -o helper_auroctie
```

#### AUPRC

**_helper_**

```bash
c++ -std=gnu++17 -pthread -W -O3 examples/aupr/proxy.cpp core/Party.cpp core/Party.h utils/constant.h utils/parse_options.cpp utils/parse_options.h utils/connection.h utils/flib.h examples/aupr/llib.h -o proxy_aupr
```

**_proxy_**

```bash
c++ -std=gnu++17 -pthread -W -O3 apps/helper.cpp core/Party.cpp core/Party.h utils/constant.h utils/parse_options.cpp utils/parse_options.h utils/connection.h utils/flib.h -o helper_aupr
```

## Usage

```bash
./helper_auroc <ip of helper> <port of helper>
./proxy_auroc Role <port of proxy 1> <ip of proxy 1> <port of helper> <ip of helper> <delta> <input>
./proxy_auroc Role <port of proxy 1> <ip of proxy 1> <port of helper> <ip of helper> <delta> <input>
```

- input = #input parties,#samples of the first input party,#samples of the second input party,...,#samples of the last input party
- delta = delta is a number that specifies how many selections are made after shuffling

```bash
./helper_auroc "172.31.43.235" 7777
./proxy_auroc 0 8888 "127.0.0.1" 7777 "127.0.0.1" 10 "8,1000,1000,1000,1000,1000,1000,1000,1000"
./proxy_auroc 1 8888 "127.0.0.1" 7777 "127.0.0.1" 10 "8,1000,1000,1000,1000,1000,1000,1000,1000"
```

## Building Blocks

`CECILIA` has several primitives implemented. The most general ones e.g. Multiplication can be found in `core`. The operations performing on XOR-shares and conversion functions are in `booleancore`. The application spesific functions are usually in their corresponding header files ()
### Core

Functions implemented in `/core` will be explained here.

### Boolean Core

Functions implemented in `/booleancore` will be explained here.

## Applications

### CNN

CNN model and data files are in `CECILIA_data` repo. Download/Clone the folders to CECILIA/apps/cnn/`
```bash
./helper_cnn <ip of helper> <port of helper> <model>
./proxy_cnn Role <port of proxy 1> <ip of proxy 1> <port of helper> <ip of helper> <model> 
./proxy_cnn Role <port of proxy 1> <ip of proxy 1> <port of helper> <ip of helper> <model> 
```

- model = model to be used: [0, 1, 2] with 0 - Chameleon, 1 - MiniONN (MiniONNs model parameter are not consistent), 2 - LeNet5, other value or none - a 4-layer CNN with random weights

```bash
./helper_cnn 127.0.0.1 7777 0
./proxy_cnn 0 8888 127.0.0.1 7777 127.0.0.1 0
./proxy_cnn 1 8888 127.0.0.1 7777 127.0.0.1 0
```

### RKN

Privacy preserving inference on a pre-trained RKN

```bash
./helper <ip of helper> <port of helper>
./proxy_rkn Role <port of proxy 1> <ip of proxy 1> <port of helper> <ip of helper> <random flag> <number of anchor points> <length of kmers> <lambda> <sigma> 
./proxy_rkn Role <port of proxy 1> <ip of proxy 1> <port of helper> <ip of helper> <random flag> <number of anchor points> <length of kmers> <lambda> <sigma> 
```

### Sorting

For the latest development version of Sort function, you need to checkout to the `Sort` branch. There are 3 main `Sort` function currently:

1. `Sort`: The first implementation of the sorting algorithm. It simply sorts the given 1D array and uses `MostSignificantBit` function for bit decomposition.. It is slow and not optimized.
2. `SortNarrow`: The latest version of sorting algorithm. In terms of permutation application the algorithm is different than the original implementation. `ComposePermutation` part is added to avoid the redundant permutation of bits. For the bit decomposition part **XOR-shares** are utilized. `AritmeticToXOR` and `XorToArithmetic` are used for the converisons. Lastly, to reduce the communication cost **ringbits** is introduced and some functions' **_narrow_** versions (they are using ringbits to reduce the number of bits sent) are added.
3. `VectorizedSort`: It gets an n-dimensional array as an input. It sorts/permutes all arrays based on the sorting of pivot array. It is currently implemented using `Sort`. It will be modified to utilze the latest and efficient version of sorting algorithm.

### Record Linkage

The relevant code is under the `recordLinkage` branch.

## Applications to be added soon

1. KNN Inference
2. Heavy Hitters
3. CNN Training

## License

[MIT](https://choosealicense.com/licenses/mit/)
