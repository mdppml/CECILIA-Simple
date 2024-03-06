# CECILIA

CECILIA is a three-party computational framework that offers a variety of building blocks to facilitate more complex algorithms in a privacy preserving manner. It is implemented in C++.

<!-- TOC -->
* [CECILIA](#cecilia)
  * [Links to the Related Papers](#links-to-the-related-papers)
  * [Installation](#installation)
  * [Building](#building)
    * [Requirements](#requirements)
    * [Building with CMake](#building-with-cmake)
  * [Usage](#usage)
  * [Building Blocks](#building-blocks)
    * [Core: Function Descriptions](#core--function-descriptions)
  * [License](#license)
<!-- TOC -->

## Links to the Related Papers

1. [CECILIA: Comprehensive Secure Machine Learning Framework](https://arxiv.org/abs/2202.03023)

2. [ppAURORA: Privacy Preserving Area Under Receiver Operating Characteristic and Precision-Recall Curves](https://arxiv.org/abs/2102.08788)

## Installation

No installation is required.

Make sure to clone the repository using the "--recurse-submodules" or "--recurse" flag to initialise the submodules as well.
```bash
git clone --recurse-submodules <repo-link>
```

If you already have a local version of this repository without submodules, use the command "git submodule update --init --recursive" to initialise the submodules.

## Building
### Requirements 

In order to build the framework you need cmake, gcc, g++, and make should be installed on your system.
Install them before continuing to building. 

### Building with CMake

After cloning the repo into directory `CECILIA-Simple`, you can build the library `CECILIA` by executing the following commands.
First, go to the directory of the library:
```bash
cd CECILIA #or CECILIA-Simple
```
Create a directory for the executables and move to that directory.
```bash
mkdir build
cd build
```
Build the library using cmake
```bash
cmake -S ../ -DCMAKE_BUILD_TYPE=Release
```
Build the executables
```bash
make
```

After the build completes, the output binaries can be found in `CECILIA/build/` directory.

## Usage
There 3 parties: Helper, Proxy0 and Proxy1
You have to run their executables separately.
```bash
./helper <ip of helper> <port of helper>
./proxy Role <port of proxy 1> <ip of proxy 1> <port of helper> <ip of helper> 
./proxy Role <port of proxy 1> <ip of proxy 1> <port of helper> <ip of helper> 
```
If you are working on your local machine you need these commands to run demo code
```bash
./helper "127.0.0.1" 7777
./demo 0 8888 "127.0.0.1" 7777 "127.0.0.1" 
./demo 1 8888 "127.0.0.1" 7777 "127.0.0.1" 
```

You can run the exercise code with the following commands
If you are working on your local machine you need these commands to run demo code
```bash
./helper "127.0.0.1" 7777
./exercise 0 8888 "127.0.0.1" 7777 "127.0.0.1" 
./exercise 1 8888 "127.0.0.1" 7777 "127.0.0.1" 
```

## Building Blocks

`CECILIA` has several primitives implemented. The most general ones e.g. Multiplication can be found in `core`. The operations performing on XOR-shares and conversion functions are in `booleancore`. The application spesific functions are usually in their corresponding header files ()
### Core: Function Descriptions

- **`Reconstruct`**: There are several overloaded versions of this function, each designed to reconstruct different data types (e.g., single values, arrays, 2D arrays, and 3D arrays) from their secret-shared representations. Reconstruction is a critical step in MPC, allowing parties to obtain the actual computation result from its secret-shared form without compromising data privacy.
  
    Parameters:
  - `a`: The secret-shared data to be reconstructed, can be a single value, a pointer to an array, or pointers to 2D/3D arrays.
  - `sz`: Size of the array (for array overloads).

- **`GenerateMultiplicationTriple`**: Used in the pre-processing phase of MPC to generate multiplication triples (Beaver triples), which are key to performing secure multiplications. Multiplication triples consist of three numbers (a, b, c) such that a * b = c, and they are distributed among parties in secret-shared form.

- **`Multiply`**: Facilitates the multiplication of secret-shared values using pre-generated multiplication triples, following the technique introduced by Beaver for secure multiplications in MPC frameworks.
  Parameters:
    - `a`: A secret-shared multiplicant, can be a single value, a pointer to an array.
    - `b`: Another secret-shared multiplicant, can be a single value, a pointer to an array.
    - `sz`: Size of the array (for array overloads).

- **`Exp`**: Computes the exponential function over secret-shared inputs. This function is particularly challenging in MPC due to the non-linear nature of the exponential operation, requiring careful handling to preserve data privacy.

- **`DotProduct`**: Calculates the dot product of two vectors represented in secret-shared form. This operation is fundamental in many machine learning algorithms, including support vector machines and neural networks.

- **`MostSignificantBit` (MSB)**: Determines the most significant bit of a secret-shared value, a crucial operation in many cryptographic protocols and algorithms, such as comparisons and division.
    - `a`: A secret-shared data, can be a single value, a pointer to an array.
    - `sz`: Size of the array (for array overloads).

- **`Compare`**: A secure comparison function that determines whether one secret-shared number is greater than another. 
    - `a`: A secret-shared data to be compared, can be a single value, a pointer to an array.
    - `b`: Another secret-shared data to be compared, can be a single value, a pointer to an array.
    - `sz`: Size of the array (for array overloads).

## License

[MIT](https://choosealicense.com/licenses/mit/)
