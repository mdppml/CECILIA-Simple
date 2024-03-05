//
// Created by Seyma Selcan on 16.10.23.
//
#include <cstdlib>
#include <iostream>
#include <chrono>
#include <tuple>
#include <iomanip>
#include <bitset>
#include "../../core/core.h"

constexpr uint64_t sz = 4;


int main(int argc, char* argv[]) {
    uint8_t role = atoi(argv[1]);
    uint16_t cport = atoi(argv[2]);
    string caddress(argv[3]);
    uint16_t hport = atoi(argv[4]);
    string haddress(argv[5]);

    /** BEFORE STARTING THIS PART OF THE TUTORIAL:
     * 1. Check the demo file in apps.demo.proxy.cpp
     * Some todos here are already written in the demo file. You can copy paste those parts.
     * For the rest, it is not exactly the same but demo file will give some insight to you.
     * /


    /*
     * Proxies are the Party instances that interact with the external world. For instance, a data owner outsources its
     * data to the proxies in secret shared form, which we will discuss and see soon.
     */
    // TODO: Create the proxies




    /* TODO: Create the random number x:
     */
    int x;

    /*
     * Secret sharing is a fundamental component of multi-party computation where multiple parties work collaboratively
     * to compute a desired function using their secret shares without revealing their secrets to each other and the
     * result of the computation. In literature, there are different secret sharing techniques. In CECILIA2, we use
     * 2-out-of-2 additive secret sharing. In this secret sharing scheme, we generate two values such that their
     * summation over a ring N gives the secret value. What we meant by this is this: Assume that the ring is 16.
     * We sum two values and take their modulus 16. We use unsigned 64-bit integers to hold the shares. To create the
     * shares, we can call CreateShare function of Party instances.
     */
    //TODO: Create secret shares of our sample data:
    uint64_t x_share;


     //TODO: Now calculate the function 5x+18



    //TODO: Let's reconstruct the resulting share


    /**
     * TODO: Now let us now lets calculate the square of our variable x
     * Note that this is not a local operation we also need to
     * inform the helper, the third computing party of CECILIA2, about the operation that we would like to perform.
     * */


    /*
     * Great, we can do operations on a single variable.
     * But usually we have more than one data point.
     * TODO: Create 2 random arrays of length 4 each
     * TODO: And create their shares
    */


    /*
     * TODO: Calculate the element-wise multiplication of these two arrays using Multiply function
     * Try the vectorized version, do not multiply them one by one
     * */

    /*
     * If we can multiply these two arrays element-wise, we can also calculate the dot product of them!
     * TODO: Calculate the dot product of these two arrays.
     * Reconstruct and check if it is correct
     */



    //TODO: Create another array of length $5$. This time add some negative integers into the list. Use the MostSignificantBit tp detect if those numbers are negative or not.

    /*
     * Once the proxies are done with the computation, they need to let the Helper know that the process is over, and
     * it can also terminate. For this, the proxies need to send coreEnd signal to the Helper.
     * Remove the comment from the following line
     */
    //proxy->SendBytes(coreEnd);

    return 0;
}