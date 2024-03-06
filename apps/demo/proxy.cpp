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

    // BEFORE STARTING THIS PART OF THE TUTORIAL:
    // 1. Make sure that FRACTIONAL_BITS is set to 0 in constant.h

    // Create the proxies:
    /*
     * Proxies are the Party instances that interact with the external world. For instance, a data owner outsources its
     * data to the proxies in secret shared form, which we will discuss and see soon.
     */
    Party *proxy;
    if (role == 0)
        proxy = new Party(proxy1, hport, haddress, cport, caddress);
    else
        proxy = new Party(proxy2, hport, haddress, cport, caddress);

    if(FRACTIONAL_BITS == 0) {
        /* Create the data:
         * Let us have a sample data that we will use in our illustrative example of CECILIA2
         */
        int d = 7;

        /* Create secret shares of our sample data:
         * Secret sharing is a fundamental component of multi-party computation where multiple parties work collaboratively
         * to compute a desired function using their secret shares without revealing their secrets to each other and the
         * result of the computation. In literature, there are different secret sharing techniques. In CECILIA2, we use
         * 2-out-of-2 additive secret sharing. In this secret sharing scheme, we generate two values such that their
         * summation over a ring N gives the secret value. What we meant by this is this: Assume that the ring is 16.
         * We sum two values and take their modulus 16. We use unsigned 64-bit integers to hold the shares. To create the
         * shares, we can call CreateShare function of Party instances.
         */
        uint64_t d_secret_share = proxy->CreateShare((double) d);

        // Let us take a look at the share of d in this proxy
        cout << "d_secret_share: " << d_secret_share << endl;

        /*
         * It looks weird and completely unrelated to the actual value, doesn't it? That's actually the whole purpose.
         * The individual shares should not make any sense to the proxy that has this share. However, they should reveal
         * the actual value once they brought together by summing them up. In order to bring the shares together and reveal
         * the actual value, we can call Reconstruct function in core.h. NOTE THAT the proxies should not call Reconstruct
         * over the shares of values as long as they are not masked. This violates the privacy!
         */
        double reconstructed_d = Reconstruct(proxy, d_secret_share);
        cout << "reconstructed_d: " << reconstructed_d << endl;

        /*
         * Now, it seems more meaningful. When the shares are summed up over their special ring, which is 2^64 in this case,
         * we can obtain the actual value.
         */

        /*
         * Now let us see an example of operations that we can perform in CECILIA2 using these shares. For this purpose, let
         * us create another value first and create its shares in proxies.
         */
        int e = 9;
        uint64_t e_secret_shared = proxy->CreateShare((double) e);
        cout << "reconstructed_e: " << Reconstruct(proxy, e_secret_shared) << endl;

        /*
         * Let's start with a simple addition operation that we can perform in local.
         * To calculate the shares of (d+e) we just need to sum them up. No need to call a function
        */
        uint64_t e_plus_d_share =  d_secret_share + e_secret_shared;

        cout << "e_plus_d_share: " << e_plus_d_share << endl;

        /*
         * We can reconstruct the shares and check whether the summation is correct or not
         * */
        double reconstructed_e_plus_d = Reconstruct(proxy, e_plus_d_share);
        cout << "reconstructed_e_plus_d: " << reconstructed_e_plus_d << endl;

        /*
         * Another local operation is scalar multiplication!
         */
        uint64_t two_times_d_share =  2 * d_secret_share ;

        cout << "two_times_d_share: " << two_times_d_share << endl;

        /*
         * Let's check the result and see if it is 2*d
         * */
        double reconstructed_two_times_d = Reconstruct(proxy, two_times_d_share);
        cout << "reconstructed_two_times_d: " << reconstructed_two_times_d << endl;

        /*
         * Let us see how we can multiply these two values using Multiply function in core.h. Note that the function returns
         * a share of the resulting multiplication - not the multiplication itself! For this purpose, we also need to
         * inform the helper, the third computing party of CECILIA2, about the operation that we would like to perform.
         * We will do it by sending the corresponding constant of the operation to the helper.
         */
        proxy->SendBytes(coreMultiply);
        uint64_t ed_secret_shared = Multiply(proxy, e_secret_shared, d_secret_share);


        cout << "ed_secret_shared: " << ed_secret_shared << endl;
        cout << "ed: " << ConvertToDouble(Reconstruct(proxy, ed_secret_shared)) << endl;

        /*
         *Up until now we only work with single values. But usually we have more than a single data point.
         *
         */

        // Compare two secret shared values
        cout << "===================== Compare two secret shared values ===========================" << endl;
        proxy->SendBytes(coreCompare);
        uint64_t de_cmp_secret_shared = Compare(proxy, d_secret_share, e_secret_shared);
        cout << "Compare(proxy, v1, v2) is called" << endl;
        cout << "v1 >=? v2: ";
        if (d >= e)
            cout << "Yes" << endl;
        else
            cout << "No" << endl;
        cout << "de_cmp_secret_shared: " << de_cmp_secret_shared << endl;
        cout << "Reconstructed de_cmp_secret_shared: " << ConvertToDouble(
                Reconstruct(proxy, de_cmp_secret_shared)) << endl;
        cout << "d >= e if the comparison result is 1. 0 otherwise." << endl << endl;

        /*
         * As you may have noticed, the functions that we call perform only a single operation. This means that the
         * multiplication operation, for instance, perform the multiplication of v1 and v2. If I have N number pairs,
         * then I need to call Multiply() function N-times for each pair separately. However, this might be costly.
         * Especially when the proxies and the helper are connected via WAN in which the round trip time is large,
         * these separate calls will require so many communication rounds, i.e. sending and receiving data between
         * parties, and cost too much time due to large round trip time between the parties. One way to handle this is
         * to gather these calls and send/receive their data as bulk. This approach is called "vectorization". CECILIA2
         * offers such optimization. We have the vectorized versions of the functions that we have seen previously.
         * Let us see how we can call these functions' vectorized versions.
         */

        // Let us create a vector of data
        int size = 3;
        double vec1[3] = {4, 2, -2};
        double vec2[3] = {-12, 2, 6};

        uint64_t* vec1_secret_shared = proxy->CreateShare(vec1, 3);
        uint64_t* vec2_secret_shared = proxy->CreateShare(vec2, 3);

        /*
         * Let us notify the Helper about the function that we want to perform.
         * We need to send the size information specifying the size of the vectors of which we want to multiply elements.
         * In this case, in addition to the name of the operation, we send param containing the size information of
         * the vectors and the size of the param vector.
         */
        cout << "===================== Vectorized multiplication of two secret shared vectors ===========================" << endl;
        uint32_t params[1] = {3};
        proxy->SendBytes(coreVectorisedMultiply, params, 1);
        uint64_t* vec1vec2_mul_secret_shared = Multiply(proxy, vec1_secret_shared, vec2_secret_shared, 3);
        double* rec_vec1vec2_mul = ConvertToDouble(Reconstruct(proxy, vec1vec2_mul_secret_shared, 3), 3);
        for (int i = 0; i < 3; i++) {
            cout << "vec1[" << i << "] * " << "vec2[" << i << "]: " << (vec1[i] * vec2[i]) << "\t";
            cout << "rec_vec1vec2_mul[" << i << "]: " << rec_vec1vec2_mul[i] << endl << endl;
        }
        cout << "===========================================================================================\n" << endl;
        // Vectorized Multiplexer
        cout << "===================== Vectorized multiplexer between two secret shared vectors using a secret shared vector of selection bits ===========================" << endl;
        params[0] = 3;
        proxy->SendBytes(coreVectorisedMultiplex, params, 1);
        double double_selection_bits[3] = {0, 1, 1};
        uint64_t* selection_bits = proxy->CreateShare(double_selection_bits, 3);
        uint64_t* vec1vec2_mux_secret_shared = Multiplex(proxy, vec1_secret_shared, vec2_secret_shared, selection_bits, 3);
        double* rec_vec1vec2_mux = ConvertToDouble(Reconstruct(proxy, vec1vec2_mux_secret_shared, 3), 3);
        for (int i = 0; i < 3; i++) {
            cout << "vec1[" << i << "] (" << vec1[i] << ") or vec2[" << i << "] (" << vec2[i] << ") for double_selection_bits[" << i << "] being "
                 << double_selection_bits[i] << ": " << (vec1[i] - double_selection_bits[i] * (vec1[i] - vec2[i])) << "\t";
            cout << "rec_vec1vec2_mux[" << i << "]: " << rec_vec1vec2_mux[i] << endl << endl;
        }
        cout << "===========================================================================================\n" << endl;


    }
        /*
         * As we expected, the secret shares look like garbage. However, the reconstructed value seems not correct! The
         * reason is that our current setting only allows us to represent integers. Now, for the rest of the tutorial,
         * please set FRACTIONAL_BITS to 10 and recompile and rerun the program. We'll wait for you. :)
         */

        // PLEASE MAKE THE REQUEST CHANGE IN THE SETTING FOR THE REST OF THE TUTORIAL
    else {
        /*
         * Now you are back, let us take a look at how our values look like.
         */
        int d = 7;
        uint64_t d_secret_share = proxy->CreateShare((double) d);

        cout << "d: " << d << endl;
        cout << "d_secret_share: " << d_secret_share << endl;
        cout << "d_reconstructed: " << Reconstruct(proxy, d_secret_share) << endl;

        /*
         * Wait? Why does the reconstructed value still look different from the actual value? Let us see their bitwise
         * representations to see if we can find a clue there.
         */
        cout << "bitwise d: " << bitset<64>(d) << endl;
        cout << "bitwise d_reconstructed: " << bitset<64>(Reconstruct(proxy, d_secret_share)) << endl;

        /*
         * The reconstructed value is actually the same as the original one, but only shifted 10 bits to left. This is
         * because we set the FRACTIONAL_BITS to 10. In our number format, we use FRACTIONAL_BITS-many least significant
         * bits to represent the decimal part of the value. The first integer-related bit is the 11th from starting from
         * the least significant bit and the least significant bit being the 1st bit. Since this value is integer,
         * fractional bits are full of 0s. Let us see how they will look like when we work with doubles.
         */
        double f = 7.625;
        uint64_t f_secret_shared = proxy->CreateShare(f);

        cout << "f: " << f << endl;
        cout << "f_secret_shared: " << f_secret_shared << endl;
        cout << "f_reconstructed: " << Reconstruct(proxy, f_secret_shared) << endl;

        /*
         * In bitwise representation of the values, when we go towards the most significant bit, that is towards left,
         * we multiply the value by 2. For instance, let the value be 1111. The least significant one is 1, the second
         * one is 2, the third one is 4 and the most significant one is 8. If we have started from the most significant
         * bit, the pattern is basically dividing by 2. This applies our number format, called fixed-point arithmetic,
         * as well. So, if the 11th bit represents 1, then 12th bit should represent 2, the 13th bit should represent
         * 4 and so on. How about the 9th bit? If we follow the pattern that we describe above, we find out that it
         * should represent 0.5. The 8th bit should represent 0.25, the 7th bit should represent 0.125 and so on. Taking
         * this into account, let us take a look at the reconstructed value again focus on the first 10 bits.
         */
        cout << "f_reconstructed: " << bitset<64>(Reconstruct(proxy, f_secret_shared)) << endl;

        /*
         * As you can see, the 10th bit is 1 and the 8th bit is 1. The rest is 0. So, the overall contribution of the
         * fractional bits to the value is 0.625, which is exactly the decimal part of the value.
         */

        /*
         * How about negative values? CECILIA2 is capable of representing and working with them too. Let's see how they
         * look like.
         */
        int g = -1 * d;
        uint64_t g_secret_share = proxy->CreateShare((double) g);

        cout << "g: " << g << endl;
        cout << "g_secret_share: " << g_secret_share << endl;
        cout << "g_reconstructed: " << Reconstruct(proxy, g_secret_share) << endl;
        cout << "g_reconstructed: " << bitset<64>(Reconstruct(proxy, g_secret_share)) << endl;

        /*
         * When we print out the value itself, it looks too big, but we can see that there are some resemblances with
         * its positive counterpart. You just need to switch all 0s to 1s, all 1s to 0s, and add 1 to the result would
         * give you the positive counterpart. You do not really need to know this conversion, but IT IS IMPORTANT TO
         * KNOW THAT the negative values' most significant bit is 1 and the positive values' most significant bit is 0.
         * In order to convert a value from fixed-point arithmetic number format to real number representation, we can
         * use
         */
        cout << "When we print out the value itself, it looks too big, but we can see that there are some resemblances with\n"
                "its positive counterpart. You just need to switch all 0s to 1s, all 1s to 0s, and add 1 to the result would\n"
                "give you the positive counterpart. You do not really need to know this conversion, but IT IS IMPORTANT TO\n"
                "KNOW THAT the negative values' most significant bit is 1 and the positive values' most significant bit is 0.\n"
                "In order to convert a value from fixed-point arithmetic number format to real number representation, we can\n"
                "use ConvertToDouble function:" << endl;
        cout << "Converted to real number g_reconstructed: " << ConvertToDouble(Reconstruct(proxy, g_secret_share)) << endl;

        /*
         * We also have a function doing the opposite, that is converting from real number format to fixed-point
         * arithmetic number format. NOTE THAT THIS IS CONVERTING THE NUMBER TO FIXED-POINT ARITHMETIC FORMAT WITHOUT
         * CREATING SHARES!
         */
        cout << "Converted to fixed point arithmetic g: " << ConvertToUint64((double) g) << endl;

        /*
         * Now, we know how the values are represented and can be converted to/from our number format, let us dive into
         * the functions provided in CECILIA2! For this purpose, let us create two values.
         */
        double v1 = 2.5;
        double v2 = 6.3;

        uint64_t v1_secret_shared = proxy->CreateShare(v1);
        uint64_t v2_secret_shared = proxy->CreateShare(v2);

        // Addition of two values
        cout << "===================== Add/Subract a scalar from secret shared value via Add function ===========================" << endl;
        uint64_t v1v2_add_secret_shared = Add(proxy, v1_secret_shared, v2_secret_shared);
        cout << "Add(proxy, v1_secret_shared, v2_secret_shared) is called" << endl;
        cout << "v1 + v2: " << (v1 + v2) << endl;
        cout << "v1v2_add_secret_shared: " << v1v2_add_secret_shared << endl;
        cout << "Reconstructed v1v2_add_secret_shared: " << ConvertToDouble(Reconstruct(proxy, v1v2_add_secret_shared)) << endl;
        cout << "===========================================================================================\n" << endl;

        // Addition/Subtraction can also be performed without calling any function
        cout << "===================== Add/Subract two secret shared values without Add function ===========================" << endl;
        cout << "Local v1v2_add_secret_shared: " << (v1_secret_shared + v2_secret_shared) << endl;
        cout << "Reconstructed local v1v2_add_secret_shared: " << ConvertToDouble(
                Reconstruct(proxy, v1_secret_shared + v2_secret_shared)) << endl << endl;

        cout << "v1 - v2: " << (v1 - v2) << endl;
        cout << "Local v1v2_sub_secret_shared: " << (v1_secret_shared - v2_secret_shared) << endl;
        cout << "Reconstructed local v1v2_sub_secret_shared: " << ConvertToDouble(
                Reconstruct(proxy, v1_secret_shared - v2_secret_shared)) << endl;
        cout << "===========================================================================================\n" << endl;

        // Add/Subtract a scalar to/from a secret shared value
        cout << "===================== Add/Subract a scalar from secret shared value ===========================" << endl;
        double add_scalar_double = 4.3;
        uint64_t add_scalar = ConvertToUint64(add_scalar_double);
        uint64_t new_v1_secret_shared = v1_secret_shared;
        if (proxy->GetPRole() == proxy1) {
            new_v1_secret_shared += add_scalar;
        }
        cout << "v1 + scalar: " << (v1 + add_scalar_double) << endl;
        cout << "v1_plus_scalar_add_secret_shared: " << new_v1_secret_shared << endl;
        cout << "Reconstructed v1_plus_scalar_add_secret_shared: " << ConvertToDouble(
                Reconstruct(proxy, new_v1_secret_shared)) << endl << endl;

        new_v1_secret_shared = v1_secret_shared;
        if (proxy->GetPRole() == proxy1) {
            new_v1_secret_shared -= add_scalar;
        }
        cout << "v1 - scalar: " << (v1 - add_scalar_double) << endl;
        cout << "v1_minus_scalar_add_secret_shared: " << new_v1_secret_shared << endl;
        cout << "Reconstructed v1_minus_scalar_add_secret_shared: " << ConvertToDouble(
                Reconstruct(proxy, new_v1_secret_shared)) << endl;
        cout << "===========================================================================================\n" << endl;

        // Multiplication of two values
        cout << "===================== Multiply two secret shared values ===========================" << endl;
        proxy->SendBytes(coreMultiply);
        uint64_t v1v2_mul_secret_shared = Multiply(proxy, v1_secret_shared, v2_secret_shared);
        cout << "Multiply(proxy, v1_secret_shared, v2_secret_shared) is called" << endl;
        cout << "v1 * v2: " << (v1 * v2) << endl;
        cout << "v1v2_mul_secret_shared: " << v1v2_mul_secret_shared << endl;
        cout << "Reconstructed v1v2_mul_secret_shared: " << ConvertToDouble(Reconstruct(proxy, v1v2_mul_secret_shared)) << endl;
        cout << "===========================================================================================\n" << endl;

        // Multiply a secret shared value with a scalar
        cout << "===================== Multiply a secret shared value by a scalar ===========================" << endl;
        double mul_scalar_double = 2.4;
        uint64_t mul_scalar = ConvertToUint64(mul_scalar_double);
        uint64_t v1_times_scalar_mul_secret_shared = LocalMultiply(v1_secret_shared, mul_scalar);
        cout << "LocalMultiply(v1_secret_shared, mul_scalar) is called" << endl;
        cout << "v1 * v2: " << (v1 * mul_scalar_double) << endl;
        cout << "v1_times_scalar_mul_secret_shared: " << v1_times_scalar_mul_secret_shared << endl;
        cout << "Reconstructed v1_times_scalar_mul_secret_shared: " << ConvertToDouble(
                Reconstruct(proxy, v1_times_scalar_mul_secret_shared)) << endl;
        cout << "===========================================================================================\n" << endl;

        // Exponentiation with a known base (i.e. e) and secret shared power
        cout << "===================== Exponentiation of secret shared value with a known base ===========================" << endl;
        proxy->SendBytes(coreExp);
        uint64_t v1_exp_secret_shared = Exp(proxy, v1_secret_shared);
        cout << "Exp(proxy, v1_secret_shared) is called" << endl;
        cout << "exp(v1): " << exp(v1) << endl;
        cout << "v1_exp_secret_shared: " << v1_exp_secret_shared << endl;
        cout << "Reconstructed v1_exp_secret_shared: " << ConvertToDouble(
                Reconstruct(proxy, v1_exp_secret_shared)) << endl;
        cout << "===========================================================================================\n" << endl;

        // Divide two secret shared value
        cout << "===================== Divide two secret shared values ===========================" << endl;
        proxy->SendBytes(coreDivide);
        uint64_t v1v2_div_secret_shared = Divide(proxy, v1, v2);
        cout << "Divide(proxy, v1, v2) is called" << endl;
        cout << "v1 / v2: " << (v1 / v2) << endl;
        cout << "v1v2_div_secret_shared: " << v1v2_div_secret_shared << endl;
        cout << "Reconstructed v1v2_div_secret_shared: " << ConvertToDouble(
                Reconstruct(proxy, v1v2_div_secret_shared)) << endl;
        cout << "===========================================================================================\n" << endl;

        // Compare two secret shared values
        cout << "===================== Compare two secret shared values ===========================" << endl;
        proxy->SendBytes(coreCompare);
        uint64_t v1v2_cmp_secret_shared = Compare(proxy, v1, v2);
        cout << "Compare(proxy, v1, v2) is called" << endl;
        cout << "v1 >=? v2: ";
        if (v1 >= v2)
            cout << "Yes" << endl;
        else
            cout << "No" << endl;
        cout << "v1v2_cmp_secret_shared: " << v1v2_cmp_secret_shared << endl;
        cout << "Reconstructed v1v2_cmp_secret_shared: " << ConvertToDouble(
                Reconstruct(proxy, v1v2_cmp_secret_shared)) << endl;
        cout << "v1 >= v2 if the comparison result is 1. 0 otherwise." << endl << endl;

        proxy->SendBytes(coreCompare);
        uint64_t v2v1_cmp_secret_shared = Compare(proxy, v2, v1);
        cout << "Compare(proxy, v2, v1) is called" << endl;
        cout << "v2 >=? v1: ";
        if (v2 >= v1)
            cout << "Yes" << endl;
        else
            cout << "No" << endl;
        cout << "v2v1_cmp_secret_shared: " << v2v1_cmp_secret_shared << endl;
        cout << "Reconstructed v2v1_cmp_secret_shared: " << ConvertToDouble(
                Reconstruct(proxy, v2v1_cmp_secret_shared)) << endl;
        cout << "v2 >= v1 if the comparison result is 1. 0 otherwise." << endl << endl;

        cout << "!!! The order of the inputs are important in Compare!" << endl;
        cout << "===========================================================================================\n" << endl;

        // Check if two secret shared values are equal

        // Select among two secret shared values
        cout << "===================== Choose between two secret shared values using a secret shared selection bit ===========================" << endl;
        proxy->SendBytes(coreMultiplex);
        uint64_t selection_bit = proxy->CreateShare((double) 1);
        uint64_t v1v2_mux_secret_shared = Multiplex(proxy, v1_secret_shared, v2_secret_shared, selection_bit);
        cout << "Multiplex(proxy, v1_secret_shared, v2_secret_shared, selection_bit) is called" << endl;
        cout << "v1v2_mux_secret_shared: " << v1v2_mux_secret_shared << endl;
        cout << "Reconstructed v1v2_mux_secret_shared: " << ConvertToDouble(
                Reconstruct(proxy, v1v2_mux_secret_shared)) << endl;
        cout << "v1 is selected if the selection bit is 0. v2 otherwise." << endl << endl;

        proxy->SendBytes(coreMultiplex);
        uint64_t v2v1_mux_secret_shared = Multiplex(proxy, v2_secret_shared, v1_secret_shared, selection_bit);
        cout << "Multiplex(proxy, v2_secret_shared, v1_secret_shared, selection_bit) is called" << endl;
        cout << "v2v1_mux_secret_shared: " << v2v1_mux_secret_shared << endl;
        cout << "Reconstructed v2v1_mux_secret_shared: " << ConvertToDouble(
                Reconstruct(proxy, v2v1_mux_secret_shared)) << endl;
        cout << "v2 is selected if the selection bit is 0. v1 otherwise." << endl << endl;

        cout << "!!! Once again: the order of the inputs is important!" << endl;
        cout << "===========================================================================================\n" << endl;

        // Determine the most significant bit of a secret shared value
        cout << "===================== Determine the most significant bit of a secret shared value ===========================" << endl;
        proxy->SendBytes(coreMostSignificantBit);
        uint64_t v1_msb_secret_shared = MostSignificantBit(proxy, v1_secret_shared);
        cout << "MostSignificantBit(proxy, v1_secret_shared) is called" << endl;
        cout << "v1: " << v1 << endl;
        cout << "v1_msb_secret_shared: " << v1_msb_secret_shared << endl;
        cout << "Reconstructed v1_msb_secret_shared: " << ConvertToDouble(
                Reconstruct(proxy, v1_msb_secret_shared)) << endl << endl;

        double opposite_sign_v1 = -1 * v1;
        uint64_t opp_v1_secret_shared = proxy->CreateShare(opposite_sign_v1);
        proxy->SendBytes(coreMostSignificantBit);
        cout << "Opposite sign version of v1: " << opposite_sign_v1 << endl;
        uint64_t opp_v1_msb_secret_shared = MostSignificantBit(proxy, opp_v1_secret_shared);
        cout << "MostSignificantBit(proxy, opp_v1_secret_shared) is called" << endl;
        cout << "opp_v1_msb_secret_shared: " << opp_v1_msb_secret_shared << endl;
        cout << "Reconstructed opp_v1_msb_secret_shared: " << ConvertToDouble(
                Reconstruct(proxy, opp_v1_msb_secret_shared)) << endl;
        cout << "===========================================================================================\n" << endl;

        /*
         * As you may have noticed, the functions that we call perform only a single operation. This means that the
         * multiplication operation, for instance, perform the multiplication of v1 and v2. If I have N number pairs,
         * then I need to call Multiply() function N-times for each pair separately. However, this might be costly.
         * Especially when the proxies and the helper are connected via WAN in which the round trip time is large,
         * these separate calls will require so many communication rounds, i.e. sending and receiving data between
         * parties, and cost too much time due to large round trip time between the parties. One way to handle this is
         * to gather these calls and send/receive their data as bulk. This approach is called "vectorization". CECILIA2
         * offers such optimization. We have the vectorized versions of the functions that we have seen previously.
         * Let us see how we can call these functions' vectorized versions.
         */

        // Let us create a vector of data
        int size = 3;
        double vec1[3] = {4.3, 2.1, -0.8};
        double vec2[3] = {-1.2, 2.0, 5.6};

        uint64_t* vec1_secret_shared = proxy->CreateShare(vec1, 3);
        uint64_t* vec2_secret_shared = proxy->CreateShare(vec2, 3);

        /*
         * Let us notify the Helper about the function that we want to perform.
         * We need to send the size information specifying the size of the vectors of which we want to multiply elements.
         * In this case, in addition to the name of the operation, we send param containing the size information of
         * the vectors and the size of the param vector.
         */
        cout << "===================== Vectorized multiplication of two secret shared vectors ===========================" << endl;
        uint32_t params[1] = {3};
        proxy->SendBytes(coreVectorisedMultiply, params, 1);
        uint64_t* vec1vec2_mul_secret_shared = Multiply(proxy, vec1_secret_shared, vec2_secret_shared, 3);
        double* rec_vec1vec2_mul = ConvertToDouble(Reconstruct(proxy, vec1vec2_mul_secret_shared, 3), 3);
        for (int i = 0; i < 3; i++) {
            cout << "vec1[" << i << "] * " << "vec2[" << i << "]: " << (vec1[i] * vec2[i]) << endl;
            cout << "vec1vec2_mul_secret_shared[" << i << "]: " << vec1vec2_mul_secret_shared[i] << endl;
            cout << "rec_vec1vec2_mul[" << i << "]: " << rec_vec1vec2_mul[i] << endl << endl;
        }
        cout << "===========================================================================================\n" << endl;

        /*
         * The addition/subtraction can be performed using either the Add() function or adding/subtracting the shares
         * of each element of the vectors locally. In case we use Add() function, we do not need to help of the Helper.
         * Therefore, we do not need to send any information to the Helper.
         */
        cout << "===================== Vectorized addition of two secret shared vectors ===========================" << endl;
        uint64_t* vec1vec2_add_secret_shared = Add(proxy, vec1_secret_shared, vec2_secret_shared, 3);
        double* rec_vec1vec2_add = ConvertToDouble(Reconstruct(proxy, vec1vec2_add_secret_shared, 3), 3);
        for (int i = 0; i < 3; i++) {
            cout << "vec1[" << i << "] + " << "vec2[" << i << "]: " << (vec1[i] + vec2[i]) << endl;
            cout << "vec1vec2_add_secret_shared[" << i << "]: " << vec1vec2_add_secret_shared[i] << endl;
            cout << "rec_vec1vec2_add[" << i << "] via Add(): " << rec_vec1vec2_add[i] << endl;
            cout << "rec_vec1vec2_add[" << i << "] via local addition/subtraction: " <<
                 ConvertToDouble(Reconstruct(proxy, vec1_secret_shared[i] + vec2_secret_shared[i])) << endl << endl;
        }
        cout << "===========================================================================================\n" << endl;

        // Vectorized Exponentiation
        cout << "===================== Vectorized exponentiation of a secret shared vector ===========================" << endl;
        params[0] = 3;
        proxy->SendBytes(coreVectorisedExp, params, 1);
        uint64_t* vec1_exp_secret_shared = Exp(proxy, vec1_secret_shared, 3);
        double* rec_vec1_exp = ConvertToDouble(Reconstruct(proxy, vec1_exp_secret_shared, 3), 3);
        for (int i = 0; i < 3; i++) {
            cout << "exp(vec1[" << i << "]): " << exp(vec1[i]) << endl;
            cout << "vec1_exp_secret_shared[" << i << "]: " << vec1_exp_secret_shared[i] << endl;
            cout << "rec_vec1_exp[" << i << "]: " << rec_vec1_exp[i] << endl << endl;
        }
        cout << "===========================================================================================\n" << endl;

        // Vectorized Division
        cout << "===================== Vectorized division of two secret shared vectors ===========================" << endl;
        params[0] = 3;
        proxy->SendBytes(coreVectorisedDivide, params, 1);
        uint64_t* vec1vec2_div_secret_shared = Divide(proxy, vec1_secret_shared, vec2_secret_shared, 3);
        double* rec_vec1vec2_div = ConvertToDouble(Reconstruct(proxy, vec1vec2_div_secret_shared, 3), 3);
        for (int i = 0; i < 3; i++) {
            cout << "vec1[" << i << "] / " << "vec2[" << i << "]: " << (vec1[i] / vec2[i]) << endl;
            cout << "vec1vec2_div_secret_shared[" << i << "]: " << vec1vec2_div_secret_shared[i] << endl;
            cout << "rec_vec1vec2_exp[" << i << "]: " << rec_vec1vec2_div[i] << endl << endl;
        }
        cout << "===========================================================================================\n" << endl;

        // Vectorized Comparison
        cout << "===================== Vectorized comparison of two secret shared vectors ===========================" << endl;
        params[0] = 3;
        proxy->SendBytes(coreVectorisedCompare, params, 1);
        uint64_t* vec1vec2_cmp_secret_shared = Compare(proxy, vec1_secret_shared, vec2_secret_shared, 3);
        double* rec_vec1vec2_cmp = ConvertToDouble(Reconstruct(proxy, vec1vec2_cmp_secret_shared, 3), 3);
        for (int i = 0; i < 3; i++) {
            cout << "vec1[" << i << "] (" << vec1[i] << ") >=? vec2[" << i << "] (" << vec2[i] << "): ";
            if (vec1[i] >= vec2[i])
                cout << "Yes" << endl;
            else
                cout << "No" << endl;
            cout << "vec1vec2_cmp_secret_shared[" << i << "]: " << vec1vec2_cmp_secret_shared[i] << endl;
            cout << "rec_vec1vec2_cmp[" << i << "]: " << rec_vec1vec2_cmp[i] << endl << endl;
        }
        cout << "===========================================================================================\n" << endl;

        // Vectorized Multiplexer
        cout << "===================== Vectorized multiplexer between two secret shared vectors using a secret shared vector of selection bits ===========================" << endl;
        params[0] = 3;
        proxy->SendBytes(coreVectorisedMultiplex, params, 1);
        double double_selection_bits[3] = {0, 1, 1};
        uint64_t* selection_bits = proxy->CreateShare(double_selection_bits, 3);
        uint64_t* vec1vec2_mux_secret_shared = Multiplex(proxy, vec1_secret_shared, vec2_secret_shared, selection_bits, 3);
        double* rec_vec1vec2_mux = ConvertToDouble(Reconstruct(proxy, vec1vec2_mux_secret_shared, 3), 3);
        for (int i = 0; i < 3; i++) {
            cout << "vec1[" << i << "] (" << vec1[i] << ") or vec2[" << i << "] (" << vec2[i] << ") for double_selection_bits[" << i << "] being "
                 << double_selection_bits[i] << ": " << (vec1[i] - double_selection_bits[i] * (vec1[i] - vec2[i])) << endl;
            cout << "vec1vec2_mux_secret_shared[" << i << "]: " << vec1vec2_mux_secret_shared[i] << endl;
            cout << "rec_vec1vec2_mux[" << i << "]: " << rec_vec1vec2_mux[i] << endl << endl;
        }
        cout << "===========================================================================================\n" << endl;

        // Vectorized Most Significant Bit
        cout << "===================== Vectorized most significant bit of a secret shared vector ===========================" << endl;
        params[0] = 3;
        proxy->SendBytes(coreVectorisedMostSignificantBit, params, 1);
        uint64_t* vec1_msb_secret_shared = MostSignificantBit(proxy, vec1_secret_shared, 3);
        double* rec_vec1_msb = ConvertToDouble(Reconstruct(proxy, vec1_msb_secret_shared, 3), 3);
        for (int i = 0; i < 3; i++) {
            cout << "vec1[" << i << "]: " << vec1[i] << endl;
            cout << "vec1_msb_secret_shared[" << i << "]: " << vec1_msb_secret_shared[i] << endl;
            cout << "rec_vec1_msb[" << i << "]: " << rec_vec1_msb[i] << endl << endl;
        }
        cout << "===========================================================================================\n" << endl;
    }

    /*
     * Once the proxies are done with the computation, they need to let the Helper know that the process is over, and
     * it can also terminate. For this, the proxies need to send coreEnd signal to the Helper.
     */
    proxy->SendBytes(coreEnd);

    return 0;
}