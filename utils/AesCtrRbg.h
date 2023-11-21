#ifndef AES_CTR_RBG_H
#define AES_CTR_RBG_H

#include <cryptopp/secblock.h>
using CryptoPP::AlignedSecByteBlock;
using CryptoPP::FixedSizeSecBlock;

#include <cryptopp/smartptr.h>
using CryptoPP::member_ptr;

#include <cryptopp/osrng.h>
using CryptoPP::OS_GenerateRandomBlock;
using CryptoPP::RandomNumberGenerator;

#include <cryptopp/aes.h>
using CryptoPP::AES;

#include <cryptopp/ccm.h>
using CryptoPP::CTR_Mode;

#include <cryptopp/sha.h>
using CryptoPP::SHA512;

#include <cryptopp/misc.h>
using CryptoPP::NotCopyable;

#include <cryptopp/config_int.h>
using CryptoPP::lword;

#include <thread>

static const long long kReseedInterval = 1LL << 48;
static const int kRandomBufferSize = 4000000;

class AesCtrRbg : public RandomNumberGenerator, public NotCopyable
{
public:
    explicit AesCtrRbg(const CryptoPP::byte *seed = nullptr, size_t length = 0)
    : m_pCipher64(new CTR_Mode<AES>::Encryption), m_pCipher8(new CTR_Mode<AES>::Encryption){
        initialised = false;
        EntropyHelper(seed, length, true);
        Initialise();
    }

    ~AesCtrRbg() override {
        delete[] current_buffer64;
        delete[] unused_buffer64;
        delete[] current_buffer8;
        delete[] unused_buffer8;
    }

    [[nodiscard]] bool CanIncorporateEntropy() const override {
        return true;
    }

    /**\brief Reseed the generator
     * @param input provided seed
     * @param length should be at least 32 for AES-128
     */
    void IncorporateEntropy(const CryptoPP::byte *input, size_t length) override {
        EntropyHelper(input, length, false);
    }

    // Does not keep track of whether the cipher has to be reseeded.
    //   Therefore, this must be done outside of this class.
    void GenerateBlock64(CryptoPP::byte *output, size_t size)
    {
        m_pCipher64->GenerateBlock(output,size*8);
    }

    void GenerateBlock8(CryptoPP::byte *output, size_t size)
    {
        m_pCipher8->GenerateBlock(output,size);
    }

    uint8_t GenerateByte() {
        if (buffer_position8 == kRandomBufferSize) {
            ReplenishBuffer8();
        }
        return current_buffer8[buffer_position8++];
    }
    uint64_t GenerateUint64() {
       if (buffer_position64 == kRandomBufferSize) {
            ReplenishBuffer64();
       }
       return current_buffer64[buffer_position64++];
    }

    /**\brief makes sure that everything is initialised
     *
     */
    void Initialise() {
        if (!initialised) {
            m_pCipher64->SetKeyWithIV(m_key, m_key.size(), m_iv, m_iv.size());
            m_pCipher8->SetKeyWithIV(m_key, m_key.size(), m_iv, m_iv.size());
            current_buffer64 = new uint64_t[kRandomBufferSize];
            unused_buffer64 = new uint64_t[kRandomBufferSize];
            FillUnusedBuffer64();
            ReplenishBuffer64();
            current_buffer8 = new uint8_t[kRandomBufferSize];
            unused_buffer8 = new uint8_t[kRandomBufferSize];
            FillUnusedBuffer8();
            ReplenishBuffer8();
            initialised = true;
        }
    }

protected:
    // Sets up to use the cipher. It's a helper to allow a throw
    //   in the constructor during initialization.
    void EntropyHelper(const CryptoPP::byte* input, size_t length, bool ctor = false) {
        if(ctor)
        {
            memset(m_key, 0x00, m_key.size());
            memset(m_iv, 0x00, m_iv.size());
        }
        // 16-byte key, 16-byte nonce
        AlignedSecByteBlock seed(16 + 16);
        SHA512 hash;
        if(input && length)
        {
            // Use the user supplied seed.
            hash.Update(input, length);
        }
        else
        {
            // No seed or size. Use the OS to gather entropy.
            OS_GenerateRandomBlock(false, seed, seed.size());
            hash.Update(seed, seed.size());
        }
        hash.Update(m_key.data(), m_key.size());
        hash.Update(m_iv.data(), m_iv.size());
        hash.TruncatedFinal(seed.data(), seed.size());
        memcpy(m_key.data(), seed.data() + 0, 16);
        memcpy(m_iv.data(), seed.data() + 16, 16);
        initialised = false;
    }

private:
    std::thread buffer_thread64{};
    std::thread buffer_thread8{};
    uint64_t * current_buffer64;
    uint64_t* unused_buffer64;
    uint8_t * current_buffer8;
    uint8_t* unused_buffer8;
    size_t buffer_position64 = 0;
    size_t buffer_position8 = 0;
    FixedSizeSecBlock<CryptoPP::byte, 16> m_key;
    FixedSizeSecBlock<CryptoPP::byte, 16> m_iv;
    member_ptr<CTR_Mode<AES>::Encryption> m_pCipher64;
    member_ptr<CTR_Mode<AES>::Encryption> m_pCipher8;
    bool initialised;

    void ReplenishBuffer64() {
        if (buffer_thread64.joinable()) {
            buffer_thread64.join();
        }
        uint64_t * swap = current_buffer64;
        current_buffer64 = unused_buffer64;
        unused_buffer64 = swap;
        buffer_position64 = 0;
        buffer_thread64 = thread(&AesCtrRbg::FillUnusedBuffer64, this);
    }

    void FillUnusedBuffer64() {
        GenerateBlock64((uint8_t*)unused_buffer64, kRandomBufferSize);
    }

    void ReplenishBuffer8() {
        if (buffer_thread8.joinable()) {
            buffer_thread8.join();
        }
        uint8_t * swap = current_buffer8;
        current_buffer8 = unused_buffer8;
        unused_buffer8 = swap;
        buffer_position8 = 0;
        buffer_thread8 = thread(&AesCtrRbg::FillUnusedBuffer8, this);
    }

    void FillUnusedBuffer8() {
        GenerateBlock8(unused_buffer8, kRandomBufferSize);
    }
};

#endif // AES_CTR_RBG_H
