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

static long long reseed_interval = 1LL << 48;

class AES_CTR_RBG : public RandomNumberGenerator, public NotCopyable
{
public:
    explicit AES_CTR_RBG(const CryptoPP::byte *seed = nullptr, size_t length = 0)
    : m_pCipher(new CTR_Mode<AES>::Encryption)
    {
        m_keyed = false;
        EntropyHelper(seed, length, true);
    }
    
    [[nodiscard]] bool CanIncorporateEntropy() const override
    {
        return true;
    }

    /**\brief Reseed the generator
     * @param input provided seed
     * @param length should be at least 32 for AES-128
     */
    void IncorporateEntropy(const CryptoPP::byte *input, size_t length) override
    {
        EntropyHelper(input, length, false);
    }

    // Does not keep track of whether the cipher has to be reseeded.
    //   Therefore, this must be done outside of this class.
    void GenerateBlock(CryptoPP::byte *output, size_t size) override
    {
        if (!m_keyed) {
            m_pCipher->SetKeyWithIV(m_key, m_key.size(), m_iv, m_iv.size());
            m_keyed = true;
        }
        m_pCipher->GenerateBlock(output, size);
    }

protected:
    // Sets up to use the cipher. It's a helper to allow a throw
    //   in the constructor during initialization.
    void EntropyHelper(const CryptoPP::byte* input, size_t length, bool ctor = false)
    {
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
        m_keyed = false;
    }
    
private:
    FixedSizeSecBlock<CryptoPP::byte, 16> m_key;
    FixedSizeSecBlock<CryptoPP::byte, 16> m_iv;
    member_ptr<CTR_Mode<AES>::Encryption> m_pCipher;
    bool m_keyed;
};
