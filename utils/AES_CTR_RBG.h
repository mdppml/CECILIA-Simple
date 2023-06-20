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

static long long RESEED_INTERVAL = 1LL << 48;
static int RANDOM_BUFFER_SIZE = BUFFER_SIZE;

class AES_CTR_RBG : public RandomNumberGenerator, public NotCopyable
{
public:
    explicit AES_CTR_RBG(const CryptoPP::byte *seed = nullptr, size_t length = 0)
    : m_pCipher(new CTR_Mode<AES>::Encryption)
    {
        initialised = false;
        EntropyHelper(seed, length, true);
        _current_buffer = new CryptoPP::byte[RANDOM_BUFFER_SIZE];
        _unused_buffer = new CryptoPP::byte[RANDOM_BUFFER_SIZE];
    }

    ~AES_CTR_RBG() override {
        delete[] _current_buffer;
        delete[] _unused_buffer;
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

    /**Does not keep track of whether the cipher has to be reseeded. Therefore, this must be done outside of this class.
     * @param output the buffer where the block should be stored
     * @param size the number of bytes to generate
     */
    void GenerateBlock(CryptoPP::byte *output, size_t size) override
    {
        initialise();
        size_t remaining_bytes = size;
        size_t transferred_bytes = 0;
        size_t bytes_in_buffer = RANDOM_BUFFER_SIZE - _buffer_position;
        while (bytes_in_buffer < remaining_bytes) {
            ::memcpy(output + transferred_bytes, _current_buffer + _buffer_position, bytes_in_buffer);
            transferred_bytes += bytes_in_buffer;
            remaining_bytes -= bytes_in_buffer;
            _replenishBuffer();
            bytes_in_buffer = RANDOM_BUFFER_SIZE;
        }
        ::memcpy(output + transferred_bytes, _current_buffer + _buffer_position, remaining_bytes);
        _buffer_position += remaining_bytes;
    }

    CryptoPP::byte GenerateByte() override
    {
        initialise();
        if (_buffer_position == RANDOM_BUFFER_SIZE) {
            _replenishBuffer();
        }
        CryptoPP::byte byte = _current_buffer[_buffer_position];
        _buffer_position += 1;
        return byte;
    }

    uint64_t GenerateLongLong() {
        initialise();
        if (_buffer_position+8 > RANDOM_BUFFER_SIZE) {
            _replenishBuffer();
        }
        uint64_t random = *(uint64_t *)(_current_buffer +_buffer_position);
        _buffer_position += 8;
        return random;
    }

    /**\brief makes sure that everything is initialised
     *
     */
    void initialise() {
        if (!initialised) {
            m_pCipher->SetKeyWithIV(m_key, m_key.size(), m_iv, m_iv.size());
            _current_buffer = new CryptoPP::byte[BUFFER_SIZE];
            _unused_buffer = new CryptoPP::byte[BUFFER_SIZE];
            _fillUnusedBuffer();
            _replenishBuffer();
            initialised = true;
        }
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
        initialised = false;
    }
    
private:
    std::thread _buffer_thread{};
    CryptoPP::byte* _current_buffer;
    CryptoPP::byte* _unused_buffer;
    size_t _buffer_position = 0;
    FixedSizeSecBlock<CryptoPP::byte, 16> m_key;
    FixedSizeSecBlock<CryptoPP::byte, 16> m_iv;
    member_ptr<CTR_Mode<AES>::Encryption> m_pCipher;
    bool initialised;

    void _replenishBuffer() {
        if (_buffer_thread.joinable()) {
            _buffer_thread.join();
        }
        CryptoPP::byte* swap = _current_buffer;
        _current_buffer = _unused_buffer;
        _unused_buffer = swap;
        _buffer_position = 0;
        _buffer_thread = thread(&AES_CTR_RBG::_fillUnusedBuffer, this);
    }

    void _fillUnusedBuffer() {
        m_pCipher->GenerateBlock(_unused_buffer, RANDOM_BUFFER_SIZE);
    }
};
