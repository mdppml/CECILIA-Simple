#ifndef PARTY_RL
#define PARTY_RL

#include "../../core/Party.h"
#include "../../utils/connection.h"
#include "../../booleancore/core.h"
#include<numeric>
#include "../../core/core.h"

#define ARRAY_SIZE 900

struct FuzzyField{
  uint64_t hasValue;
  uint64_t hammingWeight;
  uint8_t* booleanArray;

  FuzzyField(uint64_t hasValue, uint64_t hammingWeight, uint8_t* booleanArray):
  hasValue(hasValue),
  hammingWeight(hammingWeight),
  booleanArray(booleanArray){}

  ~FuzzyField() {
    delete[] booleanArray;
  }
};

struct ExactField{
  uint64_t hasValue;
  uint64_t value;
};

struct Match{
  uint64_t hasMatch;
  uint64_t matchIndex;
};

struct Record{
  FuzzyField* fuzzyFields;
  ExactField* exactFields;
  ~Record() {
    delete[] fuzzyFields;
    delete[] exactFields;
  }
};

struct Score{
  uint64_t numerator;
  uint64_t denominator;
};

class PartyRL: public Party {
public:
  /**
    *
    * @param role
    * @param helper_port
    * @param helper_ip
    * @param p0_port
    * @param p0_ip
    */
  explicit PartyRL(
    role role,
    uint16_t helperPort,
    const string &helperIP,
    uint16_t p0Port,
    const string &p0IP,
    int exactFieldCount,
    int fuzzyFieldCount,
    uint64_t* exactFieldWeights,
    uint64_t* fuzzyFieldWeights
    ) : Party(role, helperPort, helperIP, p0Port, p0IP) {
      this->exactFieldCount = exactFieldCount;
      this->fuzzyFieldCount = fuzzyFieldCount;
      this->exactFieldWeights = exactFieldWeights;
      this->fuzzyFieldWeights = fuzzyFieldWeights;
      two = createShare((uint64_t) 2);
    }

  int mapBigramToInt(char first_char,  char second_char) {
    return mapCharToInt(first_char) * 30 + mapCharToInt(second_char);
  }

  std::pair<uint8_t*, int> convertStringToBooleanArray(std::string string) {
    auto array = new uint8_t[ARRAY_SIZE]();
    int hammingWeight = 1;
    array[mapBigramToInt(' ', string[0])] = 1;
    int mappedBigram;
    for (int i = 1; i<string.length(); i++) {
      mappedBigram = mapBigramToInt(string[i-1],  string[i]);
      if (array[mappedBigram] == 0) {
        hammingWeight += 1;
        array[mappedBigram] = 1;
      }
    }
    mappedBigram = mapBigramToInt(string[string.length()-1],  ' ');
    if (array[mappedBigram] == 0) {
        hammingWeight += 1;
        array[mappedBigram] = 1;
      }
    return std::make_pair(array,  hammingWeight);
  }

  FuzzyField shareFuzzyField(uint8_t* array, uint64_t hammingWeight, bool hasValue) {
    uint8_t* booleanArray = new uint8_t[ARRAY_SIZE];
    for (int i = 0; i < ARRAY_SIZE; i++) {
      booleanArray[i] = array[i] ^ generateCommonRandomByte();
    }
    FuzzyField field(
      hasValue - generateCommonRandom(),
      hammingWeight - generateCommonRandom(),
      booleanArray
    );
    return field;
  }

  FuzzyField receiveFuzzyField() {
    uint8_t* array = new uint8_t[ARRAY_SIZE];
    for (int i = 0; i < ARRAY_SIZE; i++) {
      array[i] = this->generateCommonRandomByte();
    }
    FuzzyField field(
      generateCommonRandom(),
      generateCommonRandom(),
      array
    );
    return field;
  }

  Score computeDice(FuzzyField field1, FuzzyField field2) {
    Score score;
    uint64_t hasValue = CMP(this, two, field1.hasValue + field2.hasValue);
    score.denominator = MUL(
      this,
      field1.hammingWeight + field2.hammingWeight,
      hasValue
    );
    uint8_t* sharedArray = AND2(
      this,
      field1.booleanArray,
      field2.booleanArray,
      ARRAY_SIZE
    );
    uint64_t* arithmeticArray = XOR2Arithmetic2(this,  sharedArray, ARRAY_SIZE);
    score.numerator = MUL(
      this,
      2 * std::accumulate(arithmeticArray,  arithmeticArray + ARRAY_SIZE, 0),
      hasValue
    );
    return score;
  }

  Score computeExactFieldSimilarity(ExactField field1, ExactField field2) {
    Score score;
    uint64_t hasValue = CMP(this, two, field1.hasValue + field2.hasValue);
    score.numerator = MUL(
      this,
      EQU(this, field1.value, field2.value),
      hasValue
    );
    score.denominator = hasValue;
    return score;
  }

  Score computeRecordSimilarity(Record record1, Record record2) {
    Score totalScore;
    totalScore.numerator = createShare((uint64_t) 0);
    totalScore.denominator = createShare((uint64_t) 0);
    Score fieldScore;
    for (int i = 0; i < exactFieldCount; i++) {
      fieldScore = computeExactFieldSimilarity(record1.exactFields[i], record2.exactFields[i]);
      totalScore.numerator += fieldScore.numerator;
      totalScore.denominator += exactFieldWeights[i]*fieldScore.denominator;
    }
    for (int i = 0; i < fuzzyFieldCount; i++) {
      fieldScore = computeDice(record1.fuzzyFields[i], record2.fuzzyFields[i]);
      totalScore.numerator += fieldScore.numerator;
      totalScore.denominator += fuzzyFieldWeights[i]*fieldScore.denominator;
    }
    return totalScore;
  }

private:
  int exactFieldCount, fuzzyFieldCount;
  uint64_t two;
  uint64_t *exactFieldWeights, *fuzzyFieldWeights;
  int mapCharToInt(char character) {
    if ((97 <= character) & (character <= 122)) {
      return character - 97;
    }
    switch(character) {
      case ' ' :
        return 26;
      case '-' :
        return 27;
      case '\'' :
        return 28;
      default:
        return 29;
    }
  }
};

#endif                                                      // PART_RL
