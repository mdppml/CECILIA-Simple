#ifndef PARTY_RL
#define PARTY_RL

#include "../../core/Party.h"
#include "../../utils/connection.h"
#include "../../booleancore/core.h"
#include<numeric>
#include "../../core/core.h"
#include "../../core/cnn.h"

#define ARRAY_SIZE 900

/**
 * @brief A secret shared IDAT field that allows partial matches.
 *
 * String-encoded IDAT fields are converted to a boolean array of bigrams, similar to a Bloom filter.
 * Since we don't want to leak whether a field exists in a record, an empty boolean array is passed for non-existing
 * fields, which is why there is a value {0, 1} to denote whether the field exists.
 * Lastly, the Hamming Weight (i.e. the number of set bits in the array) is also stored to be used for Dice coefficient
 * computation.
 */
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

/**
 * @brief A secret shared IDAT field that only allows exact matches.
 *
 * Any field where we only allow exact matches, e.g. date of birth, is stored as an ExactField.
 * Since we don't want to leak whether a field exists in a record, an empty value is passed for non-existing
 * fields, which is why there is also a value {0, 1} to denote whether the field exists.
 */
struct ExactField{
  uint64_t hasValue;
  uint64_t value;
};

/**
 * @brief A simple pair of secret shared values that is returned by the computeMatch function.
 *
 * To avoid information leaking we always need to return an index even if no match was found. Therefore, there is a
 * second value {0, 1} to denote whether a match was found.
 */
struct Match{
  uint64_t hasMatch;
  uint64_t matchIndex;
};

/**
 * @brief The secret shared IDAT fields of a record. Consists of fields allowing inexact matches and fields only allowing exact matches.
 *
 */
struct Record{
  FuzzyField* fuzzyFields;
  ExactField* exactFields;
  ~Record() {
    delete[] fuzzyFields;
    delete[] exactFields;
  }
};

/**
 * @brief A secret shared match score between two records.
 *
 * To avoid having to compute too many expensive divisions, the numerator and denominator are stored separately.
 * The only time the actual score has to be computed is to determine whether the match threshold was passed.
 */
struct Score{
  uint64_t numerator;
  uint64_t denominator;
};

/**
 * @brief The secret shared maximum score and its (secret shared) index.
 *
 */
struct MaxScore{
  Score score;
  uint64_t index;
};

/**
 * @brief An extension of Party implementing the required methods to compute record linkage.
 *
 * Since most of the functions required for record linkage are very specific, they are simply member methods.
 */
class PartyRL: public Party {
public:
  /**
   * @brief Construct a Party for record linkage.
   *
   * @param role p_role: The role of the party (P1, P2 or helper)
   * @param helperPort p_helperPort: The port on which the helper is listening
   * @param helperIP p_helperIP: The IP address of the helper
   * @param p1Port p_p1Port: The port on which P1 is listening
   * @param p1IP p_p1IP: The IP address of P1
   * @param exactFieldCount p_exactFieldCount: The number of ExactFields each record has
   * @param fuzzyFieldCount p_fuzzyFieldCount: The number of FuzzyFields each record has
   * @param exactFieldWeights p_exactFieldWeights: The weights of each ExactField for computing the similarity score
   * @param fuzzyFieldWeights p_fuzzyFieldWeights: The weights of each FuzzyField for computing the similarity score
   */
  explicit PartyRL(
    role role,
    uint16_t helperPort,
    const string &helperIP,
    uint16_t p1Port,
    const string &p1IP,
    int exactFieldCount,
    int fuzzyFieldCount,
    double* exactFieldWeights,
    double* fuzzyFieldWeights
    ) : Party(role, helperPort, helperIP, p1Port, p1IP) {
      this->exactFieldCount = exactFieldCount;
      this->fuzzyFieldCount = fuzzyFieldCount;
      if (role == P1 || role == P2) {
        this->exactFieldWeights = new uint64_t[exactFieldCount];
        for (int i = 0; i < exactFieldCount; i++) {
          this->exactFieldWeights[i] = convert2uint64(exactFieldWeights[i]);
        }
        this->fuzzyFieldWeights = new uint64_t[fuzzyFieldCount];
        for (int i = 0; i < fuzzyFieldCount; i++) {
          this->fuzzyFieldWeights[i] = convert2uint64(fuzzyFieldWeights[i]);
        }
        two = createShare((uint64_t) 2);
      }
    }

    /**
     * @brief Returns a unique integer ranging from 0 to 899 for each possible bigram of the alphabet considered.
     *
     * We only consider the letters a-z, "-", " " and "'". Every other symbol is considered identical.
     * Unicode characters should be converted to ASCII using unidecode in advance.
     * @param first_char p_first_char: The first character of the bigram
     * @param second_char p_second_char: The second character of the bigram
     * @return int
     */
    int mapBigramToInt(char first_char,  char second_char) {
    return mapCharToInt(first_char) * 30 + mapCharToInt(second_char);
  }

  /**
   * @brief Computes a boolean array similar to a Bloom filter from the given string.
   *
   * A uint8_t array of size ARRAY_SIZE is zero-initialised.
   * Then, a predefined integer i is computed for each bigram contained in the string and the value at the ith index of the array is set to 1.
   * The bigrams also include space with the first character and the last character with space.
   *
   * @param string p_string: The string to be converted
   * @return std::pair< uint8_t*, int > A pointer to the array and its hamming weight (number of indices set to 1).
   */
  std::pair<uint8_t*, int> convertStringToBooleanArray(std::string string) {
    auto array = new uint8_t[ARRAY_SIZE]();
    int hammingWeight = 1; // this is set to one so the first bigram does not need to increment it
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

  /**
   * @brief Share a fuzzy field with the other proxy.
   *
   * @param array p_array: the boolean array of the fuzzy field
   * @param hammingWeight p_hammingWeight: The hamming weight of the boolean array
   * @param hasValue p_hasValue: A value in {0, 1} denoting whether the field actually exists
   * @return FuzzyField
   */
  FuzzyField shareFuzzyField(uint8_t* array, uint64_t hammingWeight, bool hasValue) {
    if (getPRole() == P1 || getPRole() == P2) {
      auto* booleanArray = new uint8_t[ARRAY_SIZE];
    for (int i = 0; i < ARRAY_SIZE; i++) {
      booleanArray[i] = array[i] ^ generateCommonRandomByte();
    }
    FuzzyField field(
      hasValue - generateCommonRandom(),
      hammingWeight - generateCommonRandom(),
      booleanArray
    );
    return field;
    } else {
      FuzzyField field(0, 0, nullptr);
      return field;
    }
  }

  /**
   * @brief Receive a fuzzy field from the other proxy.
   *
   * @return FuzzyField
   */
  FuzzyField receiveFuzzyField() {
    if (getPRole() == P1 || getPRole() == P2) {
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
    } else {
      FuzzyField field(0, 0, nullptr);
      return field;
    }
  }

  /**
   * @brief Compute the Dice coefficient (i.e. similarity score) between two FuzzyField s.
   *
   * @param field1 p_field1: The first FuzzyField
   * @param field2 p_field2: The second FuzzyField
   * @return Score
   */
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

  /**
   * @brief Compute the similarity score between two ExactField s. The score is either 1 in case of an exact match or 0.
   *
   * @param field1 p_field1: The first ExactField
   * @param field2 p_field2: The second ExactField
   * @return Score
   */
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

  /**
   * @brief Compute the similarity score between two records.
   *
   * @param record1 p_record1: The first record
   * @param record2 p_record2: The second record
   * @return Score
   */
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
    char lowerCase = tolower(character);
    if ((97 <= lowerCase) & (lowerCase <= 122)) {
      return lowerCase - 97;
    }
    switch(lowerCase) {
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

  /**
   * @brief Returns 0 if score1 is smaller than score2, 1 otherwise
   *
   * @param score1 p_score1:...
   * @param score2 p_score2:...
   * @return uint64_t
   */
  uint64_t compareScores(Score score1, Score score2) {
    uint64_t compare1 = MUL(this, score1.numerator, score2.denominator);
    uint64_t compare2 = MUL(this, score2.numerator, score1.denominator);
    return CMP(this, compare1, compare2);
  }

  uint64_t* compareScores(
    uint64_t* numerator1,
    uint64_t* numerator2,
    uint64_t* denominator1,
    uint64_t* denominator2,
    int count) {
    uint64_t* compare1 = MUL(this, numerator1, denominator2, count);
    uint64_t* compare2 = MUL(this, numerator2, denominator1, count);
    return CMP(this, compare1, compare2, count);
  }

  /**
  * Selects the maximum score from an array of shared Score s.
  * @param scores an array of shared scores
  * @param scoreCount the size of the score array
  * @return The maximum score which was found in scores.
  */
  MaxScore computeMaxScore(Score *scores, uint32_t scoreCount){
    /** MAIN IDEA:
     Repeatedly compare the first half of the remaining scores to the second
     half, using MUX to select the larger of the two scores.
     Have an array contain secretly shared indices that are selected
     with the same selection vector to also get the argmax.
     */
    bool remainderExists = false;
    int arrayLength = scoreCount;
    int halfLength;

    if ((getPRole() == P1) || getPRole() == P2) {
      // generate the initial data structures and values:
      Score remainder;
      uint64_t remainderIndex;
      uint64_t* numerators = new uint64_t[scoreCount];
      uint64_t* denominators = new uint64_t[scoreCount];
      uint64_t* indices = new uint64_t[scoreCount];
      uint64_t* selection;
      // the tmp pointers are used to avoid memory leaks at various points
      uint64_t *tmpNumerators, *tmpDenominators, *tmpIndices;
      for (int i = 0; i < scoreCount; i++) {
        numerators[i] = scores[i].numerator;
        denominators[i] = scores[i].denominator;
        indices[i] = createShare((uint64_t) i);
      }
      //iteratively halve the array by comparing the first with the second half:
      while (arrayLength > 1) {
        if (arrayLength % 2 == 1) {
          if (remainderExists) {
            /* we create new arrays with arrayLength+1 and store the remainder
               in the last index:*/
            tmpNumerators = numerators;
            tmpDenominators = denominators;
            tmpIndices = indices;
            numerators = new uint64_t[arrayLength+1];
            denominators = new uint64_t[arrayLength+1];
            indices = new uint64_t[arrayLength+1];
            memcpy(numerators, tmpNumerators, arrayLength);
            memcpy(denominators, tmpDenominators, arrayLength);
            memcpy(tmpIndices, indices, arrayLength);
            delete[] tmpDenominators;
            delete[] tmpNumerators;
            delete[] tmpIndices;
            numerators[arrayLength] = remainder.numerator;
            denominators[arrayLength] = remainder.denominator;
            indices[arrayLength] = remainderIndex;
            remainderExists = false;
            arrayLength++;
          } else { // the remainder is stored for later:
            remainder.denominator = denominators[arrayLength-1];
            remainder.numerator = numerators[arrayLength-1];
            remainderIndex = indices[arrayLength-1];
            remainderExists = true;
          }
        }
        halfLength = arrayLength / 2;
        selection = compareScores(
          numerators,
          numerators + halfLength,
          denominators,
          denominators+halfLength,
          halfLength
        );
        tmpNumerators = numerators;
        tmpDenominators = denominators;
        tmpIndices = indices;
        numerators = MUX(this, numerators, numerators+halfLength, selection, halfLength);
        denominators = MUX(this, denominators, denominators+halfLength, selection, halfLength);
        indices = MUX(this, indices, indices+halfLength, selection, halfLength);
        delete[] tmpNumerators;
        delete[] tmpDenominators;
        delete[] tmpIndices;
        delete[] selection;
        arrayLength = halfLength;
      }
      MaxScore max;
      max.score.numerator = numerators[0];
      max.score.denominator = denominators[0];
      max.index = indices[0];
      if (remainderExists) {
        uint64_t selection = compareScores(max.score, remainder);
        max.score.numerator = MUX(this, max.score.numerator, remainder.numerator, selection);
        max.score.denominator = MUX(this, max.score.denominator, remainder.denominator, selection);
        max.index = MUX(this, indices[0], remainderIndex, selection);
      }
      return max;
    } else { // HELPER
      //iteratively halve the array by comparing the first with the second half:
      while (arrayLength > 1) {
        if (arrayLength % 2 == 1) {
          if (remainderExists) {
            /* we create new arrays with arrayLength+1 and store the remainder
               in the last index:*/
            remainderExists = false;
            arrayLength++;
          } else { // the remainder is stored for later:
            remainderExists = true;
          }
        }
        halfLength = arrayLength / 2;
        compareScores(nullptr, nullptr, nullptr, nullptr, halfLength);
        MUX(this, nullptr, nullptr, nullptr, halfLength);
        MUX(this, nullptr, nullptr, nullptr, halfLength);
        MUX(this, nullptr, nullptr, nullptr, halfLength);
        arrayLength = halfLength;
      }
      if (remainderExists) {
        Score score;
        compareScores(score, score);
        MUX(this, 0, 0, 0);
        MUX(this, 0, 0, 0);
        MUX(this, 0, 0, 0);
      }
      MaxScore max;
      return max;
    }


  }
};

#endif                                                      // PART_RL
