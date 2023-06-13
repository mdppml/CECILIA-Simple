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
