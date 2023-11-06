#ifndef TEST_H
#define TEST_H
#include <memory>
#include <utility>
#include <vector>
#include "../../core/Party.h"
#include "../../core/core.h"


/* This is used as a base class for testing. Each component of CECILIA should have its own test class inheriting from this.
 * Children should construct the methods vector in their constructor, containing each test function with the name of the
 * function being tested. The constructor and private method GetGenericMethod in TestCore can be used as an example for
 * how to easily implement this. Each test function has to be of the signature bool () const; the return value denotes
 * if the tested function is broken. By following these steps, it is trivial to test all functions in the class by adding
 * two lines of code in mainTest.cpp
 */
class Test {
public:
  struct Method{
    std::string name;
    bool(Test::*reference)() const;
    Method(std::string name, bool(Test::*reference)() const) {
      this->name = std::move(name);
      this->reference = reference;
    }
  };

  const std::string name_;
  const size_t bit_precision_;

  explicit Test(std::shared_ptr<Party> proxy, std::string name, size_t bit_precision) : name_(name), bit_precision_(bit_precision) {
    this->proxy = proxy;
    if (bit_precision_ == 0 and proxy->GetPRole() == proxy2) {
      std::cout << "Bit precision is set to 0. Make sure to also test the functions with another value." << std::endl;
    }
  }

  const std::vector<Test::Method> *GetMethods() const {
    return &methods;
  };

protected:
  std::shared_ptr<Party> proxy;
  std::vector<Test::Method> methods{};

  /**
   * @brief Allows for reconstruction of 3D matrices directly to double without any memory leaks. Uses the default bit_precision_.
   */
  double*** ReconstructDouble(
    const uint64_t *const *const *const values, size_t matrix_count, size_t row_count, size_t column_count
  ) const {
    if (proxy->GetPRole() != helper) {
      uint64_t*** reconstructed = Reconstruct(proxy.get(), values, matrix_count, row_count, column_count);
      double*** converted = ConvertToDouble(reconstructed, matrix_count, row_count, column_count, bit_precision_);
      for (int i = 0; i < matrix_count; i++) {
        for (int ii = 0; ii < row_count; ii++) {
          delete [] reconstructed[i][ii];
        }
        delete [] reconstructed[i];
      }
      delete [] reconstructed;
      return converted;
    } else {
      return nullptr;
    }
  }

  /**
   * @brief Allows for reconstruction of 2D matrices directly to double without any memory leaks. Uses the default bit_precision_.
   */
  double** ReconstructDouble(const uint64_t *const *const values, size_t row_count, size_t column_count) const {
    uint64_t** reconstructed = Reconstruct(proxy.get(), values, row_count, column_count);
    double** converted = ConvertToDouble(reconstructed, row_count, column_count, bit_precision_);
    for (int i = 0; i < row_count; i++) {
      delete [] reconstructed[i];
    }
    delete [] reconstructed;
    return converted;
  }

  /**
   * @brief Allows for reconstruction of vectors directly to double without any memory leaks. Uses the default bit_precision_.
   */
  std::unique_ptr<double[]> ReconstructDouble(const uint64_t *const values, size_t size) const {
    std::unique_ptr<uint64_t[]> reconstructed(Reconstruct(proxy.get(), values, size));
    return std::unique_ptr<double[]>(ConvertToDouble(reconstructed.get(), size, bit_precision_));
  }

  /**
   * @brief Allows for reconstruction of values directly to double without any memory leaks. Uses the default bit_precision_.
   */
  double ReconstructDouble(uint64_t value) const {
    std::unique_ptr<double[]> result(ReconstructDouble(&value, 1));
    return result[0];
  }

  /**
   * @brief Converts values to their equivalent in fixed point representation.
   *
   * @param bits p_bits: The number of bits used for the decimal precision.
   */
  static int GetPreciseRepresentation(double value, int bits) {
    return (int) (value * std::pow(2, bits));
  }

  /**
   * @brief Converts values to their equivalent in fixed point representation.
   *
   * @param bits p_bits: The number of bits used for the decimal precision.
   */
  static int GetPreciseRepresentation(int value, int bits) {
    return (int) (value * std::pow(2, bits));
  }


  /**
   * @brief Checks if two values are within a margin of error of each other.
   *
   * @param bit_precision p_bit_precision: The number of bits used for the decimal precision.
   * @param margin p_margin: The tolerated error margin (in fixed point representation)
   */
  static bool IsWithinMargin(double value1, double value2, int bit_precision, int margin) {
    int rounded1 = GetPreciseRepresentation(value1, bit_precision);
    int rounded2 = GetPreciseRepresentation(value2, bit_precision);
    return ((rounded2 >= rounded1-margin) and (rounded2 <= rounded1+margin));
  }
};
#endif // TEST_H
