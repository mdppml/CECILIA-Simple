#ifndef TEST_H
#define TEST_H
#include <memory>
#include <utility>
#include <vector>
#include "../../core/Party.h"
#include "../../core/core.h"

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

  std::unique_ptr<double[]> ReconstructDouble(const uint64_t *const values, size_t size) const {
    std::unique_ptr<uint64_t[]> reconstructed(Reconstruct(proxy.get(), values, size));
    return std::unique_ptr<double[]>(ConvertToDouble(reconstructed.get(), size, bit_precision_));
  }

  double ReconstructDouble(uint64_t value) const {
    std::unique_ptr<double[]> result(ReconstructDouble(&value, 1));
    return result[0];
  }

  static int GetPreciseRepresentation(double value, int bits) {
    return (int) (value * std::pow(2, bits));
  }

  static int GetPreciseRepresentation(int value, int bits) {
    return (int) (value * std::pow(2, bits));
  }

  static bool IsWithinMargin(double value1, double value2, int bit_precision, int margin) {
    int rounded1 = GetPreciseRepresentation(value1, bit_precision);
    int rounded2 = GetPreciseRepresentation(value2, bit_precision);
    return ((rounded2 >= rounded1-margin) and (rounded2 <= rounded1+margin));
  }
};
#endif // TEST_H
