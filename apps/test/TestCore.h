#ifndef TEST_CORE_H
#define TEST_CORE_H
#include "Test.h"
#include "../../core/Party.h"
#include "../../utils/test_functions.h"

class TestCore : public Test {
public:
  explicit TestCore(std::shared_ptr<Party> proxy, size_t bit_precision): Test(proxy, "core", bit_precision) {
    methods = std::vector<Method>{
      GetGenericMethod("Compare", &TestCore::TestCompare),
      GetGenericMethod("MostSignificantBit", &TestCore::TestMostSignificantBit),
      GetGenericMethod("Equals", &TestCore::TestEquals),
      GetGenericMethod("Multiply", &TestCore::TestMultiply)
    };
  }

  bool TestCompare() const {
    bool is_broken = false;
    size_t count = 100;
    if (proxy->GetPRole() != helper) {
      std::unique_ptr<double[]> values1(Random1dData(nullptr, count, -100.0, 100.0));
      std::unique_ptr<double[]> values2(Random1dData(nullptr, count, -100.0, 100.0));
      std::unique_ptr<uint64_t[]> shared1(proxy->CreateShare(values1.get(), count, bit_precision_));
      std::unique_ptr<uint64_t[]> shared2(proxy->CreateShare(values2.get(), count, bit_precision_));
      std::unique_ptr<uint64_t[]> compare(Compare(proxy.get(), shared1.get(), shared2.get(), count, bit_precision_));
      std::unique_ptr<double[]> result(ReconstructDouble(compare.get(), count));
      for (int i = 0; i < count; i++) {
        if ((GetPreciseRepresentation(values1[i], bit_precision_) > GetPreciseRepresentation(values2[i], bit_precision_)) != result[i]) {
          is_broken = true;
        }
      }
    } else {
      Compare(proxy.get(), nullptr, nullptr, count, bit_precision_);
    }
    return is_broken;
  }

  bool TestMultiply() const {
    bool is_broken = false;
    size_t count = 100;
    if (proxy->GetPRole() != helper) {
      std::unique_ptr<double[]> values1(Random1dData(nullptr, (size_t) count, -100.0, 100.0));
      std::unique_ptr<double[]> values2(Random1dData(nullptr, (size_t) count, -100.0, 100.0));
      std::unique_ptr<uint64_t[]> shared1(proxy->CreateShare(values1.get(), count, bit_precision_));
      std::unique_ptr<uint64_t[]> shared2(proxy->CreateShare(values2.get(), count, bit_precision_));
      std::unique_ptr<uint64_t[]> multiplied(Multiply(proxy.get(), shared1.get(), shared2.get(), count, bit_precision_));
      std::unique_ptr<double[]> result(ReconstructDouble(multiplied.get(), count));
      int precision_result;
      if (bit_precision_ == 0) {
        precision_result = 0;
      } else {
        precision_result = bit_precision_/2-3;
      }
      for (int i = 0; i < count; i++) {
        if (!IsWithinMargin(values1[i] * values2[i], result[i], precision_result, 3)) {
          is_broken = true;
          std::cout<< std::format("{}*{}={}!={}", values1[i], values2[i], values1[i] * values2[i],  result[i]) << std::endl;
          std::cout << std::format("{}!={}\n", GetPreciseRepresentation(values1[i] * values2[i], precision_result), GetPreciseRepresentation(result[i], precision_result));
        }
      }
    } else {
      Multiply(proxy.get(), nullptr, nullptr, count, bit_precision_);
    }
    return is_broken;
  }

  bool TestMostSignificantBit() const {
    bool is_broken = false;
    size_t score_count = 100;
    if (proxy->GetPRole() != helper) {
      std::unique_ptr<double[]> values(Random1dData(nullptr, score_count, -100.0, 100.0));
      std::unique_ptr<uint64_t[]> shared(proxy->CreateShare(values.get(), score_count, bit_precision_));
      std::unique_ptr<uint64_t[]> msb_shared(MostSignificantBit(proxy.get(), shared.get(), score_count, bit_precision_));
      std::unique_ptr<double[]> results (ReconstructDouble(msb_shared.get(), score_count));
      for (int i = 0; i < score_count; i++) {
        if ((values[i] < 0) != results[i]) {
          is_broken = true;
        }
      }
    } else { // helper
      MostSignificantBit(proxy.get(), nullptr, score_count, bit_precision_);
    }
    return is_broken;
  }

  bool TestEquals() const {
    bool is_broken = false;
    if (proxy->GetPRole() != helper) {
      std::unique_ptr<uint64_t[]> test1 = std::make_unique_for_overwrite<uint64_t[]>(4);
      std::unique_ptr<uint64_t[]> test2 = std::make_unique_for_overwrite<uint64_t[]>(4);
      test1[0] = proxy->CreateShare(0.0, bit_precision_);
      test1[1] = proxy->CreateShare(13.0, bit_precision_);
      test1[2] = proxy->CreateShare(13.0, bit_precision_);
      test1[3] = proxy->CreateShare(2008.0, bit_precision_);
      test2[0] = proxy->CreateShare(0.0, bit_precision_);
      test2[1] = proxy->CreateShare(13.0, bit_precision_);
      test2[2] = proxy->CreateShare(5.0, bit_precision_);
      test2[3] = proxy->CreateShare(2.5, bit_precision_);
      std::unique_ptr<uint64_t[]> result(Equals(proxy.get(), test1.get(), test2.get(), 4, bit_precision_));
      std::unique_ptr<double[]> reconstructed(ReconstructDouble(result.get(), 4));
      if (
        reconstructed[0] != 1.0
        or reconstructed[1] != 1.0
        or reconstructed[2] != 0.0
        or reconstructed[3] != 0.0
      ) {
        is_broken = true;
      }
    } else {
      Equals(proxy.get(), nullptr, nullptr, 4, bit_precision_);
    }
    return is_broken;
  }

private:
  static Method GetGenericMethod(std::string method_name, bool(TestCore::*method_reference)() const) {
    return Method(method_name, static_cast<bool(Test::*)() const>(method_reference));
  }
};

#endif // TEST_CORE_H
