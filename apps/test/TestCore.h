#ifndef TEST_CORE_H
#define TEST_CORE_H
#include "Test.h"
#include "../../core/Party.h"
#include "../../utils/test_functions.h"

/* This class should test all functions of core.h. Still missing (or tests broken) are:
 * Exp
 * MatrixMatrixMultiply
 * MatrixVectorMultiply
 * ModularInverse
 * Divide
 * Normalise
 */
class TestCore : public Test {
public:
  explicit TestCore(std::shared_ptr<Party> proxy, size_t bit_precision): Test(proxy, "core", bit_precision) {
    methods = std::vector<Method>{
      GetGenericMethod("CreateShare or Reconstruct",              &TestCore::TestShareAndReconstruct),
      GetGenericMethod("Matrix CreateShare or Reconstruct",       &TestCore::TestMatrixReconstruct),
      GetGenericMethod("3D Matrix CreateShare or Reconstruct",    &TestCore::Test3dReconstruct),
      GetGenericMethod("GenerateMultiplicationTriple",            &TestCore::TestGenerateMultiplicationTriple),
      GetGenericMethod("ModularConversion or PrivateCompareBool", &TestCore::TestModularConversion),
      GetGenericMethod("Multiplex",                               &TestCore::TestMultiplex),
      GetGenericMethod("Compare",                                 &TestCore::TestCompare),
      GetGenericMethod("MostSignificantBit",                      &TestCore::TestMostSignificantBit),
      GetGenericMethod("Equals",                                  &TestCore::TestEquals),
      GetGenericMethod("Multiply",                                &TestCore::TestMultiply),
      GetGenericMethod("DotProduct",                              &TestCore::TestDotProduct)
    };
  }

  bool TestShareAndReconstruct() const {
    if (proxy->GetPRole() != helper) {
      size_t count = 100;
      std::unique_ptr<double[]> values(Random1dData(nullptr, count, (double) -100.0, (double) 100.0));
      std::unique_ptr<uint64_t[]> shared(proxy->CreateShare(values.get(), count, bit_precision_));
      std::unique_ptr<double[]> reconstructed(ReconstructDouble(shared.get(), count));
      for (int i = 0; i < count; i++) {
        if (
          GetPreciseRepresentation(values[i], bit_precision_)
          != GetPreciseRepresentation(reconstructed[i], bit_precision_)
        ) {
          // debug info can be printed here
          return true;
        }
      }
    }
    return false;
  }

  bool TestMatrixReconstruct() const {
    bool is_broken = false;
    if (proxy->GetPRole() != helper) {
      size_t columns[]{50, 100};
      size_t rows[]{100, 50};
      for (int i = 0; i < 2; i++) {
        double** values = Random2dData(proxy.get(), rows[i], columns[i], -100.0, 100.0);
        uint64_t** shared = proxy->CreateShare(values, rows[i], columns[i], bit_precision_);
        std::unique_ptr<double[]> reconstructed;
        for (int row = 0; row < rows[i]; row++) {
          reconstructed = ReconstructDouble(shared[row], columns[i]);
          for (int column = 0; column < columns[i]; column++) {
            if (
              GetPreciseRepresentation(values[row][column], bit_precision_)
              != GetPreciseRepresentation(reconstructed[column], bit_precision_)
            ) {
              // debug info can be printed here
              is_broken = true;
            }
          }
          delete[] values[row];
          delete[] shared[row];
        }
        delete[] values;
        delete[] shared;
      }
    }
    return is_broken;
  }

  bool Test3dReconstruct() const {
    bool is_broken = false;
    if (proxy->GetPRole() != helper) {
      size_t matrices[]{10, 30, 10, 20, 20};
      size_t columns[]{20, 20, 30, 30, 10};
      size_t rows[]{30, 10, 20, 10, 30};
      for (int i = 0; i < 5; i++) {
        double*** values = Random3dData(proxy.get(), matrices[i], rows[i], columns[i], -100.0, 100.0);
        uint64_t*** shared = new uint64_t**[matrices[i]];
        for (int matrix = 0; matrix < matrices[i]; matrix++) {
          shared[matrix] = proxy->CreateShare(values[matrix], rows[i], columns[i], bit_precision_);
        }
        double*** reconstructed = ReconstructDouble(shared, matrices[i], rows[i], columns[i]);
        for (int matrix = 0; matrix < matrices[i]; matrix++) {
          for (int row = 0; row < rows[i]; row++) {
            for (int column = 0; column < columns[i]; column++) {
              if (
                GetPreciseRepresentation(values[matrix][row][column], bit_precision_)
                != GetPreciseRepresentation(reconstructed[matrix][row][column], bit_precision_)
              ) {
                // debug info can be printed here
                is_broken = true;
              }
            }
            delete[] values[matrix][row];
            delete[] reconstructed[matrix][row];
            delete[] shared[matrix][row];
          }
          delete[] values[matrix];
          delete[] reconstructed[matrix];
          delete[] shared[matrix];
        }
        delete[] values;
        delete[] reconstructed;
        delete[] shared;
      }
    }
    return is_broken;
  }

  bool TestGenerateMultiplicationTriple() const {
    if (proxy->GetPRole() == proxy2) {
      size_t count = 100;
      uint64_t **share1 = new uint64_t*[3];
      uint64_t **share2 = new uint64_t*[3];
      for (int i = 0; i < 3; i++) {
        share1[i] = new uint64_t[count];
        share2[i] = new uint64_t[count];
      }
      uint64_t result1, result2;
      GenerateMultiplicationTriple(proxy.get(), share1, share2, count);
      for (int i = 0; i < count; i++) {
        result1 = (share1[0][i] + share2[0][i]) * (share1[1][i] + share2[1][i]);
        result2 = share1[2][i] + share2[2][i];
        if (result1 != result2) {
          // debug info can be printed here
          return true;
        }
      }
    }
    return false;
  }

  bool TestMultiplex() const {
    size_t count = 100;
    if (proxy->GetPRole() != helper) {
      std::unique_ptr<double[]> values1(Random1dData(proxy.get(), count));
      std::unique_ptr<double[]> values2(Random1dData(proxy.get(), count));
      std::unique_ptr<double[]> selection = std::make_unique_for_overwrite<double[]>(count);
      for (int i = 0; i < count; i++) {
        selection[i] = rand() % 2;
      }
      std::unique_ptr<uint64_t[]> shared1(proxy->CreateShare(values1.get(), count, bit_precision_));
      std::unique_ptr<uint64_t[]> shared2(proxy->CreateShare(values2.get(), count, bit_precision_));
      std::unique_ptr<uint64_t[]> shared_selection(proxy->CreateShare(selection.get(), count, bit_precision_));
      std::unique_ptr<uint64_t[]> shared_selected(
        Multiplex(proxy.get(), shared1.get(), shared2.get(), shared_selection.get(), count, bit_precision_)
      );
      std::unique_ptr<double[]> selected(ReconstructDouble(shared_selected.get(), count));
      if (proxy->GetPRole() == proxy2) {
        double selected_value;
        for (int i = 0; i < count; i++) {
          if (selection[i] == 0) {
            selected_value = values1[i];
          } else {
            selected_value = values2[i];
          }
          if (
            GetPreciseRepresentation(selected[i], bit_precision_)
            != GetPreciseRepresentation(selected_value, bit_precision_)
          ) {
            // debug info can be printed here
            return true;
          }
        }
      }
    } else {
      Multiplex(proxy.get(), nullptr, nullptr, nullptr, count, bit_precision_);
    }
    return false;
  }

  bool TestModularConversion() const {
    size_t count = 100;
    if (proxy->GetPRole() != helper) {
      std::unique_ptr<uint64_t[]> shared = std::make_unique_for_overwrite<uint64_t[]>(count);
      for (int i = 0; i < count; i++) {
        shared[i] = proxy->GenerateRandom() & N1_MASK;
      }
      std::unique_ptr<uint64_t[]> shared_converted(ModularConversion(proxy.get(), shared.get(), count));
      std::unique_ptr<uint64_t[]> values(Reconstruct(proxy.get(), shared.get(), count, N1_MASK));
      std::unique_ptr<uint64_t[]> converted(Reconstruct(proxy.get(), shared_converted.get(), count));
      for (int i = 0; i < count; i++) {
        if (values[i] != converted[i]) {
          // debug info can be printed here
          return true;
        }
      }
    } else {
      ModularConversion(proxy.get(), nullptr, count);
    }
    return false;
  }

  bool TestCompare() const {
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
          // debug info can be printed here
          return true;
        }
      }
    } else {
      Compare(proxy.get(), nullptr, nullptr, count, bit_precision_);
    }
    return false;
  }

  bool TestMultiply() const {
    size_t count = 100;
    if (proxy->GetPRole() != helper) {
      std::unique_ptr<double[]> values1(Random1dData(nullptr, count, -100.0, 100.0));
      std::unique_ptr<double[]> values2(Random1dData(nullptr, count, -100.0, 100.0));
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
          // debug info can be printed here
          return true;
        }
      }
    } else {
      Multiply(proxy.get(), nullptr, nullptr, count, bit_precision_);
    }
    return false;
  }

  bool TestMostSignificantBit() const {
    size_t score_count = 100;
    bool is_broken = false;
    if (proxy->GetPRole() != helper) {
      std::unique_ptr<double[]> values(Random1dData(nullptr, score_count, -100.0, 100.0));
      std::unique_ptr<uint64_t[]> shared(proxy->CreateShare(values.get(), score_count, 0));
      std::unique_ptr<uint64_t[]> msb_shared(MostSignificantBit(proxy.get(), shared.get(), score_count, 0));
      std::unique_ptr<double[]> results (ReconstructDouble(msb_shared.get(), score_count));
      for (int i = 0; i < score_count; i++) {
        if ((GetPreciseRepresentation(values[i], bit_precision_) < 0) != results[i]) {
          // debug info can be printed here
          is_broken = true;
        }
      }
    } else { // helper
      MostSignificantBit(proxy.get(), nullptr, score_count, bit_precision_);
    }
    return is_broken;
  }

  bool TestEquals() const {
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
        return true;
      }
    } else {
      Equals(proxy.get(), nullptr, nullptr, 4, bit_precision_);
    }
    return false;
  }

  bool TestDotProduct() const {
    size_t vector_size = 10;
    size_t vector_count = 10;
    size_t total_size = vector_size * vector_count;
    int test_precision = bit_precision_/2-4;
    if (proxy->GetPRole() != helper) {
      std::unique_ptr<double[]> vector1(Random1dData(proxy.get(), total_size));
      std::unique_ptr<uint64_t[]> shared1(proxy->CreateShare(vector1.get(), total_size, bit_precision_));
      std::unique_ptr<double[]> vector2(Random1dData(proxy.get(), total_size));
      std::unique_ptr<uint64_t[]> shared2(proxy->CreateShare(vector2.get(), total_size, bit_precision_));
      std::unique_ptr<uint64_t[]> shared_result(
        DotProduct(proxy.get(), shared1.get(), shared2.get(), total_size, vector_size, bit_precision_)
      );
      std::unique_ptr<double[]> result(ReconstructDouble(shared_result.get(), vector_count));
      if (proxy->GetPRole() == proxy2) {
        double sum;
        for (int i = 0; i < vector_count; i++) {
          sum = 0;
          for (int ii = 0; ii < vector_size; ii++) {
            sum += vector1[i*vector_size+ii] * vector2[i*vector_size+ii];
          }
          if (GetPreciseRepresentation(sum, test_precision) != GetPreciseRepresentation(result[i], test_precision)) {
            // debug info can be printed here
            return true;
          }
        }
      }
    } else {
      DotProduct(proxy.get(), nullptr, nullptr, total_size, vector_size, bit_precision_);
    }
    return false;
  }

  // TODO fix and add to the methods vector
  bool TestExp() const {
    size_t count = 100;
    bool is_broken = false;
    if (proxy->GetPRole() != helper) {
      std::unique_ptr<double[]> values(
        Random1dData(proxy.get(), count, proxy->GetMinPower() + 10, proxy->GetMaxPower() - 10)
      );
      std::unique_ptr<uint64_t[]> shared(proxy->CreateShare(values.get(), count, bit_precision_));
      std::unique_ptr<uint64_t[]> shared_result(Exp(proxy.get(), shared.get(), count, bit_precision_));
      std::unique_ptr<double[]> result(ReconstructDouble(shared_result.get(), count));
      if (proxy->GetPRole() == proxy2) {
        double correct_result;
        double difference;
        for (int i = 0; i < count; i++) {
          correct_result = exp(values[i]);
          difference = abs(correct_result - result[i]);
          if ((difference * 100.0 / abs(correct_result)) >= 1) {
            // debug info can be printed here
            is_broken = true;
          }
        }
      }
    } else {
      Exp(proxy.get(), nullptr, count, bit_precision_);
    }
    return is_broken;
  }

  // TODO fix and add to the methods vector
  bool TestMatrixMatrixMultiply() const {
    int matrix_count = 3;
    int row_count_a = 50;
    int column_count_a = 60;
    int column_count_b = 30;
    bool is_broken = false;
    if (proxy->GetPRole() != helper) {
      double min_val = -100;
      double max_val = 100;
      int precision_result;
      if (bit_precision_ == 0) {
        precision_result = 0;
      } else {
        precision_result = bit_precision_/2-3;
      }
      double ***matrices1 = new double **[matrix_count];
      double ***matrices2 = new double **[matrix_count];
      uint64_t ***shared1 = new uint64_t **[matrix_count];
      uint64_t ***shared2 = new uint64_t **[matrix_count];
      for (int i = 0; i < matrix_count; i++) {
        matrices1[i] = Random2dData(proxy.get(), row_count_a, column_count_a, min_val, max_val);
        shared1[i] = proxy->CreateShare(matrices1[i], row_count_a, column_count_a, bit_precision_);
        matrices2[i] = Random2dData(proxy.get(), column_count_a, column_count_b, min_val, max_val);
        shared2[i] = proxy->CreateShare(matrices2[i], column_count_a, column_count_b, bit_precision_);
      }
      uint64_t ***shared_result = MatrixMatrixMultiply(
        proxy.get(), shared1, shared2, matrix_count, row_count_a, column_count_a, column_count_b, bit_precision_
      );
      double **result;
      double **correct_result;
      for (int i = 0; i < matrix_count; i++) {
        result = ReconstructDouble(shared_result[i], row_count_a, column_count_b);
        correct_result = MultiplyMatrices(matrices1[i], matrices2[i], row_count_a, column_count_a, column_count_b);
        for (int row = 0; row < row_count_a; row++) {
          for (int column = 0; column < column_count_b; column++) {
            if (!IsWithinMargin(correct_result[row][column], result[row][column], precision_result, 3)) {
              // debug info can be printed here
              is_broken = true;
            }
          }
        }
        Delete2dMatrix(result, row_count_a);
        Delete2dMatrix(correct_result, row_count_a);
      }
      Delete3dMatrix(matrices1, matrix_count, row_count_a);
      Delete3dMatrix(shared1, matrix_count, row_count_a);
      Delete3dMatrix(matrices2, matrix_count, column_count_a);
      Delete3dMatrix(shared2, matrix_count, column_count_a);
      Delete3dMatrix(shared_result, matrix_count, row_count_a);
    } else {
      MatrixMatrixMultiply(
        proxy.get(), nullptr, nullptr, matrix_count, row_count_a, column_count_a, column_count_b, bit_precision_
      );
    }
    return is_broken;
  }

  // TODO fix and add to the methods vector
  bool TestMatrixVectorMultiply() const {
    bool is_broken = false;
    size_t row_count = 30;
    size_t column_count = 15;
    size_t matrix_count = 5;
    if (proxy->GetPRole() != helper) {
      double min_val = -0.784378;
      double max_val = 120.76;
      int precision_result;
      if (bit_precision_ == 0) {
        precision_result = 0;
      } else {
        precision_result = bit_precision_/2-3;
      }
      double ***matrices = new double **[matrix_count];
      double **vectors = Random2dData(proxy.get(), matrix_count, column_count, min_val, max_val);
      uint64_t ***shared_matrices = new uint64_t **[matrix_count];
      uint64_t **shared_vectors = proxy->CreateShare(vectors, matrix_count, column_count);
      for (int i = 0; i < matrix_count; i++) {
        matrices[i] = Random2dData(proxy.get(), row_count, column_count, min_val, max_val);
        shared_matrices[i] = proxy->CreateShare(matrices[i], row_count, column_count);
      }
      uint64_t **shared_result = MatrixVectorMultiply(
        proxy.get(), shared_matrices, shared_vectors, matrix_count, row_count, column_count, bit_precision_
      );
      Delete3dMatrix(shared_matrices, matrix_count, row_count);
      Delete2dMatrix(shared_vectors, matrix_count);
      double **result = ReconstructDouble(shared_result, matrix_count, row_count);
      Delete2dMatrix(shared_result, matrix_count);
      std::unique_ptr<double[]> correct_result;
      for (int i = 0; i < matrix_count; i++) {
        correct_result = std::unique_ptr<double[]>(MultiplyMatrixVector(matrices[i], vectors[i], row_count, column_count));
        for (int row = 0; row < row_count; row++) {
          if (!IsWithinMargin(correct_result[row], result[i][row], precision_result, 3)) {
            // debug info can be printed here
            is_broken = true;
          }
        }
      }
      Delete2dMatrix(result, matrix_count);
      Delete3dMatrix(matrices, matrix_count, row_count);
      Delete2dMatrix(vectors, matrix_count);
    } else { // helper
      MatrixVectorMultiply(proxy.get(), nullptr, nullptr, matrix_count, row_count, column_count, bit_precision_);
    }
    return is_broken;
  }

private:
  static Method GetGenericMethod(std::string method_name, bool(TestCore::*method_reference)() const) {
    return Method(method_name, static_cast<bool(Test::*)() const>(method_reference));
  }
};

#endif // TEST_CORE_H
