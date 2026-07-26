#include "Matrix.h"
#include <iostream>
#include <ranges>

#include <cstdlib>
#include <ctime>

void Print(auto const &matrix) {
  std::cout << "Matrix dimensions: row=" << matrix.GetRow()
            << ", col=" << matrix.GetCol() << '\n';
  auto rows = matrix.elems | std::views::chunk(matrix.GetCol());

  for (auto const &view : rows) {
    for (auto const val : view) {
      std::cout << val << " ";
    }
    std::cout << '\n';
  }
  std::cout << '\n';
}

int main() {
  std::size_t const row = 4, col = 4;
  Math::Matrix<float, row, col> m;
  Math::Matrix<int, row, col> n;
  for (std::size_t idx = 0; idx < row * col; idx++) {
    m.elems[idx] = 1; // static_cast<float>(idx + 1);
    n.elems[idx] = 1; // static_cast<int>(idx + 1);
  }

  // std::cout << "Det (m): " << Math::Determinant(m) << '\n';
  // std::cout << "Det (n): " << Math::Determinant(n) << '\n';
  std::cout << "Symmetric (m): " << Math::IsSymmetric(m) << '\n';
  std::cout << "Symmetric (n): " << Math::IsSymmetric(n) << '\n';
}
