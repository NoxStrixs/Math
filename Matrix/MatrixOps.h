#pragma once

#include "Matrix.h"
#include <limits>

namespace Math {

/*
   Matrix Dot Product
*/
template <Numeric T, std::size_t ROW1, std::size_t COL1, std::size_t COL2>
[[nodiscard]] constexpr auto
operator*(Matrix<T, ROW1, COL1> const &lhs,
          Matrix<T, COL1, COL2> const &rhs) noexcept {
  Matrix<T, ROW1, COL2> result{};
  for (std::size_t i = 0; i < ROW1; i++) {
    for (std::size_t j = 0; j < COL2; j++) {
      T sum = static_cast<T>(0);
      for (std::size_t k = 0; k < COL1; k++) {
        sum += lhs.elems[i * COL1 + k] * rhs.elems[k * COL2 + j];
      }
      result.elems[i * COL2 + j] = sum;
    }
  }
  return result;
}

/*
   Hadamard
*/
namespace Hadamard {
template <Numeric T, std::size_t ROW, std::size_t COL>
[[nodiscard]] constexpr auto Product(Matrix<T, ROW, COL> const &lhs,
                                     Matrix<T, ROW, COL> const &rhs) noexcept {
  Matrix<T, ROW, COL> result{};
  for (std::size_t idx = 0; idx < ROW * COL; idx++) {
    result.elems[idx] = lhs.elems[idx] * rhs.elems[idx];
  }
  return result;
}

template <Integral T, std::size_t ROW, std::size_t COL>
[[nodiscard]] constexpr auto Division(Matrix<T, ROW, COL> const &lhs,
                                      Matrix<T, ROW, COL> const &rhs) {
  Matrix<T, ROW, COL> result{};
  for (std::size_t idx = 0; idx < ROW * COL; idx++) {
      if (rhs.elems[idx] == static_cast<T>(0)
      {
      /// Add policy to set the maximum value or throw????
      }
    result.elems[idx] = lhs.elems[idx] / rhs.elems[idx];
  }
  return result;
}

template <Numeric T, std::size_t ROW, std::size_t COL>
[[nodiscard]] constexpr auto Division(Matrix<T, ROW, COL> const &lhs,
                                      Matrix<T, ROW, COL> const &rhs) noexcept {
  Matrix<T, ROW, COL> result{};
  for (std::size_t idx = 0; idx < ROW * COL; idx++) {
    result.elems[idx] = lhs.elems[idx] / rhs.elems[idx];
  }
  return result;
}

template <Numeric T, std::size_t ROW, std::size_t COL>
[[nodiscard]] constexpr auto
Remainder(Matrix<T, ROW, COL> const &lhs,
          Matrix<T, ROW, COL> const &rhs) noexcept {
  Matrix<T, ROW, COL> result{};
  for (std::size_t idx = 0; idx < ROW * COL; idx++) {
    result.elems[idx] = lhs.elems[idx] % rhs.elems[idx];
  }
  return result;
}
} // namespace Hadamard

namespace Schur {
template <Numeric T, std::size_t ROW, std::size_t COL>
[[nodiscard]] constexpr auto Product(Matrix<T, ROW, COL> const &lhs,
                                     Matrix<T, ROW, COL> const &rhs) noexcept {
  return Hadamard::Product(lhs, rhs);
}

template <Numeric T, std::size_t ROW, std::size_t COL>
[[nodiscard]] constexpr auto Division(Matrix<T, ROW, COL> const &lhs,
                                      Matrix<T, ROW, COL> const &rhs) noexcept {
  return Hadamard::Division(lhs, rhs);
}

template <Numeric T, std::size_t ROW, std::size_t COL>
[[nodiscard]] constexpr auto
Remainder(Matrix<T, ROW, COL> const &lhs,
          Matrix<T, ROW, COL> const &rhs) noexcept {
  return Hadamard::Division(lhs, rhs);
}

} // namespace Schur

template <Numeric T, std::size_t ROW, std::size_t COL>
[[nodiscard]] constexpr bool
IsZeroMatrix(Matrix<T, ROW, COL> const &matrix) noexcept {
  if constexpr (ROW == 1 && COL == 1) {
    return matrix.elems[0] == static_cast<T>(0);
  }
  for (auto const elem : matrix.elems) {
    if (elem != static_cast<T>(0))
      return false;
  }
  return true;
}

template <Numeric T, std::size_t ROW, std::size_t COL>
[[nodiscard]] constexpr bool
IsSymmetric(Matrix<T, ROW, COL> const &matrix) noexcept {
  if constexpr (ROW != COL)
    return false;
  if constexpr (ROW <= 1)
    return true;

  auto const &elems = matrix.elems;
  for (std::size_t row = 0; row < ROW - 1; row++) {
    for (std::size_t col = row + 1; col < ROW; col++) {
      if (elems[COL * row + col] != elems[COL * col + row])
        return false;
    }
  }
  return true;
}

template <Numeric T, std::size_t ROW, std::size_t COL>
[[nodiscard]] constexpr Matrix<T, COL, ROW>
Transpose(Matrix<T, ROW, COL> const &matrix) noexcept {
  Matrix<T, COL, ROW> result;
  if constexpr (ROW == 1 || COL == 1) {
    result.elems = matrix.elems;
    return result;
  }
  for (std::size_t row = 0; row < ROW; row++) {
    for (std::size_t col = 0; col < COL; col++) {
      result.elems[ROW * col + row] = matrix.elems[COL * row + col];
    }
  }
  return result;
}

template <Numeric T, std::size_t ROW, std::size_t COL>
[[nodiscard]] constexpr T
Determinant(Matrix<T, ROW, COL> const &matrix) noexcept(ROW == COL) {
  if constexpr (ROW != COL) {
    throw std::invalid_argument("Must be a square matrix.");
  }

  auto const &elems = matrix.elems;
  if constexpr (ROW == 1)
    return elems[0];
  if constexpr (ROW == 2)
    return (elems[0] * elems[3] - elems[1] * elems[2]);

  if constexpr (ROW == 3) {
    return (elems[0] * (elems[4] * elems[8] - elems[5] * elems[7]) -
            elems[1] * (elems[3] * elems[8] - elems[5] * elems[6]) +
            elems[2] * (elems[3] * elems[7] - elems[4] * elems[6]));
  }

  if constexpr (ROW == 4) {
    auto GetSubMatrix = [&elems](std::size_t excludedCol) {
      Matrix<T, 3, 3> subMatrix;
      std::size_t subIdx = 0;
      for (std::size_t row = 1; row < 4; row++) {
        for (std::size_t col = 0; col < 4; col++) {
          if (col == excludedCol)
            continue;
          subMatrix.elems[subIdx++] = elems[4 * row + col];
        }
      }
      return subMatrix;
    };
    return elems[0] * Determinant(GetSubMatrix(0)) -
           elems[1] * Determinant(GetSubMatrix(1)) +
           elems[2] * Determinant(GetSubMatrix(2)) -
           elems[3] * Determinant(GetSubMatrix(3));
  }

  constexpr T k_epsilonThresh = static_cast<T>(1e-9);
  auto workingElems = elems;
  T det = static_cast<T>(1);
  std::size_t swaps = 0;
  [[maybe_unused]] T previousPivot = static_cast<T>(1);

  for (std::size_t idx = 0; idx < ROW; idx++) {
    std::size_t pivot = idx;
    T max = std::abs(workingElems[COL * idx + idx]);

    for (std::size_t row = idx + 1; row < ROW; row++) {
      T val = std::abs(workingElems[COL * row + idx]);
      if (val > max) {
        max = val;
        pivot = row;
      }
    }

    if constexpr (Real<T>) {
      if (max < k_epsilonThresh)
        return static_cast<T>(0);
    } else {
      if (max == static_cast<T>(0))
        return static_cast<T>(0);
    }

    if (pivot != idx) {
      for (std::size_t col = 0; col < COL; col++) {
        std::swap(workingElems[COL * idx + col],
                  workingElems[COL * pivot + col]);
      }
      swaps++;
    }

    if constexpr (Real<T>) {
      for (std::size_t row = idx + 1; row < ROW; row++) {
        T factor =
            workingElems[COL * row + idx] / workingElems[COL * idx + idx];
        for (std::size_t col = idx; col < COL; col++) {
          workingElems[COL * row + col] -=
              factor * workingElems[COL * idx + col];
        }
      }
    } else {
      for (std::size_t row = idx + 1; row < ROW; row++) {
        for (std::size_t col = idx + 1; col < COL; col++) {
          if (previousPivot == static_cast<T>(0)) {
            workingElems[COL * row + col] = std::numeric_limits<T>::max();
          } else {
            T numerator =
                (workingElems[COL * row + col] *
                 workingElems[COL * idx + idx]) -
                (workingElems[COL * row + idx] * workingElems[COL * idx + col]);
            workingElems[COL * row + col] = numerator / previousPivot;
          }
        }
      }
      previousPivot = workingElems[COL * idx + idx];
    }
  }

  if constexpr (Real<T>) {
    for (std::size_t idx = 0; idx < ROW; idx++)
      det *= workingElems[COL * idx + idx];
    if (std::abs(det) < k_epsilonThresh)
      return static_cast<T>(0);
  } else {
    det = workingElems[(ROW - 1) * COL + (COL - 1)];
  }

  return (swaps & 1) ? -det : det;
}

template <Numeric T, std::size_t ROW, std::size_t COL = ROW>
[[nodiscard]] constexpr auto Identity() noexcept -> Matrix<T, ROW, COL> {
  static_assert(ROW == COL, "Identity matrix must be square.");
  Matrix<T, ROW, COL> result{};
  for (std::size_t row = 0; row < ROW; row++) {
    result.elems[ROW * row + row] = static_cast<T>(1);
  }
  return result;
}

template <Numeric T, std::size_t ROW, std::size_t COL>
[[nodiscard]] constexpr auto Identity(Matrix<T, ROW, COL> const &) noexcept {
  return Identity<T, ROW, COL>();
}

template <Numeric T, std::size_t ROW, std::size_t COL>
[[nodiscard]] constexpr auto Inverse(Matrix<T, ROW, COL> const &matrix) {
  static_assert(ROW == COL, "Matrix must be square to be inverted.");
  if constexpr (!Real<T>) {
    throw std::runtime_error(
        "Matrix inversion is only supported for real/floating-point types.");
  }

  if (Determinant(matrix) == static_cast<T>(0)) {
    throw std::runtime_error("Matrix is singular and cannot be inverted.");
  }

  Matrix<T, ROW, COL> augmented = matrix;
  Matrix<T, ROW, COL> inverse = Identity<T, ROW, COL>();

  for (std::size_t idx = 0; idx < ROW; idx++) {
    std::size_t pivot = idx;
    T max = std::abs(augmented.elems[COL * idx + idx]);
    for (std::size_t row = idx + 1; row < ROW; row++) {
      T val = std::abs(augmented.elems[COL * row + idx]);
      if (val > max) {
        max = val;
        pivot = row;
      }
    }

    if (pivot != idx) {
      for (std::size_t col = 0; col < COL; col++) {
        std::swap(augmented.elems[COL * idx + col],
                  augmented.elems[COL * pivot + col]);
        std::swap(inverse.elems[COL * idx + col],
                  inverse.elems[COL * pivot + col]);
      }
    }

    T pivotVal = augmented.elems[COL * idx + idx];
    for (std::size_t col = 0; col < COL; col++) {
      augmented.elems[COL * idx + col] /= pivotVal;
      inverse.elems[COL * idx + col] /= pivotVal;
    }

    for (std::size_t row = 0; row < ROW; row++) {
      if (row == idx)
        continue;
      T factor = augmented.elems[COL * row + idx];
      for (std::size_t col = 0; col < COL; col++) {
        augmented.elems[COL * row + col] -=
            factor * augmented.elems[COL * idx + col];
        inverse.elems[COL * row + col] -=
            factor * inverse.elems[COL * idx + col];
      }
    }
  }

  return inverse;
}

} // namespace Math
