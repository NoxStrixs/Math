#pragma once

#include <array>
#include <concepts>
#include <stdexcept>
#include <type_traits>
#include <utility>

template <typename T>
concept Numeric = std::is_arithmetic_v<T>;

template <typename T>
concept Integral = std::is_integral<T>::value;

template <typename T>
concept SignedIntegral = Integral<T> && std::is_signed_v<T>;

template <typename T>
concept UnsignedIntegral = Integral<T> && !std::is_signed_v<T>;

template <typename T>
concept Real = std::is_floating_point_v<T>;

namespace Math {

template <Numeric T, std::size_t ROW, std::size_t COL> struct Matrix {
  static_assert(ROW > 0 && COL > 0, "Matrix row or column cannot be 0.");

  std::array<T, ROW * COL> elems{};

  /*
     Helpers & Accessors
   */
  [[nodiscard]] constexpr std::size_t GetRow() const noexcept { return ROW; }
  [[nodiscard]] constexpr std::size_t GetCol() const noexcept { return COL; }
  [[nodiscard]] constexpr auto GetDimensions() const noexcept {
    return std::pair{ROW, COL};
  }

  // Element access via operator()(row, col)
  [[nodiscard]] constexpr T &operator()(std::size_t row,
                                        std::size_t col) noexcept {
    return elems[col + row * COL];
  }

  [[nodiscard]] constexpr T const &operator()(std::size_t row,
                                              std::size_t col) const noexcept {
    return elems[col + row * COL];
  }

  /*
     Private Internal Helpers
   */
  template <typename F>
  constexpr Matrix &ApplyScalarInPlace_(T const scalar, F &&op) noexcept {
    for (std::size_t idx = 0; idx < elems.size(); idx++) {
      elems[idx] = op(elems[idx], scalar);
    }
    return *this;
  }

  template <typename F>
  constexpr Matrix &ApplyMatrixInPlace_(Matrix const &other, F &&op) noexcept {
    for (std::size_t idx = 0; idx < elems.size(); idx++) {
      elems[idx] = op(elems[idx], other.elems[idx]);
    }
    return *this;
  }

  /*
     Scalar Operations
   */
  [[nodiscard]] constexpr auto operator+(T const scalar) const noexcept {
    Matrix res = *this;
    return res += scalar;
  }
  [[nodiscard]] constexpr auto operator-(T const scalar) const noexcept {
    Matrix res = *this;
    return res -= scalar;
  }
  [[nodiscard]] constexpr auto operator*(T const scalar) const noexcept {
    Matrix res = *this;
    return res *= scalar;
  }
  [[nodiscard]] constexpr auto operator/(T const scalar) const noexcept {
    Matrix res = *this;
    return res /= scalar;
  }
  [[nodiscard]] constexpr auto operator%(T const scalar) const noexcept {
    Matrix res = *this;
    return res %= scalar;
  }

  /*
     Scalar Assignment Operations
   */
  constexpr auto operator+=(T const scalar) noexcept -> Matrix & {
    return ApplyScalarInPlace_(scalar, [](T a, T b) { return a + b; });
  }
  constexpr auto operator-=(T const scalar) noexcept -> Matrix & {
    return ApplyScalarInPlace_(scalar, [](T a, T b) { return a - b; });
  }
  constexpr auto operator*=(T const scalar) noexcept -> Matrix & {
    return ApplyScalarInPlace_(scalar, [](T a, T b) { return a * b; });
  }
  constexpr auto operator/=(T const scalar) noexcept -> Matrix & {
    return ApplyScalarInPlace_(scalar, [](T a, T b) { return a / b; });
  }
  constexpr auto operator%=(T const scalar) noexcept -> Matrix & {
    return ApplyScalarInPlace_(scalar, [](T a, T b) { return a % b; });
  }

  /*
     Matrix Element-Wise Operations
   */
  [[nodiscard]] constexpr auto operator+(Matrix const &other) const noexcept {
    Matrix res = *this;
    return res += other;
  }
  [[nodiscard]] constexpr auto operator-(Matrix const &other) const noexcept {
    Matrix res = *this;
    return res -= other;
  }

  /*
     Matrix Element-Wise Assignment Operations
   */
  constexpr auto operator+=(Matrix const &other) noexcept -> Matrix & {
    return ApplyMatrixInPlace_(other, [](T a, T b) { return a + b; });
  }
  constexpr auto operator-=(Matrix const &other) noexcept -> Matrix & {
    return ApplyMatrixInPlace_(other, [](T a, T b) { return a - b; });
  }
};

} // namespace Math
