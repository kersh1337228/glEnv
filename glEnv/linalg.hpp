#ifndef GL_ENV_LINALG_HPP
#define GL_ENV_LINALG_HPP

#include <iostream>
#include <execution>
#include <cassert>
#include <concepts>
#include "utils.hpp"

#define PI std::acos(-1.0f)

template<typename T>
[[nodiscard]] inline T radians(T&& degrees) requires std::is_arithmetic_v<T> {
    return PI / 180. * std::forward<T>(degrees);
}

template <typename T, std::size_t n>
requires std::is_arithmetic_v<T> && (n > 1)
struct vec final {
private:
    T components[n] {};
public:
    vec() = default;

    template<is_all_convertible<T>... U>
    constexpr vec(U&&... a)
    requires (sizeof...(U) == n)
        : components {std::forward<T>(a)...} {}

    template<is_all_convertible<T> U>
    vec(const U(&init)[n]) {
        unpack<n>(init, components);
    }

    template <std::size_t m>
    explicit vec(const vec<T, m>& v) requires (n <= m) {
        unpack<n>(v, components);
    }
private:
    template <size_t... Idxs, is_all_convertible<T>... U>
    vec(std::index_sequence<Idxs...>, U&&... a): components {Idxs..., std::forward<T>(a)...} {}
public:
    template <std::size_t m, is_all_convertible<T>... U>
    vec(const vec<T, m>& v, U&&... a)
    requires  (n > m) && (sizeof...(U) == n - m)
    : vec(std::make_index_sequence<m>{}, std::forward<T>(a)...) {
        unpack<m>(v, components);
    }

    constexpr vec(T value): vec(subs_n_times<n, vec, T>(value)) {}

    [[nodiscard]] constexpr static inline vec<T, n> null() noexcept { return {}; }

    [[nodiscard]] constexpr static inline vec<T, n> identity() { return 1; }
public:
    [[nodiscard]] inline T* begin() const noexcept {
        return const_cast<T*>(components);
    }

    [[nodiscard]] inline T* end() const noexcept {
        return const_cast<T*>(components + n);
    }

    [[nodiscard]] inline T& x() const noexcept {
        return const_cast<T&>(components[0]);
    }

    [[nodiscard]] inline T& y() const noexcept {
        return const_cast<T&>(components[1]);
    }

    [[nodiscard]] inline T& z() const noexcept requires (n > 2) {
        return const_cast<T&>(components[2]);
    }

    [[nodiscard]] inline T& w() const noexcept requires (n > 3) {
        return const_cast<T&>(components[3]);
    }

    [[nodiscard]] inline T& operator [] (int i) const {
        return const_cast<T&>(components[i]);
    }
public:
    [[nodiscard]] inline T len2() const noexcept {
        return std::transform_reduce(
            std::execution::par_unseq,
            components, components + n, T(0),
            std::plus{}, [](auto&& a) { return a * a; }
        );
    }

    [[nodiscard]] inline T len() const noexcept {
        return std::sqrt(this->len2());
    }

    [[nodiscard]] inline vec<T, n> normalize() const noexcept {
        return *this / this->len();
    }

    [[nodiscard]] inline T dot(const vec<T, n>& other) const noexcept {
        return std::transform_reduce(
            std::execution::par_unseq,
            components, components + n, other.components, T(0),
            std::plus{}, std::multiplies{}
        );
    };

    [[nodiscard]] constexpr inline vec<T, 3> cross(const vec<T, 3>& other) const noexcept {
        return {
            components[1] * other.components[2] - components[2] * other.components[1],
            components[2] * other.components[0] - components[0] * other.components[2],
            components[0] * other.components[1] - components[1] * other.components[0]
        };
    };

    [[nodiscard]] inline vec<T, n> operator + (const vec<T, n>& other) const noexcept {
        vec<T, n> result;
        std::transform(
            std::execution::par_unseq,
            components, components + n, other.components,
            result.components, std::plus{}
        );
        return result;
    }

    [[nodiscard]] inline vec<T, n> operator - (const vec<T, n>& other) const noexcept {
        vec<T, n> result;
        std::transform(
            std::execution::par_unseq,
            components, components + n, other.components,
            result.components, std::minus{}
        );
        return result;
    }

    [[nodiscard]] inline vec<T, n> operator * (const vec<T, n>& other) const noexcept {
        vec<T, n> result;
        std::transform(
            std::execution::par_unseq,
            components, components + n, other.components,
            result.components, std::multiplies{}
        );
        return result;
    }

    [[nodiscard]] inline vec<T, n> operator / (const vec<T, n>& other) const noexcept {
        vec<T, n> result;
        std::transform(
            std::execution::par_unseq,
            components, components + n, other.components,
            result.components, std::divides{}
        );
        return result;
    }

    [[nodiscard]] inline vec<T, n> operator - () const noexcept {
        vec<T, n> result;
        std::transform(
            std::execution::par_unseq,
            components, components + n,
            result.components, std::negate{}
        );
        return result;
    }

    template <typename U>
    [[nodiscard]] inline vec<T, n> operator * (U lam) const noexcept requires std::is_arithmetic_v<U> {
        vec<T, n> result;
        std::transform(
            std::execution::par_unseq,
            components, components + n,
            result.components, [&lam](auto&& a) { return a * lam; }
        );
        return result;
    };

    template <typename U>
    [[nodiscard]] inline vec<T, n> operator / (U lam) const noexcept requires std::is_arithmetic_v<U> {
        vec<T, n> result;
        std::transform(
            std::execution::par_unseq,
            components, components + n,
            result.components, [&lam](auto&& a) { return a / lam; }
        );
        return result;
    };

    inline vec<T, n>& operator += (const vec<T, n>& other) noexcept {
        *this = *this + other;
        return *this;
    };

    inline vec<T, n>& operator -= (const vec<T, n>& other) noexcept {
        *this = *this - other;
        return *this;
    };

    inline vec<T, n>& operator *= (const vec<T, n>& other) noexcept {
        *this = *this * other;
        return *this;
    };

    inline vec<T, n>& operator /= (const vec<T, n>& other) noexcept {
        *this = *this / other;
        return *this;
    };

    template <typename U>
    inline vec<T, n>& operator *= (U lam) noexcept requires std::is_arithmetic_v<U> {
        *this = *this * lam;
        return *this;
    }

    template <typename U>
    inline vec<T, n>& operator /= (U lam) noexcept requires std::is_arithmetic_v<U> {
        *this = *this / lam;
        return *this;
    }

    friend std::ostream& operator << (std::ostream& out, const vec<T, n>& v) {
        out << "vec" << n << "{ ";
        std::copy_n(
            std::execution::par_unseq, v.components,
            n, std::ostream_iterator<T>(out, ", ")
        );
        out << "\b\b }";
        return out;
    }
};

using vec2 = vec<float, 2>;
using vec3 = vec<float, 3>;
using vec4 = vec<float, 4>;

using dvec2 = vec<double, 2>;
using dvec3 = vec<double, 3>;
using dvec4 = vec<double, 4>;

using ivec2 = vec<int, 2>;
using ivec3 = vec<int, 3>;
using ivec4 = vec<int, 4>;

using lvec2 = vec<long, 2>;
using lvec3 = vec<long, 3>;
using lvec4 = vec<long, 4>;

template<typename T, std::size_t n, std::size_t m>
requires std::is_arithmetic_v<T> && (n > 1)
struct mat final {
    using diag_iterator = n_step_iterator<T, m + 1, std::bidirectional_iterator_tag>;
    using col_iterator = n_step_iterator<T, m, std::bidirectional_iterator_tag>;
    using row_t = vec<T, m>;
    using col_t = vec<T, n>;
private:
    row_t rows[n] {};
public:
    constexpr mat() = default;

    template<is_all_convertible<T>... U>
    constexpr mat(const T(&a)[m], const U(&...b)[m])
    requires (1 + sizeof...(U) == n)
    : rows {row_t(a), row_t(b)...} {}

    mat(const vec<T, std::min(n, m)>& diag) {
        std::copy_n(
            std::execution::par_unseq,
            reinterpret_cast<const T*>(&diag), std::min(n, m),
            diag_iterator(reinterpret_cast<T*>(rows))
        );
    }

    mat(vec<T, std::min(n, m)>&& diag) {
        std::copy_n(
            std::execution::par_unseq,
            reinterpret_cast<const T*>(std::move(diag).begin()), std::min(n, m),
            diag_iterator(reinterpret_cast<T*>(rows))
        );
    }

    mat(T fill_diag) {
        std::fill_n(
            std::execution::par_unseq,
            diag_iterator(reinterpret_cast<T*>(rows)),
            std::min(n, m), fill_diag
        );
    }

    template <std::size_t n_, std::size_t m_>
    mat(const mat<T, n_, m_>& other) requires (n <= n_) && (m <= m_) {
        std::for_each_n(
            std::execution::par_unseq,
            parallel_iterator(rows, reinterpret_cast<const std::remove_cvref_t<decltype(other)>::row_t*>(&other)), n,
            [&other](auto&& row) {
                std::copy_n(
                    std::execution::par_unseq,
                    reinterpret_cast<const T*>(&std::get<1>(row)),
                    m, reinterpret_cast<T*>(&std::get<0>(row))
                );
            }
        );
    }

    [[nodiscard]] constexpr static inline mat<T, n, m> null() noexcept { return {}; }

    [[nodiscard]] constexpr static inline mat<T, n, m> identity() noexcept { return 1; }

    [[nodiscard]] inline T* begin() const noexcept {
        return rows;
    }

    [[nodiscard]] inline T* end() const noexcept {
        return rows + n;
    }

    inline row_t& operator [] (int i) const {
        return const_cast<row_t&>(rows[i]);
    }

    [[nodiscard]] inline vec<T, n> col (int j) const {
        vec<T, n> result;
        std::copy_n(
            std::execution::par_unseq,
            col_iterator(const_cast<T*>(reinterpret_cast<const T*>(rows)) + j),
            n, reinterpret_cast<T*>(&result)
        );
        return result;
    }

    template<std::size_t l>
    [[nodiscard]] mat<T, n, l> matmul (const mat<T, m, l>& other) const {
        mat<T, n, l> result;
        auto* cols = other.transpose().rows;
        std::transform(
            std::execution::par_unseq, rows, rows + n,
            result.rows,
            [&cols](auto&& row) {
                vec<T, l> temp;
                std::transform(
                    std::execution::par_unseq, cols, cols + l,
                    reinterpret_cast<T*>(&temp),
                    [&row](auto&& col) { return row.dot(col); }
                );
                return temp;
            }
        );
        return result;
    };

    [[nodiscard]] inline vec<T, n> matmul(const row_t& v) const noexcept {
        vec<T, n> result;
        std::transform(
            std::execution::par_unseq, rows, rows + n,
            reinterpret_cast<T*>(&result),
            [&v](auto&& row) { return row.dot(v); }
        );
        return result;
    }

    [[nodiscard]] inline mat<T, m, n> transpose() const {
        mat<T, m, n> result;
        std::for_each_n(
            std::execution::par_unseq,
            enum_iterator(&result.rows[0]), n, [this](auto&& row) {
            std::copy_n(
                std::execution::par_unseq,
                col_iterator(const_cast<T*>(reinterpret_cast<const T*>(rows)) + std::get<1>(row)),
                n, reinterpret_cast<T*>(&std::get<0>(row))
            );
        });
        return result;
    }

    [[nodiscard]] inline mat<T, n - 1, m - 1> submatrix(int i_, int j_) const {
        mat<T, n - 1, m - 1> result;
        std::for_each_n(
            std::execution::par_unseq,
            enum_iterator(rows), n,
            [i_, j_, &result](auto&& row) {
                const auto index = std::get<1>(row);
                if (index != i_) {
                    auto data = std::get<0>(row);
                    std::copy_if(
                        std::execution::par_unseq,
                        enum_iterator(data.begin()),
                        enum_iterator(data.end()),
                        enum_iterator(reinterpret_cast<T*>(&result) + (index < i_ ? index : index - 1) * (m - 1)),
                        [j_](auto&& val) { return std::get<1>(val) != j_; }
                    );
                }
            }
        );
        return result;
    }

    [[nodiscard]] inline T minor(int i_, int j_) const requires (n == m) {
        return submatrix(i_, j_).det();
    }

    [[nodiscard]] inline T det() const noexcept requires (n == 2) && (n == m) {
        return rows[0][0] * rows[1][1] - rows[0][1] * rows[1][0];
    }

    [[nodiscard]] inline T det() const noexcept requires (n == 3) && (n == m) {
        return (
            rows[0][0] * rows[1][1] * rows[2][2] +
            rows[0][1] * rows[1][2] * rows[2][0] +
            rows[1][0] * rows[0][2] * rows[2][1] -
            rows[0][2] * rows[1][1] * rows[2][0] -
            rows[1][0] * rows[0][1] * rows[2][2] -
            rows[0][0] * rows[2][1] * rows[1][2]
        );
    }

    [[nodiscard]] inline T det() const requires (n == m) && (n > 3) {
        return std::transform_reduce(
            std::execution::par_unseq,
            enum_iterator(rows[0].begin()),
            enum_iterator(rows[0].end()),
            T(0), std::plus{}, [this](auto&& a) {
                const auto j = std::get<1>(a);
                return std::get<0>(a) * (1 - 2 * (j % 2)) * minor(0, j);
            }
        );
    }

    [[nodiscard]] mat<T, n, n> inverse() const requires (n == m) {
        mat<T, n, n> adj;
        std::transform(
            std::execution::par_unseq,
            enum_iterator(adj.rows),
            enum_iterator(adj.rows + n),
            adj.rows, [this](auto&& row) {
                auto data = std::get<0>(row);
                std::transform(
                    std::execution::par_unseq,
                    enum_iterator(data.begin()),
                    enum_iterator(data.end()),
                    data.begin(),
                    [j = std::get<1>(row), this](auto&& val) {
                        const auto i = std::get<1>(val);
                        return (1 - 2 * ((i + j) % 2)) * minor(i, j);
                    }
                );
                return data;
            }
        );
        return adj / this->det();
    }

    [[nodiscard]] T track() const {
        return std::reduce(
            std::execution::par_unseq,
            diag_iterator(const_cast<T*>(reinterpret_cast<const T*>(rows))),
            diag_iterator(const_cast<T*>(reinterpret_cast<const T*>(rows))) + std::min(n, m),
            T(0), std::plus{}
        );
    }

    [[nodiscard]] T* raw() const {
        return const_cast<T*>(&rows[0][0]);
    }

    [[nodiscard]] constexpr static mat<T, 4, 4> perspective(float FoV, float ar, float near, float far) noexcept requires (n == 4) && (n == m) {
        return {
            {static_cast<float>(1. / ar / tan(FoV / 2.)), T(0), T(0), T(0)},
            {T(0), static_cast<float>(1. / tan(FoV / 2.)), T(0), T(0)},
            {T(0), T(0), (near + far) / (near - far), 2 * near * far / (near - far)},
            {T(0), T(0), T(-1), T(0)}
        };
    }
    [[nodiscard]] constexpr static mat<T, 4, 4> ortho(float w, float h, float near, float far) noexcept requires (n == 4) && (n == m) {
        return {
            {2. / w, T(0), T(0), T(0)},
            {T(0), 2. / h, T(0), T(0)},
            {T(0), T(0), 2. / (near - far), (near + far) / (near - far)},
            {T(0), T(0), T(0), T(1)}
        };
    }

    [[nodiscard]] constexpr static inline mat<T, 4, 4> translate(const vec<T, 3>& d) noexcept requires (n == 4) && (n == m) {
        return {
            {T(1), T(0), T(0), d.x()},
            {T(0), T(1), T(0), d.y()},
            {T(0), T(0), T(1), d.z()},
            {T(0), T(0), T(0), T(1)}
        };
    }

    [[nodiscard]] constexpr static inline mat<T, 3, 3> scale(const vec<T, 3>& s) noexcept requires (n == 3) && (n == m) {
        return {
            {s.x(), T(0), T(0)},
            {T(0), s.y(), T(0)},
            {T(0), T(0), s.z()},
        };
    }

    [[nodiscard]] constexpr static inline mat<T, 4, 4> scale(const vec<T, 3>& s) noexcept requires (n == 4) && (n == m) {
        return {
            {s.x(), T(0), T(0), T(0)},
            {T(0), s.y(), T(0), T(0)},
            {T(0), T(0), s.z(), T(0)},
            {T(0), T(0), T(0), T(1)},
        };
    }

    [[nodiscard]] constexpr static inline mat<T, 2, 2> rotate(T angle) noexcept {
        const T c = T(std::cos(angle));
        const T s = T(std::sin(angle));
        return {
            {c, -s},
            {s, c}
        };
    }

    [[nodiscard]] constexpr static inline mat<T, 3, 3> rotate(T angle, const vec<T, 3>& axis) noexcept requires (n == m == 3) {
        const vec3 a = axis.normalize();
        const T c = T(std::cos(angle));
        const T s = T(std::sin(angle));
        return {
            {
                c + (1 - c) * a.x() * a.x(),
                (1 - c) * a.x() * a.y() - s * a.z(),
                (1 - c) * a.x() * a.z() + s * a.y()
            }, {
                (1 - c) * a.x() * a.y() + s * a.z(),
                c + (1 - c) * a.y() * a.y(),
                (1 - c) * a.y() * a.z() - s * a.x()
            }, {
                (1 - c) * a.x() * a.z() - s * a.y(),
                (1 - c) * a.y() * a.z() + s * a.x(),
                c + (1 - c) * a.z() * a.z()
            }
        };
    }

    [[nodiscard]] constexpr static inline mat<T, 4, 4> rotate(T angle, const vec3& axis) noexcept requires (n == m == 4) {
        const vec3 a = axis.normalize();
        const T c = T(std::cos(angle));
        const T s = T(std::sin(angle));
        return {
            {
                c + (1 - c) * a.x() * a.x(),
                (1 - c) * a.x() * a.y() - s * a.z(),
                (1 - c) * a.x() * a.z() + s * a.y(), T(0)
            }, {
                (1 - c) * a.x() * a.y() + s * a.z(),
                c + (1 - c) * a.y() * a.y(),
                (1 - c) * a.y() * a.z() - s * a.x(), T(0)
            }, {
                (1 - c) * a.x() * a.z() - s * a.y(),
                (1 - c) * a.y() * a.z() + s * a.x(),
                c + (1 - c) * a.z() * a.z(), T(0)
            }, {T(0), T(0), T(0), T(1)}
        };
    }

    [[nodiscard]] inline mat<T, n, m> operator + (const mat<T, n, m>& other) const noexcept {
        mat<T, n, m> result;
        std::transform(
            std::execution::par_unseq,
            rows, rows + n, other.rows,
            result.rows, std::plus{}
        );
        return result;
    }

    [[nodiscard]] inline mat<T, n, m> operator - (const mat<T, n, m>& other) const noexcept {
        mat<T, n, m> result;
        std::transform(
            std::execution::par_unseq,
            rows, rows + n, other.rows,
            result.rows, std::minus{}
        );
        return result;
    }

    [[nodiscard]] inline mat<T, n, m> operator * (const mat<T, n, m>& other) const noexcept {
        mat<T, n, m> result;
        std::transform(
            std::execution::par_unseq,
            rows, rows + n, other.rows,
            result.rows, std::multiplies{}
        );
        return result;
    }

    [[nodiscard]] inline mat<T, n, m> operator / (const mat<T, n, m>& other) const noexcept {
        mat<T, n, m> result;
        std::transform(
            std::execution::par_unseq,
            rows, rows + n, other.rows,
            result.rows, std::divides{}
        );
        return result;
    }

    [[nodiscard]] inline mat<T, n, m> operator - () const noexcept {
        mat<T, n, m> result;
        std::transform(
            std::execution::par_unseq,
            rows, rows + n,
            result.rows, std::negate{}
        );
        return result;
    }

    template <typename U>
    [[nodiscard]] inline mat<T, n, m> operator * (U lam) const noexcept requires std::is_arithmetic_v<U> {
        mat<T, n, m> result;
        std::transform(
            std::execution::par_unseq,
            rows, rows + n,
            result.rows, [&lam](auto&& row) { return row * lam; }
        );
        return result;
    };

    template <typename U>
    [[nodiscard]] inline mat<T, n, m> operator / (U lam) const noexcept requires std::is_arithmetic_v<U> {
        mat<T, n, m> result;
        std::transform(
            std::execution::par_unseq,
            rows, rows + n,
            result.rows, [&lam](auto&& row) { return row / lam; }
        );
        return result;
    };

    inline mat<T, n, m>& operator += (const vec<T, n>& other) noexcept {
        *this = *this + other;
        return *this;
    };

    inline mat<T, n, m>& operator -= (const vec<T, n>& other) noexcept {
        *this = *this - other;
        return *this;
    };

    inline mat<T, n, m>& operator *= (const vec<T, n>& other) noexcept {
        *this = *this * other;
        return *this;
    };

    inline mat<T, n, m>& operator /= (const vec<T, n>& other) noexcept {
        *this = *this / other;
        return *this;
    };

    template <typename U>
    inline mat<T, n, m>& operator *= (U lam) noexcept requires std::is_arithmetic_v<U> {
        *this = *this * lam;
        return *this;
    }

    template <typename U>
    inline mat<T, n, m>& operator /= (U lam) noexcept requires std::is_arithmetic_v<U> {
        *this = *this / lam;
        return *this;
    }

    friend std::ostream& operator << (std::ostream& out, const mat<T, n, m>& matrix) {
        out << "mat" << n << "x" << m << "{\n";
        std::for_each_n(
            std::execution::par_unseq,
            matrix.rows, n, [&out](auto&& row) {
                out << "    { ";
                std::copy_n(
                    std::execution::par_unseq,
                    reinterpret_cast<const T*>(&row), m,
                    std::ostream_iterator<T>(out, ", ")
                );
                out << "\b\b },\n";
            }
        );
        out << "\b\b}";
        return out;
    }
};

using mat2 = mat<float, 2, 2>;
using mat3 = mat<float, 3, 3>;
using mat4 = mat<float, 4, 4>;

using dmat2 = mat<double, 2, 2>;
using dmat3 = mat<double, 3, 3>;
using dmat4 = mat<double, 4, 4>;

using imat2 = mat<int, 2, 2>;
using imat3 = mat<int, 3, 3>;
using imat4 = mat<int, 4, 4>;

using lmat2 = mat<long, 2, 2>;
using lmat3 = mat<long, 3, 3>;
using lmat4 = mat<long, 4, 4>;

struct quat final {
public:
    vec3 v;
    float w;
public:
    constexpr quat(float w = 0, float x = 0, float y = 0, float z = 0): w(w), v(x, y, z) {}
    constexpr quat(float w = 0, const vec3& v = {}): w(w), v(v) {}

    [[nodiscard]] consteval static inline quat null () noexcept {
        return {0.0f, 0.0f, 0.0f, 0.0f};
    }

    [[nodiscard]] consteval static inline quat identity () noexcept {
        return {1.0f, 0.0f, 0.0f, 0.0f};
    }

    [[nodiscard]] constexpr float len() const noexcept {
        return std::sqrt(w * w + v.x() * v.x() + v.y() * v.y() + v.z() * v.z());
    }

    [[nodiscard]] constexpr quat norm() const noexcept {
        return *this / this->len();
    }

    [[nodiscard]] constexpr inline quat operator + (const quat& other) const noexcept {
        return { w + other.w, v + other.v };
    };

    [[nodiscard]] constexpr inline quat operator - (const quat& other) const noexcept {
        return { w - other.w, v - other.v };
    };

    template <typename T>
    [[nodiscard]] constexpr inline quat operator * (T lam) const noexcept {
        return { w * lam, v * lam };
    };

    [[nodiscard]] constexpr inline quat operator * (const quat& other) const noexcept {
        return { w * other.w, v * other.v };
    };

    template <typename T>
    [[nodiscard]] constexpr inline quat operator / (T lam) const noexcept {
        return { w / lam, v / lam };
    };

    [[nodiscard]] constexpr inline quat operator / (const quat& other) const noexcept {
        return { w / other.w, v / other.v };
    };

    [[nodiscard]] constexpr inline static quat quatmul(const quat& a, const quat& b) noexcept {
        return {
            a.w * b.w - a.v.dot(b.v),
            b.v * a.w + a.v * b.w + a.v.cross(b.v)
        };
    }

    [[nodiscard]] constexpr inline static quat rotate(float angle, const vec3& axis) noexcept {
        return {static_cast<float>(cos(angle / 2.)), axis.normalize() * std::sin(angle / 2.)};
    }

    [[nodiscard]] constexpr inline quat quatmul(const quat& other) const noexcept {
        return {
            w * other.w - v.dot(other.v),
            other.v * w + v * other.w + v.cross(other.v)
        };
    };
    [[nodiscard]] constexpr inline quat conj() const noexcept {
        return {w, -v};
    }
    [[nodiscard]] constexpr inline mat4 to_mat4() const noexcept {
        return {
            {1 - 2 * (v.y() * v.y() + v.z() * v.z()), 2 * (v.x() * v.y() - v.z() * w), 2 * (v.x() * v.z() + v.y() * w), 0.0f},
            {2 * (v.x() * v.y() + v.z() * w), 1 - 2 * (v.x() * v.x() + v.z() * v.z()), 2 * (v.y() * v.z() - v.x() * w), 0.0f},
            {2 * (v.x() * v.z() - v.y() * w), 2 * (v.y() * v.z() + v.x() * w), 1 - 2 * (v.x() * v.x() + v.y() * v.y()), 0.0f},
            {0.0f, 0.0f, 0.0f, 1.0f}
        };
    }
    [[nodiscard]] constexpr inline mat3 to_mat3() const noexcept {
        return {
            {1 - 2 * (v.y() * v.y() + v.z() * v.z()), 2 * (v.x() * v.y() - v.z() * w), 2 * (v.x() * v.z() + v.y() * w)},
            {2 * (v.x() * v.y() + v.z() * w), 1 - 2 * (v.x() * v.x() + v.z() * v.z()), 2 * (v.y() * v.z() - v.x() * w)},
            {2 * (v.x() * v.z() - v.y() * w), 2 * (v.y() * v.z() + v.x() * w), 1 - 2 * (v.x() * v.x() + v.y() * v.y())},
        };
    }
};

#endif
