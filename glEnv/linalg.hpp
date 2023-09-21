#ifndef GL_TEST_LINALG_HPP
#define GL_TEST_LINALG_HPP

#include <iostream>
#include <execution>

#define PI std::acos(-1.0f)

template<typename T>
[[nodiscard]] inline T radians(T&& degrees) {
    static_assert(std::is_arithmetic_v<T>);
    return PI / 180. * std::forward<T>(degrees);
}

template <typename T, std::size_t n>
struct vec final {
    static_assert(std::is_arithmetic_v<T>, "Only arithmetic types are allowed.");
    static_assert(n > 1, "Vector must have at least 2 components.");
private:
    T components[n] {};
public:
    vec() = default;

    template<typename T_, typename... U>
    vec(T_&& a, U&&... b)
    requires (
        std::is_nothrow_convertible_v<std::remove_reference_t<T_>, T> &&
        (std::is_nothrow_convertible_v<std::remove_reference_t<U>, T> && ...)
    ) && (1 + sizeof...(U) == n) : components {std::forward<T>(a), std::forward<T>(b)...} {}

    template<typename T_>
    vec(const T_(&init)[n]) {
        std::move(init, init + n, components);
    }

    template <std::size_t m>
    explicit vec(const vec<T, m>& v) requires (n <= m) {
        const T* vp = reinterpret_cast<const T*>(&v);
        std::copy_n(vp, n, components);
    }

    template <std::size_t m, typename T_, typename... U>
    vec(const vec<T, m>& v, T_&& a, U&&... b)
    requires (n > m) && (
        std::is_nothrow_convertible_v<std::remove_reference_t<T_>, T> &&
        (std::is_nothrow_convertible_v<std::remove_reference_t<U>, T> && ...)
    ) && (1 + sizeof...(U) == n - m) : components {std::forward<T>(a), std::forward<T>(b)...}  {
        std::move(components, components + n - m, components + m);
        std::copy_n(reinterpret_cast<const T*>(&v), m, components);
    }

    [[nodiscard]] constexpr static inline vec<T, n> null() noexcept {
        return {};
    }

    [[nodiscard]] constexpr static inline vec<T, n> identity() {
        vec<T, n> result;
        std::fill_n(reinterpret_cast<T*>(&result), n, T(1));
        return result;
    }

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

    [[nodiscard]] inline T len2() const noexcept {
        return std::transform_reduce(
            std::execution::par_unseq,
            components, components + n, T(0),
            std::plus{}, [](const T& a) { return a * a; }
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

    [[nodiscard]] inline vec<T, 3> cross(const vec<T, 3>& other) const noexcept {
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
            result.components, [&lam](const T& a) { return a * lam; }
        );
        return result;
    };

    template <typename U>
    [[nodiscard]] inline vec<T, n> operator / (U lam) const noexcept requires std::is_arithmetic_v<U> {
        vec<T, n> result;
        std::transform(
            std::execution::par_unseq,
            components, components + n,
            result.components, [&lam](const T& a) { return a / lam; }
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
        std::copy_n(v.components, n, std::ostream_iterator<T>(out, ", "));
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
struct mat final {
    static_assert(std::is_arithmetic_v<T>, "Only arithmetic types are allowed.");
    static_assert(n > 1, "Matrix must have at least 2 rows.");
private:
    vec<T, m> rows[n] {};
public:
    mat() = default;

    template<typename... U>
    mat(const T(&a)[m], const U(&...b)[m])
    requires (std::is_nothrow_convertible_v<T, U> && ...) && (1 + sizeof...(U) == n)
    : rows {vec<T, m>(a), vec<T, m>(b)...} {}

    mat(const vec<T, std::min(n, m)>& diag) {
        for (std::size_t i = 0, end = std::min(n, m); i < end; ++i)
            rows[i][i] = diag[i];
    }

    mat(vec<T, std::min(n, m)>&& diag) {
        for (std::size_t i = 0, end = std::min(n, m); i < end; ++i)
            rows[i][i] = diag[i];
    }

    explicit mat(T fill_diag) {
        for (std::size_t i = 0, end = std::min(n, m); i < end; ++i)
            rows[i][i] = fill_diag;
    }

    template <std::size_t n_, std::size_t m_>
    mat(mat<T, n_, m_>& other) requires (n <= n_) && (m <= m_) {
        for (std::size_t i = 0; i < n; ++i)
            std::copy_n(reinterpret_cast<T*>(&other) + i * m_, m, reinterpret_cast<T*>(rows) + i * m);
    }

    [[nodiscard]] constexpr static inline mat<T, n, m> null() noexcept {
        return {};
    }

    [[nodiscard]] constexpr static inline mat<T, n, m> identity() noexcept {
        return mat<T, n, m>(T(1));
    }

    [[nodiscard]] inline T* begin() const noexcept {
        return rows;
    }

    [[nodiscard]] inline T* end() const noexcept {
        return rows + n;
    }

    inline vec<T, m>& operator [] (int i) const {
        return const_cast<vec<T, m>&>(rows[i]);
    }

    [[nodiscard]] inline vec<T, n> col (int j) const {
        vec<T, n> result;
        for (std::size_t i = 0; i < n; ++i)
            result[i] = rows[i][j];
        return result;
    }

    template<std::size_t l>
    [[nodiscard]] mat<T, n, m> matmul (const mat<T, m, l>& other) const {
        mat<T, n, m> result;
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < l; ++j)
                result.rows[i][j] = rows[i].dot(other.col(j));
        return result;
    };

    [[nodiscard]] inline vec<T, n> matmul(const vec<T, m>& v) const noexcept {
        vec<T, n> result;
        std::transform(
            std::execution::par_unseq, rows, rows + n,
            reinterpret_cast<T*>(&result),
            [&v](const vec<T, m>& row) { return row.dot(v); }
        );
        return result;
    }

    [[nodiscard]] mat<T, m, n> transpose() const {
        mat<T, m, n> result;
        for (std::size_t j = 0; j < m; ++j)
            result.rows[j] = this->col(j);
        return result;
    }

    [[nodiscard]] inline T minor(int i_, int j_) const requires (n == m) {
        mat<T, n - 1, m - 1> submatrix;
        int _i(0);
        for (std::size_t i = 0; i < n; ++i) {
            int _j(0);
            if (i == i_)
                continue;
            for (std::size_t j = 0; j < m; ++j) {
                if (j == j_)
                    continue;
                submatrix[_i][_j] = rows[i][j];
                ++_j;
            }
            ++_i;
        }
        return submatrix.det();
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
        float result(0.0f);
        for (int j = 0; j < m; ++j)
            result += rows[0][j] * (1 - 2 * (j % 2)) * minor(0, j);
        return result;
    }

    [[nodiscard]] mat<T, n, n> inverse() const requires (n == m) {
        mat<T, n, n> adj;
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < m; ++j)
                adj[j][i] = (1 - 2 * ((i + j) % 2)) * minor(i, j);
        return adj / this->det();
    }

    [[nodiscard]] float track() const {
        float result(0.0f);
        for (std::size_t i = 0; i < n; ++i)
            result += rows[i][i];
        return result;
    }

    [[nodiscard]] T* raw() const {
        return const_cast<T*>(&rows[0][0]);
    }

    [[nodiscard]] static mat<T, 4, 4> perspective(float FoV, float ar, float near, float far) noexcept requires (n == 4) && (n == m) {
        return {
            {static_cast<float>(1. / ar / tan(FoV / 2.)), T(0), T(0), T(0)},
            {T(0), static_cast<float>(1. / tan(FoV / 2.)), T(0), T(0)},
            {T(0), T(0), (near + far) / (near - far), 2 * near * far / (near - far)},
            {T(0), T(0), T(-1), T(0)}
        };
    }
    [[nodiscard]] static mat<T, 4, 4> ortho(float w, float h, float near, float far) noexcept requires (n == 4) && (n == m) {
        return {
            {2. / w, T(0), T(0), T(0)},
            {T(0), 2. / h, T(0), T(0)},
            {T(0), T(0), 2. / (near - far), (near + far) / (near - far)},
            {T(0), T(0), T(0), T(1)}
        };
    }

    [[nodiscard]] static inline mat<T, 4, 4> translate(const vec<T, 3>& d) noexcept requires (n == 4) && (n == m) {
        return {
            {T(1), T(0), T(0), d.x()},
            {T(0), T(1), T(0), d.y()},
            {T(0), T(0), T(1), d.z()},
            {T(0), T(0), T(0), T(1)}
        };
    }

    [[nodiscard]] static inline mat<T, 3, 3> scale(const vec<T, 3>& s) noexcept requires (n == 3) && (n == m) {
        return {
            {s.x(), T(0), T(0)},
            {T(0), s.y(), T(0)},
            {T(0), T(0), s.z()},
        };
    }

    template<typename ...U>
    [[nodiscard]] static inline mat<T, 4, 4> scale(const vec<T, 3>& s) noexcept requires (n == 4) && (n == m) {
        return {
            {s.x(), T(0), T(0), T(0)},
            {T(0), s.y(), T(0), T(0)},
            {T(0), T(0), s.z(), T(0)},
            {T(0), T(0), T(0), T(1)},
        };
    }

    [[nodiscard]] static inline mat<T, 2, 2> rotate(T angle) noexcept {
        const T c = T(std::cos(angle));
        const T s = T(std::sin(angle));
        return {
            {c, -s},
            {s, c}
        };
    }

    [[nodiscard]] static inline mat<T, 3, 3> rotate(T angle, const vec<T, 3>& axis) noexcept requires (n == m == 3) {
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

    [[nodiscard]] static inline mat<T, 4, 4> rotate(T angle, const vec3& axis) noexcept requires (n == m == 4) {
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
            result.rows, [&lam](const vec<T, m>& row) { return row * lam; }
        );
        return result;
    };

    template <typename U>
    [[nodiscard]] inline mat<T, n, m> operator / (U lam) const noexcept requires std::is_arithmetic_v<U> {
        mat<T, n, m> result;
        std::transform(
            std::execution::par_unseq,
            rows, rows + n,
            result.rows, [&lam](const vec<T, m>& row) { return row / lam; }
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
        for (const vec<T, m>& row : matrix.rows) {
            out << "    { ";
            std::copy_n(reinterpret_cast<const T*>(&row), m, std::ostream_iterator<T>(out, ", "));
            out << "\b\b },\n";
        }
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
    quat(float w = 0, float x = 0, float y = 0, float z = 0): w(w), v(x, y, z) {}
    quat(float w = 0, const vec3& v = {}): w(w), v(v) {}

    [[nodiscard]] static inline quat null () noexcept {
        return {0.0f, 0.0f, 0.0f, 0.0f};
    }

    [[nodiscard]] static inline quat identity () noexcept {
        return {1.0f, 0.0f, 0.0f, 0.0f};
    }

    [[nodiscard]] float len() const noexcept {
        return std::sqrt(w * w + v.x() * v.x() + v.y() * v.y() + v.z() * v.z());
    }

    [[nodiscard]] quat norm() const noexcept {
        return *this / this->len();
    }

    [[nodiscard]] inline quat operator + (const quat& other) const noexcept {
        return { w + other.w, v + other.v };
    };

    [[nodiscard]] inline quat operator - (const quat& other) const noexcept {
        return { w - other.w, v - other.v };
    };

    template <typename T>
    [[nodiscard]] inline quat operator * (T lam) const noexcept {
        return { w * lam, v * lam };
    };

    [[nodiscard]] inline quat operator * (const quat& other) const noexcept {
        return { w * other.w, v * other.v };
    };

    template <typename T>
    [[nodiscard]] inline quat operator / (T lam) const noexcept {
        return { w / lam, v / lam };
    };

    [[nodiscard]] inline quat operator / (const quat& other) const noexcept {
        return { w / other.w, v / other.v };
    };

    [[nodiscard]] inline static quat quatmul(const quat& a, const quat& b) noexcept {
        return {
            a.w * b.w - a.v.dot(b.v),
            b.v * a.w + a.v * b.w + a.v.cross(b.v)
        };
    }

    [[nodiscard]] inline static quat rotate(float angle, const vec3& axis) noexcept {
        return {static_cast<float>(cos(angle / 2.)), axis.normalize() * std::sin(angle / 2.)};
    }

    [[nodiscard]] inline quat quatmul(const quat& other) const noexcept {
        return {
            w * other.w - v.dot(other.v),
            other.v * w + v * other.w + v.cross(other.v)
        };
    };
    [[nodiscard]] inline quat conj() const noexcept {
        return {w, -v};
    }
    [[nodiscard]] inline mat4 to_mat4() const noexcept {
        return {
            {1 - 2 * (v.y() * v.y() + v.z() * v.z()), 2 * (v.x() * v.y() - v.z() * w), 2 * (v.x() * v.z() + v.y() * w), 0.0f},
            {2 * (v.x() * v.y() + v.z() * w), 1 - 2 * (v.x() * v.x() + v.z() * v.z()), 2 * (v.y() * v.z() - v.x() * w), 0.0f},
            {2 * (v.x() * v.z() - v.y() * w), 2 * (v.y() * v.z() + v.x() * w), 1 - 2 * (v.x() * v.x() + v.y() * v.y()), 0.0f},
            {0.0f, 0.0f, 0.0f, 1.0f}
        };
    }
    [[nodiscard]] inline mat3 to_mat3() const noexcept {
        return {
            {1 - 2 * (v.y() * v.y() + v.z() * v.z()), 2 * (v.x() * v.y() - v.z() * w), 2 * (v.x() * v.z() + v.y() * w)},
            {2 * (v.x() * v.y() + v.z() * w), 1 - 2 * (v.x() * v.x() + v.z() * v.z()), 2 * (v.y() * v.z() - v.x() * w)},
            {2 * (v.x() * v.z() - v.y() * w), 2 * (v.y() * v.z() + v.x() * w), 1 - 2 * (v.x() * v.x() + v.y() * v.y())},
        };
    }
};

#endif
