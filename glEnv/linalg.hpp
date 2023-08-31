#ifndef GL_TEST_LINALG_HPP
#define GL_TEST_LINALG_HPP

#include "../../../../Software/IDE/CLion/bin/mingw/lib/gcc/x86_64-w64-mingw32/13.1.0/include/c++/iostream"
#include "../../../../Software/IDE/CLion/bin/mingw/lib/gcc/x86_64-w64-mingw32/13.1.0/include/c++/cmath"
#include "../../../../Software/IDE/CLion/bin/mingw/lib/gcc/x86_64-w64-mingw32/13.1.0/include/c++/numeric"
#include "../../../../Software/IDE/CLion/bin/mingw/lib/gcc/x86_64-w64-mingw32/13.1.0/include/c++/memory"

#define PI std::acos(-1.0f)

template<typename T>
[[nodiscard]] inline T radians(const T degrees) {
    static_assert(std::is_floating_point_v<T>);
    return PI / 180. * degrees;
}

template <typename T, std::size_t n>
struct vec final {
    static_assert(std::is_arithmetic_v<T>, "Only arithmetic types are allowed.");
    static_assert(n > 1, "Vector must have at least 2 components.");
private:
    T components[n] {};
public:
    vec() = default;

    template<typename... U>
    vec(T a, U... b)
    requires (std::is_nothrow_convertible_v<T, U> && ...) && (1 + sizeof...(U) == n)
    : components {T(a), T(b)...} {}

    template<typename T_>
    vec(const T_(&init)[n]) {
        std::move(init, init + n, components);
    }

    template <std::size_t m>
    explicit vec(const vec<T, m>& v) requires (n <= m) {
        const T* vp = reinterpret_cast<const T*>(&v);
        std::copy_n(vp, n, components);
    }

    template <std::size_t m, typename... U>
    vec(const vec<T, m>& v, T a, U... b)
    requires (n > m) && (std::is_nothrow_convertible_v<T, U> && ...) && (1 + sizeof...(U) == n - m)
    : components {T(a), T(b)...}  {
        std::move(components, components + n - m, components + m);
        const T* vp = reinterpret_cast<const T*>(&v);
        std::copy_n(vp, m, components);
    }

    [[nodiscard]] constexpr static inline vec<T, n> null() {
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

    [[nodiscard]] inline T len() const {
        return std::sqrt(std::accumulate(components, components + n, T(0), [](T sum, T val){
            return sum + val * val;
        }));
    }

    [[nodiscard]] inline vec<T, n> normalize() const {
        return *this / this->len();
    }

    [[nodiscard]] T dot(const vec<T, n>& other) const {
        T result(0);
        for (std::size_t i = 0; i < n; ++i)
            result += components[i] * other.components[i];
        return result;
    };

    [[nodiscard]] inline vec<T, 3> cross(const vec<T, 3>& other) const {
        return {
            components[1] * other.components[2] - components[2] * other.components[1],
            components[2] * other.components[0] - components[0] * other.components[2],
            components[0] * other.components[1] - components[1] * other.components[0]
        };
    };

    vec<T, n> operator + (const vec<T, n>& other) const {
        vec<T, n> result;
        for (std::size_t i = 0; i < n; ++i)
            result.components[i] = components[i] + other.components[i];
        return result;
    }

    vec<T, n> operator - (const vec<T, n>& other) const {
        vec<T, n> result;
        for (std::size_t i = 0; i < n; ++i)
            result.components[i] = components[i] - other.components[i];
        return result;
    }

    vec<T, n> operator * (const vec<T, n>& other) const {
        vec<T, n> result;
        for (std::size_t i = 0; i < n; ++i)
            result.components[i] = components[i] * other.components[i];
        return result;
    }

    vec<T, n> operator / (const vec<T, n>& other) const {
        vec<T, n> result;
        for (std::size_t i = 0; i < n; ++i)
            result.components[i] = components[i] / other.components[i];
        return result;
    }

    vec<T, n> operator - () const {
        vec<T, n> result;
        for (std::size_t i = 0; i < n; ++i)
            result.components[i] = -components[i];
        return result;
    }

    template <typename U>
    vec<T, n> operator * (U lam) const requires std::is_arithmetic_v<U> {
        vec<T, n> result;
        for (std::size_t i = 0; i < n; ++i)
            result.components[i] = components[i] * lam;
        return result;
    };

    template <typename U>
    vec<T, n> operator / (U lam) const requires std::is_arithmetic_v<U> {
        vec<T, n> result;
        for (std::size_t i = 0; i < n; ++i)
            result.components[i] = components[i] / lam;
        return result;
    };

    vec<T, n>& operator += (const vec<T, n>& other) {
        *this = *this + other;
        return *this;
    };

    vec<T, n>& operator -= (const vec<T, n>& other) {
        *this = *this - other;
        return *this;
    };

    vec<T, n>& operator *= (const vec<T, n>& other) {
        *this = *this * other;
        return *this;
    };

    vec<T, n>& operator /= (const vec<T, n>& other) {
        *this = *this / other;
        return *this;
    };

    template <typename U>
    vec<T, n>& operator *= (U lam) const requires std::is_arithmetic_v<U> {
        *this = *this * lam;
        return *this;
    }

    template <typename U>
    vec<T, n>& operator /= (U lam) const requires std::is_arithmetic_v<U> {
        *this = *this / lam;
        return *this;
    }

    friend std::ostream& operator << (std::ostream& out, const vec<T, n>& v) {
        std::string output("vec" + std::to_string(n) + "{ ");
        for (const T& component : v.components)
            output.append(std::to_string(component)).append(", ");
        out << output << "\b\b }";
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

    explicit mat(T fill_diag) {
        for (std::size_t i = 0, end = std::min(n, m); i < end; ++i)
            rows[i][i] = fill_diag;
    }

    template <std::size_t n_, std::size_t m_>
    mat(mat<T, n_, m_>& other) requires (n <= n_) && (m <= m_) {
        for (std::size_t i = 0; i < n; ++i)
            std::copy_n(reinterpret_cast<T*>(&other) + i * m_, m, reinterpret_cast<T*>(rows) + i * m);
    }

    [[nodiscard]] constexpr static inline mat<T, n, m> null() {
        return {};
    }

    [[nodiscard]] constexpr static inline mat<T, n, m> identity() {
        return mat<T, n, m>(T(1));
    }

    [[nodiscard]] inline T* begin() const {
        return rows;
    }

    [[nodiscard]] inline T* end() const {
        return rows + n;
    }

    inline vec<T, m>& operator [] (int i) const {
        return const_cast<vec<T, m>&>(rows[i]);
    }

    inline vec<T, n> col (int j) const {
        vec<T, n> result;
        for (std::size_t i = 0; i < n; ++i)
            result[i] = rows[i][j];
        return result;
    }

    template<std::size_t l>
    [[nodiscard]] mat<T, n, m> matmul (const mat<T, m, l>& other) const {
        mat<T, n, m> result;
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = 0; j < l; ++j) {
                result.rows[i][j] = rows[i].dot(other.col(j));
            }
        }
        return result;
    };

    [[nodiscard]] vec<T, n> matmul(const vec<T, m>& v) const {
        vec<T, n> result;
        for (std::size_t i = 0; i < n; ++i)
            result[i] = rows[i].dot(v);
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

    [[nodiscard]] inline T det() const requires (n == 2) && (n == m) {
        return rows[0][0] * rows[1][1] - rows[0][1] * rows[1][0];
    }

    [[nodiscard]] inline T det() const requires (n == 3) && (n == m) {
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

    [[nodiscard]] static mat<T, 4, 4> perspective(float FoV, float ar, float near, float far) requires (n == 4) && (n == m) {
        return {
            {static_cast<float>(1. / ar / tan(FoV / 2.)), T(0), T(0), T(0)},
            {T(0), static_cast<float>(1. / tan(FoV / 2.)), T(0), T(0)},
            {T(0), T(0), (near + far) / (near - far), 2 * near * far / (near - far)},
            {T(0), T(0), T(-1), T(0)}
        };
    }
    [[nodiscard]] static mat<T, 4, 4> ortho(float w, float h, float near, float far) requires (n == 4) && (n == m) {
        return {
            {2. / w, T(0), T(0), T(0)},
            {T(0), 2. / h, T(0), T(0)},
            {T(0), T(0), 2. / (near - far), (near + far) / (near - far)},
            {T(0), T(0), T(0), T(1)}
        };
    }

    [[nodiscard]] static inline mat<T, 4, 4> translate(const vec<T, 3>& d) requires (n == 4) && (n == m) {
        return {
            {T(1), T(0), T(0), d.x()},
            {T(0), T(1), T(0), d.y()},
            {T(0), T(0), T(1), d.z()},
            {T(0), T(0), T(0), T(1)}
        };
    }

    [[nodiscard]] static inline mat<T, 3, 3> scale(const vec<T, 3>& s) requires (n == 3) && (n == m) {
        return {
            {s.x(), T(0), T(0)},
            {T(0), s.y(), T(0)},
            {T(0), T(0), s.z()},
        };
    }

    template<typename ...U>
    [[nodiscard]] static inline mat<T, 4, 4> scale(const vec<T, 3>& s) requires (n == 4) && (n == m) {
        return {
            {s.x(), T(0), T(0), T(0)},
            {T(0), s.y(), T(0), T(0)},
            {T(0), T(0), s.z(), T(0)},
            {T(0), T(0), T(0), T(1)},
        };
    }

    [[nodiscard]] static inline mat<T, 2, 2> rotate(T angle) {
        const T c = T(std::cos(angle));
        const T s = T(std::sin(angle));
        return {
            {c, -s},
            {s, c}
        };
    }

    [[nodiscard]] static inline mat<T, 3, 3> rotate(T angle, const vec<T, 3>& axis) requires (n == m == 3) {
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

    [[nodiscard]] static inline mat<T, 4, 4> rotate(T angle, const vec3& axis) requires (n == m == 4) {
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

    mat<T, n, m> operator + (const mat<T, n, m>& other) const {
        mat<T, n, m> result;
        for (std::size_t i = 0; i < n; ++i)
            result.rows[i] = rows[i] + other.rows[i];
        return result;
    }

    mat<T, n, m> operator - (const mat<T, n, m>& other) const {
        mat<T, n, m> result;
        for (std::size_t i = 0; i < n; ++i)
            result.rows[i] = rows[i] - other.rows[i];
        return result;
    }

    mat<T, n, m> operator * (const mat<T, n, m>& other) const {
        mat<T, n, m> result;
        for (std::size_t i = 0; i < n; ++i)
            result.rows[i] = rows[i] * other.rows[i];
        return result;
    }

    mat<T, n, m> operator / (const mat<T, n, m>& other) const {
        mat<T, n, m> result;
        for (std::size_t i = 0; i < n; ++i)
            result.rows[i] = rows[i] / other.rows[i];
        return result;
    }

    mat<T, n, m> operator - () const {
        mat<T, n, m> result;
        for (std::size_t i = 0; i < n; ++i)
            result.rows[i] = -rows[i];
        return result;
    }

    template <typename U>
    mat<T, n, m> operator * (U lam) const requires std::is_arithmetic_v<U> {
        mat<T, n, m> result;
        for (std::size_t i = 0; i < n; ++i)
            result.rows[i] = rows[i] * lam;
        return result;
    };

    template <typename U>
    mat<T, n, m> operator / (U lam) const requires std::is_arithmetic_v<U> {
        mat<T, n, m> result;
        for (std::size_t i = 0; i < n; ++i)
            result.rows[i] = rows[i] / lam;
        return result;
    };

    mat<T, n, m>& operator += (const vec<T, n>& other) {
        *this = *this + other;
        return *this;
    };

    mat<T, n, m>& operator -= (const vec<T, n>& other) {
        *this = *this - other;
        return *this;
    };

    mat<T, n, m>& operator *= (const vec<T, n>& other) {
        *this = *this * other;
        return *this;
    };

    mat<T, n, m>& operator /= (const vec<T, n>& other) {
        *this = *this / other;
        return *this;
    };

    template <typename U>
    mat<T, n, m>& operator *= (U lam) const requires std::is_arithmetic_v<U> {
        *this = *this * lam;
        return *this;
    }

    template <typename U>
    mat<T, n, m>& operator /= (U lam) const requires std::is_arithmetic_v<U> {
        *this = *this / lam;
        return *this;
    }

    friend std::ostream& operator << (std::ostream& out, const mat<T, n, m>& matrix) {
        std::string output("mat" + std::to_string(n) + "x" + std::to_string(m) + "{\n");
        for (const vec<T, m>& row : matrix.rows) {
            output.append("    { ");
            for (const T& component : row)
                output.append(std::to_string(component)).append(", ");
            output.append("\b\b },\n");
        }
        out << output << "}";
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

struct quat {
public:
    vec3 v;
    float w;
public:
    quat(float w = 0, float x = 0, float y = 0, float z = 0): w(w), v(x, y, z) {}
    quat(float w = 0, const vec3& v = {}): w(w), v(v) {}

    [[nodiscard]] static inline quat null () {
        return {0.0f, 0.0f, 0.0f, 0.0f};
    }

    [[nodiscard]] static inline quat identity () {
        return {1.0f, 0.0f, 0.0f, 0.0f};
    }

    [[nodiscard]] float len() const {
        return std::sqrt(w * w + v.x() * v.x() + v.y() * v.y() + v.z() * v.z());
    }

    [[nodiscard]] quat norm() const {
        return *this / this->len();
    }

    quat operator + (const quat& other) const {
        return { w + other.w, v + other.v };
    };

    quat operator - (const quat& other) const {
        return { w - other.w, v - other.v };
    };

    template <typename T>
    quat operator * (T lam) const {
        return { w * lam, v * lam };
    };

    quat operator * (const quat& other) const {
        return { w * other.w, v * other.v };
    };

    template <typename T>
    quat operator / (T lam) const {
        return { w / lam, v / lam };
    };

    quat operator / (const quat& other) const {
        return { w / other.w, v / other.v };
    };

    static quat quatmul(const quat& a, const quat& b) {
        const vec3 v = b.v * a.w + a.v * b.w + a.v.cross(b.v);
        return {
            a.w * b.w - a.v.dot(b.v),
            v.x(), v.y(), v.z()
        };
    }

    static quat rotate(float angle, const vec3& axis) {
        const vec3 u = axis.normalize() * std::sin(angle / 2.);
        return {static_cast<float>(cos(angle / 2.)), u.x(), u.y(), u.z()};
    }

    [[nodiscard]] quat quatmul(const quat& other) const {
        const vec3 vector = other.v * w + v * other.w + v.cross(other.v);
        return {
            w * other.w - v.dot(other.v),
            vector.x(), vector.y(), vector.z()
        };
    };
    [[nodiscard]] quat conj() const {
        return {w, -v};
    }
    [[nodiscard]] mat4 to_mat4() const {
        return {
            {1 - 2 * (v.y() * v.y() + v.z() * v.z()), 2 * (v.x() * v.y() - v.z() * w), 2 * (v.x() * v.z() + v.y() * w), 0.0f},
            {2 * (v.x() * v.y() + v.z() * w), 1 - 2 * (v.x() * v.x() + v.z() * v.z()), 2 * (v.y() * v.z() - v.x() * w), 0.0f},
            {2 * (v.x() * v.z() - v.y() * w), 2 * (v.y() * v.z() + v.x() * w), 1 - 2 * (v.x() * v.x() + v.y() * v.y()), 0.0f},
            {0.0f, 0.0f, 0.0f, 1.0f}
        };
    }
    [[nodiscard]] mat3 to_mat3() const {
        return {
            {1 - 2 * (v.y() * v.y() + v.z() * v.z()), 2 * (v.x() * v.y() - v.z() * w), 2 * (v.x() * v.z() + v.y() * w)},
            {2 * (v.x() * v.y() + v.z() * w), 1 - 2 * (v.x() * v.x() + v.z() * v.z()), 2 * (v.y() * v.z() - v.x() * w)},
            {2 * (v.x() * v.z() - v.y() * w), 2 * (v.y() * v.z() + v.x() * w), 1 - 2 * (v.x() * v.x() + v.y() * v.y())},
        };
    }
};

#endif
