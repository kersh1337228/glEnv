#ifndef GL_ENV_UTILS_HPP
#define GL_ENV_UTILS_HPP

template <typename... T>
constexpr void do_nothing(T&&...) {}

template<typename To, typename... From>
concept is_all_convertible = requires(To to, From... from) {
    (std::is_convertible_v<From, To> && ...);
};

template<typename T, typename Category = std::contiguous_iterator_tag>
class enum_iterator final {
public:
    using iterator_category = Category;
    using difference_type = ptrdiff_t;
    using value_type = std::tuple<T, difference_type>;
    using pointer = std::tuple<T*, difference_type>;
    using reference = std::tuple<T&, difference_type>;
private:
    T* ptr;
    difference_type i;
public:
    explicit enum_iterator(T* ptr = nullptr, difference_type i = 0): ptr(ptr), i(i) {}
    enum_iterator(const enum_iterator& other): ptr(other.ptr), i(other.i) {}
public:
    inline reference operator*() const { return { *ptr, i }; }
    inline pointer operator->() const { return { ptr, i }; }
public:
    inline enum_iterator& operator++() { ++ptr; ++i; return *this; }
    inline enum_iterator& operator--() { --ptr; --i; return *this; }
public:
    inline enum_iterator operator++(int) { auto old = *this; ++ptr; ++i; return old; }
    inline enum_iterator operator--(int) { auto old = *this; --ptr; --i; return old; }
public:
    inline enum_iterator& operator+=(difference_type d) { ptr += d; i += d; return *this; }
    inline enum_iterator& operator-=(difference_type d) { ptr -= d; i -= d; return *this; }
    inline enum_iterator operator+(difference_type d) const { return enum_iterator(ptr + d); }
    inline enum_iterator operator-(difference_type d) const { return enum_iterator(ptr - d); }
    friend inline enum_iterator operator+(difference_type d, const enum_iterator& it) { return enum_iterator(d + it.ptr); }
    friend inline enum_iterator operator-(difference_type d, const enum_iterator& it) { return enum_iterator(d - it.ptr); }
public:
    inline reference operator[](difference_type d) const { return { ptr[d], i + d }; }
    inline difference_type operator-(const enum_iterator& other) const { return i - other.i; }
public:
    inline auto operator<=>(const enum_iterator& other) const { return ptr <=> other.ptr; }
    inline bool operator==(const enum_iterator& other) const { return ptr == other.ptr; }
    inline bool operator!=(const enum_iterator& other) const { return ptr != other.ptr; }
};

template<
    typename Category = std::contiguous_iterator_tag,
    typename... T
> class parallel_iterator final {
public:
    using iterator_category = Category;
    using value_type = std::tuple<T...>;
    using difference_type = ptrdiff_t;
    using pointer = std::tuple<T*...>;
    using reference = std::tuple<T&...>;
private:
    pointer ptrs;
public:
    explicit parallel_iterator(T*... ptrs): ptrs(ptrs...) {}
    parallel_iterator(const parallel_iterator& other): ptrs(other.ptrs) {}
public:
    inline reference operator*() const {
        return std::apply([](auto&& ...ptr){ return std::make_tuple(std::ref(*ptr)...); }, ptrs);
    }
    inline pointer operator->() const { return ptrs; }
public:
    inline parallel_iterator& operator++() {
        std::apply([](auto&& ...ptr){ do_nothing(++ptr...); }, ptrs); return *this;
    }
    inline parallel_iterator& operator--() {
        std::apply([](auto&& ...ptr){ do_nothing(--ptr...); }, ptrs); return *this;
    }
public:
    inline parallel_iterator operator++(int) {
        auto old = *this; std::apply([](auto&& ...ptr){ do_nothing(++ptr...); }, ptrs); return old;
    }
    inline parallel_iterator operator--(int) {
        auto old = *this; std::apply([](auto&& ...ptr){ do_nothing(--ptr...); }, ptrs); return old;
    }
public:
    inline parallel_iterator& operator+=(difference_type d) {
        std::apply([&d](auto&& ...ptr){ do_nothing(ptr+=d...); }, ptrs); return *this;
    }
    inline parallel_iterator& operator-=(difference_type d) {
        std::apply([&d](auto&& ...ptr){ do_nothing(ptr-=d...); }, ptrs); return *this;
    }
    inline parallel_iterator operator+(difference_type d) const {
        return parallel_iterator(std::apply([&d](auto&& ...ptr){ return std::make_tuple(ptr + d...); }, ptrs));
    }
    inline parallel_iterator operator-(difference_type d) const {
        return parallel_iterator(std::apply([&d](auto&& ...ptr){ return std::make_tuple(ptr - d...); }, ptrs));
    }
    friend inline parallel_iterator operator+(
        difference_type d, const parallel_iterator& it
    ) { return parallel_iterator(std::apply([&d](auto&& ...ptr){ return std::make_tuple(d + ptr...); }, it.ptrs)); }
    friend inline parallel_iterator operator-(
        difference_type d, const parallel_iterator& it
    ) { return parallel_iterator(std::apply([&d](auto&& ...ptr){ return std::make_tuple(d - ptr...); }, it.ptrs)); }
public:
    inline reference operator[](difference_type d) const {
        return std::apply([&d](auto&& ...ptr){ return std::make_tuple(ptr[d]...); }, ptrs);
    }
    inline difference_type operator-(const parallel_iterator& other) { return std::get<0>(ptrs) - std::get<0>(other); }
public:
    inline friend auto operator<=>(const parallel_iterator& a, const parallel_iterator& b) { return a.ptrs <=> b.ptrs; }
    inline friend bool operator==(const parallel_iterator& a, const parallel_iterator& b) { return a.ptrs == b.ptrs; }
    inline friend bool operator!=(const parallel_iterator& a, const parallel_iterator& b) { return a.ptrs != b.ptrs; }
};

template<
    typename T, std::size_t n,
    typename Category = std::contiguous_iterator_tag
> class n_step_iterator final {
public:
    using iterator_category = Category;
    using value_type = T;
    using difference_type = ptrdiff_t;
    using pointer = T*;
    using reference = T&;
private:
    pointer ptr;
public:
    explicit n_step_iterator(T* ptr = nullptr): ptr(ptr) {}
    n_step_iterator(const n_step_iterator& other): ptr(other.ptr) {}
public:
    inline reference operator*() const { return *ptr; }
    inline pointer operator->() const { return ptr; }
public:
    inline n_step_iterator& operator++() { ptr += n; return *this; }
    inline n_step_iterator& operator--() { ptr -= n; return *this; }
public:
    inline n_step_iterator operator++(int) { auto old = *this; ptr += n; return old; }
    inline n_step_iterator operator--(int) { auto old = *this; ptr -= n; return old; }
public:
    inline n_step_iterator& operator+=(difference_type d) { ptr += n * d; return *this; }
    inline n_step_iterator& operator-=(difference_type d) { ptr -= n * d; return *this; }
    inline n_step_iterator operator+(difference_type d) const { return enum_iterator(ptr + n * d); }
    inline n_step_iterator operator-(difference_type d) const { return enum_iterator(ptr - n * d); }
    friend inline n_step_iterator operator+(difference_type d, const n_step_iterator& it) { return n_step_iterator(n * d + it.ptr); }
    friend inline n_step_iterator operator-(difference_type d, const n_step_iterator& it) { return n_step_iterator(n * d - it.ptr); }
public:
    inline reference operator[](difference_type d) const { return ptr[d]; }
    inline difference_type operator-(const n_step_iterator& other) { return (ptr - other.ptr) / n; }
public:
    inline friend auto operator<=>(const n_step_iterator& a, const n_step_iterator& b) { return a.ptr <=> b.ptr; }
    inline friend bool operator==(const n_step_iterator& a, const n_step_iterator& b) { return a.ptr == b.ptr; }
    inline friend bool operator!=(const n_step_iterator& a, const n_step_iterator& b) { return a.ptr != b.ptr; }
};

template<typename T, typename... Skip>
constexpr T take_first(T val, Skip...){
    return val;
}

template<typename T, typename Arg, size_t... Idxs>
constexpr T subs(Arg arg, std::index_sequence<Idxs...>){
    return T(take_first(arg, Idxs)...);
}

template<std::size_t n, typename T, typename Arg>
constexpr T subs_n_times(Arg arg) {
    return subs<T, Arg>(arg, std::make_index_sequence<n>{});
}

template<typename T, typename U, size_t... Idxs>
constexpr void do_unpack(T from, U to, std::index_sequence<Idxs...>){
    do_nothing(to[Idxs] = from[Idxs]...);
}

template<std::size_t n, typename T, typename U>
constexpr void unpack(T from, U to) {
    do_unpack<T, U>(from, to, std::make_index_sequence<n>{});
}

#endif
