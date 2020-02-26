#pragma once
#include "mandoline/coord_mask.hpp"

namespace mandoline {

    // bind only a single coordinate (represents a plane)
   template <int D, typename T>
    coord_mask<D,T>::coord_mask(int idx, int coord) {
        reset(idx,coord);
    }

    // bind everything (this is for grid vertices)
   template <int D, typename T>
    coord_mask<D,T>::coord_mask(const coord_type &o) {
        reset(o);
    }


template<int D, typename T>
void coord_mask<D, T>::reset() { *this = {}; }
template<int D, typename T>
void coord_mask<D, T>::reset(int idx, int coord) {
    reset();
    (*this)[idx] = coord;
}
template<int D, typename T>
void coord_mask<D, T>::reset(const coord_type &o) {
    std::copy(o.begin(), o.end(), this->begin());
}

template<int D, typename T>
bool coord_mask<D, T>::is_bound(size_t idx) const {
    assert(idx < D);
    return bool((*this)[idx]);
}

template<int D, typename T>
int coord_mask<D, T>::unbound_axis() const {
    assert(count() == D - 1);
    for (size_t i = 0; i < D; ++i) {
        if (!is_bound(i)) {
            return i;
        }
    }
    return -1;
}


template<int D, typename T>
int coord_mask<D, T>::bound_axis() const {
    assert(count() == 1);
    for (size_t i = 0; i < D; ++i) {
        if (is_bound(i)) {
            return i;
        }
    }
    return -1;
}
template<int D, typename T>
auto coord_mask<D, T>::as_array() const -> coord_type {
    assert(count() == D);
    coord_type m;
    std::transform(this->begin(), this->end(), m.begin(), [](auto &&a) {
        return *a;
    });
    return m;
}
template<int D, typename T>
size_t coord_mask<D, T>::count() const {
    return std::count_if(this->begin(), this->end(), [](const std::optional<T> &v) -> bool { return bool(v); });
}

template<int D, typename T>
bool coord_mask<D, T>::all() const {
    return std::all_of(this->begin(), this->end(), [](const std::optional<T> &v) -> bool { return bool(v); });
}
template<int D, typename T>
bool coord_mask<D, T>::active() const {
    return std::any_of(this->begin(), this->end(), [](const std::optional<T> &v) -> bool { return bool(v); });
}
template<int D, typename T>
auto coord_mask<D, T>::operator&=(const coord_mask &o) -> coord_mask & {
    for (auto &&[a, b] : mtao::iterator::zip(*this, o)) {
        if (a && b && (*a == *b)) {
        } else {
            a = {};
        }
    }
    return *this;
}

template<int D, typename T>
auto coord_mask<D, T>::operator&(const coord_mask &o) const -> coord_mask {
    coord_mask ret = *this;
    ret &= o;
    return ret;
}
template<int D, typename T>
auto coord_mask<D, T>::operator|=(const coord_mask &o) -> coord_mask & {
    for (auto &&[a, b] : mtao::iterator::zip(*this, o)) {
        if (a && b && (*a == *b)) {
        } else {
            a = {};
        }
    }
    return *this;
}

template<int D, typename T>
auto coord_mask<D, T>::operator|(const coord_mask &o) const -> coord_mask {
    coord_mask ret = *this;
    ret |= o;
    return ret;
}
template<int D, typename T>
auto coord_mask<D, T>::operator-=(const coord_mask &o) -> coord_mask & {
    for (auto &&[a, b] : mtao::iterator::zip(*this, o)) {
        if (b) {
            a = {};
        }
    }
    return *this;
}
template<int D, typename T>
auto coord_mask<D, T>::operator-(const coord_mask &o) const -> coord_mask {
    coord_mask ret = *this;
    ret -= o;
    return ret;
}

//integer lattice partial ordering
template<int D, typename T>
auto coord_mask<D, T>::partial_ordering(const coord_mask &other) const -> PartialOrdering {
    //*this > other
    PartialOrdering po = PartialOrdering::Equal;
    for (auto &&[a, b] : mtao::iterator::zip(*this, other)) {
        if (a && b) {
            if (*a != *b) {
                return PartialOrdering::Unknown;
            }
        } else if (a) {//but ont b  i.e (a > b)
            switch (po) {
            case PartialOrdering::Greater:
                continue;
            case PartialOrdering::Less:
                return PartialOrdering::Unknown;
            case PartialOrdering::Equal:
                po = PartialOrdering::Greater;
                continue;
            case PartialOrdering::Unknown:
                return PartialOrdering::Unknown;
            }
        } else if (b) {//but ont b  i.e (a < b)
            switch (po) {
            case PartialOrdering::Less:
                continue;
            case PartialOrdering::Greater:
                return PartialOrdering::Unknown;
            case PartialOrdering::Equal:
                po = PartialOrdering::Less;
                continue;
            case PartialOrdering::Unknown:
                return PartialOrdering::Unknown;
            }
        }
    }
    return po;
}
template<int D, typename T>
bool coord_mask<D, T>::subsumes(const coord_mask &other) const {
    PartialOrdering po = partial_ordering(other);
    return po == PartialOrdering::Greater || po == PartialOrdering::Equal;
}
template<int D, typename T>
bool coord_mask<D, T>::strict_subsumes(const coord_mask &other) const {
    PartialOrdering po = partial_ordering(other);
    return po == PartialOrdering::Greater;
}

template<int D, typename T>
coord_mask<D, T>::operator std::string() const {
    std::stringstream ss;
    ss << "(";
    for (int i = 0; i < D - 1; ++i) {
        if ((*this)[i]) {
            ss << *(*this)[i] << ",";
        } else {
            ss << "_.";
        }
    }
    if ((*this)[D - 1]) {
        ss << *(*this)[D - 1];
    } else {
        ss << "_";
    }
    ss << ")";
    return ss.str();
}
template<int D, typename T>
void coord_mask<D, T>::clamp(coord_type &vec) const {
    for (auto &&[v, t] : mtao::iterator::zip(vec, *this)) {
        if (t) {
            v = *t;
        }
    }
}
template<int D, typename T>
void coord_mask<D, T>::clamp(Vertex<D> &vec) const {
    for (auto &&[v, q, t] :
         mtao::iterator::zip(vec.coord, mtao::eigen::iterable(vec.quot), *this)) {
        if (t) {
            v = *t;
            q = 0;
        }
    }
    vec.clamped_indices |= as_bitset();
}
template<int D, typename T>
std::bitset<D> coord_mask<D, T>::as_bitset() const {
    std::bitset<D> bs;
    for (auto &&[i, t] : mtao::iterator::enumerate(*this)) {
        bs[i] = bool(t);
    }
    return bs;
}

template<int D, typename T>
auto coord_mask<D, T>::possible_cells(const coord_type &coord) const -> std::set<coord_type> {
    std::bitset<D> clamped_indices = as_bitset();
#if defined(_DEBUG)
    for (int i = 0; i < D; ++i) {
        if (clamped_indices[i]) {
            assert((*this)[i] == coord[i]);
        }
    }
#endif
    std::set<coord_type> ret;
    for (int i = 0; i < (2 << D); ++i) {
        std::bitset<D> bs(i);
        if ((bs & clamped_indices) == bs) {

            coord_type cc = coord;
            for (int j = 0; j < D; ++j) {
                cc[j] -= bs[j] ? 1 : 0;
            }
            ret.emplace(cc);
        }
    }
    return ret;
}
}// namespace mandoline
