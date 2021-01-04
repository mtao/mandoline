#pragma once
#include <mtao/eigen/stl2eigen.hpp>
#include <mtao/eigen/iterable.hpp>
#include <mtao/iterator/zip.hpp>
#include <mtao/iterator/interval.hpp>
#include "mandoline/construction/vertex_types.hpp"
#include <type_traits>
namespace mandoline {
template<int D>
Vertex<D>::Vertex(const coord_type &c) : coord(c), clamped_indices((1 << D) - 1) {}

template<int D>
template<typename Derived>
Vertex<D>::Vertex(const coord_type &c, const Eigen::MatrixBase<Derived> &q) : coord(c), quot(q), clamped_indices(bs_from_quot(q)) {}


template<int D>
template<typename Derived>
Vertex<D>::Vertex(const coord_type &c, const Eigen::MatrixBase<Derived> &q, const std::bitset<D> &clamped_indices) : coord(c), quot(q), clamped_indices(clamped_indices) {}

template<int D>
template<typename Derived>
Vertex<D> Vertex<D>::from_vertex(const Eigen::MatrixBase<Derived> &p) {
    Vertex ret{ {}, p };
    ret.repair();
    ret.clamped_indices = ret.bs_from_quot(ret.quot);
    return ret;
}

template<int D>
template<typename Derived>
std::bitset<D> Vertex<D>::bs_from_quot(const Eigen::MatrixBase<Derived> &q) {
    std::bitset<D> bs;
    for (int i = 0; i < D; ++i) {
        bs[i] = (q(i) == 0);
    }
    return bs;
}


template<int D>
bool Vertex<D>::operator==(const Vertex &o) const {
    return quot == o.quot && coord == o.coord;
}

template<int D>
bool Vertex<D>::operator!=(const Vertex &o) const {
    return quot != o.quot || coord != o.coord;
}

template<int D>
bool Vertex<D>::operator<(const Vertex &o) const {
    //If on the same line we can truely compare them
    if (coord < o.coord) {
        return true;
    } else if (coord == o.coord) {
        for (int i = 0; i < D; ++i) {
            auto &a = quot(i);
            auto &b = o.quot(i);
            if (a < b) {
                return true;
            } else if (a == b) {
                continue;
            } else {
                return false;
            }
        }
    }
    return false;
}

template<int D>
bool Vertex<D>::approx(const Vertex &o) const {
    for (int i = 0; i < D; ++i) {
        auto &q = quot(i);
        auto &c = coord[i];
        auto &oq = o.quot(i);
        auto &oc = o.coord[i];
        if (q < .5) {
            if (oc == c) {
                if (std::abs(q - oq) > threshold_epsilon) {
                    return false;
                }
            } else if (oc + 1 == c) {
                if (std::abs(q - (oq + 1)) > threshold_epsilon) {
                    return false;
                }
            } else {
                return false;
            }
        } else {
            if (oc == c) {
                if (std::abs(q - oq) > threshold_epsilon) {
                    return false;
                }
            } else if (oc - 1 == c) {
                if (std::abs(oq - (q + 1)) > threshold_epsilon) {
                    return false;
                }
            } else {
                return false;
            }
        }
    }
    return true;
}


template<int D>
auto Vertex<D>::p() const -> Vec {
    return Eigen::Map<const mtao::Vector<int, D>>(&*coord.begin()).template cast<double>() + quot;
}

template<int D>
auto Vertex<D>::optional_index() const -> OptInd {
    OptInd ret;
    for (int i = 0; i < D; ++i) {
        if (clamped_indices[i]) {
            ret[i] = coord[i];
        }
    }
    return ret;
}

template<int D>
auto Vertex<D>::mask() const -> MaskType {
    MaskType m;
    for (int i = 0; i < D; ++i) {
        if (quot(i) == 0) {
            m[i] = coord[i];
        }
    }
    return m;
}

template<int D>
Vertex<D>::operator std::string() const {
    std::stringstream ss;
    ss << "(" << mtao::eigen::stl2eigen(coord).transpose() << "{" << quot.transpose() << "}"
       << "[";
    for (int i = 0; i < clamped_indices.size(); ++i) {
        ss << clamped_indices[i];
    }
    ss << "])";
    return ss.str();
}


template <int D>
template <typename IndexerType>
bool Vertex<D>::is_in_grid(const IndexerType& g) const {

    const auto& vertex_shape = g.shape();
    for(int j = 0; j < D; ++j) {
        auto&& c = coord[j];
        auto&& vs = vertex_shape[j];
        if(c >= 0 && c < vs) {
            // its on the interior of this axis, move on
            continue;
        } else if(c == vs && clamped_indices[j]) {
            // it's on the boundary of this axis, move on
            continue;
        } else {
            return false;
        }
    }
    return true;

}

template<int D>
bool Vertex<D>::clamped(int index) const {
    return clamped_indices[index];
}

template<int D>
size_t Vertex<D>::clamped_count() const {
    return clamped_indices.count();
}

template<int D>
bool Vertex<D>::is_grid_vertex() const {
    return clamped_indices.all();
}

template<int D>
auto Vertex<D>::possible_cells() const -> std::set<coord_type> {
    return mask().possible_cells(coord);
    /*
            std::set<coord_type> ret;
            for(int i = 0; i < (2 << D); ++i) {
                std::bitset<D> bs(i);
                if((bs&clamped_indices) == bs) {

                    coord_type cc = coord;
                    for(int j = 0; j < D; ++j) {
                        cc[j] -= bs[j]?1:0;
                    }
                    ret.emplace(cc);

                }
            }
            return ret;
            */
}


template<int D>
void Vertex<D>::repair() {
    auto rcoord = mtao::eigen::stl2eigen(coord);
    Vec qfloor = quot.array().floor();
    rcoord += qfloor.template cast<int>();
    quot = quot - qfloor;
}

template<int D>
void Vertex<D>::apply_thresholding() {
    apply_thresholding(threshold_epsilon);
}

template<int D>
void Vertex<D>::apply_thresholding(double thresh) {
    for (int i = 0; i < D; ++i) {
        auto &q = quot(i);
        auto &c = coord[i];
        if (q < .5) {
            if (q < thresh) {
                clamped_indices[i] = true;
                q = 0;
            }
        } else {
            if (q - 1 > -thresh * (std::abs(q + 1))) {
                clamped_indices[i] = true;
                q = 0;
                c++;
            }
        }
    }
}


template<int D>
auto Vertex<D>::operator+(const Vertex &o) const -> Vertex<D> {
    using namespace mtao::eigen;
    Vertex ret;
    auto mcm = stl2eigen(coord);
    auto ocm = stl2eigen(o.coord);
    auto rcm = stl2eigen(ret.coord);
    rcm = mcm + ocm;
    ret.clamped_indices = clamped_indices & o.clamped_indices;
    //ret.clamped_indices = 0;
    ret.quot = o.quot + quot;
    ret.repair();
    return ret;
}

template<int D>
auto Vertex<D>::operator-(const Vertex &o) const -> Vertex<D> {
    using namespace mtao::eigen;
    Vertex ret;
    auto mcm = stl2eigen(coord);
    auto ocm = stl2eigen(o.coord);
    auto rcm = stl2eigen(ret.coord);
    rcm = mcm - ocm;
    ret.clamped_indices = clamped_indices & o.clamped_indices;
    //ret.clamped_indices = 0;
    ret.quot = quot - o.quot;
    ret.repair();
    return ret;
}

template<int D>
auto Vertex<D>::operator*(double val) const -> Vertex<D> {
    using namespace mtao::eigen;
    auto mcm = stl2eigen(coord);
    Vertex ret = from_vertex(val * mcm.template cast<double>());
    ret.quot += val * quot;
    ret.repair();
    ret.clamped_indices = 0;
    return ret;
}

template<int D>
auto Vertex<D>::lerp(const Vertex &other, double t) const -> Vertex<D> {

    Vertex ret = (*this) * (1 - t) + other * t;
    for (int i = 0; i < D; ++i) {
        if (coord[i] == other.coord[i] && clamped_indices[i] && other.clamped_indices[i]) {
            ret.clamped_indices[i] = true;
        }
    }

    return ret;
}
}// namespace mandoline
