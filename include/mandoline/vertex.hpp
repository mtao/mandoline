#pragma once
#include <set>
#include <mtao/types.hpp>
#include <variant>
#include "mandoline/coord_mask.hpp"

#define __STATIC_TEST__


namespace mandoline {
#if defined(_MSC_VER) || defined(__APPLE__) || defined(__STATIC_TEST__)
extern const double threshold_epsilon;
#else
//TODO: sqrt isn't constepxr on mac
constexpr static double threshold_epsilon = std::sqrt(std::numeric_limits<double>::epsilon());
#endif
//Base vertex class, holds points in grid space.
//Vertex is at p() = coord+quot and clamped_indices stores what type of grid entry the vertex belongs to
template<int D>
struct Vertex {


    //Definitions
    using Vec = mtao::Vector<double, D>;
    using coord_type = std::array<int, D>;
    using MaskType = coord_mask<D>;
    using OptInd = std::array<std::optional<int>, D>;

    //Members
    coord_type coord = {};
    Vec quot = Vec::Zero();
    std::bitset<D> clamped_indices = {};


    //Constructors
    Vertex() = default;
    Vertex(const coord_type &);
    template<typename Derived>
    Vertex(const coord_type &c, const Eigen::MatrixBase<Derived> &q);
    template<typename Derived>
    Vertex(const coord_type &c, const Eigen::MatrixBase<Derived> &q, const std::bitset<D> &bs);
    Vertex(const Vertex &) = default;
    Vertex(Vertex &&) = default;
    Vertex &operator=(const Vertex &) = default;
    Vertex &operator=(Vertex &&) = default;

    //Construction helpers
    template<typename Derived>
    static std::bitset<D> bs_from_quot(const Eigen::MatrixBase<Derived> &q);
    template<typename Derived>
    static Vertex from_vertex(const Eigen::MatrixBase<Derived> &p);

    //Comparison
    bool operator==(const Vertex &o) const;
    bool operator!=(const Vertex &o) const;
    bool operator<(const Vertex &o) const;
    bool approx(const Vertex &o) const;

    //Transformations
    Vec p() const;
    OptInd optional_index() const;
    MaskType mask() const;
    operator std::string() const;

    //Bitmask convenience functions
    bool clamped(int index) const;
    size_t clamped_count() const;
    bool is_grid_vertex() const;
    bool is_in_cell(const coord_type &c) const;//closed cell concept
    // takes in a structure for the vertices of the input grid
    template <typename IndexerType>
    bool is_in_grid(const IndexerType& g) const;
    std::set<coord_type> possible_cells() const;
    //std::set<coord_type> possible_faces() const;

    //Thresholding
    void repair();
    void apply_thresholding();
    void apply_thresholding(double thresh);

    //Vector field operations
    Vertex operator+(const Vertex &o) const;
    Vertex operator-(const Vertex &o) const;
    Vertex operator*(double val) const;
    Vertex lerp(const Vertex &other, double t) const;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};
}// namespace mandoline
#include "vertex_impl.hpp"
