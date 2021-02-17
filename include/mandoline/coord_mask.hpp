#pragma once
#include <algorithm>
#include <array>
#include <bitset>
#include <mtao/eigen/iterable.hpp>
#include <mtao/iterator/enumerate.hpp>
#include <mtao/iterator/zip.hpp>
#include <optional>
#include <set>
#include <sstream>

namespace mandoline {
template <int D>
struct Vertex;
template <int D, typename T = int>
struct coord_mask : public std::array<std::optional<T>, D> {
    // coord mask partial ordering has all possible PO cases
    enum class PartialOrdering { Less, Greater, Equal, Unknown };
    using coord_type = std::array<T, D>;

    bool operator==(const coord_mask &other) const = default;

    // rule of 5
    coord_mask() = default;
    coord_mask(const coord_mask &) = default;
    coord_mask(coord_mask &&) = default;
    coord_mask &operator=(const coord_mask &) = default;
    coord_mask &operator=(coord_mask &&) = default;

    // bind only a single coordinate (represents a plane)
    coord_mask(int idx, int coord);

    // bind everything (this is for grid vertices)
    coord_mask(const coord_type &o);

    void reset();
    // set to be a plane
    void reset(int idx, int coord);
    // set to be a grid vertex
    void reset(const coord_type &o);

    // returns true if the axis is bound
    bool is_bound(size_t idx) const;

    // if we're sure we have D-1 elements bound we can pick out the one unbound
    // one this selects the axis that an axis-aligned edge lies on
    int unbound_axis() const;

    // if we're sure we have one bound axis (i.e we have a plane)
    // this selects the axis the plane lies on
    int bound_axis() const;

    // for a grid vertex we may want to pull the vertex coord out
    coord_type as_array() const;

    // the number of bound coordinates (0 = arbitrary position, D = grid vertex)
    size_t count() const;

    // returns if every is grid-aligned; i.e this is a vertex
    bool all() const;

    // this lies on some grid plane
    bool active() const;

    // return the planes shared by two masks
    coord_mask &operator&=(const coord_mask &o);
    coord_mask operator&(const coord_mask &o) const;

    // return the planes that either plane has (not sure why this would
    // geometrically make sense)
    coord_mask &operator|=(const coord_mask &o);
    coord_mask operator|(const coord_mask &o) const;

    // get unshared planes; useful when comparing partial ordered elements when
    // we dont care about some axes
    coord_mask &operator-=(const coord_mask &o);
    coord_mask operator-(const coord_mask &o) const;

    // check the integer lattice partial ordering
    // checks *this < other
    PartialOrdering partial_ordering(const coord_mask &other) const;
    // checks *this <= other
    bool subsumes(const coord_mask &other) const;
    // checks *this < other
    bool strict_subsumes(const coord_mask &other) const;

    // string for visualization
    operator std::string() const;

    // if we know this mask applies to vec, we may want to make sure vec is
    // assigned properly this sets the bound entries in vec to match up with the
    // mask (i.e set q to 0, make sure integer parts are set right) integer part
    // might be wrong if something is set to (N-1,.99999999999) -> (N,0)
    void clamp(Vertex<D> &vec) const;
    // mostly a helper for the prior
    void clamp(coord_type &vec) const;

    // returns which entries are masked
    std::bitset<D> as_bitset() const;

    // Given a fully specified coordinate from a vertex (or
    // coord-masked polygon), returns the grid cells that the object can be
    // considered to be part of
    std::set<coord_type> possible_cells(const coord_type &coord) const;
};
}  // namespace mandoline

#include "mandoline/coord_mask_impl.hpp"
