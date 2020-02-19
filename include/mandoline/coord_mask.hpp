#pragma once
#include <bitset>
#include <array>
#include <sstream>
#include <algorithm>
#include <mtao/iterator/zip.hpp>
#include <mtao/iterator/enumerate.hpp>
#include <mtao/eigen/iterable.hpp>
#include <optional>


namespace mandoline {
    template <int D>
       struct Vertex;
    template <int D, typename T=int>
        struct coord_mask: public std::array<std::optional<T>,D> {
            coord_mask() = default;
            coord_mask(const coord_mask&) = default;
            coord_mask(coord_mask&&) = default;
            coord_mask& operator=(const coord_mask&) = default;
            coord_mask& operator=(coord_mask&&) = default;

            coord_mask(int idx, int coord) { (*this)[idx] = coord; }
            coord_mask(const std::array<T,D>& o) { reset(o); }


            void reset() { *this = {}; }
            void reset(int idx, int coord) {
                reset();
                (*this)[idx] = coord;
            }
            void reset(const std::array<T,D>& o) {
                std::copy(o.begin(),o.end(),this->begin());
            }

            int unbound_axis() const {
                assert(count() == D-1);
                for(int i = 0; i < D; ++i) {
                    if(!(*this)[i]) {
                        return i;
                    }
                }
                return -1;
            }


            int bound_axis() const {
                assert(count() == 1);
                for(int i = 0; i < D; ++i) {
                    if((*this)[i]) {
                        return i;
                    }
                }
                return -1;
            }
            std::array<int,D> as_array() const {
                std::array<int,D> m;
                std::transform(this->begin(),this->end(),m.begin(),[](auto&& a) {
                        return *a;
                        });
                return m;
            }
            size_t count() const {
                return std::count_if(this->begin(),this->end(),[](const std::optional<T>& v) -> bool {return bool(v);});
            }

            bool all() const {
                return std::all_of(this->begin(),this->end(),[](const std::optional<T>& v) -> bool {return bool(v);});
            }
            bool active() const {
                return std::any_of(this->begin(),this->end(),[](const std::optional<T>& v) -> bool {return bool(v);});
            }
            coord_mask& operator&=(const coord_mask& o) {
                for(auto&& [a,b]: mtao::iterator::zip(*this,o)) {
                    if(a && b && (*a == *b)) {
                    } else {
                        a = {};
                    }
                }
                return *this;
            }

            coord_mask operator&(const coord_mask& o) const {
                coord_mask ret = *this;
                ret &= o;
                return ret;
            }
            coord_mask& operator|=(const coord_mask& o) {
                for(auto&& [a,b]: mtao::iterator::zip(*this,o)) {
                    if(a && b && (*a == *b)) {
                    } else {
                        a = {};
                    }
                }
                return *this;
            }

            coord_mask operator|(const coord_mask& o) const {
                coord_mask ret = *this;
                ret |= o;
                return ret;
            }
            coord_mask operator-=(const coord_mask& o){
                for(auto&& [a,b]: mtao::iterator::zip(*this,o)) {
                    if(b) {
                        a = {};
                    }
                }
                return *this;
            }
            coord_mask operator-(const coord_mask& o) const {
                coord_mask ret = *this;
                ret -= o;
                return ret;
            }

            enum class PartialOrdering { Less, Greater, Equal, Unknown };
            //integer lattice partial ordering
            PartialOrdering partial_ordering(const coord_mask& other) const {
                //*this > other
                PartialOrdering po = PartialOrdering::Equal;
                for(auto&& [a,b]: mtao::iterator::zip(*this,other)) {
                    if(a && b) {
                        if(*a != *b) {
                            return PartialOrdering::Unknown;
                        }
                    } else if(a) {//but ont b  i.e (a > b)
                        switch(po) {
                            case PartialOrdering::Greater: continue;
                            case PartialOrdering::Less: return PartialOrdering::Unknown;
                            case PartialOrdering::Equal: po = PartialOrdering::Greater; continue;
                            case PartialOrdering::Unknown: return PartialOrdering::Unknown;
                        }
                    } else if(b) {//but ont b  i.e (a < b)
                        switch(po) {
                            case PartialOrdering::Less: continue;
                            case PartialOrdering::Greater: return PartialOrdering::Unknown;
                            case PartialOrdering::Equal:po = PartialOrdering::Less; continue;
                            case PartialOrdering::Unknown: return PartialOrdering::Unknown;
                        }
                    }
                }
                return po;
            }
            bool subsumes(const coord_mask& other) const {
                PartialOrdering po = partial_ordering(other);
                return po == PartialOrdering::Greater || po == PartialOrdering::Equal;
            }
            bool strict_subsumes(const coord_mask& other) const {
                PartialOrdering po = partial_ordering(other);
                return po == PartialOrdering::Greater;
            }

            operator std::string() const {
                std::stringstream ss;
                ss << "(";
                for(int i = 0; i < D-1; ++i) {
                    if((*this)[i]) {
                        ss << *(*this)[i] <<",";
                    } else {
                        ss << "_.";
                    }
                }
                if((*this)[D-1]) {
                    ss << *(*this)[D-1];
                } else {
                    ss << "_";
                }
                ss << ")";
                return ss.str();

            }
            void clamp(std::array<T,D>& vec) const {
                for(auto&& [v,t]: mtao::iterator::zip(vec,*this)) {
                    if(t) {
                        v = *t;
                    }
                }
            }
            void clamp(Vertex<D>& vec) const {
                for(auto&& [v,q,t]:
                        mtao::iterator::zip(vec.coord,mtao::eigen::iterable(vec.quot),*this)) {
                    if(t) {
                        v = *t;
                        q = 0;
                    }
                }
                vec.clamped_indices |= as_bitset();
            }
            std::bitset<D> as_bitset() const {
                std::bitset<D> bs;
                for(auto&& [i,t]: mtao::iterator::enumerate(*this)) {
                    bs[i] = bool(t);
                }
                return bs;
            }
        };
}
