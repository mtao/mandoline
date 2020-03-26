#include "mandoline/line.hpp"


class PLCurve2 {
    public:
        using LineType = mandoline::Line<2>;
        using Vec = mtao::Vector<double,2>;
        using ColVecs = mtao::ColVectors<double,2>;

        PLCurve2() = default;
        PLCurve2(const PLCurve2&) = default;
        PLCurve2(PLCurve2&&) = default;
        PLCurve2& operator=(const PLCurve2&) = default;
        PLCurve2& operator=(PLCurve2&&) = default;

        PLCurve2(const std::string& filename);

        mtao::vector<LineType> lines() const;
        void add_point(const Vec& p);
        mtao::ColVectors<unsigned int,2> edges() const;
        ColVecs points() const;
        mtao::vector<Vec> stl_points() const;
        void clear() { curves.clear(); }
        void pop_back() { curves.pop_back(); }
        mtao::vector<Vec>& current_curve() { return std::get<0>(curves.back()); }
        const mtao::vector<Vec>& current_curve() const { return std::get<0>(curves.back()); }
        void create_curve() { curves.emplace_back(); }
        void clear_curve() { if(!curves.empty()) current_curve().clear(); }
        void pop_back_point() { if(!curves.empty()) {auto&& c = current_curve(); if(!c.empty()) c.pop_back();} }
        void toggle_closed() { if(!curves.empty()) { bool& v = std::get<1>(curves.back()); v = !v;} }
        bool is_closed() const { if(!curves.empty()) { return std::get<1>(curves.back());} else { return false; } }
        size_t size() const { return curves.size(); }
        size_t current_curve_size() const { if(!curves.empty()) { return current_curve().size(); } else { return 0; } }
        bool empty() const { return curves.empty(); }

        void load(const std::string& filename);
        void save(const std::string& filename) const;

    private:
        std::vector<std::tuple<mtao::vector<Vec>,bool>> curves;
};
