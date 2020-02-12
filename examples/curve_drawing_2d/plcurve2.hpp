#include "mandoline/line.hpp"


class PLCurve2 {
    public:
        using LineType = mandoline::Line<2>;
        using Vec = mtao::Vector<double,2>;
        using ColVecs = mtao::ColVectors<double,2>;

        mtao::vector<LineType> lines() const;
        void add_point(const Vec& p);
        mtao::ColVectors<unsigned int,2> edges() const;
        ColVecs points() const;
        mtao::vector<Vec> stl_points() const;
        void clear() { pts.clear(); }
        void pop_back() { pts.pop_back(); }
        void toggle_closed() { closed = !closed; }
        size_t size() const { return pts.size(); }
        bool empty() const { return pts.empty(); }

    private:
        mtao::vector<Vec> pts;
        bool closed = false;
};
