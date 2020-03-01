

mandoline {
    struct TriangulatedPolygon {
        std::optional<ColVecs> triangulated_vertices;
        std::optional<mtao::ColVecs3i> triangulation;

        template<typename Derived>
            mtao::Vec3d brep_centroid(const Eigen::MatrixBase<Derived> &V, bool use_triangulation = false) const;

        template<typename Derived>
            double brep_volume(const Eigen::MatrixBase<Derived> &V, bool use_triangulation = false) const;

        mtao::ColVecs3i triangulate_fan() const;
        mtao::ColVecs3i triangulate_earclipping(const mtao::ColVecs2d &V) const;
        std::tuple<ColVecs, mtao::ColVecs3i> triangulate_triangle(const mtao::ColVecs2d &V, bool add_vertices = false) const;
        mtao::ColVecs3i triangulate(const std::array<mtao::ColVecs2d, 3> &V) const;
        std::tuple<ColVecs, mtao::ColVecs3i> triangulate(const std::array<mtao::ColVecs2d, 3> &V, bool add_vertices) const;
        std::tuple<ColVecs, mtao::ColVecs3i> triangulate(const mtao::ColVecs2d &V, bool add_vertices) const;
        void cache_triangulation(const std::array<mtao::ColVecs2d, 3> &V, bool add_verts = true);
        void cache_triangulation(const mtao::ColVecs3i &F);
        void cache_triangulation(const ColVecs &V, const mtao::ColVecs3i &F);
    };
}
