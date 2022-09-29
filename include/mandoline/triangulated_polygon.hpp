

mandoline {
    struct TriangulatedPolygon {
        std::optional<ColVecs> triangulated_vertices;
        std::optional<balsa::eigen::ColVecs3i> triangulation;

        template<typename Derived>
            balsa::eigen::Vec3d brep_centroid(const Eigen::MatrixBase<Derived> &V, bool use_triangulation = false) const;

        template<typename Derived>
            double brep_volume(const Eigen::MatrixBase<Derived> &V, bool use_triangulation = false) const;

        balsa::eigen::ColVecs3i triangulate_fan() const;
        balsa::eigen::ColVecs3i triangulate_earclipping(const balsa::eigen::ColVecs2d &V) const;
        std::tuple<ColVecs, balsa::eigen::ColVecs3i> triangulate_triangle(const balsa::eigen::ColVecs2d &V, bool add_vertices = false) const;
        balsa::eigen::ColVecs3i triangulate(const std::array<balsa::eigen::ColVecs2d, 3> &V) const;
        std::tuple<ColVecs, balsa::eigen::ColVecs3i> triangulate(const std::array<balsa::eigen::ColVecs2d, 3> &V, bool add_vertices) const;
        std::tuple<ColVecs, balsa::eigen::ColVecs3i> triangulate(const balsa::eigen::ColVecs2d &V, bool add_vertices) const;
        void cache_triangulation(const std::array<balsa::eigen::ColVecs2d, 3> &V, bool add_verts = true);
        void cache_triangulation(const balsa::eigen::ColVecs3i &F);
        void cache_triangulation(const ColVecs &V, const balsa::eigen::ColVecs3i &F);
    };
}
