project('mandoline', 'cpp',
  version : '0.1',
  default_options : ['warning_level=3', 'cpp_std=c++20'])






balsa = subproject('balsa', default_options: ['visualization=false'])



spdlog_dep = balsa.get_variable('spdlog_dep')
fmt_dep = balsa.get_variable('fmt_dep')
eigen_dep = balsa.get_variable('eigen_dep')
range_v3_dep = balsa.get_variable('range_v3_dep')
protobuf_dep = balsa.get_variable('protobuf_dep')
nlohmann_json_dep = balsa.get_variable('nlohmann_json_dep')

balsa_core_dep = balsa.get_variable('core_dep')
balsa_geometry_dep = balsa.get_variable('geometry_dep')
#balsa_geometry_dep = balsa.get_variable('geometry_lib')

tbb_dep = dependency('tbb')

cmake = import('cmake')

mtao_core_opts = cmake.subproject_options()
mtao_core_opts.add_cmake_defines({'MTAO_USE_OPENGL': false})
mtao_core_proj = cmake.subproject('mtao_core', options: mtao_core_opts)


igl_opts = cmake.subproject_options()
igl_opts.add_cmake_defines({
  'LIBIGL_BUILD_TUTORIALS': false,
  'LIBIGL_BUILD_TESTS': false,
  'LIBIGL_NOT_TOPLEVEL_PROJECT': true,
  })
igl_proj = cmake.subproject('libigl', options: igl_opts)
igl_core_dep = igl_proj.get_variable('igl_core_dep')


mtao_headers_dep = mtao_core_proj.get_variable('mtao_headers_dep')
mtao_geometry_dep = mtao_core_proj.get_variable('mtao_geometry_dep')


protoc = find_program('protoc')
protobuf_dep = dependency('protobuf')

gen = generator(protoc, \
  output    : ['@BASENAME@.pb.cc', '@BASENAME@.pb.h'],
  arguments : ['--proto_path=@CURRENT_SOURCE_DIR@/proto', '--cpp_out=@BUILD_DIR@', '@INPUT@'])

proto_srcs = ['proto/cutmesh.proto', 'proto/common.proto']
proto_genned_srcs = []
foreach proto: proto_srcs
  proto_genned_srcs += gen.process(proto)
endforeach

#mandoline_deps = [mtao_headers_dep, mtao_geometry_dep, igl_core_dep, protobuf_dep, tbb_dep)
mandoline_deps = [mtao_headers_dep, mtao_geometry_dep, protobuf_dep, tbb_dep, balsa_core_dep, igl_core_dep]


common_srcs = [
    #src/diffgeo_utils.cpp'
    'src/vertex.cpp'
    ,'src/domain_boundary.cpp'
    ,'src/operators/boundary.cpp'
    ,'src/interpolated_edge.cpp'
    ]


cutmesh2_srcs = [
    'src/mesh2.cpp'
    ,'src/operators/region_boundaries2.cpp'
    ,'src/operators/boundary2.cpp'
    ,'src/operators/interpolation2.cpp'
    ,'src/operators/volume2.cpp'
    ,'src/operators/diffgeo2.cpp'
    ,'src/cutface2.cpp'
    ]
cutmesh3_srcs = [
    'src/mesh3.cpp'
    ,'src/cutcell.cpp'
    ,'src/proto_util.cpp'
    ,'src/barycentric_triangle_face.cpp'
    ,'src/cutface3.cpp'
    ,'src/adaptive_grid.cpp'
    ,'src/operators/boundary3.cpp'
    ,'src/operators/diffgeo3.cpp'
    ,'src/operators/interpolation3.cpp'
    ,'src/operators/masks.cpp'
    ,'src/operators/volume3.cpp'
    ,'src/operators/centroids3.cpp'
    ,'src/operators/nearest_facet3.cpp'
    ,'src/operators/cell_indices3.cpp'
    ,'src/operators/region_boundaries3.cpp'
    ,'src/operators/grid_cell_ownership.cpp'
    ]

construction_tool_srcs = [
    'src/construction/tools/remesh_self_intersections.cpp'
    ,'src/construction/tools/preprocess_mesh.cpp'
    ,'src/construction/tools/read_mesh.cpp'
    ]

include_dir = [include_directories('include')]

mandoline_build_args = [
  '-DMANDOLINE_USE_ADAPTIVE_GRID'
  ]

mandoline = static_library('mandoline', proto_genned_srcs,
  common_srcs,
  cutmesh2_srcs,
  cutmesh3_srcs,
  construction_tool_srcs,
  include_directories: include_dir,
  dependencies: mandoline_deps, pic: true,
  cpp_args: mandoline_build_args
  )


#if get_option('protobuf')
#  protobuf_dep = dependency('protobuf', required: false)
#  if not protobuf_dep.found()
#    protobuf_proj = subproject('protobuf')
#    protobuf_dep = dependency('protobuf')
#  endif
#endif
#
#
#if get_option('imgui')
#  imgui_dep = dependency('imgui', required: false)
#  if not imgui_dep.found()
#    imgui_proj = subproject('imgui')
#    imgui_dep = dependency('imgui')
#  endif
#endif
#
#glm_dep = dependency('glm', required: false)
#if not glm_dep.found()
#  glm_proj = subproject('glm')
#  glm_dep = dependency('glm')
#endif
#
#
#TBB_dep = dependency('TBB')






#===========================================
# Core 
#===========================================

include_dir = [include_directories('include')]


core_deps = [spdlog_dep, eigen_dep,range_v3_dep, fmt_dep]






#===========================================
# Geometry
#===========================================



#subdir('tests')
#subdir('tools')
