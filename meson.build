project('balsa', 'cpp',
  version : '0.1',
  default_options : ['warning_level=3', 'cpp_std=c++20'])






balsa_dep = dependency('balsa')




if get_option('protobuf')
  protobuf_dep = dependency('protobuf', required: false)
  if not protobuf_dep.found()
    protobuf_proj = subproject('protobuf')
    protobuf_dep = dependency('protobuf')
  endif
endif


if get_option('imgui')
  imgui_dep = dependency('imgui', required: false)
  if not imgui_dep.found()
    imgui_proj = subproject('imgui')
    imgui_dep = dependency('imgui')
  endif
endif

glm_dep = dependency('glm', required: false)
if not glm_dep.found()
  glm_proj = subproject('glm')
  glm_dep = dependency('glm')
endif


TBB_dep = dependency('TBB')






#===========================================
# Core 
#===========================================

include_dir = [include_directories('include')]


core_deps = [spdlog_dep, eigen_dep,range_v3_dep, fmt_dep]

if get_option('json')
  nlohmann_json_dep = dependency('nlohmann_json', required: false)
  if not nlohmann_json_dep.found()
    nlohmann_json_proj = subproject('nlohmann_json')
    nlohmann_json_dep = dependency('nlohmann_json')
  endif
  core_deps += nlohmann_json_dep
endif






#===========================================
# Geometry
#===========================================



#subdir('tests')
#subdir('tools')
