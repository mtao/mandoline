
catch2_dep = dependency('catch2', required: false)
if not catch2_dep.found()
  catch2_proj = subproject('catch2')
  catch2_dep = dependency('catch2')
endif

mandoline_test_sources = [
  'face_collapser_test.cpp',
  'exterior_domain_test.cpp',
  '2d_boundary_test.cpp',
  'diffgeo_tests.cpp',
  'interp_tests.cpp',
  'simple_cutmesh3_tests.cpp',
  'cell_collapser_test.cpp',
  'cell_identification_test.cpp',
  'facet_test.cpp'
]
test_mandoline = executable('test_mandoline', 'main.cpp', mandoline_test_sources , dependencies: [mandoline_dep,catch2_dep])
test('mandoline test', test_mandoline)

