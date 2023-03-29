project('nbody', 'c', 'cpp',
        version: '0.1.0',
        default_options: [
          'buildtype=release',
          'warning_level=3',
          'cpp_std=c++11',
          'default_library=static',
        ],
        meson_version: '>=0.54.0',
)

cc = meson.get_compiler('c')
cxx = meson.get_compiler('cpp')

mpi_cxx = find_program('mpicxx')
mpirun_args = ['-n', '4']
mpi_dep = dependency('mpi', language: 'cpp')

# mpi_lib = find_library('mpi', required: true)

subdir('src')
