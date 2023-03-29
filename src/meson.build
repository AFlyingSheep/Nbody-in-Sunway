sources = files([
    'main.cc'
])

executable('nbody',
    sources,
    dependencies: [mpi_dep],
#   cpp_compiler: mpi_cxx,
    link_args: ['-lmpi']
)
