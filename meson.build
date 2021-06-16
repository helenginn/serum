project('something', 'cpp', 'c')
qt5 = import('qt5')

qt5_dep = dependency('qt5', modules: ['Core', 'Gui', 'Widgets'], required : true)
dep_gl = dependency('gl', required : true)
png_dep = dependency('libpng')
c4x_dep = dependency('libcluster4x')
helen3d = dependency('helen3d')
helencore = dependency('helencore')

cc = meson.get_compiler('c')
m_dep = cc.find_library('m', required : false)

# Run the version and hash generator
hash_script = find_program('get_hash.sh')

gen_src = custom_target('commit-output',
                        input : ['src/commit.h.pre'],
                        output : ['commit.h'],
                        command : [hash_script, '@INPUT@', '@OUTPUT@'],
                        build_by_default: true,
                        build_always_stale: true)

# Graphical user interface
# (including all the Qt and OpenGL stuff)

# Pre-processing by Qt5

moc_files = qt5.preprocess(moc_headers : [
'src/Loader.h', 
'src/SerumView.h', 
],
		moc_extra_arguments: ['-DMAKES_MY_MOC_HEADER_COMPILE'])

executable('something', gen_src, moc_files,
'src/Loader.cpp', 
'src/MyDictator.cpp', 
'src/SerumView.cpp', 
'src/Strain.cpp', 
'src/Serum.cpp', 
'src/main.cpp', 
cpp_args: ['-std=c++11'], dependencies : [ 
qt5_dep, 
c4x_dep, 
dep_gl, 
png_dep,
helen3d,
helencore
], install: true)
