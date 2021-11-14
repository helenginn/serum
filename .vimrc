set path=./,src/**,~/vagabond/c4xsrc/**,~/vagabond/subprojects/**,

command! Tags !ctags -R src/**
command! Ninja :wa|!ninja -C build/current
command! Dinja :wa|!ninja -C build/debug

command! Doxy !doxygen Doxyfile



