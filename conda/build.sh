mkdir -p build
cd build

export LIBRARY_PATH=$PREFIX/lib

echo "USER: $(id)"

cmake -G "Ninja" \
      -D CMAKE_INSTALL_PREFIX=${PREFIX} \
      -D CMAKE_BUILD_TYPE=Release \
      -D BUILD_TESTS=ON \
      -D BUILD_VTK_WRITER=ON \
      ..

ninja install