cd ../src
make clean
cd ../
aclocal -I ./tools/m4 --install
autoreconf -fi
chmod a+x configure
