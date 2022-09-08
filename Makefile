setup: Dependencies/htslib-1.9/libhts.a
	python setup.py install
Dependencies/htslib-1.9/libhts.a: Dependencies/htslib-1.9
	cd Dependencies/htslib-1.9; ./configure --disable-bz2 --disable-lzma CFLAGS=-fpic && make
