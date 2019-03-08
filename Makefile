.PHONY:all python clean mostlyclean
CXX?=g++
CC?=gcc
ifndef DBG
DBG=-DNDEBUG
else
DBG=
endif
WARNINGS=-Wall -Wextra -Wno-char-subscripts \
		 -Wpointer-arith -Wwrite-strings -Wdisabled-optimization \
		 -Wformat -Wcast-align -Wno-unused-function -Wno-unused-parameter \
		 -pedantic -Wunused-variable -Wno-attributes -Wno-ignored-attributes
FLAGS=-O3 -funroll-loops -pipe -march=native -msse2 -mavx2 -Ivec/blaze -Ivec -I. -fpic -Wall $(WARNINGS) \
     -fno-strict-aliasing \
    -Wreorder -DXXH_INLINE_ALL  \
	-Wno-attributes -Wno-pragmas -Ibagminhash/c++/ \
    # -fsanitize=address -fsanitize=undefined # -Wsuggest-attribute=malloc

ifeq ($(shell uname),Darwin)
    UNDEFSTR=-undefined dynamic_lookup
else
    UNDEFSTR=
endif


EX=$(patsubst src/%.cpp,%,$(wildcard src/*.cpp))
all: $(EX)
run_tests: $(EX) lztest
	for i in $(EX) lztest; do ./$$i; done

STD?=-std=c++14

INCLUDES=-I`python3-config --includes` -Ipybind11/include
SUF=`python3-config --extension-suffix`
OBJS=$(patsubst %.cpp,%$(SUF),$(wildcard *.cpp))
HEADERS=$(wildcard *.h)

python: _hll.cpython.so
	python -c "import subprocess;import site; subprocess.check_call('cp hll.py "*`python3-config --extension-suffix`" %s' % site.getsitepackages()[0], shell=True)"

%.cpython.so: %.cpp
	$(CXX) $(UNDEFSTR) $(INCLUDES) -O3 -Wall $(FLAGS) -shared $(STD) -fPIC `python3 -m pybind11 --includes` $< -o $*$(SUF) -lz && \
    ln -fs $*$(SUF) $@

%.o: %.cpp
	$(CXX) -c $(FLAGS) $(STD)	$< -o $@

%.o: %.c
	$(CC) -c $(FLAGS)	$< -o $@

%: src/%.cpp kthread.o $(HEADERS)
	$(CXX) $(FLAGS)	$(STD) -Wno-unused-parameter -pthread kthread.o $< -o $@ -lz

%_d: src/%.cpp kthread.o $(HEADERS)
	$(CXX) $(FLAGS)	$(STD) -fsanitize=leak -fsanitize=undefined -Wno-unused-parameter -pthread kthread.o $< -o $@ -lz

lztest: src/test.cpp kthread.o $(HEADERS)
	$(CXX) $(FLAGS)	$(STD) -Wno-unused-parameter -pthread kthread.o -DLZ_COUNTER $< -o $@ -lz

dev_test_p: dev_test.cpp kthread.o hll.h
	$(CXX) $(FLAGS)	$(STD) -Wno-unused-parameter -pthread kthread.o -static-libstdc++ -static-libgcc $< -o $@ -lz

clean:
	rm -f test.o test hll.o kthread.o *hll*cpython*so $(EX)

#mctest: mctest.cpp ccm.h
#mctest_d: mctest.cpp ccm.h

mostlyclean: clean
