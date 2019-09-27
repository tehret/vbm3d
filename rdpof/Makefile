CC=gcc
CXX=g++
CFLAGS=-Wall -Wextra -O3 -Werror -I../iio/
LFLAGS=iio.o -lpng -ljpeg -ltiff -fopenmp

main: main.cpp robust_expo_methods.h mask.h bicubic_interpolation.h smoothness.h zoom.h gaussian.h iio.o
	$(CXX) -std=c++98 $(CFLAGS) -o main main.cpp  $(LFLAGS)

iio.o: ../iio/iio.c
	$(CC) -std=c99 $(CFLAGS) -Wno-unused -pedantic -DNDEBUG -D_GNU_SOURCE -c ../iio/iio.c

clean:
	rm -f iio.o main.o main
