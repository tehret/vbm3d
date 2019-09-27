CC=gcc
CXX=g++
CFLAGS=-Wall -Wextra -Werror -O4 -I../iio/
LFLAGS=iio.o -lpng -ljpeg -ltiff -fopenmp

rof: main.cpp brox_optic_flow.h mask.h bicubic_interpolation.h  zoom.h gaussian.h iio.o
	$(CXX) -std=c++98 $(CFLAGS) -o rof main.cpp  $(LFLAGS)


iio.o: ../iio/iio.c
	$(CC) -std=c99 $(CFLAGS) -Wno-unused -pedantic -DNDEBUG -D_GNU_SOURCE -c ../iio/iio.c

clean:
	rm -f iio.o main.o rof
