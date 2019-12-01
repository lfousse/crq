#GMP = $(HOME)/local
GMP = /usr
MPFR = $(GMP)
CC = gcc
CPPFLAGS = -DDBASE=10 -DVERBOSE=1 -DNDEBUG
CFLAGS = -ggdb -O2 -Wall -ansi
#DFLAGS = -pg -O2 -DDEBUG
