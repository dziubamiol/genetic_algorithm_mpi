#!/bin/sh
rm main main.o
mpicc -c main.cpp && mpicc -o main main.o -lm -lstdc++ 
ls
