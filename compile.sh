#!/bin/sh
g++ -O4 -fopenmp -g  svdDynamic.c RayTracer.c utils.c -lm -o RayTracer
