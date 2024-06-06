gcc -c main.c svd.c image.c
gcc *.o -o app -L cblas.a -lblas -llapack -llapacke -lm