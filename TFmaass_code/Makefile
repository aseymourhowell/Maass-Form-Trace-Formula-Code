INC=/data/local/include
LIB=/data/local/lib

maass: maasscomp.c
	gcc -O2 -I$(INC) -I$(INC)/flint -o maass maasscomp.c $(LIB)/libarb.a $(LIB)/libflint.a $(LIB)/libmpfr.a $(LIB)/libgmp.a -lm ./factor64-master/factor64.o
