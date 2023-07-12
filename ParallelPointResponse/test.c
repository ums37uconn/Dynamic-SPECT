#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include<string.h>


FILE *file_ptr;
main()
{

file_ptr= (double*)  malloc(24);

file_ptr=fopen("test.txt","r");

fprintf(stderr, "Expected %d paramaters in file %s, only read %d paramaters\n", file_ptr);

}
