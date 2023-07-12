INCLUDE = -I.

CC	  = gcc
CFLAGS	  = -g -O2 -m32 -D_FILE_OFFSET_BITS=64 $(INCLUDE)

#use_sm=SMFILE TOMOFILE RECONFILE R NUMBEROFITERATIONS

use_sm: use_sm.o $(OBJS)
	$(CC) $(CFLAGS) -o $@ $@.o $(OBJS) 



