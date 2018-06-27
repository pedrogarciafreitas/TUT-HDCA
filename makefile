CC=gcc
CFLAGS=-I.
DEPS = include/gen_types.hh include/warpingFunctions.hh
OBJ_ENCODER = tut-hdca-encoder.o 
OBJ_DECODER = tut-hdca-decoder.o 

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

tut-hdca-encoder: $(OBJ_ENCODER)
	$(CC) -o $@ $^ $(CFLAGS)
	
tut-hdca-decoder: $(OBJ_DECODER)
	$(CC) -o $@ $^ $(CFLAGS)
	
all: tut-hdca-encoder tut-hdca-decoder