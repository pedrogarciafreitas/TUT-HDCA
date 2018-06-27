CC=g++
CFLAGS=-I. -std=c++11
DEPS = include/gen_types.hh include/warpingFunctions.hh 
OBJ_ENCODER = TUT-HDCA-Encoder/tut-hdca-encoder.o
OBJ_DECODER = TUT-HDCA-Decoder/tut-hdca-decoder.o 

%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

tut-hdca-encoder: $(OBJ_ENCODER)
	$(CC) -o $@ $^ $(CFLAGS)
	
tut-hdca-decoder: $(OBJ_DECODER)
	$(CC) -o $@ $^ $(CFLAGS)
	
all: tut-hdca-encoder tut-hdca-decoder

clean:
	rm TUT-HDCA-Decoder/tut-hdca-decoder.o TUT-HDCA-Encoder/tut-hdca-encoder.o
