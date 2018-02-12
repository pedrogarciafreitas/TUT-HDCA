CC=gcc -g -O3 
CFLAGS=-I. 
DEPS = 

%.o: %.cc $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

cerv.a: adaptive_model.o arithmetic_decode.o arithmetic_encode.o bit_input.o  bit_output.o  cerv.o 
	ar rsc cerv.a adaptive_model.o arithmetic_decode.o arithmetic_encode.o bit_input.o  bit_output.o  cerv.o

clean:
	rm *.o *.a
