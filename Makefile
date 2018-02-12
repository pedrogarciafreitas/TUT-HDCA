CC=g++ -O3 -std=c++0x
CFLAGS=-I. 
DEPS = 

%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

encoder: GolombCoder/bit_input.o GolombCoder/bit_output.o GolombCoder/golomb_coder.o TUT-HDCA-Encoder/tut-hdca-encoder.o
	$(CC) -o TUT-HDCA-Encoder_bin GolombCoder/bit_input.o GolombCoder/bit_output.o GolombCoder/golomb_coder.o TUT-HDCA-Encoder/tut-hdca-encoder.o cerv/cerv.a $(CFLAGS) 
decoder: GolombCoder/bit_input.o GolombCoder/bit_output.o GolombCoder/golomb_coder.o TUT-HDCA-Decoder/tut-hdca-decoder.o
	$(CC) -o TUT-HDCA-Decoder_bin GolombCoder/bit_input.o GolombCoder/bit_output.o GolombCoder/golomb_coder.o TUT-HDCA-Decoder/tut-hdca-decoder.o cerv/cerv.a $(CFLAGS) 
all: encoder decoder

clean:
	rm GolombCoder/bit_input.o
	rm GolombCoder/bit_output.o
	rm TUT-HDCA-Encoder/tut-hdca-encoder.o
	rm TUT-HDCA-Decoder/tut-hdca-decoder.o