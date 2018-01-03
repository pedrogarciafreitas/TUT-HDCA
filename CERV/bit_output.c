/* BIT OUTPUT ROUTINES. */

#include <stdio.h>


/* THE BIT BUFFER. */

static int buffer;		/* Bits buffered for output                 */
static int bits_to_go;		/* Number of bits free in buffer            */
static int buffer2;		/* Bits waiting to be input                 */
static int bits_to_go2;		/* Number of bits still in buffer           */

/* INITIALIZE FOR BIT OUTPUT. */

start_outputing_bits()
{   buffer = 0;					/* Buffer is empty to start */
    bits_to_go= 8;				/* with.                    */
}


/* OUTPUT A BIT. */

output_bit(int bit,FILE* outputFile)
{   buffer >>= 1; if (bit) buffer |= 0x80;	/* Put bit in top of buffer.*/
    bits_to_go -= 1;
    if (bits_to_go==0) {			/* Output buffer if it is   */
        //putc(buffer,stdout);			/* now full.                */
		fwrite(&buffer,1,1,outputFile);
        bits_to_go = 8;
    }
}


/* FLUSH OUT THE LAST BITS. */

done_outputing_bits(FILE* outputFile)
{   
	//putc(buffer>>bits_to_go,stdout);
	buffer = (buffer>>bits_to_go);
	fwrite(&buffer,1,1,outputFile);
}
savePos()
{
	buffer2 = buffer;
	bits_to_go2 = bits_to_go;
	saveReg();
}
loadPos()
{
	buffer = buffer2;
	bits_to_go = bits_to_go2;
	loadReg();
}

