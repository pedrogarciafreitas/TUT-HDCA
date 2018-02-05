/* BIT INPUT ROUTINES. */

#include "bit_input.h"
#include <stdio.h>


/* THE BIT BUFFER. */

static int buffer ;		/* Bits waiting to be input                 */
static int bits_to_go;		/* Number of bits still in buffer           */


/* INITIALIZE BIT INPUT. */

void start_inputing_bits()
{   
	bits_to_go = 0;				/* Buffer starts out with   */
}						/* no bits in it.           */

/* SAVE THE POSITION */
/*int savePos(FILE* outputFile)
{
	FILE * fid;
	long int pos;
	fid = fopen("positions.txt","wt");
	if (fid == NULL)
		return 1;
	pos = ftell(outputFile);
	fprintf(fid, "%d %d %ld ", bits_to_go, buffer, pos );
	fclose(fid);
	saveReg();
	return 0;
}*/

/* LOAD THE POSITION*/
/*int loadPos(FILE* outputFile)
{
	FILE * fid;
	int val, lo, hi;
	long int pos = 0;
	fid = fopen("positions.txt","rt");
	if (fid == NULL)
		return 1;
	fscanf(fid, "%d %d %ld %d %d %d ", &bits_to_go, &buffer, &pos, &val, &lo, &hi);
	loadReg(val, lo, hi);
	fseek(outputFile, pos, SEEK_SET);
	fclose(fid);
	return 0;
}
*/

/* INPUT A BIT. */

int input_bit(FILE* outputFile)
{   int t;
    bits_to_go -= 1;
    if (bits_to_go<0) {				/* Read the next byte if no */
        if(fread(&buffer,1,1,outputFile) != 1) {
          printf("Can not read!\n");
        }
		//buffer = getc(stdin);			/* bits are left in the     */
        bits_to_go = 7;				/* buffer. Return anything  */
    }			   			/* after end-of-file.       */
    t = buffer&1;
    buffer >>= 1;				/* Return the next bit from */
    return t;					/* the bottom of the byte.  */
}
