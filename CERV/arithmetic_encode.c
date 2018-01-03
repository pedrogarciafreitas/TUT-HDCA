/* ARITHMETIC ENCODING ALGORITHM. */

#include "arithmetic_coding.h"
#include <stdio.h>

static void bit_plus_follow();	/* Routine that follows                     */


/* CURRENT STATE OF THE ENCODING. */

static code_value low, high;	/* Ends of the current code region          */
static long bits_to_follow;	/* Number of opposite bits to output after  */
				/* the next bit.                            */
static code_value bits_to_follow2;	/* Currently-seen code value                */
static code_value low2, high2;	/* Ends of current code region              */

void saveReg(void)
{
	bits_to_follow2 = bits_to_follow;
	low2 = low;
	high2 = high;

}
void loadReg(void)
{
	bits_to_follow = bits_to_follow2;
	low = low2;
	high = high2;
}

/* START ENCODING A STREAM OF SYMBOLS. */

start_encoding()
{   low = 0;					/* Full code range.         */
    high = Top_value;
    bits_to_follow = 0;				/* No bits to follow next.  */
}


/* ENCODE A SYMBOL. */

encode_symbol(symbol, cum_freq, outputFile)
    int symbol;			/* Symbol to encode                         */
    int cum_freq[];		/* Cumulative symbol frequencies            */
	FILE* outputFile;
{   long range;			/* Size of the current code region          */
    range = (long)(high-low)+1;
    high = low +				/* Narrow the code region   */
      (range*cum_freq[symbol-1])/cum_freq[0]-1;	/* to that allotted to this */
    low = low + 				/* symbol.                  */
      (range*cum_freq[symbol])/cum_freq[0];
    for (;;) {					/* Loop to output bits.     */
        if (high<Half) {
            bit_plus_follow(0, outputFile);			/* Output 0 if in low half. */
        } 
        else if (low>=Half) {			/* Output 1 if in high half.*/
            bit_plus_follow(1, outputFile);
            low -= Half;
            high -= Half;			/* Subtract offset to top.  */
        }
        else if (low>=First_qtr			/* Output an opposite bit   */
              && high<Third_qtr) {		/* later if in middle half. */
            bits_to_follow += 1;
            low -= First_qtr;			/* Subtract offset to middle*/
            high -= First_qtr;
        }
        else break;				/* Otherwise exit loop.     */
        low = 2*low;
        high = 2*high+1;			/* Scale up code range.     */
    }
}


/* FINISH ENCODING THE STREAM. */

done_encoding(FILE* outputFile)
{   bits_to_follow += 1;			/* Output two bits that     */
    if (low<First_qtr) 
		bit_plus_follow(0, outputFile);	/* select the quarter that  */
    else 
		bit_plus_follow(1, outputFile);			/* the current code range   */
}						/* contains.                */


/* OUTPUT BITS PLUS FOLLOWING OPPOSITE BITS. */
long int nrBiti;
int estim =0;
static void bit_plus_follow(int bit,FILE* outputFile)
{   
	if (estim)
		nrBiti++;
	else
		output_bit(bit, outputFile);				/* Output the bit.          */
    while (bits_to_follow>0) {
		if (estim)
			nrBiti++;
		else
			output_bit(!bit, outputFile);			/* Output bits_to_follow    */
        bits_to_follow -= 1;			/* opposite bits. Set       */
    }						/* bits_to_follow to zero.  */
}

void startEstim(int estimation)
{
	nrBiti=0;
	estim = estimation;
}
long int getEstim(void)
{
	return nrBiti;
}