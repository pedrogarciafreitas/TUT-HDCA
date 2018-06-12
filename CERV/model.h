/* INTERFACE TO THE MODEL. */

#define NCON  18

#define NS_NCON 1024*256
#define  LargetsCount 1080*1920*2*2 // for a 4K IMAGE
// #define NCON  5
// #define NS_NCON 1024
#define NS 2
#define MAX_NICE 5*4


/* THE SET OF SYMBOLS THAT MAY BE ENCODED. */
int My_no_of_symbols;
#define EOF_symbol (My_no_of_symbols+1)	/* Index of EOF symbol              */

#define No_of_symbols (My_no_of_symbols+1)	/* Total number of symbols          */
//#define MAX_NICE 6
int My_no_of_rows;

int MODEL;

/* CUMULATIVE FREQUENCY TABLE. */

#define Max_frequency 16383		/* Maximum allowed frequency count  */

int cum_freq[Max_frequency];		    /* Cumulative symbol frequencies */
int freq[Max_frequency];	/* Symbol frequencies            */

int **cum_freq2,**cum_freq3;
int **freq2,**freq3;
