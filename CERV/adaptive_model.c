/* THE ADAPTIVE SOURCE MODEL */

#include "model.h"
#include <stdlib.h>

/* INITIALIZE THE MODEL. */

void start_model()
{   
	int i;
	cum_freq[My_no_of_symbols] = 1;
	for (i=0;i<16383;i++)
	{
		freq[i] = 0;
		cum_freq[i] = 0;
	}
    for (i = 0; i<=My_no_of_symbols; i++)
	{	
		freq[i] = 1;
        cum_freq[i] = No_of_symbols-i;
    }
	freq[0] = 0;
}




/* UPDATE THE MODEL TO ACCOUNT FOR A NEW SYMBOL. */

update_model(symbol)
int symbol;
{
	int i;			                     /* New index for symbol     */
    if (cum_freq[0]>=Max_frequency) 
	{		                             /* See if frequency counts  */
        int cum;				         /* are at their maximum.    */
        cum = 0;
        for (i = No_of_symbols; i>=0; i--) 
		{	                            /* If so, halve all the     */
            freq[i] = (freq[i]+1)/2;	/* counts (keeping them     */
            cum_freq[i] = cum; 			/* non-zero).               */
            cum += freq[i];
        }
    }
	
	if (MODEL==1)
		freq[symbol] += 1;				/* Increment the frequency  */
	else
		freq[symbol] += No_of_symbols-1;

    while (symbol>0) {				/* count for the symbol and */
        symbol -= 1;					/* update the cumulative    */
        //cum_freq[symbol] += 1;		/* frequencies.             */
		if (MODEL==1)
			cum_freq[symbol] += 1;				/* Increment the frequency  */
		else
			cum_freq[symbol] += No_of_symbols-1;
    }

}

int** alocaMatrice2(int m, int n)
{
	int **matrix, i, j;
	matrix = (int**)malloc(m*sizeof(int*));
	if (matrix == NULL) 
		exit(1);
	for (i=0; i<m; i++)
	{
		matrix[i] = (int*)malloc(n*sizeof(int));
		if (matrix[i] == NULL) 
			exit(1);
		for (j=0; j<n; j++)
			matrix[i][j] = 0;
	}
	return matrix;
}


/* INITIALIZE THE MODEL. */
//*
void start_model2()
{   
	int i,j;
	
	//aloca
	cum_freq2 = alocaMatrice2(My_no_of_rows, NS+1);
	freq2 = alocaMatrice2(My_no_of_rows, NS+1);
	// initializeaza
	for (i=0; i<My_no_of_rows; i++)
	{
		for (j=0; j<=My_no_of_symbols; j++)
		{	
			freq2[i][j] = 1;
			cum_freq2[i][j] = NS-j;  // cum_freq2[i][j] = My_no_of_symbols-j;
		}
		freq2[i][0] = 0;
	}
}



update_model2(symbol, row)
int symbol;
int row;
{
	int i;			                    
    if ( ((row==0) & (cum_freq2[row][0]>=Max_frequency)) | ((row>0) & (cum_freq2[row][0]>=400)) )  // 16383  100
	{		                            
        int cum;				         
        cum = 0;
        for (i = NS; i>=0; i--) 
		{	                           
            freq2[row][i] = (freq2[row][i]+1)/2;	
            cum_freq2[row][i] = cum; 			
            cum += freq2[row][i];
        }
    }
	
	if (MODEL==1)
		freq2[row][symbol] += 1;				/* Increment the frequency  */
	else
		freq2[row][symbol] += NS;				
    while (symbol>0) {				
        symbol -= 1;					
        //cum_freq2[row][symbol] += 1;	
		if (MODEL==1)
			cum_freq2[row][symbol] += 1;				/* Increment the frequency  */
		else
			cum_freq2[row][symbol] += NS;
    }

}

/* INITIALIZE THE MODEL. */
//*
void start_model3()
{   
	int i,j;
	// number of symbols in each of the 15 contexts
	int nsf3[15] = {256,256,256,256,256,2*MAX_NICE+1,2*MAX_NICE+1,2*MAX_NICE+1,2*MAX_NICE+1,2*MAX_NICE+1,2,2,2,2,2};
	
	//aloca
	cum_freq3 = alocaMatrice2(15, 256+1);  // conservative, the largest symbol is 256
	freq3 = alocaMatrice2(15, 256+1);
	// initializeaza
	for (i=0; i<15; i++)
	{
		for (j=0; j<=nsf3[i]; j++)
		{	
			freq3[i][j] = 1;
			cum_freq3[i][j] = nsf3[i]-j;  // cum_freq2[i][j] = My_no_of_symbols-j;
		}
		freq3[i][0] = 0;
	}
}



update_model3(symbol, row)
int symbol;
int row;
{
	int i;			
	// number of symbols in each of the 15 contexts
	int nsf3[15] = {256,256,256,256,256,2*MAX_NICE+1,2*MAX_NICE+1,2*MAX_NICE+1,2*MAX_NICE+1,2*MAX_NICE+1,2,2,2,2,2};

    if (cum_freq3[row][0]>=Max_frequency) 
	{		                            
        int cum;				         
        cum = 0;
        for (i = nsf3[row]; i>=0; i--) 
		{	                           
            freq3[row][i] = (freq3[row][i]+1)/2;	
            cum_freq3[row][i] = cum; 			
            cum += freq3[row][i];
        }
    }
	
	if (MODEL==1)
		freq3[row][symbol] += 1;				/* Increment the frequency  */
	else
		freq3[row][symbol] += nsf3[row];				
    while (symbol>0) {				
        symbol -= 1;					
        //cum_freq2[row][symbol] += 1;	
		if (MODEL==1)
			cum_freq3[row][symbol] += 1;				/* Increment the frequency  */
		else
			cum_freq3[row][symbol] += nsf3[row];
    }

}

 

 

free_model()
{
	free(cum_freq2);
	free(freq2);
}
free_model2()
{
	int i;
	for (i=0;i<My_no_of_rows;i++)
	{
		free(cum_freq2[i]);
		free(freq2[i]);
	}
}
free_model3()
{
	int i;
	for (i=0;i<15;i++)
	{
		free(cum_freq3[i]);
		free(freq3[i]);
	}
}
//*/
