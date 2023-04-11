#include        <stdio.h>
#include        <fcntl.h>
#include        <math.h>
#include    <unistd.h>
#include <stdlib.h>
/*
 *                         Filter to compute logarithm of input data
 *
 *    call 
 *
 *    power compute power of input data. 
 *
 *  data input is r*4
 *
 *
*/
#define STDINP     0
#define STDOUT     1
#define ISAPIPE    1
#define PRGID      "stokes"
#define TRUE       1
#define FALSE      0
#define MAX(a,b)  ((a>b)?a:b)
#define MIN(a,b)  ((a<b)?a:b)

int main(int argc,char **argv)
{
        float    inbuf1[4096];
        float    inbuf2[4096];
        float    outbuf[4096];
        int      dsize;   
        int      wordsinput;
        int      numi1;                              /* read in*/
        int      numi2;                              /* read in*/
        int      num_out;
        int     i,j,numout;
        FILE    *f1;
        FILE    *f2;

        dsize = 8;
        if (argc < 3) {
	  fprintf(stderr, "Usage: stokes file1 file2\n");
          exit(1);
        }
        f1 = fopen(argv[1], "r");
        if (!f1) perror;
        f2 = fopen(argv[2], "r");
        if (!f2) perror;
        if (! (f1 && f2)) exit(1);

        for (;;){
/*
 *      input a record 
*/
                numi1=fread((void *)inbuf1, dsize, sizeof(inbuf1)/dsize, f1);
                numi2=fread((void *)inbuf2, dsize, sizeof(inbuf2)/dsize, f2);
                   if (!(numi1 && numi2)) goto done;  
                if (numi1 != numi2) fprintf(stderr, "Got %d and %d bytes on rad\n");
/*  A B* -> a c + b d, b c - a d */
                 for (i = 0; i < numi1*2; i+=2) {
                   outbuf[i]=inbuf1[i]*inbuf2[i] + inbuf1[i+1]*inbuf2[i+1];
                   outbuf[i+1]=inbuf1[i+1]*inbuf2[i] - inbuf1[i]*inbuf2[i+1];
                 }
                num_out=write(STDOUT,(char *)outbuf,dsize*numi1);
                if (num_out != dsize*numi1){
                   fprintf(stderr,"stokes: output bytes != input bytes\n");
                   exit(-1);
                }
        }
done:   exit(0);
        /*NOTREACHED*/
}    
