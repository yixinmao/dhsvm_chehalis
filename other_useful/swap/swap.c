#include <stdio.h>
#include <stdlib.h>
#include <string.h>

main(int argc, char *argv[])
{
  char   h,k,l,m;
  int    a;
  FILE  *infile, *outfile;

  if(argc!=3) {
    fprintf(stdout,"Use: %s  <input file> <output file>\n",argv[0]);
    exit(0);
  }

  if((infile=fopen(argv[1], "rb"))==NULL) {
    fprintf(stderr,"%s: ERROR: Unable to open %s\n",argv[0],argv[1]);
  }

  if((outfile=fopen(argv[2],"wb"))==NULL) {
    fprintf(stderr,"%s: ERROR: Unable to open %s\n",argv[0],argv[2]);
  }

  a=0;
  while(fread(&h,1,1,infile)==1) {
    fread(&k,1,1,infile);
    fread(&l,1,1,infile);
    fread(&m,1,1,infile);
//    printf("%f\n",h);
    fwrite(&m,1,1,outfile);
    fwrite(&l,1,1,outfile);
    fwrite(&k,1,1,outfile);
    fwrite(&h,1,1,outfile);
    a++;
  }
  fprintf(stdout,"swap: count is %d\n",a);
  fclose(infile);
  fclose(outfile);
}
