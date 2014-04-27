#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Reads the mm5 fields (extracted by extract_surface_fields.f) and reformats
them to be meteorological maps for the defined domain. These maps
are input directly into DHSVM */

#define nclosest 9/*maximum number of points to use in interpolation */
#define nlook 2    /*number of cells away from centercell to look valid interp points*/
#define nholder 25 /*size of (nlook*2+1)^2  */
/*nlook = 1 nholder = 9
		    nlook = 2 nholder 25 */

/* to compile cc -lm read_field_weather.c reproject.sub.c byteswap.c -o read_field_weather */

/*modified by Ed M December 2002 to add the aggregation factor */
/*This can be used for example to aggregate hourly mm5 output to 3-hourly */

struct subhead
{
  int junk1,year,month,day,hour,min,sec;
  float fcsttime;
  int field,dim1,dim2,dim3,dim4;
  int junk2;
};

struct fieldweather
{
  float metdata[20];
};

struct word
{
  char ch[80];
};

struct weights
{
  int cell[nclosest];
  float weight[nclosest];
};

struct mainhead
/*note that the headers will read transposed and start from 0 in C*/
{
  int junk1;
  int bhi[50][20];
  float bhr[20][20];
  struct word bhic[50][20];
  struct word bhrc[20][20];
  int junk2;
};

void byte_swap_short(short *buffer, int number_of_swaps);
void byte_swap_long(long *buffer, int number_of_swaps);
int GetNumber(char *numberStr);
float GetFloat(char *numberStr);
int   GetAllIJElev();
int   prepareforecast();
int   aggreg_outputforecast(int *stepcounter);
int   outputforecast();
int   GetInterpWeights();
int   ll2ne(float lat, float lon, float *pparam, float *north, float *east);
float CalcDistance(float lat1, float lon1, float lat2, float lon2);
int isclose(float val1, float val2);
FILE  *fieldsfile;
char elevfilename[255],latlonfilename[255],fieldsfilename[255];
char headerfilename[255],outputsuffix[50],outputfilename[255];
char outputpath[255];
static char *outputprefix[255] = {"elev","temp","hum","wsp","sw","lw","prec","lpse"}; 
unsigned long OffSet;
int dumpsize;   
int i,j;
int nread;
int aggreg_factor,stepcounter;
struct subhead head;
long int fileend;
int ilat,ilon;
int gotdata;
float *data;
/*8,9,10,11,12,13,14,15,16 = sfp,sw,lw,lapse,precip,u10,v10,t2,q2*/
float lasttime,lastoutputtime,inittime;
float starthour, endhour;
int icount,ncells,nrows,ncols;
FILE *elevfile,*latlonfile;
float *targetlat,*targetlon,*allelev,*iposition,*jposition,*lastprecip;
float *mm5elev,*mm5lat,*mm5lon,*localelev;
struct weights *icells;
struct fieldweather  *rawdata;
float *outt,*outhum,*outwsp,*outsw,*outlw,*outprec,*outlpse,*outmm5elev;
float *outtag,*outhumag,*outwspag,*outswag,*outlwag,*outprecag,*outlpseag,*outmm5elevag;
FILE *outputfile[8];
int byteswapin,byteswapout;
int main(int argc, char **argv)
{
  
  if(argc<12) {
    printf("usage is:\n\tread_field_weather byte_swap_input(0/1) byte_swap_output(0/1) nrows ncols elevfile\n\tlatlonfile fieldsfile headerfile outputpath outputsuffix aggregation_factor\n");
    printf("output files will named xxxx.outputsuffix where xxxx is elev,temp,hum,wsp,sw,lw,prec,lpse\n");
    printf("make sure that byte order of fields and elev file match and byteswap is set correctly\n");
    printf("aggregation factoris the number of time steps that are aggregated into one (set to 1 for no aggregation)\n");
    exit(-1);
  }
  byteswapin = GetNumber(argv[1]);
  byteswapout = GetNumber(argv[2]);
  nrows = GetNumber(argv[3]);
  ncols = GetNumber(argv[4]);
  strcpy(elevfilename, argv[5]);
  strcpy(latlonfilename, argv[6]);
  strcpy(fieldsfilename, argv[7]);
  strcpy(headerfilename, argv[8]);
  strcpy(outputpath,   argv[9]);
  strcpy(outputsuffix, argv[10]);
  aggreg_factor = GetNumber(argv[11]);
  ncells = nrows*ncols;
  targetlat = calloc(ncells, sizeof(float));
  targetlon = calloc(ncells, sizeof(float));
  icells    = calloc(ncells, sizeof(struct weights));
  allelev   = calloc(ncells, sizeof(float));
  iposition = calloc(ncells, sizeof(float));
  jposition = calloc(ncells, sizeof(float));
  lastprecip= calloc(ncells, sizeof(float));
  mm5elev   = calloc(ncells, sizeof(float));
  mm5lat    = calloc(ncells, sizeof(float));
  mm5lon    = calloc(ncells, sizeof(float));
  localelev = calloc(ncells, sizeof(float));
  outt = calloc(ncells, sizeof(float));
  outmm5elev = calloc(ncells,sizeof(float));
  outhum = calloc(ncells, sizeof(float));
  outwsp = calloc(ncells, sizeof(float));
  outsw = calloc(ncells, sizeof(float));
  outlw = calloc(ncells, sizeof(float));
  outprec = calloc(ncells, sizeof(float));
  outlpse = calloc(ncells, sizeof(float));
  outtag = calloc(ncells, sizeof(float));
  outmm5elevag = calloc(ncells,sizeof(float));
  outhumag = calloc(ncells, sizeof(float));
  outwspag = calloc(ncells, sizeof(float));
  outswag = calloc(ncells, sizeof(float));
  outlwag = calloc(ncells, sizeof(float));
  outprecag = calloc(ncells, sizeof(float));
  outlpseag = calloc(ncells, sizeof(float));
  rawdata   = calloc(ncells, sizeof(struct fieldweather));
  
  stepcounter = 0;
  lasttime = 0.0;
  lastoutputtime = 0.0;
  inittime = -10000.0;
  for(i=0;i<ncells;i++)
  lastprecip[i] = 0.0;

  if (!(elevfile = fopen(elevfilename, "rb")))
    {
      printf("could not open %s \n",elevfilename);
      exit(-1);
    }
  if (fread(localelev,sizeof(float),ncells,elevfile)!=ncells) {
    printf("expected %d cells in %s file\n",ncells,elevfilename);
    exit(-1);
}
  if(byteswapin) byte_swap_long(localelev,ncells);

  if (!(latlonfile = fopen(latlonfilename, "rt")))
    {
      printf("could not open %s \n",latlonfilename);
      exit(-1);
    }
  if (!(fieldsfile = fopen(fieldsfilename, "rb")))
    {
      printf("could not open %s \n",fieldsfilename);
      exit(-1);
    }
  
  for(i=0;i<8;i++) {
     sprintf(outputfilename,"%s/%s.%s",outputpath,outputprefix[i],outputsuffix);
     printf("%s\n",outputfilename);
     outputfile[i] = fopen(outputfilename,"wb");
  }
 
  for(i=0;i<ncells;i++)
    fscanf(latlonfile,"%f %f\n",&targetlon[i],&targetlat[i]);
  
  printf("getting mm5 domain info and interpolation weights\n");
  GetAllIJElev();
  printf("done with that\n");
  
  fseek(fieldsfile,0,SEEK_END);
  fileend = ftell(fieldsfile);
  fseek(fieldsfile,0,SEEK_SET);
  /*  printf("fileend= %d bytes\n",fileend);*/

  while(ftell(fieldsfile)<fileend) { 
    fread(&head,sizeof(struct subhead),1,fieldsfile);
    if(byteswapin) byte_swap_long(&head,14);
    /* printf("head.fcsttime= %f inittime= %f\n",head.fcsttime,inittime);*/
    if(inittime < -100.0 && head.fcsttime>inittime) {
	inittime = head.fcsttime;
	lastoutputtime = inittime;
    }
    if(head.fcsttime>lastoutputtime) {
      stepcounter++;
      printf("preparing field forecast\n");
      prepareforecast();
      if(aggreg_factor != 1) {
	printf("aggregating hour %0.2f\n\n",lastoutputtime/60);
	aggreg_outputforecast(&stepcounter);
      }
      else {
	printf("outputting hour %0.2f\n\n",lastoutputtime/60);
	outputforecast();
      }
      gotdata=0;
      lastoutputtime = head.fcsttime;
    }

    /* a total of nine fields are extracted */
    if(head.field >= 8 && head.field <= 16) {
      gotdata++;
      if(gotdata==1)printf("reading fields for hour %0.2f\n",head.fcsttime/60);
	data = calloc(head.dim1*head.dim2, sizeof(float));
	fread(&dumpsize,sizeof(int),1,fieldsfile);
	fread(data, sizeof(float), head.dim1*head.dim2, fieldsfile);
	if(byteswapin) byte_swap_long(data,head.dim1*head.dim2);
	fread(&dumpsize,sizeof(int),1,fieldsfile);
	for(i=0;i<ncells;i++) {
	  rawdata[i].metdata[head.field] = 0;
	  for(j=0;j<nclosest;j++) {
	    rawdata[i].metdata[head.field] += data[icells[i].cell[j]]*icells[i].weight[j];
	  }
	}
	/*	printf("mid sample of this field %f\n",rawdata[(int)(ncells/2)].metdata[head.field]); */
	free(data);
      }
    else
      {
	fseek(fieldsfile,(head.dim1*head.dim2+2)*4,SEEK_CUR);    
      }
  } /*end of file reached*/

  /* print out last time step to output file, if record complete */
  if(gotdata==9) {
    stepcounter++;
    prepareforecast();
    if(aggreg_factor != 1) {
      printf("aggregating time %0.2f\n",lastoutputtime/60);
      aggreg_outputforecast(&stepcounter);
      }
    else {
      printf("outputting hour %0.2f\n",lastoutputtime/60);
      outputforecast();
    }
  }
  

}


int prepareforecast() {
  
  int i; 
  int k;
  float sfp,sw,lw,lps,accp,u,v,t,q,ea,es,wsp,wdir,rh;
  float curp;

  for(i=0;i<ncells;i++) {
    sfp = rawdata[i].metdata[8];
    outsw[i]  = rawdata[i].metdata[9];
    outlw[i]  = rawdata[i].metdata[10];
    outlpse[i] = rawdata[i].metdata[11];
    accp= rawdata[i].metdata[12];
    u   = rawdata[i].metdata[13];
    v   = rawdata[i].metdata[14];
    outt[i]   = rawdata[i].metdata[15]-273.15; /*convert from K to C*/
    q   = rawdata[i].metdata[16];
    
    
    /* calculate RH before lapse temp to avoid
       having elevation artifacts in surface humidity*/
    es = 611.2*exp(17.67*outt[i]/(outt[i]+243.5));
    ea = q*sfp/0.622;
    outhum[i] = ea/es*100.0;
    if(outhum[i] > 100.0) outhum[i] = 100.0;
    if(outhum[i] < 10.0) outhum[i] = 10.0;

    /*lapse temperature */
    if(localelev[i] > -999) {
      outt[i] += (localelev[i]-mm5elev[i])*outlpse[i];
    }
   
    
    /* calculate precip */
    curp = (accp - lastprecip[i])/100.; /* convert from cm to m*/
    outprec[i] = curp;
    lastprecip[i] = accp;
    if(curp < 0) curp = 0;
    
    /* calculate wsp */
    outwsp[i] = sqrt(u*u+v*v);
    outmm5elev[i] = mm5elev[i];
  }
/*	k=(int)(ncells/2);
	printf("midpoint sample %f %f %f %f %f %f %f %f\n",
		mm5elev[k],outt[k],outhum[k],outwsp[k],outsw[k],outlw[k],
		outprec[k],outlpse[k]);*/
	
}

int outputforecast() {

  if(byteswapout) {
    byte_swap_long(outmm5elev,ncells);
    byte_swap_long(outwsp,ncells);
    byte_swap_long(outt,ncells);
    byte_swap_long(outhum,ncells);
    byte_swap_long(outsw,ncells);
    byte_swap_long(outlw,ncells);
    byte_swap_long(outprec,ncells);
    byte_swap_long(outlpse,ncells);
  }
  fwrite(outmm5elev,sizeof(float),ncells,outputfile[0]);
  fwrite(outt,sizeof(float),ncells,outputfile[1]);
  fwrite(outhum,sizeof(float),ncells,outputfile[2]);
  fwrite(outwsp,sizeof(float),ncells,outputfile[3]);
  fwrite(outsw,sizeof(float),ncells,outputfile[4]);
  fwrite(outlw,sizeof(float),ncells,outputfile[5]);
  fwrite(outprec,sizeof(float),ncells,outputfile[6]);
  fwrite(outlpse,sizeof(float),ncells,outputfile[7]);
}

int aggreg_outputforecast(int *stepcounter) {

  int i; 

  /* aggregate output */
  for(i=0;i<ncells;i++) {
    outmm5elevag[i] = outmm5elev[i];
    outwspag[i] += outwsp[i]/aggreg_factor;
    outtag[i] += outt[i]/aggreg_factor;
    outhumag[i] += outhum[i]/aggreg_factor;
    outswag[i] += outsw[i]/aggreg_factor;
    outlwag[i] += outlw[i]/aggreg_factor;
    outprecag[i] += outprec[i];
    outlpseag[i] += outlpse[i]/aggreg_factor;
  }

  /* output and reset aggregated values to zero */
  if(*stepcounter == aggreg_factor) {
    printf("writing aggregated output\n");
    if(byteswapout) {
      byte_swap_long(outmm5elevag,ncells);
      byte_swap_long(outwspag,ncells);
      byte_swap_long(outtag,ncells);
      byte_swap_long(outhumag,ncells);
      byte_swap_long(outswag,ncells);
      byte_swap_long(outlwag,ncells);
      byte_swap_long(outprecag,ncells);
      byte_swap_long(outlpseag,ncells);
    }
    fwrite(outmm5elevag,sizeof(float),ncells,outputfile[0]);
    fwrite(outtag,sizeof(float),ncells,outputfile[1]);
    fwrite(outhumag,sizeof(float),ncells,outputfile[2]);
    fwrite(outwspag,sizeof(float),ncells,outputfile[3]);
    fwrite(outswag,sizeof(float),ncells,outputfile[4]);
    fwrite(outlwag,sizeof(float),ncells,outputfile[5]);
    fwrite(outprecag,sizeof(float),ncells,outputfile[6]);
    fwrite(outlpseag,sizeof(float),ncells,outputfile[7]);
    for(i=0;i<ncells;i++) {
      outwspag[i] = 0;
      outtag[i] = 0;
      outhumag[i] = 0; 
      outswag[i] = 0; 
      outlwag[i] = 0;
      outprecag[i] = 0; 
      outlpseag[i] = 0;
      *stepcounter=0;
    }
  }
}


int GetAllIJElev()
{
  FILE  *headerfile;
  unsigned long OffSet;
  int dumpsize;   
  int i,j;
  int nread;
  struct subhead head;
  struct mainhead mhead;
  long int fileend;
  float *elev;
  float *lat,*lon;
  int icoord,jcoord;
  float centerlat,centerlon,truelat1,truelat2,gridspacing;
  float north,east;
  float ioffset,joffset;
  int   icell,icoarse,jcoarse,scale;
  float distcoarse;
  int   gridi,gridj;
  float pparam[4],sumdist;
  int centercell,currentcell,icount,ii,jj;
  int exact;
  float mindist;
  int localcell[nholder];
  float localdist[nholder];
  float temp;
  float itemp;
  int lookn,looks,looke,lookw;

  printf("byteswapin is %d\n",byteswapin);

  if(!(headerfile = fopen(headerfilename,"rb"))) {
    printf("error opening file header file %s\n",headerfilename);
    exit(-1);
  }
  /* ok, so how big is this thing */
  fseek(headerfile,0,SEEK_END);
  fileend = ftell(headerfile);
  /*got the size, now move back to the beginning*/
  fseek(headerfile,0,SEEK_SET);
  fread(&mhead,sizeof(struct mainhead),1,headerfile);
  fread(&head,sizeof(struct subhead),1,headerfile);
  if(byteswapin) byte_swap_long(&head,14);
  printf("size is %d %d\n",head.dim1,head.dim2);
   
  elev = calloc(head.dim1*head.dim2, sizeof(float));
  fread(&dumpsize,sizeof(int),1,headerfile);
  fread(elev, sizeof(float), head.dim1*head.dim2, headerfile);
  if(byteswapin) byte_swap_long(elev,head.dim1*head.dim2); 
  fread(&dumpsize,sizeof(int),1,headerfile);
  printf("got elev\n");

  fread(&head,sizeof(struct subhead),1,headerfile);
  if(byteswapin) byte_swap_long(&head,14);
  lat = calloc(head.dim1*head.dim2, sizeof(float));
  fread(&dumpsize,sizeof(int),1,headerfile);
  fread(lat, sizeof(float), head.dim1*head.dim2, headerfile);
  if(byteswapin) byte_swap_long(lat,head.dim1*head.dim2);  
  fread(&dumpsize,sizeof(int),1,headerfile);
  printf("got lat\n");

  fread(&head,sizeof(struct subhead),1,headerfile);
  if(byteswapin) byte_swap_long(&head,14);
  lon = calloc(head.dim1*head.dim2, sizeof(float));
  fread(&dumpsize,sizeof(int),1,headerfile);
  fread(lon, sizeof(float), head.dim1*head.dim2, headerfile); 
  if(byteswapin) byte_swap_long(lon,head.dim1*head.dim2); 
  fread(&dumpsize,sizeof(int),1,headerfile);
  printf("got lon\n");
  fclose(headerfile);

  centerlat = mhead.bhr[0][1];
  centerlon = mhead.bhr[0][2];
  truelat1  = mhead.bhr[0][4];
  truelat2  = mhead.bhr[0][5];
  gridspacing = mhead.bhr[0][8];
  icoarse     = mhead.bhi[0][4];
  jcoarse     = mhead.bhi[0][5];
  distcoarse  = mhead.bhr[0][0];
  ioffset     = mhead.bhr[0][9];
  joffset     = mhead.bhr[0][10];
  scale       = mhead.bhi[0][19];
  gridi       = mhead.bhi[0][15];
  gridj       = mhead.bhi[0][16];
  pparam[0]   = mhead.bhr[0][4];
  pparam[1]   = mhead.bhr[0][5];
  pparam[2]   = mhead.bhr[0][2];
  pparam[3]   = mhead.bhr[0][1];
  if(byteswapin) {
  byte_swap_long(&centerlat,1);
  byte_swap_long(&centerlon,1);
  byte_swap_long(&truelat1,1);
  byte_swap_long(&truelat2,1);
  byte_swap_long(&gridspacing,1);
  byte_swap_long(&icoarse,1);
  byte_swap_long(&jcoarse,1);
  byte_swap_long(&distcoarse,1);
  byte_swap_long(&ioffset,1);
  byte_swap_long(&joffset,1);
  byte_swap_long(&scale,1);
  byte_swap_long(&gridi, 1);
  byte_swap_long(&gridj, 1);
  byte_swap_long(pparam, 4);
  }
  printf("%f %f %f %f\n",centerlat,centerlon,truelat1,truelat2);
  for(i=0;i<ncells;i++) {
    /*    printf("doing %f %f \n",targetlat[i],targetlon[i]);*/
    ll2ne(targetlat[i], targetlon[i], pparam, &north, &east);
    iposition[i] = ((float)icoarse)/2+north/distcoarse;
    jposition[i] = ((float)jcoarse)/2+east/distcoarse;
    /* in dot point space of this grid */
    iposition[i] = ((float)iposition[i]-ioffset)*(float)scale+(float)scale/2;
    jposition[i] = ((float)jposition[i]-joffset)*(float)scale+(float)scale/2;
    if(iposition[i] < 0. || iposition[i]>(float)gridi || jposition[i]<0. || jposition[i] > (float)gridj) {
      printf("location %f %f is not in current mm5grid\n",targetlat[i],targetlon[i]);
    }
    icoord = (int)iposition[i];
    jcoord = (int)jposition[i];
    /*   printf("%f %f\n",iposition[i],jposition[i]);*/
    centercell = icoord + jcoord*gridi;
    icount = 0;
    sumdist = 0;
    exact = -1;

     for(ii=-nlook;ii<=nlook;ii++)
      for(jj=-nlook;jj<=nlook;jj++) {
	currentcell = (icoord + ii) + (jcoord + jj)*gridi;
	localcell[icount] = currentcell;
	localdist[icount] = CalcDistance(lat[currentcell],lon[currentcell],targetlat[i],targetlon[i]);
	if(localdist[icount] <= 0) exact = icount;
	icount++;
      }

     for(ii=0;ii<nholder;ii++) 
       for(jj=ii+1;jj<nholder;jj++) {
	 if(localdist[jj]<localdist[ii]) {
	   temp = localdist[ii];
	   localdist[ii] = localdist[jj];
	   localdist[jj] = temp;
	   itemp = localcell[ii];
	   localcell[ii] = localcell[jj];
	   localcell[jj] = itemp;
	 }
       }

     if(exact>=0) {
       for(ii=0;ii<nclosest;ii++) {
	 icells[i].cell[ii] = localcell[ii];
	 if(ii==exact) 
	   icells[i].weight[ii]=1.0;
	 else
	   icells[i].weight[ii]=0.0;
       }
       sumdist = 1.0;
     }
     else{
     sumdist = 0;
     for(ii=0;ii<nclosest;ii++) {
       icells[i].cell[ii] = localcell[ii];
       icells[i].weight[ii] = 1/(localdist[ii]*localdist[ii]);
       sumdist += icells[i].weight[ii];
     }
     }

     for(ii=0;ii<nclosest;ii++) 
       icells[i].weight[ii] /= sumdist;
 
     mm5elev[i]=0;
          for(ii=0;ii<nclosest;ii++){
     mm5elev[i] += elev[icells[i].cell[ii]]*icells[i].weight[ii];
     }
     mm5lat[i]  = lat[centercell];
     mm5lon[i]  = lon[centercell];
  }  
}

float CalcDistance(float lat1, float lon1, float lat2, float lon2) {
  float distance;
  distance = sqrt((lat1-lat2)*(lat1-lat2)+(lon1-lon2)*(lon1-lon2));
  return distance;
  
}

/*****************************************************************************
  GetNumber()
*****************************************************************************/
int GetNumber(char *numberStr) 
{
  char *endPtr;
  int number = 0;
  
  number = (int) strtol(numberStr, &endPtr, 0);
  if (*endPtr != '\0')
    exit(-1);
  
  return number;
}

float GetFloat(char *numberstr) 
{
  char *endptr;
  float number = 0;
  
  number = (float) strtod(numberstr, &endptr);
  if (*endptr != '\0')
    exit(-1);
  
  return number;
}

int isclose(float val1, float val2) {
  if(fabs(val1-val2)<5.) {
    return 1;
  }
  else
    {
      return 0;
    }

}








