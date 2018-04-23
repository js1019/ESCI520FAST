#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sac.h"
#include "sacio.c"
#include "ccFreqRM.c"
#define USAGE "%s [-f template_sac_file] [-s long_sac_file] [-g high frequency ][-l low frequency][-o outFile] -h help\n"
static  char    *dataFile1 = NULL;
static  char    *dataFile2 = NULL;
static  char    *outFile = NULL;
void _rmean_(float *aa,int n);
void xapiir( float data[],int nsamps,char *aproto,double trbndw,double a,int iord,char *type,double flo,double fhi,double ts,int  passes);
//double* ccFreq_(const double *x, int nx, const double *y, int ny);
int main(int argc, char *argv[]) 
{

  extern   int getopt();
  extern    char *optarg;
  extern    int optind;
  char          *poro,*type;
  double        *tmp, *cont,*ccout;
  double        trbndw,a,flo,fhi;
  float	        *wavedataM,*wavedataC,*output,dt_m,dt_c,ts,tb;     
  int 		j,nn_m,nn_c,n,n_o,c;
  int           iord,passes;
  size_t	temp_size,cont_size,n_size,no_size;
  SACHEAD       hdm,hdc,hd_o;
  //BP parameters
  poro=strdup("BU");
  type=strdup("BP");
  trbndw=0.0;
  a=0.0;
  iord=4;
  flo=0.0;
  fhi=0.0;
  passes=2;
  ts=0.0;tb=0.0;
//  printf("%lf",a);
  /* input parameters */
  while( (c=getopt( argc, argv, "f:s:g:l:o:h" )) != (-1) )
  {
      switch( c ) 
      {
	  case 'f':
	      dataFile1 = optarg;
	      break;
	  case 's':
	      dataFile2 = optarg;
	      break;
	  case 'g':
	      fhi = atof(optarg);
	      break;
	  case 'l':
	      flo = atof(optarg);
	      break;
	  case 'o':
	      outFile = optarg;
	      break;
	  case 'h':
	      fprintf( stderr, USAGE, argv[0] );
	      exit(1);
	  default:
	      fprintf(stderr, USAGE, argv[0] );
	      exit(1);
      }
  }
  if(argc == 1 || dataFile1 == NULL || dataFile2 == NULL) 
  {
      fprintf(stderr, USAGE, argv[0] );
      exit(1);
  }
//  printf("%s  %s \n ",dataFile1,dataFile2);
  /* input template trace  and contwave */
 if( (wavedataM = read_sac(dataFile1,&hdm)) == NULL ){ return -1;}
 if( (wavedataC = read_sac(dataFile2,&hdc)) == NULL ){ return -1;}
  nn_m = hdm.npts;
  dt_m = hdm.delta;
  nn_c = hdc.npts;
  dt_c = hdc.delta;
    tb = hdc.b;  
  ts = dt_m;
/*  compute the ccoefient */
  //allow the memory
  temp_size = nn_m*sizeof(double);
  cont_size = nn_c*sizeof(double);
  tmp = (double *)malloc(temp_size);
  cont = (double *)malloc(cont_size);
  n = nn_m + nn_c;
  n_size = n*sizeof(double);
  ccout = (double *)malloc(n_size); 
  n_o = nn_c - nn_m + 1;
  no_size = n_o*sizeof(float);
  output = (float *)malloc(no_size); 
  ////////////////////////////////////
//printf("%d  %d\n",nn_m,nn_c);
      for( j = 0 ; j < nn_c ; j++ ) 
      {
	  if(fabs(wavedataC[j])>=pow(10,9))
	  {
	      if(j==0) wavedataC[j] = 0.0;
	      else wavedataC[j] = wavedataC[j-1];
	  }
       }
      for( j = 0 ; j < nn_m ; j++ ) 
      {
	  if(fabs(wavedataM[j])>=pow(10,9))
	  {
	      if(j==0) wavedataM[j] = 0.0;
	      else wavedataM[j] = wavedataM[j-1];
	  }
       }
      _rmean_(wavedataM,nn_m);

      if(flo != 0.0 && fhi != 0.0)
      {
//	  xapiir( wavedataM,nn_m,poro, trbndw, a, iord,type, flo, fhi, ts, passes);
	  xapiir( wavedataC,nn_c,poro, trbndw, a, iord,type, flo, fhi, ts, passes);
      }
      for( j = 0 ; j < nn_m ; j++ )
      {
	  tmp[j] = wavedataM[j];
      }
      for( j = 0 ; j < nn_c ; j++ ) 
      {
	  cont[j] = wavedataC[j];
      }
  ccout = ccFreq_(tmp, nn_m, cont, nn_c);
      for ( j = nn_m ; j < n-nn_m+1 ; j++) 
      {
	  output[j-nn_m] = ccout[j];
	  if(fabs(output[j-nn_m])>1|| isnan(output[j-nn_m]) != 0 ) 
	  {
	      output[j-nn_m]=0;
	  }
      }
      hd_o = sachdr(ts,n_o,tb);
      write_sac(outFile,hd_o,output);
//  fscanf(fp,"%s",outfinal);
//  Stack_shift(output,btime,staNum_m,n_o,ts,outfinal);
//  fclose(fp3);
  free(tmp);
  free(cont);
  free(ccout);
  free(output);
  return 0;
}

/*void	_taper_(double *aa, int n, double taper)
{
  int i, m;
  double tt, pi1;
  m = taper*n;
  pi1 = 3.1415926/m;
  for (i=0; i<m; i++) {
    tt = 0.5*(1.-cos(i*pi1));
    aa[i] *= tt;
    aa[n-i-1] *= tt;
  }
}
*/
void _rmean_(float *aa,int n)
{
    int i;
    float mean=0;
    for(i=0;i<n;i++)
    {
	mean+=aa[i];

    }
    mean=mean/(float)n;
    for(i=0;i<n;i++)
    {
	aa[i]=aa[i]-mean;
    }

}

