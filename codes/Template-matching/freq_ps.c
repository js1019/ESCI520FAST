#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include "comThreshold.c"
//#include "median.c"
#include "sac.h"
#include "sacio.c"
#include "MaxValue.c"
#define USAGE "%s [-p pcc_file][-s scc_file][-m pdt_file] [-n sdt_file] [-q tcc ][-o outFile] -h help\n"
#define test 1
static  char    *pccFile = NULL;
static  char    *sccFile = NULL;
static  char    *outFile = NULL;
void _rmean_(float *aa,int n);
float maxf_(float *aa,int n,int *s);
int min_(int *aa,int n);
int main(int argc, char *argv[]) 
{

   printf("start\n");
   extern   int getopt();
   extern    char *optarg;
   extern    int optind;
   char     sac[300];
   char     *PdtFile=NULL;
   char     *SdtFile = NULL;
   float    *data,**wfccp,**wfccs,**stack,**grid,**shift;
   float    *maxSta;
   float    b[100],ts = 0.0,delta=0.0;
   float    threshold=0.0;
   float    outT=0.0;
   float     tstack;
   int   **shiftNP,**shiftNS;
   int   *shNmin;
   int      *eve;
   int 	    c=0,i=0,j=0,k=0,l=0,staNum=0,npts=0,npts0=0,evet=0;
   int      shNm=0,shNp=0,shNs=0;
   SACHEAD  hd;
  //BP parameters
  /* input parameters */
  while( (c=getopt( argc, argv, "p:s:m:n:q:o:h" )) != (-1) )
  {
      switch( c ) 
      {
	  case 'p':
	      pccFile = optarg;
	      break;
	  case 's':
	      sccFile = optarg;
	      break;
	  case 'm':
	      PdtFile = optarg;
	      break;
	  case 'n':
	      SdtFile = optarg;
	      break;
	  case 'q':
	      outT= atof(optarg);
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
  if(argc == 1 || pccFile == NULL ) 
  {
      fprintf(stderr, USAGE, argv[0] );
      exit(1);
  }
//  printf("%s  %s  %s\n",dataFile1,dataFile2,dtFile);
  FILE *fp1;
  fp1=fopen(pccFile,"r");
  int ch;
  ch=0;
  while((ch = fgetc(fp1))!=EOF)
  {
      if(ch == '\n') staNum++;
  }
  printf("staNum=%d\n",staNum);
  fclose(fp1);
  wfccp = (float**)malloc(sizeof(float*) * staNum);
  wfccs = (float**)malloc(sizeof(float*) * staNum);
  FILE *fp2;
  fp2=fopen(pccFile,"r");
  for( i = 0 ; i < staNum ; i++ )
  {
      b[i] = 0.0;
      fscanf(fp2,"%s",sac);
      //printf("%s\n",sac);
      if(( data = read_sac(sac,&hd) ) == NULL ){printf("error");}
      delta = hd.delta;
       npts = hd.npts;
       b[i] = hd.b;
      ts = hd.delta;
//      printf("%d\n",npts);
      wfccp[i] = (float *)malloc(npts*sizeof(float));
      memset(wfccp[i],0x0, sizeof(float) * npts);
      for( j = 0 ; j < npts ; j++)
      {
	  wfccp[i][j] = data[j];
      }
      free(data);
      data = NULL;
  }
  npts0 = npts;
  fp2=fopen(sccFile,"r");
  for( i = 0 ; i < staNum ; i++ )
  {
      fscanf(fp2,"%s",sac);
      //printf("%s\n",sac);
      if(( data = read_sac(sac,&hd) ) == NULL ){printf("error");}
      delta = hd.delta;
      npts = hd.npts;
      ts = hd.delta;
  //    printf("%d\n",npts);
      wfccs[i] = (float *)malloc(npts*sizeof(float));
      memset(wfccs[i],0x0, sizeof(float) * npts);
      for( j = 0 ; j < npts ; j++)
      {
	  wfccs[i][j] = data[j];
      }
      free(data);
      data = NULL;
  }
  fclose(fp2);
  if(npts0>npts) npts0=npts;
  int m = 0;
  printf("m=%d",m);
  FILE *fp;
  fp = fopen(PdtFile,"r");
  ch = 0;
  while(( ch = fgetc(fp)) != EOF )
  {
      if( ch == '\n' ) m++;
  }
  fclose(fp);
  printf("m=%d\n",m);
  stack = (float**)malloc(sizeof(float*)*m);
  grid = (float**)malloc(sizeof(float*)*m);
  shift = (float**)malloc(sizeof(float*)*m);
  maxSta = (float*)malloc(sizeof(float*)*npts0);
  shiftNP = (int**)malloc(sizeof(int*)*m);
  shiftNS = (int**)malloc(sizeof(int*)*m);
  shNmin = (int*)malloc(sizeof(int)*m);
  memset(shNmin,0x0, sizeof(int) * m);
  memset(maxSta,0x0, sizeof(float) * npts0);
  for(i = 0 ; i < m ; i ++)
  {
      grid[i] = (float*)malloc(sizeof(float)*3);
      shift[i] = (float*)malloc(sizeof(float)*staNum);
      shiftNP[i] = (int*)malloc(sizeof(int)*staNum);
      shiftNS[i] = (int*)malloc(sizeof(int)*staNum);
      grid[i][0] = 0.0; grid[i][1] = 0.0;  grid[i][2] = 0.0;
      memset(shiftNP[i],0x0, sizeof(int) * staNum);
      memset(shiftNS[i],0x0, sizeof(int) * staNum);
      memset(shift[i],0x0, sizeof(int) * staNum);
  }
  FILE  *fp3;
  fp3=fopen(PdtFile,"r");
  for (i=0;i<m;i++){
      fscanf(fp3,"%f\t%f\t%f\t",&grid[i][0],&grid[i][1],&grid[i][2]);
//      printf("%f\t%f\t%f\t\n",grid[i][0],grid[i][1],grid[i][2]);
      for(j=0;j<staNum;j++)
      {
	  fscanf(fp3,"%f\t",&shift[i][j]);
	  if(shift[i][j]>0)
	  {
	      shiftNP[i][j]=(int)(0.5+((shift[i][j])/ts));
	  }
	  else
	  {
	      shiftNP[i][j]=(int)((shift[i][j])/ts);
	  }
      }
  }
  fp3=fopen(SdtFile,"r");
  for (i=0;i<m;i++){
      grid[i][0]=0.0;
      grid[i][1]=0.0;
      grid[i][2]=0.0;
      fscanf(fp3,"%f\t%f\t%f\t",&grid[i][0],&grid[i][1],&grid[i][2]);
      for(j=0;j<staNum;j++)
      {
	  shift[i][j]=0.0;
	  fscanf(fp3,"%f\t",&shift[i][j]);
	  if(shift[i][j]>0)
	  {
	      shiftNS[i][j]=(int)(0.5+(shift[i][j])/ts);
	  }
	  else
	  {
	      shiftNS[i][j]=(int)((shift[i][j])/ts);
	  }
      }
  }
  fclose(fp3);
  for(i=0;i<m;i++)
  {
    stack[i]=(float*)malloc(npts0*sizeof(float));
    memset(stack[i],0x0, sizeof(float) * npts0);
  }
  eve = (int *)malloc(npts0*sizeof(int));
  memset(eve,0x0, sizeof(int) * npts0);
  printf("m=%d\n",m);
  #pragma omp parallel for private(i,j,k,shNp,shNs,shNm) shared(npts0,staNum,m) num_threads(5)
  for(i = 0;i < m;i++)
  {
      shNmin[i]=min_(shiftNP[i],staNum);
      shNm = shNmin[i];
      for(k = 0;k < npts0;k++)
      {
	  for(j=0;j<staNum;j++)
	  {
	      shNp = shiftNP[i][j]-shNm;
	      shNs = shiftNS[i][j]-shNm;
	      if(k+shNp>=0&&k+shNs<npts0)
	      {
		  stack[i][k] += wfccp[j][k+shNp];
		  stack[i][k] += wfccs[j][k+shNs];
	      }
	  }
      }
    }
  printf("%.4f\n",ts);
#pragma omp parallel for private(i,k,tstack) num_threads(2)
  for( i = 0 ; i < m ; i ++ )
  {
      for ( k  =0 ; k < npts0 ; k++ )
      {
	  tstack = stack[i][k] / (2*(float)staNum);
	  stack[i][k] = tstack;
      }
  }
//#pragma omp parallel for private(i,k,threshold) shared(npts0,evet) num_threads(2)
  for( i = 0 ; i < m ; i ++ )
  {
      threshold = comThreshold(stack[i],npts0);
     //printf("%.4f\n",threshold);
      for ( k  =0 ; k < npts0 ; k++ )
      {
	  evet = eve[k];
	  if( (stack[i][k] >= threshold)&&(evet ==0 ) )
	  {   
	      eve[k] = k;
	      //printf("stack=%f threshold=%f\n",stack[i][k],threshold);

	  }
      }
  }

  int s=0;
  FILE *fp4;
  fp4=fopen(outFile,"wb");
#pragma omp parallel for private(i,k,s) shared(ts,npts0,m) num_threads(2)
  for(i = 0;i < npts0;i++)
  {
      maxSta[i] = 0.0;
      if(eve[i]!=0)
      {
      for(k = 0 ; k < m ; k++)
      {
	  if(maxSta[i] <= stack[k][i])
	  {
	      maxSta[i] = stack[k][i];
	      s=k;
	  }
      }
      fprintf(fp4,"%-.4f  %-.6f   %-.4f  %-.4f  %-.4f\n",i*ts,maxSta[i],grid[s][0],grid[s][1],grid[s][2]);
      }
  }
  fclose(fp4);
/*float maxAll;
  maxAll = maxf_(maxSta,npts,&s);
//  s=(int)(outT/ts);
  for(i=0;i<m;i++)
  {
      printf("%-.4f  %-.4f  %-.4f  %-.4f  %-.4f %d %d\n",s*ts,stack[i][s],grid[i][0],grid[i][1],grid[i][2],shiftNS[i][0],shiftNP[i][0]);
  }
  */
  for (i=0;i<staNum;i++)
  {
      free(wfccp[i]);
  }
  for (i=0;i<m;i++)
  {
      free(stack[i]);
  }
  free(eve);
  printf("OK,END\n");
  return 0;
}
void _rmean_(float *aa,int n)
{
    int i;
    float mean=0;
    for(i=0;i<n;i++)
    {
	mean+=aa[i];

    }
    mean=mean/n;
    for(i=0;i<n;i++)
    {
	aa[i]=aa[i]-mean;
    }

}
int min_(int *aa,int n)
{
    int minA;
    int i;
    minA=aa[0];
    for(i=0;i<n;i++)
    {
	if(minA>aa[i])
	{
	    minA=aa[i];
	}
    }
    return minA;
}
float maxf_(float *aa,int n,int *s)
{
    float maxA;
    int i;
    maxA=aa[0];
    *s=0;
    for(i=0;i<n;i++)
    {
	if(maxA<aa[i])
	{
	    maxA=aa[i];
	    *s = i;
	}
    }
    return maxA;
}

