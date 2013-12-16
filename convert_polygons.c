#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int main(int argc,char *argv[]){

  FILE *mt_input, *dv_output;
  int i=0,j=0,num_ver=0;
  char buff[1000],poly_file[40];
  double coords[1000];

  if(argc<2){
    printf("Usage: %s polygon_file\n",argv[0]);
    exit(0);
  }
  strcpy(poly_file,argv[1]); //read in polygon filename

  if((mt_input=fopen(poly_file, "r"))==NULL){
    printf("cannot open %s\n", poly_file);
    exit(-1);
  } 
  else printf("opened file %s\n",poly_file); //check that polygon file opens correctly
  
  fgets(buff,1000,mt_input); //read contents to buffer
  
  for(i=0;i<strlen(buff);i++){
    if(buff[i]=='{' || buff[i]==' '){ // coordinates are preceded by '{' or ' '
      if(buff[i+1]=='{'){ // start of the file has a double {{ so skip
        for(j=i;j<strlen(buff);j++){
          buff[j]=buff[j+1];
        }
        i--;
      }else{ //we have a new coordinate so read it from buffer
        coords[i]=atof(&buff[i+1]);
        num_ver++;
      }
    }else{ //discard anything else
      for(j=i;j<strlen(buff);j++){
        buff[j]=buff[j+1];
      }
      i--;
    }
  }

  num_ver/=2; //number of vertices is half that of the number of coordinates
  fclose(mt_input);

  for(j=0;j<i;j++){ 
    printf("%f %f\n",coords[j],coords[j+1]);
    j++;
  }
  printf("There are %i vertices in the polygon\n",num_ver);

  dv_output=fopen("data_view.polygon","w");
  printf("Opening output file 'dataview.polygon'\n");
  
  fprintf(dv_output,"%d\n",num_ver);
  for(i=0;i<num_ver*2;i++){
    fprintf(dv_output,"%f %f\n",coords[i],coords[i+1]);
    i++;
  }
  fclose(dv_output);
  printf("Conversion successful\n");

  return 0;
}
