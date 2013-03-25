#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#define BUFFER_SIZE 255


int main(int arc, char *argv[]){
    FILE *ifp, *ofp;
    int row, i;
    double r, r_f;
    char *mode_1 = "r", *tmp, *num;
    char buf[BUFFER_SIZE];
    char outputfilename[] = "out.list", inputfilename[20], suf[10];
    char r_test[]="R VALUE       ",r_f_test[]="FREE R VALUE  ";

for(i = 1; i <= 150; i++){
    strcpy(inputfilename, "out.pdb"); 
    sprintf(suf, "%d", i);
    strcat(inputfilename, suf);
    
    ifp = fopen(inputfilename, mode_1);
    
    if (ifp == NULL){
      fprintf(stderr, "Can't open input file %s!\n", inputfilename);
      exit(EXIT_FAILURE);
}

    ofp = fopen(outputfilename, "a");
       
    if (ofp == NULL) {
      fprintf(stderr, "Can't create output file %s!\n", outputfilename);
      exit(EXIT_FAILURE);
}    
/*    while(fscanf(ifp, "%i %lf %lf\n", &row, &r, &r_f) != EOF){
      fprintf(ofp, "%2i %lf %lf\n", row, r, r_f);
}*/
    while(fgets(buf, BUFFER_SIZE, ifp) != NULL ){
      tmp = (char *)malloc(15*sizeof( *tmp ));
      strncpy(tmp, buf+13, 14);
      row = strcmp(tmp, r_test);
      if( row == 0 ){ 
        num = (char *)malloc(10 * sizeof( *num ));
        strncpy(num, buf+49, 7);
        if(sscanf(num, "%lf", &r)!=1)
        {
            fprintf(stderr, "r value input err\n");
            exit(EXIT_FAILURE); 
        }
}
      else if( strcmp(tmp, r_f_test) == 0 ){
        strncpy(num, buf+49, 7);
        if(sscanf(num, "%lf", &r_f)!=1)
        {
            fprintf(stderr, "r_f value input err\n");
            exit(EXIT_FAILURE);
        }
        free(num);
        free(tmp);
        break;
}
}     
    fprintf(ofp, "%3d %.5lf %.5lf %s %s\n", i, r, r_f, argv[1], argv[2]); 
    fclose(ifp);
    fclose(ofp);
}
    return EXIT_SUCCESS;
}
