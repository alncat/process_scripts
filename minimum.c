#include<stdio.h>
#include<stdlib.h>
#include<string.h>

typedef struct r_value
{
    int row;
    double crys;
    double free;
    int number;
    char method[7];
} r_value;

int comp(const void * a, const void * b)
{
    double acrys = ((r_value*)a)->crys;
    double bcrys = ((r_value*)b)->crys;

    return acrys < bcrys ? -1 : acrys > bcrys ? 1 : 0;
}
int main(int argc, char *argv[])
{
    int tlsi=0, nmi=0;
    r_value tls[150], nm[150], tmpr;
    char buf[200], tls_test[]="tlsref", nm_test[]="nmref";
    FILE *ifp; 

    ifp = fopen(argv[1], "r");
    if(ifp == NULL)
    { 
        fprintf(stderr, "cannot open file %s\n", argv[1]);
        exit(EXIT_FAILURE);
    }
 //   ofp = fopen("out.list_sort", "a");

    
    while(fgets(buf, 200, ifp) != NULL)
    {
        if(sscanf(buf, "%d %lf %lf %d %s", &(tmpr.row), &(tmpr.crys), &(tmpr.free), &(tmpr.number), &(tmpr.method)) != 5)
        {
            fprintf(stderr, "input value error!");
            exit(EXIT_FAILURE);
        }
        if(strcmp(tmpr.method, tls_test) == 0 )
        {
            tls[tlsi]=tmpr;
            tlsi++;
        }
        else if(strcmp(tmpr.method, nm_test) == 0 )
        {
            nm[nmi]=tmpr;
            nmi++;
        }
    }
    fclose(ifp);
    
    if(tlsi != 0)
        qsort(tls, tlsi, sizeof(r_value), comp);
    if(nmi != 0)
        qsort(nm, nmi, sizeof(r_value), comp);


    printf("%3d %.5lf %.5lf %3d %s\n", tls[0].row, tls[0].crys, tls[0].free, tls[0].number, tls[0].method);
   printf("%3d %.5lf %.5lf %3d %s\n", nm[0].row, nm[0].crys, nm[0].free, nm[0].number, nm[0].method);

    return EXIT_SUCCESS;

}


