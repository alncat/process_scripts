#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#define BUFFER_SIZE 85
int sort(char (*p)[BUFFER_SIZE], int , int );
int printresidue(char (*p)[BUFFER_SIZE], int , FILE * );

int main(int arc, char *argv[])
{
    FILE *ifp, *ofp;
    int seq_no, i, seqmark, middlemark, terminate;
    char recname[7], atom[5], termiatom[5], seqno[6], atommark[] = " C  ", prefix[] = "ATOM  ", termimark[] = " OXT";
    char buf[BUFFER_SIZE], residue[100][BUFFER_SIZE];
    atom[4]='\0';
    terminate=1;

    ifp = fopen(argv[1], "r");

    if(ifp == NULL)
    {
        fprintf(stderr, "Cannot open input file %s!\n", argv[1]);
        exit(EXIT_FAILURE);
    }

    ofp = fopen(argv[2], "a");

    if(ofp == NULL)
    {
        fprintf(stderr, "Cannot create output file %s!\n", argv[2]);
        exit(EXIT_FAILURE);
    }

    while(fgets(buf, BUFFER_SIZE, ifp) != NULL)
    {
        strncpy(recname, buf, 6);
        if(strcmp(recname, prefix) != 0)
            fprintf(ofp, "%s", buf);
        else
            break;
    }
    seqmark = 1, i = 0;
    do
    {
        strncpy(seqno, buf+22, 5);
        if(sscanf(seqno, "%d", &seq_no) != 1)
        {
            fprintf(stderr, "wrong sequence number!\n");
            exit(EXIT_FAILURE);
        }
        if(seq_no == seqmark)
        {
            strcpy(residue[i], buf);
            strncpy(atom, buf+12, 4);
            if(strcmp(atom, atommark) == 0)
            {
                middlemark = i;
            }
            fgets(buf, BUFFER_SIZE, ifp);
            i++;
        }
        else
        {
            strncpy(termiatom, residue[i-1]+12, 4);
            terminate = strcmp(termiatom, termimark);
            if(terminate != 0)
            {
                sort(residue, middlemark-1, i-1);
            }
            printresidue(residue, i-1, ofp);
            seqmark = seq_no;
            i = 0;
        }
            
    } while(strcmp(termiatom, termimark) != 0);
    fprintf(ofp, "%s", buf);
    while(fgets(buf, BUFFER_SIZE, ifp) != NULL)
    {
            fprintf(ofp, "%s", buf);
    }
    fclose(ifp);
    fclose(ofp);
    return EXIT_SUCCESS;
} 

int sort(char (*residue)[BUFFER_SIZE], int middle, int boundary)
{
    int i;
    char tmp[4][BUFFER_SIZE];
   
    for(i = middle; i <= boundary; i++)
    {
        strcpy(tmp[i+3-boundary], residue[i]);
    }

    for(i = middle - 1; i >=  4; i--)
    {
        strcpy(residue[i+4], residue[i]);
    }

    for(i = 4; i <= 7; i++)
    {
       strcpy(residue[i], tmp[i-4]);
    }
    return EXIT_SUCCESS;
}

int printresidue(char (*residue)[BUFFER_SIZE], int boundary, FILE *ofp)
{
    int i;
    for(i = 0; i <= boundary; i++)
        fprintf(ofp, "%s", residue[i]);
    return EXIT_SUCCESS;
}
