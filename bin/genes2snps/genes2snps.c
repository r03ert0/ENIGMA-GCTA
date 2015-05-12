/*
roberto toro, 30 July 2013
make a list of snps from a list of genes
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

int main(int argc, char *argv[])
{
	// argv[1]	my gene set file
	// argv[2]	gene/snp correspondence list file
	
	FILE	*fmy;	// gene set file
	FILE	*fref;	// gene/snp file
	char	s1[512];
	char	s2[512];
	int		found;

	fmy=fopen(argv[1],"r");
	fref=fopen(argv[2],"r");
	
	while(!feof(fmy))
	{
		fgets(s1,512,fmy); // get a gene from my list
		s1[strlen(s1)-1]=(char)0;
		//printf("-----> %s\n",s1);
		found=0;
		do
		{
			fgets(s2,512,fref);	// get a gene from gene/snp list
			s2[strlen(s2)-1]=(char)0;
			//printf("\t\t%s\n",s2);
			if(strcmp(s1,s2)==0)
				found=1;
			do
			{
				fgets(s2,512,fref);
				s2[strlen(s2)-1]=(char)0;
				//printf("\t\t\t\t[%s][%s]:%i\n",s2,"END",strcmp(s2,"END"));
				if(found==1 && strcmp(s2,"END")!=0)
					printf("%s\n",s2);
			}
			while(strcmp(s2,"END")!=0 && !feof(fref));
			fgets(s2,512,fref); // skip empty line
			if(feof(fref))
			{
				printf("WARNING: gene %s not found\n",s1);
				fclose(fref);
				found=1;
				fref=fopen(argv[2],"r");
			}
		}
		while(found==0);
	}
	
	return 0;
}
