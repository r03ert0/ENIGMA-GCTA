/*
	mygcta
	makes it nicer to use gcta
	roberto toro, v2 30 July 2013
	roberto toro, v1 12 July 2012
*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>

int main(int argc, char *argv[])
{
	int		i=1;
	int		np=0,nq=0,nc=0;
	char	cmd[4096];
	char	cmd2[4096];
	int		rnd;
	
	struct	timeval	t;
	gettimeofday(&t, 0);
	srand(t.tv_usec);
	rnd=rand();
	
	strcpy(cmd,argv[i++]);
	
	while(i<argc)
	{
		if(strcmp(argv[i],"--pheno")==0)
		{
			if(np==0)
			{
				sprintf(cmd2,"sort %s > tmp_pheno%i",argv[++i],rnd);
				system(cmd2);
			}
			else
			{
				sprintf(cmd2,"sort %s >tmp%i;join -o $(%s%i)$(%s%i) tmp_pheno%i tmp%i > tmp2%i;rm tmp_pheno%i;rm tmp%i;mv tmp2%i tmp_pheno%i",
								argv[++i],
								rnd,
								"awk 'END{for(i=1;i<=NF;i++) printf(\"1.%i,\",i);}' tmp_pheno",rnd,
								"awk 'END{for(i=3;i<NF;i++) printf(\"2.%i,\",i);printf(\"2.%i\",i);}' tmp",rnd,
								rnd,rnd,rnd,rnd,rnd,rnd,rnd);
				system(cmd2);
			}
			printf("COMMAND 0: [%s]\n",cmd2);
			np++;
		}
		else
		if(strcmp(argv[i],"--qcovar")==0)
		{
			if(nq==0)
			{
				sprintf(cmd2,"sort %s > tmp_qcovar%i",argv[++i],rnd);
				system(cmd2);
			}
			else
			{
				sprintf(cmd2,"sort %s >tmp%i;join -o $(%s%i)$(%s%i) tmp_qcovar%i tmp%i > tmp2%i;rm tmp_qcovar%i;rm tmp%i;mv tmp2%i tmp_qcovar%i",
								argv[++i],
								rnd,
								"awk 'END{for(i=1;i<=NF;i++) printf(\"1.%i,\",i);}' tmp_qcovar",rnd,
								"awk 'END{for(i=3;i<NF;i++) printf(\"2.%i,\",i);printf(\"2.%i\",i);}' tmp",rnd,
								rnd,rnd,rnd,rnd,rnd,rnd,rnd);
				system(cmd2);
			}
			printf("COMMAND 1: [%s]\n",cmd2);
			nq++;
		}
		else
		if(strcmp(argv[i],"--covar")==0)
		{
			if(nc==0)
			{
				sprintf(cmd2,"sort %s > tmp_covar%i",argv[++i],rnd);
				system(cmd2);
			}
			else
			{
				sprintf(cmd2,"sort %s >tmp%i;join -o $(%s%i)$(%s%i) tmp_covar%i tmp%i > tmp2%i;rm tmp_covar%i;rm tmp%i;mv tmp2%i tmp_covar%i",
								argv[++i],
								rnd,
								"awk 'END{for(i=1;i<=NF;i++) printf(\"1.%i,\",i);}' tmp_covar",rnd,
								"awk 'END{for(i=3;i<NF;i++) printf(\"2.%i,\",i);printf(\"2.%i\",i);}' tmp",rnd,
								rnd,rnd,rnd,rnd,rnd,rnd,rnd);
				system(cmd2);
			}
			printf("COMMAND 2: [%s]\n",cmd2);
			nc++;
		}
		else
			sprintf(cmd+strlen(cmd)," %s",argv[i]);
		i++;
	}
	
	printf("npheno=%i, nqcovar=%i, ncovar=%i\n",np,nq,nc);
	
	if(np)
	{
		sprintf(cmd+strlen(cmd)," --pheno tmp_pheno%i",rnd);
	}
	if(nq)
	{
		sprintf(cmd+strlen(cmd)," --qcovar tmp_qcovar%i",rnd);
	}
	if(nc)
	{
		sprintf(cmd+strlen(cmd)," --covar tmp_covar%i",rnd);
	}
	printf("%s\n",cmd);
	system(cmd);

	if(np)
	{
		sprintf(cmd,"rm tmp_pheno%i",rnd);
		system(cmd);
	}
	if(nq)
	{
		sprintf(cmd,"rm tmp_qcovar%i",rnd);
		system(cmd);
	}
	if(nc)
	{
		sprintf(cmd,"rm tmp_covar%i",rnd);
		system(cmd);
	}

	return 0;
}

