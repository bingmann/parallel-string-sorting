#include "copy.h"

void vswap(char **a, char **b, int n)
{
	while (n-- > 0)
	{
		char *t = *a;
		*a++ = *b;
		*b++ = t;
	}
}

char **med3(char **a, char **b, char **c, int d)
{
	int i, j, k;

	if ((i = P2C(a, d)) == (j = P2C(b, d)))
		return a;

	if ((k = P2C(c, d)) == i || k == j)
		return c;
	return i < j ? (j < k ? b : (i < k ? c : a)) : (j > k ? b : (i < k ? a : c));
}

char **med3s(char **a, char **b, char **c)
{
	string s, t, u; s = *a; t = *b; u = *c;
	return s < t ? (t < u ? b : (s < u ? c : a)) : (t > u ? b : (s < u ? a : c));
}

void cr()
{
	printf("\n"); fflush(stdout);
}

void dot()
{
	printf(". ");
}

void say(string s)
{
	printf("%s ", s);
	fflush(stdout);
}

void sayln(string s)
{
	printf("%s\n", s);
	fflush(stdout);
}

void isay(int i)
{
	printf("%d ", i);
	fflush(stdout);
}

void dsay(double d, int n)
{
	switch (n)
	{
		case 0:
			printf("%.0f ", d);
			break;
		case 1:
			printf("%.1f ", d);
			break;
		case 2:
			printf("%.2f ", d);
			break;
		case 3:
			printf("%.3f ", d);
			break;
		case 5:
			printf("%.5f ", d);
			break;
		case 6:
			printf("%.6f ", d);
				  break;
		case 4: default:
				  printf("%.2f ", d);
				  break;
	}
}

void esay(string s)
{
	fprintf(ERRORLOG, "%s ", s);
	fflush(ERRORLOG); say(s);
}

void br()
{
	while (getchar() != '\n') ;
}

void brp(string s)
{
	esay(s); br();
}

string allof(string s) {if (*s==0) --s; while (*s) --s; return(++s);}

int tokenize(string s, string *tok, char d, char dd)
{
	int n;
	string t;
	char c; n = 1; t = s;

	while ((c = *s) != 0)
	{
		if (c == d || c == dd)
		{
			++n; *s++ = 0; *tok++ = t; t = s;
		}
		else
			++s;
	}
	*tok = t;
	return(n);
}

string fp(string dir, string fn, string ft)
{
	string s = (string) malloc(100);
	strcpy(s, dir);
	strcat(s, fn); strcat(s, "."); strcat(s, ft);
	return(s);
}

void treset(int nt, int nr)
{
	 int i, j;

	 for (i=0;i<nt;++i)
	 {
		 TR[i].n=nr;
		 for (j=0;j<nr;++j)
			  TR[i].t[j]=0;
	 }
	 TMR=TR+(TIMERPHASE=0);
}

void ton(int i)
{
	 static clock_t t1=0;
	 clock_t t2;

	 t2=clock();
	 TMR->t[REP]+=t2-t1;
	 TMR->n=REP+1;
	 TIMERPHASE=i;
	 TMR=TR+i;
	 t1=t2;
}

timer *tsort(int i)
{
	timer *nu;
	int l, r, n;
	double lt, rt;

	nu = (timer *) calloc(1, sizeof(timer));
	ALLOCATED += sizeof(timer);
	nu->n = n = TR[i].n;

	for (r = 0; r < n; ++r)
		nu->t[r] = TR[i].t[r];

	for (r = 1; r < n; ++r)
	{
		rt = nu->t[r];

		for (l = r; l > 0; --l)
		{
			lt = nu->t[l - 1];
			if (rt < lt) nu->t[l] = lt;
			else
				break;
		}

		nu->t[l] = rt;
	}
	return(nu);
}

double tmean(int i)
{
	 int j;
	 double d;

	 for (d = j = 0; j < TR[i].n; ++j)
		 d += TR[i].t[j];
	 d /= TR[i].n;
	 return(d * MSEC_PER_CLOCK);
}

double tmed(int i)
{
	timer *nu; int m, n; double d;

	n = TR[i].n;
	if (n <= 0)
		return(0);
	else if (n == 1)
		return(TR[i].t[0] * MSEC_PER_CLOCK);

	nu = tsort(i);
	m = n>>1;
	d = n % 2 ? nu->t[m] : ((nu->t[m] + nu->t[--m]) / 2);

	FREE(nu, sizeof(timer));
	return(d * MSEC_PER_CLOCK);
}

void tmin(int nr)
{
	int j, k;
	double d;
	double finaltime[MAXREPS];

	d=1000000000;

	/* add all parts for final time for each run */

	for(k = 0; k < nr; k++)
		finaltime[k] = TR[ti].t[k] + TR[tb].t[k] + TR[tg].t[k] + TR[tc].t[k];

	for (j = 0; j < nr; ++j)
	{
		if (d > finaltime[j])
			d = finaltime[j];
	}

	TMIN = d;
	TNORM = TMIN/(NKEYS * log10(NKEYS));
}

void tminsay()
{
	dsay(TMIN * MSEC_PER_CLOCK, 0);
}

void tnormsay()
{
	dsay(TNORM * MSEC_PER_CLOCK, 6);
}

void isaye(int i,string s) {printf("%s=%d ",s,i); fflush(stdout);}
