
int lenactual (char *s1, int s1_len) {

/* Local variables */
    int i;

/* ****************************************************************** */

    for (i = 1; i <= 256; ++i)
        if (s1[i - 1] == ' ') return(i-1);
    return(256);
}



void fstrcat_ (char *s1, char *s2, char *scomb, int s1_len, int s2_len, 
	int scomb_len) {

/* Local variables */
    int i, j;
    int ls1, ls2;

/* ****************************************************************** */

    ls1 = lenactual (s1, s1_len);
    ls2 = lenactual (s2, s2_len);

    for (i = 1; i <= ls1; ++i) 
	scomb[i - 1] = s1[i - 1];
    for (i = ls1 + 1; i <= ls1+ls2; ++i) {
	j = i - ls1;
	scomb[i - 1] = s2[j - 1];
    }
    for (i = ls1 + ls2 + 1; i <= 256; ++i) 
	scomb[i - 1] = ' ';
} 

