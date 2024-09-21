# LDPC_code
Project of Error control coding - LDPC code

##Consider the (1023, 781) Low-Density Parity-Check (LDCP) code with  
block length $n = 1023$, dimension $k = 781$, row weight $\rho = 32$, and column weight $\gamma = 32$.

##Block Diagram
![Block Diagram](https://github.com/H-Y-Hs/LDPC_code/blob/main/Block_Diagram.jpg?raw=true)

###Notes
1. Use the recursion $u_{l+6} = u_{l+1} ⊕ u_l$, for $l ≥ 0$ with the initial conditions $u_0 =1,\ u_1 =u_2 =u_3 =u_4 =u_5 =0$ to generate $k=781$ information bits.  
2. The generated sequence is 100000100001 . . . with period 63.

####Pseudo code
```c
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

main ()
{
  ...
  long *idum;
  idum = (long *)malloc(sizeof(long));
  *idum = SEED; //SEED must be a negative integer.
  ...
  ...
}
```

#####Use **normal()** to output two independent normal random variables, $n_1$ and $n_2$
<img src="https://github.com/H-Y-Hs/LDPC_code/blob/main/normal.jpg?raw=true" alt="normal()" width="50%">


######Use **ran1()** to generate a random variable uniformly distributed in the interval (0, 1).
```c
ran1(long *idum)
{
  int j;
  long k;
  static long iy=O;
  static long iv[NTAB];
  double temp;
  if (*idum <= 0 || !iy){
    if (-(*idum) < 1) *idum=1;
    else *idum = -(*idum);
    for (j=NTAB+7;j>=0;j--){
      k= (*idum) / IQ;
      *idum=IA*(*idum-k*IQ)-IR*k;
      if (*idum < 0) *idum += IM;
      if (j<NTAB) iv[j]=*idum;
    }
  iy=iv[0];
  }
  k= (*idum) / IQ;
  *idum=IA*(*idum-k*IQ)-IR*k;
  if (*idum < 0) *idum += IM;
  j=iy/NDIV;
  iy=iv[j];
  iv[j] = *idum;
  if ((temp=AM*iy) > RNMX) return RNMX,
  else return temp;
}
```
