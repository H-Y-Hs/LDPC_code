#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>

#define IA 16807
#define IM 2147483647
#define AM (1.0 / IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1 + (IM - 1) / NTAB)
#define EPS 1.2e-7
#define RNMX (1.0 - EPS)

long static *idum;
double ran1()
{
    int j;
    long k;
    static long iy = 0;
    static long iv[NTAB];
    double temp;
    if (*idum <= 0 || !iy)
    {
        if (-(*idum) < 1)
            *idum = 1;
        else
            *idum = -(*idum);
        for (j = NTAB + 7; j >= 0; j--)
        {
            k = (*idum) / IQ;
            *idum = IA * (*idum - k * IQ) - IR * k;
            if (*idum < 0)
                *idum += IM;
            if (j < NTAB)
                iv[j] = *idum;
        }
        iy = iv[0];
    }
    k = (*idum) / IQ;
    *idum = IA * (*idum - k * IQ) - IR * k;
    if (*idum < 0)
        *idum += IM;
    j = iy / NDIV;
    iy = iv[j];
    iv[j] = *idum;
    if ((temp = AM * iy) > RNMX)
        return RNMX;
    else
        return temp;
}
// void normal
void normal(double *n, double sigma)
{
    double x1, x2, s;
    long idum;
    do
    {
        x1 = ran1();
        x2 = ran1();
        x1 = 2.0 * x1 - 1.0;
        x2 = 2.0 * x2 - 1.0;
        s = pow(x1, 2) + pow(x2, 2);
    } while (s >= 1.0);
    n[0] = sigma * x1 * pow((-2) * log(s) / s, 0.5);
    n[1] = sigma * x2 * pow((-2) * log(s) / s, 0.5);
}

int main()
{
    double SNR;
    for (int iteration = 17; iteration < 19; iteration++)   
    {
        SNR = 0.2 * iteration;
        int N = 781;
        int Seed = -2023;
        int n = 1023;
        int k = 781;
        int iteration_block = 1;
        int iteration = 100; 
        // read file
        FILE *G_in;
        FILE *H_in;
        G_in = fopen("ldpc_G_1023.txt", "r");
        H_in = fopen("ldpc_H_1023.txt", "r");        

        int **G = (int **)malloc(k * sizeof(int *));
        for (int i = 0; i < k; i++)
        {
            G[i] = (int *)malloc(n * sizeof(int));
        }
        for (int i = 0; i < k; i++)
        {
            for (int j = 0; j < n; j++)
            {

                fscanf(G_in, "%1d", &G[i][j]);
            }
        }
        int **H = (int **)malloc(2046 * sizeof(int *));
        for (int i = 0; i < 2046; i++)
        {
            H[i] = (int *)malloc(32 * sizeof(int));
        }
        for (int i = 0; i < 2046; i++)
        {
            for (int j = 0; j < 32; j++)
            {
                fscanf(H_in, "%d", &H[i][j]);
                H[i][j]--;
            }
        }

        int t = 0;
        int total_error = 0;
        int block_error = 0;
        // seed
        idum = (long *)malloc(sizeof(long));
        *idum = Seed;
        
        while (1)
        {
            // m sequence
            int *l = (int *)malloc(N * sizeof(int));
            int m = 6;
            l[0] = 1;
            for (int i = 1; i < 6; i++)
            {
                l[i] = 0;
            }
            for (int i = 6; i < N; i++)
            {
                l[i] = (l[i - m] + l[i - m + 1]) % 2;
            } 
            int *U = (int *)malloc(N * sizeof(int));
            for (int i = 0; i < N; i++)
            {
                U[i] = 0;
            }
            for (int i = 0; i < N; i++)
            {
                U[i] = l[(0 + (i * 1)) % 63];
            }

            //Encode
            int *c = (int *)malloc(n * sizeof(int));
            for (int i = 0; i < n; i++)
            {
                c[i] = 0;
            }
            // c = U * G
            for (int i = 0; i < k; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    c[j] += U[i] * G[i][j];
                    c[j] = c[j] % 2;
                }
            }            
            int *x = (int *)malloc(n * sizeof(int));
            for (int i = 0; i < n; i++)
            {
                x[i] = 0;
            }
            // x = (-1) ^ c
            for (int i = 0; i < n; i++)
            {
                x[i] = (c[i] == 0) ? 1 : -1;
            }

            // sigma, noise
            double sigma;
            double R = 0.76344086; // 781/1023
            sigma = sqrt(1 / (2 * R * pow(10, (SNR / 10))));
            double noise[2] = {0};

            // initial y
            double *y = (double *)malloc(n * sizeof(double));
            for (int i = 0; i < n; i++)
            {
                y[i] = 0.0;
            }
            // y = x + n
            for (int i = 0; i < n; i++)
            {
                normal(noise, sigma);
                // printf("noise1: \n%lf\n ", noise[0]);
                // printf("noise2: \n%lf\n\n ", noise[1]);
                y[i] = x[i] + noise[1];
                // y[i][1] = x[i][1] + noise[1];
            } 

            //Decode  
            // initial L
            double *L = (double *)malloc(n * sizeof(double));
            for (int i = 0; i < n; i++)
            {
                L[i] = 0.0;
            }
            //L_l公式
            double L_c = 2 / pow(sigma, 2);
            for (int i = 0; i < n; i++)
            {
                L[i] = L_c * y[i];
            }
            double **q = (double **)malloc(n * sizeof(double *));
            for (int i = 0; i < n; i++)
            {
                q[i] = (double *)malloc(32 * sizeof(double));
            }
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < 32; j++)
                {
                    q[i][j] = 0;
                }
            }
            double **r = (double **)malloc(n * sizeof(double *));
            for (int i = 0; i < n; i++)
            {
                r[i] = (double *)malloc(32 * sizeof(double));
            }
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < 32; j++)
                {
                    r[i][j] = 0;
                }
            }
            double *q_l = (double *)malloc(n * sizeof(double)); 
            for (int i = 0; i < n; i++)
            {
                q_l[i] = 0;
            }
            int *x_hat = (int *)malloc(n * sizeof(int));
            for (int i = 0; i < n; i++)
            {
                x_hat[i] = 0;
            }
            int x_check = 0;
            int *check = (int *)malloc(n * sizeof(int));
            for (int i = 0; i < n; i++)
            {
                check[i] = 0;
            }
            int *temp = (int *)malloc(n * sizeof(int));
            for (int i = 0; i < n; i++)
            {
                temp[i] = 0;
            }
            int *temp2 = (int *)malloc(n * sizeof(int));
            for (int i = 0; i < n; i++)
            {
                temp2[i] = 0;
            }
            int *temp3 = (int *)malloc(n * sizeof(int));
            for (int i = 0; i < n; i++)
            {
                temp3[i] = 0;
            }
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < 32; j++)
                {
                    q[i][j] = L[i];
                }
            }            

            int s = 0;
            while (s < iteration)
            {
                for (int i = 0; i < n; i++)
                {
                    for (int j = 0; j < 32; j++)
                    {
                        //Bottom up
                        r[i][j] = 2 * atanh(tanh(q[H[i][(j + 1) % 32]][temp[H[i][(j + 1) % 32]]] / 2) * tanh(q[H[i][(j + 2) % 32]][temp[H[i][(j + 2) % 32]]] / 2));
                        //r[i][j] = log((1 + (tanh(q[H[i][(j + 1) % 32]][temp[H[i][(j + 1) % 32]]] / 2)) * (tanh(q[H[i][(j + 2) % 32]][temp[H[i][(j + 2) % 32]]] / 2)))/(1 - (tanh(q[H[i][(j + 1) % 32]][temp[H[i][(j + 1) % 32]]] / 2)) * (tanh(q[H[i][(j + 2) % 32]][temp[H[i][(j + 2) % 32]]] / 2))));
                        for (int k = 3; k < 32; k++)
                        {
                            r[i][j] = 2 * atanh(tanh((r[i][j]) / 2) * tanh((q[H[i][(j + k) % 32]][temp[H[i][(j + k) % 32]]]) / 2));
                            //r[i][j] = log((1 + (tanh((r[i][j]) / 2) * tanh((q[H[i][(j + k) % 32]][temp[H[i][(j + k) % 32]]]) / 2)))/(1 - (tanh((r[i][j]) / 2) * tanh((q[H[i][(j + k) % 32]][temp[H[i][(j + k) % 32]]]) / 2))));
                        }
                    }                 
                    for (int j = 0; j < 32; j++)
                    {
                        temp[H[i][j]] += 1;
                    }
                }
                for (int i = 0; i < n; i++)
                {
                    temp[i] = 0;
                }
                //Top down
                for (int i = 0; i < n; i++)
                {
                    for (int j = 0; j < 32; j++)
                    {
                        q[i][j] = L[i];
                        for (int k = 1; k < 32; k++)
                        {
                            q[i][j] += r[H[i + 1023][(j + k) % 32]][temp2[H[i + 1023][(j + k) % 32]]];
                        }
                    }
                    for (int j = 0; j < 32; j++)
                        temp2[H[i + 1023][j]] += 1;
                }
                for (int i = 0; i < n; i++)
                {
                    temp2[i] = 0;
                }

                for (int i = 0; i < n; i++)
                {
                    q_l[i] = L[i];
                    for (int j = 0; j < 32; j++)
                    {
                        q_l[i] = q_l[i] + r[H[i + 1023][j]][temp3[H[i + 1023][j]]];
                        temp3[H[i + 1023][j]] += 1;
                    }
                }
                for (int i = 0; i < n; i++)
                {
                    temp3[i] = 0;
                }

                for (int i = 0; i < n; i++)
                {
                    if (q_l[i] >= 0)
                    {
                        x_hat[i] = 0;
                    }
                    else
                    {
                        x_hat[i] = 1;
                    }
                }

                for (int i = 0; i < n; i++)
                {
                    for (int j = 0; j < 32; j++)
                    {
                        x_check += x_hat[H[i][j]];
                    }
                    check[i] = x_check % 2;
                    x_check = 0;
                }
                int temp4 = 1;
                for (int i = 0; i < n; i++)
                {
                    if (check[i] == 1)
                    {
                        temp4 = 0;
                        break;
                    }
                }
                if (temp4 == 1)
                    break;

                s++;
            }
            t++;

            int error = 0;
            double BER = 0.0;
            for (int i = 0; i < N; i++)
            {
                if (U[i] - x_hat[i] != 0)
                    error += 1;
            }
            //printf("Error: %d\n", error);
            BER = error * 1.0 / (N);
            //printf("Bit Error Rate: %f\n", BER);
            //printf("\n");
            total_error += error;
            if (error != 0)
            {
                block_error += 1;
				//printf("Error : %d\n", block_error);
				//printf("iteration_block: %d \n\n", t);				
            }
            
            free(U);
            free(c);
            free(x);
            free(y);
            free(L);
            free(q);
            free(r);
            free(q_l);
            free(x_hat);
            free(check);
            free(temp);
            free(temp2);
            free(temp3);
            
            if (block_error > 99)
            {
                break;
            }
            iteration_block++;
        }

        double BER = 0.0;
        double BLER = 0.0;
        printf("SNR = %lf\n", SNR);
        printf("Error: %d\n", total_error);
        printf("Decoded Bits: %d\n", (N * iteration_block));
        BER = total_error * 1.0 / (N * iteration_block);
        printf("Bit Error Rate: %f\n", BER);
        BLER = block_error * 1.0 / iteration_block;
        printf("Block: %d\n", iteration_block);
        printf("Block Error Rate: %f\n\n", BLER);

        free(G);
        free(H);
        free(idum);
        fclose(G_in);
        fclose(H_in);
    }
    return 0;
}