    #define MAXDIM 50   
    #include <cstdlib>   
    #include <math.h>
    #include <algorithm>  
    #include <iostream>
    #include <stddef.h> /* defines NULL */
    #include <sys/time.h>
    #include <cstdio>

    using namespace std;

    /*Initial variable declaration*/
    int IA[MAXDIM], N = MAXDIM;
    double AV[MAXDIM], BV[MAXDIM], CV[MAXDIM];
    double OP[MAXDIM][MAXDIM], ID[MAXDIM][MAXDIM];
    double AM[MAXDIM][MAXDIM], BM[MAXDIM][MAXDIM];
    double CM[MAXDIM][MAXDIM], DM[MAXDIM][MAXDIM];

    #ifdef NOPRECISIONLOSS
        double start, finish, sum;
    #else
        float start, finish, sum;
    #endif

    /*Function Prototypes*/
    double walltime_();
    double cputime_();
    void idcheck(int N, double *check, double AV[], double BV[], double ID[][MAXDIM]);
    double trig(int i, int j);

    /* cputime.cc */

    double cputime_()
    {
        return  (double)clock() / (double)CLOCKS_PER_SEC;
    }

    double factor = 1.0e-6;

    /* walltime.cc */
  	double walltime_() 
  	{
  		struct timeval tp;
  		int rtn;
  		double seconds;
  
  		rtn=gettimeofday(&tp, NULL);
  
  		seconds = tp.tv_sec + factor * tp.tv_usec;
  
  		return  seconds ; 
  	}
  
  	/* conrand.c */
  	double conrand(double *seed) 
  	{
  		double a, m;
  		double temp;
  		a = 16807.00000000000;
  		m = 2147483647.000000;
  		temp = a* *seed; 
  		*seed = temp - m * (int) (temp/m);
  		return *seed / m; 
  	}

    int main()
    {
		double check, BOT, TOP, HOLDA, HOLDB, TRACE3, seed, sub, jn0, jn1;
		int jn2;
		cout.precision(17);
		start = cputime_();
		seed = 1;

      /*Loop 10 Series -- Filling Arrays*/

      for (int i = 0; i < N; i++)
      {
        jn1 = conrand(&seed);
        jn2 = 10 * conrand(&seed);
        jn0 = jn1 * pow(-1, jn2 % N);
        AV[i] = jn(0, jn0);
      }

      for (int i = 0; i < N; i++)
      {
        jn1 = conrand(&seed);
        jn2 = 10 * conrand(&seed);
        jn0 = jn1 * pow(-1, jn2 % N);
        BV[i] = jn(1, jn0);
      }

      check = 0.0;

      for (int i = 0; i < N; i++)
      {
        int ival = N;
        check = check + AV[i] * BV[i];
        idcheck(ival, &check, AV, BV, ID);
      }

    /*Compute |AV><BV|*/

        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                idcheck(N, &check, AV, BV, ID);
                if (check > 0.5)
                    OP[i][j] = AV[i] * BV[j] / BV[i];
                else
                    OP[i][j] = AV[j] * BV[i] / BV[j];
            }
            IA[i] = i + 1;
        }

        for (int i = 0; i < N; i++)
        {
            for (int j = -1; j <= i; j += 8)
            {
                IA[i] = ((((i + 1) + (j + 1)) % N) % N);
            }
        }
		
    /* Loop 20 */

        for (int i = 0; i < N; i++)
        {
            idcheck(N, &check, AV, BV, ID);
            CV[IA[i]] = (AV[IA[i]] + BV[IA[i]]) / check;
		}
		
    /* Loop 30 */
 
        for (int i = 1; i < N; i++)
        {
            idcheck(N, &check, AV, BV, ID);
            AV[i] = AV[i - 1] * BV[i] + CV[i];  
        }
		
    /* Loop 40 */
        for (int i = 0; i < N; i++)
        {
            idcheck(N, &check, AV, BV, ID);
			for (int j = 0; j < N; j++)
            {
                if (check > 0.5)
                {
                    BOT = OP[i][j];
                    TOP = AV[j] * BV[j];
                    HOLDA = AV[j];
                    AV[j] = BV[j] + CV[j] / (TOP - BOT) * ID[i][i];
                    BV[j] = HOLDA + CV[j] / (TOP - BOT) * ID[j][j];
                    AM[i][j] = AV[j] * trig(IA[i], IA[j]);
                    BM[i][j] = BV[j] * trig(IA[j], IA[i]);
                }
                else
                {
                    BOT = OP[i][j];
                    TOP = AV[j] * BV[j];
                    HOLDA = AV[j];
                    AV[j] = BV[j] - CV[j] / (TOP - BOT) * ID[j][j];
                    BV[j] = HOLDA - CV[j] / (TOP - BOT) * ID[i][i];
                    AM[i][j] = AV[j] / trig(IA[i], IA[j]);
                    BM[i][j] = BV[j] / trig(IA[j], IA[i]);
                }
            }
        }

#ifdef FOUR
    /* Loop 50 */
		for (int i = 0; i < N; i++)
		{
			for (int j = i + 1; j < N; j++)
			{
				CM[i][j] = 0.0;
				for (int k = 0; k < N; k++)
				{
					CM[i][j] = CM[i][j] - AM[i][k] * BM[k][j] / check;
				}
			}
		}
		
        for (int i = 0; i < N; i++)
        {
            for (int j = i; j >= 0; j--)
            {
                CM[i][j] = 0.0;
                for (int k = 0; k < N; k++)
                {
					CM[i][j] = CM[i][j] + AM[i][k] * BM[k][j] / check;
                }
            }
        }
#else
		for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                CM[i][j] = 0.0;
                for (int k = 0; k < N; k++)
                {
                    if (i < j)
                        CM[i][j] = CM[i][j] - AM[i][k] * BM[k][j] / check;
                    else
                        CM[i][j] = CM[i][j] + AM[i][k] * BM[k][j] / check;
                }
            }
        }
#endif

    /* Loop 60 */
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                sum = 0.0;
                for (int k = 0; k < N; k++)
                    sum = sum + CM[i][k] * AM[j][k];
                DM[i][j] = sum;
            }
        }

        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
                CM[i][j] = DM[i][j];
        }

    /* Loop 70 */
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                sum = 0.0;
                for (int k = 0; k < N; k++)
				{
					sum = sum - CM[i][k] * BM[j][k];
				}          
                DM[i][j] = sum;
            }
        }


        HOLDA = fabs(AM[1][1]);
        HOLDB = fabs(BM[1][1]);

        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                HOLDA = max(HOLDA, fabs(AM[i][j]));
                HOLDB = max(HOLDB, fabs(BM[i][j]));
            }
        }

        TRACE3 = 0.0;

    /* Loop 80 */
#ifdef TWO
		sub = HOLDA * HOLDB;
#else
#endif
        for (int i = 0; i < N; i++)
        {
#ifdef TWO
            TRACE3 = TRACE3 + (AM[IA[i]][IA[i]] + BM[IA[i]][IA[i]]
                - DM[IA[i]][IA[i]]) / sub;
#else
            TRACE3 = TRACE3 + (AM[IA[i]][IA[i]] + BM[IA[i]][IA[i]]
                - DM[IA[i]][IA[i]]) / (HOLDA * HOLDB);
#endif
        }

        finish = cputime_();
		
        cout << "*" << "Final trace = " << TRACE3 << " and IDCHECK " << check << endl;
        cout << "*" << "-- RUNTIME -> " << finish - start << " second" << endl;
        return 0;
    }

    double trig(int i, int j)
    {
        double x, y, z;
#ifdef NOPRECISIONLOSS
		double pi = acos(-1.0);
#else
		float pi = acos(-1.0);
#endif
        x = double(i + 1) - double(j + 1);
        y = double(i + 1) + double(j + 1);
#ifdef ONE
        z = exp(sin(sqrt(x*x + y*y)*pi));
#else
        z = exp(sin(sqrt(pow(x, 2) + pow(y, 2))*pi));
#endif // ONE
        double trig = x + y + log10(fabs(1.0 + z + (x*y*z))) / (fabs(x) + fabs(y));
        return trig;
    }

    void idcheck(int N, double *check, double AV[], double BV[], double ID[][MAXDIM])
    {
        double l2 = 0, check2 = 0, a = 0, b = 0, c = 0, d = 0, pastA = 0, pastB = 0, pastC = 0, sub = 0, dcheck = *check;
        int aux = 0; 
		float twoPi = 2.0 * acos(-1.0);

#ifdef FOUR
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				ID[i][j] = cos(dcheck + (i + 1) * twoPi / N) + 2.0 * sin(dcheck + (j + 1) * twoPi / N);
			}
		}

        int j = 0;
        for (int i = 0; i < N; i++)
        {
            if (((AV[i] < 0) && (BV[j] > 0)) || ((AV[i] > 0) && (BV[j] < 0))) ID[i][j] = -1.0;
            else ID[i][j] = 1.0;
			j += 1;
        }
#else
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				if (i == j)
				{
					if ((AV[i] < 0) && (BV[j] < 0)) ID[i][j] = 1.0;
					else if ((AV[i] < 0) && (BV[j] > 0)) ID[i][j] = -1.0;
					else if ((AV[i] > 0) && (BV[j] < 0)) ID[i][j] = -1.0;
					else ID[i][j] = 1.0;
				}
				else if (i != j)
				{
					ID[i][j] = cos(dcheck + (i + 1) * twoPi / N) + 2.0 * sin(dcheck + (j + 1) * twoPi / N);
				}				
			}
		}
#endif
		
        l2 = 0.0;
#ifdef ONE
        for (int i = 0; i < N; i++)
            l2 = l2 + AV[i] * AV[i];
#else
        for (int i = 0; i < N; i++)
            l2 = l2 + pow(AV[i],2);
#endif

        l2 = sqrt(l2);
        for (int i = 0; i < N; i++)
            AV[i] = AV[i] / l2;

#ifdef ONE
        l2 = 0.0;
        for (int i = 0; i < N; i++)
            l2 = l2 + BV[i] * BV[i];
#else
        l2 = 0.0;
        for (int i = 0; i < N; i++)
            l2 = l2 + pow(BV[i], 2);
#endif

        l2 = sqrt(l2);
        for (int i = 0; i < N; i++)
            BV[i] = BV[i] / l2;

#ifdef FOUR	  	 
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
#ifdef TWO
				sub = AV[i] * BV[j];
#else
#endif
				for (int k = 3 - (((j + i) + 2) % 4); k < N; k += 4)
				{
				    pastA = a;
#ifdef TWO
				    a = a + sub * ID[j][k];
#else
				    a = a + AV[i] * BV[j] * ID[j][k];
#endif
				}
            }
        }
		
		
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
#ifdef TWO
				sub = AV[j] * BV[i];
#else
#endif
				for (int k = 3 - ((j + i + 1) % 4); k < N; k += 4)
				{
				    pastB = b;
#ifdef TWO
				    b = b + sub * ID[k][j];
#else
				    b = b + AV[j] * BV[i] * ID[k][j];
#endif
                }
			}
		}

        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
#ifdef TWO
				sub = AV[i] * BV[j];
#else
#endif
				for (int k = 3 - ((j + i) % 4); k < N; k += 4)
				{
				    pastC = c;
#ifdef TWO
				    c = c - sub * ID[k][j];
#else
				    c = c - AV[i] * BV[j] * ID[k][j];
#endif
				}
			}
        }

        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
#ifdef TWO
				sub = AV[j] * BV[i];
#else
#endif
				for (int k = 3 - ((j + i + 3) % 4); k < N; k += 4)
				{
#ifdef TWO
					d = d - sub * ID[j][k];
#else
				    d = d - AV[j] * BV[i] * ID[j][k];
#endif
				}
            }
        }
		
        aux = 3 * N;
        aux = aux % 4;
		switch(aux)
		{
			case 0:
		        dcheck = sqrt(b*b + c*c) + a;
		        check2 = pastA + b + c + d;
				break;
			case 1:
		        dcheck = sqrt(pastB*pastB + c*c) + a - b;
		        check2 = pastA + pastB + c + d;
				break;
			case 2:
		        dcheck = sqrt(b*b + c*c);
		        check2 = pastA + pastB + pastC + d; 
				break;
			case 3:
		        dcheck = sqrt(b*b + c*c);
		        check2 = a + b + c + d;
				break;
			default:
				break;	
		}
#else  //FOUR
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				for (int k = 0; k < N; k++)
				{
					int number = int((i + j + k + 3) % 4);
					//cout << number << endl;
					//getchar();
					if (number == 0)
					{ 
						a = a + AV[i] * BV[j] * ID[j][k];
						dcheck = dcheck + a;
					}
					else if (number == 1)
					{
						b = b + AV[j] * BV[i] * ID[k][j];
						dcheck = dcheck - b;
					} 
					else if (number == 2)
					{
						c = c - AV[i] * BV[j] * ID[k][j];
#ifdef ONE
						dcheck = sqrt(b*b + c*c);
#else
						dcheck = sqrt(pow(b, 2) + pow(c, 2));
#endif
					}
					else if (number == 3)
					{
						d = d - AV[j] * BV[i] * ID[j][k]; 
						check2 = a + b + c + d;
					}
                }
            }
        }
#endif // FOUR
		
        *check = min(fabs(check2), fabs(dcheck)) / max(fabs(check2), fabs(dcheck));
		/*cout << *check << ' ' << check2 << " " << a << " " << b << ' ' << c <<  ' ' << d << ' ' << endl;
		getchar();*/
    }