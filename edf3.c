/*****************************************************************************/
/*                                                                           */
/*   edf3.c                                                                  */
/*                                                                           */
/*   Calculation functions for edf program                                   */
/*   OBSOLETE - Replaced by edf4.c - Do not include in project               */
/*                                                                           */
/*   (c) 2000  W. Riley  Hamilton Technical Services  All Rights Reserved    */
/*                                                                           */
/*****************************************************************************/

// Includes
#include "edf.h"
#include "resource.h"

/***************************************************************************/
/*                                                                         */
/*                              Notice:                                    */
/*                                                                         */
/*  Algorithm by C.A Greenhall, JPL                                        */
/*  Code by W.J. Riley, Hamilton Technical Services                        */
/*  Rev. 1.20 5/14/99.                                                     */
/*  This code may be freely used for non-commercial evaluation purposes    */
/*                                                                         */
/***************************************************************************/

/***************************************************************************/
/*                                                                         */
/*                              HadamardEDF()                              */
/*                                                                         */
/*  Function to calc the estimated # of chi-squared degrees of freedom     */
/*  for the overlapping Hadamard variance of a power law noise process.    */
/*                                                                         */
/*      Parameters:    int N     = # phase data points                     */
/*                     int m     = averaging factor                        */
/*                     int b     = beta (-2 to -6)                         */
/*                                                                         */
/*      Return:        float     = Hadamard edf                            */
/*                                 or -1 if error                          */
/*                                                                         */
/*      Reference:     C.Greenhall, "Hadamard edf Algorithm", May 1999     */
/*                     (private communication via e-mail)                  */
/*                                                                         */
/*      Revision Record:                                                   */
/*          05/07/99   Created                                             */
/*          05/08/99   Debugged                                            */
/*			05/10/99   More debugging                                      */
/*			05/11/99   Ditto - May be OK                                   */
/*			05/12/99   Fixed error in smaller sum method                   */
/*			05/14/99   Fixed another error in smaller sum method           */
/*                                                                         */
/* (c) Copyright 1999 Hamilton Technical Services  All Rights Reserved     */
/*                                                                         */
/***************************************************************************/

float WINAPI HadamardEDF(int N, int m, int b)
{
	// Local macros
	#define JMAX 100                // Max # terms to calc as sum

	// Local variables
	int M=N-3*m;                    // M=N-3m
	int j;                          // Index
	int jmax;                       // Max j for summation
	int m_prime;					// Smaller sum m
	double a0;                      // edf upper bound 1st term
	double a1;                      // edf upper bound 2nd term
	double edfinv;                  // Inverse of Hadamard edf
	double r;                       // r(t)
	double r0;                      // r(0)
	double e;                       // edf summation term
	double t;                       // Argument of r(t)
	double p;                       // M/m

	// Check for allowed beta
	if(b<-6 || b>-2)
	{
		return(-1);                 // Error code
	}

	// Check for enough points: N must be > 3m
	if(M<=0)
	{
		return(-1);                 // Error code
	}

	switch(b)                       // Different calc for each noise type
	{
		case -2:                    // White FM
		{
			jmax=min(M,3*m);

			if(jmax<JMAX)			// Calc edf as sum
			{
				// TRACE_STR("Calc edf as sum");

				r0=12;
				r=0;

				for(j=1; j<=jmax; j++)
				{
					t=(double)j/(double)m;

					// Rx(t) = -|t|
					e=-20*fabs(t)+
						15*fabs(t+1)+
						15*fabs(t-1)-
						6*fabs(t+2)-
						6*fabs(t-2)+
						fabs(t+3)+
						fabs(t-3);

					r+=e*e*(1-((double)j/(double)M));
				}

				r*=2.0;
				r+=r0*r0;
				edfinv=r/(M*r0*r0);
			}
			else 					// Too many terms
			{
				if(M>=3*m)			// Calc edf as upper bound
				{
					// TRACE_STR("Calc edf as upper bound");
					
					a0=7.0/9.0;
					a1=0.5;
					p=(double)M/(double)m;
					edfinv=(1/p)*(a0-(a1/p));
				}
				else  				// Use smaller sum
				{
					// TRACE_STR("Use smaller sum");
					
					// Find nearest integer to JMAX/p

					#if(1) // Added 05/12/99
					p=(double)M/(double)m;
					#endif

					m_prime=(int)(((double)JMAX/p)+0.5);

					// Calc sum using this m

					#if(0) // Original code
					jmax=min(M,3*m_prime);
					#endif

					#if(1) // New code 05/14/99
					jmax=min(JMAX,3*m_prime);
					#endif

					r0=12;
					r=0;

					for(j=1; j<=jmax; j++)
					{
						t=(double)j/(double)m_prime;

						// Rx(t) = -|t|
						e=-20*fabs(t)+
							15*fabs(t+1)+
							15*fabs(t-1)-
							6*fabs(t+2)-
							6*fabs(t-2)+
							fabs(t+3)+
							fabs(t-3);

						#if(0) // original code
						r+=e*e*(1-((double)j/(double)M));
						#endif

						#if(1) // New code 05/12/99
						r+=e*e*(1-((double)j/(double)JMAX));
						#endif
					}

					r*=2.0;
					r+=r0*r0;

					#if(0) // Original code
					edfinv=r/(M*r0*r0);
					#endif

					#if(1) // New code 05/12/99
					edfinv=r/(JMAX*r0*r0);
					#endif
				}
			}
		}
		break;

		case -3:                    // Flicker FM
		{
			jmax=min(M,3*m);

			if(jmax<JMAX)			// Calc edf as sum
			{
				r0=48*log(2)-18*log(3);	// 13.49604
				r=0;

				for(j=1; j<=jmax; j++)
				{
					t=(double)j/(double)m;

					// Rx(t) = t*t*ln(|t|)
					e=20*t*t*Log(fabs(t))-
						15*(t+1)*(t+1)*Log(fabs((t+1)))-
						15*(t-1)*(t-1)*Log(fabs((t-1)))+
						6*(t+2)*(t+2)*Log(fabs((t+2)))+
						6*(t-2)*(t-2)*Log(fabs((t-2)))-
						(t+3)*(t+3)*Log(fabs((t+3)))-
						(t-3)*(t-3)*Log(fabs((t-3)));

					r+=e*e*(1-((double)j/(double)M));
				}

				r*=2.0;
				r+=r0*r0;
				edfinv=r/(M*r0*r0);
			}
			else				// Too many terms
			{
				if(M>=3*m)		// Calc edf as upper bound
				{
					a0=1.0;
					a1=0.62;
					p=(double)M/(double)m;
					edfinv=(1/p)*(a0-(a1/p));
				}
				else			// Use smaller sum
				{
					// Find nearest integer to JMAX/p

					#if(1) // Added 05/12/99
					p=(double)M/(double)m;
					#endif

					m_prime=(int)(((double)JMAX/p)+0.5);

					// Calc sum using this m

					#if(0) // Original code
					jmax=min(M,3*m_prime);
					#endif

					#if(1) // New code 05/14/99
					jmax=min(JMAX,3*m_prime);
					#endif

					r0=48*log(2)-18*log(3);	// 13.49604
					r=0;

					for(j=1; j<=jmax; j++)
					{
						t=(double)j/(double)m_prime;

						// Rx(t) = t*t*ln(|t|)
						e=20*t*t*Log(fabs(t))-
							15*(t+1)*(t+1)*Log(fabs((t+1)))-
							15*(t-1)*(t-1)*Log(fabs((t-1)))+
							6*(t+2)*(t+2)*Log(fabs((t+2)))+
							6*(t-2)*(t-2)*Log(fabs((t-2)))-
							(t+3)*(t+3)*Log(fabs((t+3)))-
							(t-3)*(t-3)*Log(fabs((t-3)));


						#if(0) // original code
						r+=e*e*(1-((double)j/(double)M));
						#endif

						#if(1) // New code 05/12/99
						r+=e*e*(1-((double)j/(double)JMAX));
						#endif
					}

					r*=2.0;
					r+=r0*r0;

					#if(0) // Original code
					edfinv=r/(M*r0*r0);
					#endif

					#if(1) // New code 05/12/99
					edfinv=r/(JMAX*r0*r0);
					#endif
				}
			}
		}
		break;

		case -4:                    // Random Walk FM
		{
			jmax=min(M,3*m);

			if(jmax<JMAX)           // Calc edf as sum
			{
				r0=12;
				r=0;

				for(j=1; j<=jmax; j++)
				{
					t=(double)j/(double)m;

					// Rx(t) = -|t|^3
					e=20*fabs(t*t*t)-
						15*fabs((t+1)*(t+1)*(t+1))-
						15*fabs((t-1)*(t-1)*(t-1))+
						6*fabs((t+2)*(t+2)*(t+2))+
						6*fabs((t-2)*(t-2)*(t-2))-
						fabs((t+3)*(t+3)*(t+3))-
						fabs((t-3)*(t-3)*(t-3));

					r+=e*e*(1-((double)j/(double)M));
				}

				r*=2.0;
				r+=r0*r0;
				edfinv=r/(M*r0*r0);
			}
			else                    // Too many terms
			{
				if(M>=3*m)			// Calc edf as upper bound
				{
					a0=(double)31/(double)30;
					a1=(double)17/(double)28;
					p=(double)M/(double)m;
					edfinv=(1/p)*(a0-(a1/p));
				}
				else				// Use smaller sum
				{
					// Find nearest integer to JMAX/p

					#if(1) // Added 05/12/99
					p=(double)M/(double)m;
					#endif

					m_prime=(int)(((double)JMAX/p)+0.5);

					// Calc sum using this m

					#if(0) // Original code
					jmax=min(M,3*m_prime);
					#endif

					#if(1) // New code 05/14/99
					jmax=min(JMAX,3*m_prime);
					#endif

					r0=12;
					r=0;

					for(j=1; j<=jmax; j++)
					{
						t=(double)j/(double)m_prime;

						// Rx(t) = -|t|^3
						e=20*fabs(t*t*t)-
							15*fabs((t+1)*(t+1)*(t+1))-
							15*fabs((t-1)*(t-1)*(t-1))+
							6*fabs((t+2)*(t+2)*(t+2))+
							6*fabs((t-2)*(t-2)*(t-2))-
							fabs((t+3)*(t+3)*(t+3))-
							fabs((t-3)*(t-3)*(t-3));


						#if(0) // original code
						r+=e*e*(1-((double)j/(double)M));
						#endif

						#if(1) // New code 05/12/99
						r+=e*e*(1-((double)j/(double)JMAX));
						#endif
					}

					r*=2.0;
					r+=r0*r0;

					#if(0) // Original code
					edfinv=r/(M*r0*r0);
					#endif

					#if(1) // New code 05/12/99
					edfinv=r/(JMAX*r0*r0);
					#endif
				}
			}
		}
		break;

		case -5:                    // Flicker Walk FM
		{
			jmax=min(M,3*m);

			if(jmax<JMAX)           // Calc edf as sum
			{
				r0=-192*log(2)+162*log(3);	// 44.89093
				r=0;

				for(j=1; j<=jmax; j++)
				{
					t=(double)j/(double)m;

					// Rx(t) = -t*t*t*t*ln(|t|)
					e=-20*t*t*t*t*Log((fabs(t)))+
						15*(t+1)*(t+1)*(t+1)*(t+1)*Log(fabs((t+1)))+
						15*(t-1)*(t-1)*(t-1)*(t-1)*Log(fabs((t-1)))-
						6*(t+2)*(t+2)*(t+2)*(t+2)*Log(fabs((t+2)))-
						6*(t-2)*(t-2)*(t-2)*(t-2)*Log(fabs((t-2)))+
						(t+3)*(t+3)*(t+3)*(t+3)*Log(fabs((t+3)))+
						(t-3)*(t-3)*(t-3)*(t-3)*Log(fabs((t-3)));

					r+=e*e*(1-((double)j/(double)M));
				}

				r*=2.0;
				r+=r0*r0;
				edfinv=r/(M*r0*r0);
			}
			else                    // Too many terms
			{
				if(M>=3*m)			// Calc edf as upper bound
				{
					a0=1.06;
					a1=0.53;
					p=(double)M/(double)m;
					edfinv=(1/p)*(a0-(a1/p));
				}
				else 				// Use smaller sum
				{
					// Find nearest integer to JMAX/p

					#if(1) // Added 05/12/99
					p=(double)M/(double)m;
					#endif

					m_prime=(int)(((double)JMAX/p)+0.5);

					// Calc sum using this m

					#if(0) // Original code
					jmax=min(M,3*m_prime);
					#endif

					#if(1) // New code 05/14/99
					jmax=min(JMAX,3*m_prime);
					#endif

					r0=-192*log(2)+162*log(3);	// 44.89093
					r=0;

					for(j=1; j<=jmax; j++)
					{
						t=(double)j/(double)m_prime;

						// Rx(t) = -t*t*t*t*ln(|t|)
						e=-20*t*t*t*t*Log((fabs(t)))+
							15*(t+1)*(t+1)*(t+1)*(t+1)*Log(fabs((t+1)))+
							15*(t-1)*(t-1)*(t-1)*(t-1)*Log(fabs((t-1)))-
							6*(t+2)*(t+2)*(t+2)*(t+2)*Log(fabs((t+2)))-
							6*(t-2)*(t-2)*(t-2)*(t-2)*Log(fabs((t-2)))+
							(t+3)*(t+3)*(t+3)*(t+3)*Log(fabs((t+3)))+
							(t-3)*(t-3)*(t-3)*(t-3)*Log(fabs((t-3)));


						#if(0) // original code
						r+=e*e*(1-((double)j/(double)M));
						#endif

						#if(1) // New code 05/12/99
						r+=e*e*(1-((double)j/(double)JMAX));
						#endif
					}

					r*=2.0;
					r+=r0*r0;

					#if(0) // Original code
					edfinv=r/(M*r0*r0);
					#endif

					#if(1) // New code 05/12/99
					edfinv=r/(JMAX*r0*r0);
					#endif
				}
			}
		}
		break;

		case -6:                    // Random Run FM
		{
			jmax=min(M,3*m);

			if(jmax<JMAX)           // Calc edf as sum
			{
				r0=132;
				r=0;

				for(j=1; j<=jmax; j++)
				{
					t=(double)j/(double)m;

					// Rx(t) = -|t|^5
					e=-20*fabs(t*t*t*t*t)+
						15*fabs((t+1)*(t+1)*(t+1)*(t+1)*(t+1))+
						15*fabs((t-1)*(t-1)*(t-1)*(t-1)*(t-1))-
						6*fabs((t+2)*(t+2)*(t+2)*(t+2)*(t+2))-
						6*fabs((t-2)*(t-2)*(t-2)*(t-2)*(t-2))+
						fabs((t+3)*(t+3)*(t+3)*(t+3)*(t+3))+
						fabs((t-3)*(t-3)*(t-3)*(t-3)*(t-3));

					r+=e*e*(1-((double)j/(double)M));
				}

				r*=2.0;
				r+=r0*r0;
				edfinv=r/(M*r0*r0);
			}
			else                    // Too many terms
			{
				if(M>=3*m)			// Calc edf as upper bound
				{
					a0=1.30;
					a1=0.54;
					p=(double)M/(double)m;
					edfinv=(1/p)*(a0-(a1/p));
				}
				else 				// Use smaller sum
				{
					// Find nearest integer to JMAX/p

					#if(1) // Added 05/12/99
					p=(double)M/(double)m;
					#endif

					m_prime=(int)(((double)JMAX/p)+0.5);

					// Calc sum using this m

					#if(0) // Original code
					jmax=min(M,3*m_prime);
					#endif

					#if(1) // New code 05/14/99
					jmax=min(JMAX,3*m_prime);
					#endif

					r0=132;
					r=0;

					for(j=1; j<=jmax; j++)
					{
						t=(double)j/(double)m_prime;

						// Rx(t) = -|t|^5
						e=-20*fabs(t*t*t*t*t)+
							15*fabs((t+1)*(t+1)*(t+1)*(t+1)*(t+1))+
							15*fabs((t-1)*(t-1)*(t-1)*(t-1)*(t-1))-
							6*fabs((t+2)*(t+2)*(t+2)*(t+2)*(t+2))-
							6*fabs((t-2)*(t-2)*(t-2)*(t-2)*(t-2))+
							fabs((t+3)*(t+3)*(t+3)*(t+3)*(t+3))+
							fabs((t-3)*(t-3)*(t-3)*(t-3)*(t-3));


						#if(0) // original code
						r+=e*e*(1-((double)j/(double)M));
						#endif

						#if(1) // New code 05/12/99
						r+=e*e*(1-((double)j/(double)JMAX));
						#endif
					}

					r*=2.0;
					r+=r0*r0;

					#if(0) // Original code
					edfinv=r/(M*r0*r0);
					#endif

					#if(1) // New code 05/12/99
					edfinv=r/(JMAX*r0*r0);
					#endif
				}
			}
		}
		break;

		default:                    // Unallowed beta
		{
			return(-1);             // Error code
		}
	}

	return((float)(1/edfinv));
}

/***************************************************************************/
/*                                                                         */
/*                                Log()                                    */
/*                                                                         */
/*  Wrapper function for log(x) modified to return 0 for x=0               */
/*                                                                         */
/*      Parameters:    double x  = argument for log(x)                     */
/*                                                                         */
/*      Return:        double    = log(x)                                  */
/*                                 or 0 if x=0                             */
/*                                                                         */
/*      Note: x can be 0 but must not be negative                          */
/*                                                                         */
/*      Revision Record:                                                   */
/*          05/10/99   Created                                             */
/*                                                                         */
/* (c) Copyright 1999 Hamilton Technical Services  All Rights Reserved     */
/*                                                                         */
/***************************************************************************/

double WINAPI Log(double x)
{
	if(x)
	{
		return log(x);
	}
	else
	{
		return 0;
	}
}

/***************************************************************************/


