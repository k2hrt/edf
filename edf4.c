/***************************************************************************/
/*                                                                         */
/*   edf4.c                                                                */
/*                                                                         */
/*   Calculation functions for hedf program                                */
/*                                                                         */
/*   (c) 2000-3  W. Riley Hamilton Technical Services All Rights Reserved  */
/*                                                                         */
/***************************************************************************/

// Includes
#include "edf.h"

/***************************************************************************/
/*                                                                         */
/*                              Notice:                                    */
/*                                                                         */
/*  Hadamard and General EDF algorithm by C.A Greenhall, JPL               */
/*	Thêo1 EDF algorithm by D.A. Howe, NIST                                 */
/*  Code by W.J. Riley, Hamilton Technical Services                        */
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
/*      Reference:     C. Greenhall, "Hadamard edf Algorithm", May 1999    */
/*                     (private communication via e-mail)                  */
/*                                                                         */
/*		Note:		   HadamardEDF() applies only to alpha<=0              */
/*					   Stable32 uses the following code to handle all      */
/*					   power-law noise types:                              */
/*                                                                         */
/*					   // Calc overlapping Hadamard edf                    */
/*					   // Use for alpha=0 to -4, beta=-2 to -6             */
/*					   if(nAlpha<1)                                        */
/*					   {                                                   */
/*					      // Args = # phase data pts, avg factor & beta    */
/*					      fDF=HadamardEDF(nNum, nAF, nAlpha-2);            */
/*					   }                                                   */
/*					   else // F PM or W PM (alpha 1 or 2)                 */
/*					   {                                                   */
/*					      // Calc overlapping Allan df                     */
/*						  // 2nd arg = # phase data points                 */
/*						  fDF=CalcDegFree(nAlpha, nNum, nAF);              */
/*					   }                                                   */
/*                                                                         */
/*      Revision Record:                                                   */
/*          05/15/99   Created from hedf3.c per C. Greenhall suggestions   */
/*                                                                         */
/* (c) Copyright 1999 Hamilton Technical Services  All Rights Reserved     */
/*                                                                         */
/***************************************************************************/

float WINAPI HadamardEDF(int N, int m, int b)
{
	// Local macros
	#define JMAX 100                // Max # terms to calc as sum
	#define	SUM 0					// Calc using sum
	#define BOUND 1					// Calc using upper bound
	#define SMALL 2					// Calc using smaller sum

	// Local variables
	int M=N-3*m;                    // M=N-3m
	int j;                          // Index
	int jmax;                       // Max j for summation
	int m_prime;					// Smaller sum m
	int c;							// Calculation type (SUM, BOUND or SMALL)
	double a0;                      // edf upper bound 1st term
	double a1;                      // edf upper bound 2nd term
	double edfinv;                  // Inverse of Hadamard edf
	double r=0;                     // r(t)
	double r0;                      // r(0)
	double e;                       // edf summation term
	double t;                       // Argument of r(t)
	double p;                       // M/m

	// Check for allowed beta
	if(b<-6 || b>-2)
	{
		return -1;                  // Error code
	}

	// Check for enough points: N must be >3m
	if(M<=0)
	{
		return -1;                  // Error code
	}

	// Find jmax
	jmax=min(M, 3*m);

	switch(b)                       // Different calc for each noise type
	{
		case -2:                    // White FM
		{
			// Set r0
			r0=12;

			// Test jmax
			if(jmax<JMAX)			// Calc edf as sum
			{
				// Set calc type
				c=SUM;

				// Find sum
				for(j=1; j<=jmax; j++)
				{
					t=(double)j/(double)m;
					e=Rx(t, b);
					r+=e*e*(1-((double)j/(double)M));
				}
			}
			else 					// Too many terms
			{
				// Test M
				if(M>=3*m)			// Calc edf as upper bound
				{
					// Set calc type
					c=BOUND;

					// Set coefficients
					a0=0.777778;  // 7/9
					a1=0.5;
				}
				else  				// Use smaller sum
				{
					// Set calc type
					c=SMALL;

					// Find nearest integer to JMAX/p
					p=(double)M/(double)m;
					m_prime=(int)(((double)JMAX/p)+0.5);

					// Calc sum using this m
					jmax=min(JMAX, 3*m_prime);

					// Find smaller sum
					for(j=1; j<=jmax; j++)
					{
						t=(double)j/(double)m_prime;
					    e=Rx(t, b);
						r+=e*e*(1-((double)j/(double)JMAX));
					}
				}
			}
		}
		break;

		case -3:                    // Flicker FM
		{
			// Set r0
			r0=13.49604;			// 48*log(2)-18*log(3)

			// Test jmax
			if(jmax<JMAX)			// Calc edf as sum
			{
				// Set calc type
				c=SUM;

				// Find sum
				for(j=1; j<=jmax; j++)
				{
					t=(double)j/(double)m;
					e=Rx(t, b);
					r+=e*e*(1-((double)j/(double)M));
				}
			}
			else				    // Too many terms
			{
				// Test M
				if(M>=3*m)		    // Calc edf as upper bound
				{
					// Set calc type
					c=BOUND;
					
					// Set coefficients
					a0=1.0;
					a1=0.62;
				}
				else			    // Use smaller sum
				{
					c=SMALL;

					// Find nearest integer to JMAX/p
					p=(double)M/(double)m;
					m_prime=(int)(((double)JMAX/p)+0.5);

					// Calc sum using this m
					jmax=min(JMAX, 3*m_prime);

					// Find smaller sum
					for(j=1; j<=jmax; j++)
					{
						t=(double)j/(double)m_prime;
					    e=Rx(t, b);
						r+=e*e*(1-((double)j/(double)JMAX));
					}
				}
			}
		}
		break;

		case -4:                    // Random Walk FM
		{
			// Set r0
			r0=12;

			// Test jmax
			if(jmax<JMAX)           // Calc edf as sum
			{
				// Set calc type
				c=SUM;

				// Find sum
				for(j=1; j<=jmax; j++)
				{
					t=(double)j/(double)m;
					e=Rx(t, b);
					r+=e*e*(1-((double)j/(double)M));
				}
			}
			else                    // Too many terms
			{
				// Test M
				if(M>=3*m)			// Calc edf as upper bound
				{
					// Set calc type
					c=BOUND;

					// Set coefficients
					a0=1.033334;  // 31/30
					a1=0.607143;  // 17/28
				}
				else				// Use smaller sum
				{
					// Set calc type
					c=SMALL;

					// Find nearest integer to JMAX/p
					p=(double)M/(double)m;
					m_prime=(int)(((double)JMAX/p)+0.5);

					// Calc sum using this m
					jmax=min(JMAX, 3*m_prime);

					// Find smaller sum
					for(j=1; j<=jmax; j++)
					{
						t=(double)j/(double)m_prime;
					    e=Rx(t, b);
						r+=e*e*(1-((double)j/(double)JMAX));
					}
				}
			}
		}
		break;

		case -5:                    // Flicker Walk FM
		{
			// Set r0
			r0=44.89093;			// -192*log(2)+162*log(3)

			// Test jmax
			if(jmax<JMAX)           // Calc edf as sum
			{
				// Set calc type
				c=SUM;

				// Find sum
				for(j=1; j<=jmax; j++)
				{
					t=(double)j/(double)m;
					e=Rx(t, b);
					r+=e*e*(1-((double)j/(double)M));
				}
			}
			else                    // Too many terms
			{
				// Test M
				if(M>=3*m)			// Calc edf as upper bound
				{
					// Set calc type
					c=BOUND;
					
					// Set coefficients
					a0=1.06;
					a1=0.53;
				}
				else 				// Use smaller sum
				{
					// Set calc type
					c=SMALL;

					// Find nearest integer to JMAX/p
					p=(double)M/(double)m;
					m_prime=(int)(((double)JMAX/p)+0.5);

					// Calc sum using this m
					jmax=min(JMAX, 3*m_prime);

					// Find smaller sum
					for(j=1; j<=jmax; j++)
					{
						t=(double)j/(double)m_prime;
					    e=Rx(t, b);
						r+=e*e*(1-((double)j/(double)JMAX));
					}
				}
			}
		}
		break;

		case -6:                    // Random Run FM
		{
			// Set r0
			r0=132;

			// Test jmax
			if(jmax<JMAX)           // Calc edf as sum
			{
				// Set calc type
				c=SUM;

				// Find sum
				for(j=1; j<=jmax; j++)
				{
					t=(double)j/(double)m;
					e=Rx(t, b);
					r+=e*e*(1-((double)j/(double)M));
				}
			}
			else                    // Too many terms
			{
				// Test M
				if(M>=3*m)			// Calc edf as upper bound
				{
					// Set calc type
					c=BOUND;
					
					// Set coefficients
					a0=1.30;
					a1=0.54;
				}
				else 				// Use smaller sum
				{
					// Set calc type
					c=SMALL;

					// Find nearest integer to JMAX/p
					p=(double)M/(double)m;
					m_prime=(int)(((double)JMAX/p)+0.5);

					// Calc sum using this m
					jmax=min(JMAX, 3*m_prime);

					// Find smaller sum
					for(j=1; j<=jmax; j++)
					{
						t=(double)j/(double)m_prime;
					    e=Rx(t, b);
						r+=e*e*(1-((double)j/(double)JMAX));
					}
				}
			}
		}
		break;
	}

	switch(c)                       // Different way for each calc type
	{
		case SUM:                   // Calc using sum
		{
			r*=2.0;
			r+=r0*r0;
			edfinv=r/(M*r0*r0);
		}
		break;

		case BOUND:                 // Calc using upper bound
		{
			p=(double)M/(double)m;
			edfinv=(1/p)*(a0-(a1/p));
		}
		break;

		case SMALL:                 // Calc using smaller sum
		{
			r*=2.0;
			r+=r0*r0;
			edfinv=r/(JMAX*r0*r0);
		}
		break;
	}

	return (float)(1/edfinv);		// Return Hadamard edf
}

/***************************************************************************/
/*                                                                         */
/*                                 Rx()                                    */
/*                                                                         */
/*  Function to calculate Rx(t,b) for Hamamard edf                         */
/*                                                                         */
/*      Parameters:    double t  = time (>0)                               */
/*					   int b     = beta (-2 to -6)                         */
/*                                                                         */
/*      Return:        double    = Rx(t,b)                                 */
/*                                                                         */
/*      Revision Record:                                                   */
/*          05/15/99   Created                                             */
/*                                                                         */
/* (c) Copyright 1999 Hamilton Technical Services  All Rights Reserved     */
/*                                                                         */
/***************************************************************************/

double WINAPI Rx(double t, int b)
{
	switch(b)       // Different expression for each noise type
	{
		case -2:    // White FM: Rx(t)=-|t|
		{
			return -20*t+15*(t+1)+15*fabs(t-1)-6*(t+2)-6*fabs(t-2)+
				(t+3)+fabs(t-3);
		}

		case -3:    // Flicker FM: Rx(t)=t*t*ln(|t|)
		{
			return 20*t*t*log(t)-15*(t+1)*(t+1)*Log(t+1)-
				15*(t-1)*(t-1)*Log(fabs((t-1)))+
				6*(t+2)*(t+2)*Log((t+2))+
				6*(t-2)*(t-2)*Log(fabs((t-2)))-
				(t+3)*(t+3)*Log((t+3))-(t-3)*(t-3)*Log(fabs((t-3)));
		}

		case -4:    // Random Walk FM: Rx(t)=-|t|^3
		{
			return 20*t*t*t-
				15*(t+1)*(t+1)*(t+1)-
				15*fabs((t-1)*(t-1)*(t-1))+
				6*(t+2)*(t+2)*(t+2)+
				6*fabs((t-2)*(t-2)*(t-2))-
				(t+3)*(t+3)*(t+3)-
				fabs((t-3)*(t-3)*(t-3));
		}

		case -5:    // Flicker Walk FM: Rx(t)=-t*t*t*t*ln(|t|)
		{
			return -20*t*t*t*t*log((t))+
				15*(t+1)*(t+1)*(t+1)*(t+1)*Log((t+1))+
				15*(t-1)*(t-1)*(t-1)*(t-1)*Log(fabs((t-1)))-
				6*(t+2)*(t+2)*(t+2)*(t+2)*Log((t+2))-
				6*(t-2)*(t-2)*(t-2)*(t-2)*Log(fabs((t-2)))+
				(t+3)*(t+3)*(t+3)*(t+3)*Log((t+3))+
				(t-3)*(t-3)*(t-3)*(t-3)*Log(fabs((t-3)));
		}

		case -6:    // Random Run FM: Rx(t)=-|t|^5
		{
			return -20*t*t*t*t*t+
				15*(t+1)*(t+1)*(t+1)*(t+1)*(t+1)+
				15*fabs((t-1)*(t-1)*(t-1)*(t-1)*(t-1))-
				6*(t+2)*(t+2)*(t+2)*(t+2)*(t+2)-
				6*fabs((t-2)*(t-2)*(t-2)*(t-2)*(t-2))+
				(t+3)*(t+3)*(t+3)*(t+3)*(t+3)+
				fabs((t-3)*(t-3)*(t-3)*(t-3)*(t-3));
		}
	}
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
/*                                                                         */
/* 		                       ModTotvarEDF()                              */
/*                                                                         */
/* 	Function to calc the estimated # of chi-squared degrees of freedom     */
/*  for the modified total variance (MTOT) of a power law noise process.   */
/*                                                                         */
/*  Params:			int   nAlpha = alpha                                   */
/*                  float fRatio = T/tau = N/m                             */
/*                                                                         */
/*	Return:         float edf = estimated degrees of freedom               */
/*								or -1 if error                             */
/*                                                                         */
/*  The edf is modeled by the expression b(T/tau) - c, where:              */
/*                                                                         */
/*  Alpha    Noise     b      c        Original values                     */
/*   2        W PM    1.90   2.10                                          */
/*   1        F PM    1.20   1.40                                          */
/*   0        W FM    1.10   1.20                                          */
/*  -1        F FM    0.85   0.50                                          */
/*  -2       RW FM    0.75   0.31                                          */
/*                                                                         */
/*  Alpha    Noise     b      c        New values per D. Howe              */
/*   2        W PM    2.16   2.28      e-mail of 05/20/00                  */
/*   1        F PM    1.73   2.99                                          */
/*   0        W FM    1.33   1.89                                          */
/*  -1        F FM    0.92   0.76                                          */
/*  -2       RW FM    0.79   0.37                                          */
/*                                                                         */
/*	Notes:			(1) This function applies only to -2 >= alpha >= 2     */
/*						(W PM, F PM, W FM, F FM and RW FM noise)           */
/*					(2) This function applies only for m > 8.0             */
/*					(3) T is defined as tau0 * # phase data points         */
/*					(4) The t/tau ratio should never be < 3                */
/*                                                                         */
/* 	References:		E-Mails from D.A. Howe/NIST, 3/13/00 & 5/20/00         */
/*                                                                         */
/*  Revision record:                                                       */
/*      03/25/00    Created                                                */
/*		03/30/00	Parameter verification revised                         */
/*		05/21/00	Coefficients changed                                   */
/*                                                                         */
/* (c) Copyright 2000  Hamilton Technical Services   All Rights Reserved   */
/*                                                                         */
/***************************************************************************/

float WINAPI ModTotvarEDF(int nAlpha, float fRatio)
{
	// Local variables
	float fEDF;

	// Verify parameters
	if(fRatio<3.0)
	{
		return -1.0; // Error
	}

	if(nAlpha<-2 || nAlpha>2)
	{
		return -1.0; // Error
	}

	switch(nAlpha)
	{
		case 2:	// W PM noise
		{
			// fEDF=(float)(1.90*fRatio - 2.10);
			fEDF=(float)(2.16*fRatio - 2.28);
		}
		break;

		case 1:	// F PM noise
		{
			// fEDF=(float)(1.20*fRatio - 1.40);
			fEDF=(float)(1.73*fRatio - 2.99);
		}
		break;

		case 0:	// W FM noise
		{
			// fEDF=(float)(1.10*fRatio - 1.20);
			fEDF=(float)(1.33*fRatio - 1.89);
		}
		break;

		case -1: // F FM noise
		{
			// fEDF=(float)(0.85*fRatio - 0.50);
			fEDF=(float)(0.92*fRatio - 0.76);
		}
		break;

		case -2: // RW FM noise
		{
			// fEDF=(float)(0.75*fRatio - 0.31);
			fEDF=(float)(0.79*fRatio - 0.37);
		}
		break;
	}

	return fEDF;
}

/***************************************************************************/
/*                                                                         */
/* 		                        EDF()                                      */
/*                                                                         */
/* 	Function to calc the estimated # of chi-squared degrees of freedom     */
/*  for the modified Allan variance (MVAR) of a power law noise process.   */
/*                                                                         */
/*  Based on edf() of STAB_F38.C from Stable/DOS                           */
/*                                                                         */
/*	Original:                                                              */
/*	Parameters:		int beta  = exponent of Sx(f) of power law noise       */
/*                              where: 0 = W  í                            */
/*                                    -1 = F  í                            */
/*                                    -2 = W  f                            */
/*                                    -3 = F  f	                           */
/*                                    -4 = RW f                            */
/*					int m     = averaging factor                           */
/*					float p   = M/m                                        */
/*							    where: M = # phase analysis points         */
/*										 = N-3m+1                          */
/*										   where: N = # phase data points  */
/*  New Params:		int nBeta = beta                                       */
/*                  int nAvgFactor = m                                     */
/*                  int nNum = M                                           */
/*                                                                         */
/*                  nNum and nAvgFactor are used to find p                 */
/*                                                                         */
/*	Return:         float edf = estimated degrees of freedom               */
/*                                                                         */
/*	Note:			This function calculates the MVAR edf for the usual    */
/*					case of fully overlapping samples (m1=1 & r=m)         */
/*                                                                         */
/*	Reference:		C. H. Greenhall, "Estimating the Modified Allan        */
/*                  Variance", 1995 IEEE International Frequency Control   */
/*                  Symposium (to be published)                            */
/*				                                                           */
/*	Revision record:                                                       */
/*			01/01/98	Changes for MS VC++ compatibility & warnings       */
/*                                                                         */
/*	(c) Copyright 1996-8 Hamilton Technical Services   All Rights Reserved */
/*                                                                         */
/***************************************************************************/

float WINAPI EDF(int nBeta, int nAvgFactor, int nNum)
{
	// Local variables
	// the a0 and a1 coefficients for estimating MVAR edf
	// are stored in the following static arrays
	// á varies from 0 to -4 ( -á varies from 0 to 4)
	// m varies from 1 to >2 (m-1 varies from 0 to 2)

	// a0 term: a0[-á][m-1]
	static float a0[5][3] =
	{
		{ 0.51429F, 0.93506F, 1.22450F },
		{ 0.57640F, 0.97339F, 1.00300F },
		{ 0.66667F, 1.01010F, 0.96774F },
		{ 0.81057F, 1.02660F, 0.94663F },
		{ 1.00000F, 0.86580F, 0.76791F }
	};

	// a1 term: a1[-á][m-1]
	static float a1[5][3] =
	{
		{ 0.00000F, 0.00000F, 0.58929F },
		{ 0.00000F, 0.00000F, 0.60163F },
		{ 0.00000F, 0.00000F, 0.57124F },
		{ 0.00000F, 0.00000F, 0.41643F },
		{ 0.00000F, 0.00000F, 0.41115F }
	};

	// Convert args to those of original edf() function
	int beta=nBeta;
	int m=nAvgFactor;
	float p=(float)(nNum)/(float)nAvgFactor;

	// truncate m
	if(m>3)
	{
		m=3;
	}

	// calc and return edf
	return (a0[-beta][m-1]*p)/(1-(a1[-beta][m-1]/p));
}

/***************************************************************************/
/*                                                                         */
/*							CalcDegFree()                                  */
/*                                                                         */
/*      Function to calculate the # of chi-square degrees of freedom       */
/*      for phase or frequency data.                                       */
/*                                                                         */
/*		Reference:     NIST Technical Note 1337, p.TN-85, Table 12-4       */
/*                     (as corrected).                                     */
/*                                                                         */
/*		Parameters:    int alpha = power-law exponent of Sy(f) noise type  */
/*					   int n     = # phase data points                     */
/*					             = # frequency data points + 1             */
/*					   int m     = fully-overlapping averaging factor      */
/*                                                                         */
/*		Return:        float     = # chi-square degrees of freedom,        */
/*                                 or -1 if error                          */
/*                                                                         */
/*		Revision record:                                                   */
/* 			12/30/91   Created                                             */
/*			12/31/91   Renamed                                             */
/*			01/04/92   Changed switch from (2-a) to a                      */
/*					   Added default case for NG alpha                     */
/*					   Changed n and m parameters to ints                  */
/*					   Do calcs with double versions dn and dm             */
/*					   Added error traps                                   */
/*			01/18/92   Edited title block                                  */
/*			02/03/96   Modified for use as Win 3.1 DLL                     */
/*			01/01/98   Changes for MS VC++ compatibility & warnings        */
/*                                                                         */
/* (c) Copyright 1991-8  Hamilton Technical Services  All Rights Reserved  */
/*                                                                         */
/***************************************************************************/

float WINAPI CalcDegFree(int alpha, int n, int m)
{
	double dn;                                      /* double version of n */
	double dm;                                      /* double version of m */
	double df;                                     /* # degrees of freedom */
	double df1,df2,df3;                         /* terms of df calculation */

	if( (n==0) || (m==0) )
	{
		return(-1.0);                                        /* error code */
	}

	dn=(double) abs(n);                         /* make sure n is positive */
	dm=(double) abs(m);                         /* make sure m is positive */

	switch(alpha)          /* alpha is exponent of power-law noise process */
	{                                    /* slope of Sy(f) on log-log plot */
			case 2:
			{
				if(dn==dm)
				{
					return(-1.0);                            /* error code */
				}
				
				// Note: If dm=dn/2 (max allowable avg factor), get df=0
				// Should get df=1.  Handle as special case
				if(m+m==n)
				{
					return(1.0);
				}

				df=((dn+1)*(dn-(2*dm)))/(2*(dn-dm));  /* white phase noise */
				break;                                          /* alpha=2 */
			}
			case 1:
			{
				if(dn==1)
				{
					return(-1.0);                           /* error code */
				}

				// Note: If dm=dn/2 (max allowable avg factor), get neg df2
				// which causes error return.  Should get df=1
				// Handle as special case
				if(m+m==n)
				{
					return(1.0);
				}

				df1=log((((2*dm)+1)*(dn-1))/(4));   /* flicker phase noise */
				df2=log((dn-1)/(2*dm));                        /* alpha=1) */

				if(df1*df2<0)
				{
					return(-1.0);                            /* error code */
				}

				df=exp(sqrt(df1*df2));
				break;
			}
			case 0:
			{
				df1=(4*dm*dm)/((4*dm*dm)+5);      /* white frequency noise */
				df2=(3*(dn-1))/(2*dm);                          /* alpha=0 */
				df3=(2*(dn-2))/(dn);
				df=(df2-df3)*df1;
				break;
			}
			case -1:
			{
				if(dm==1.0)                     /* flicker frequency noise */
				{                                              /* alpha=-1 */
					df=(2*(dn-2)*(dn-2))/((2.3*dn)-4.9);
					break;
				}                                                /* end if */
				else
				{
					df=(5*dn*dn)/((4*dm)*(dn+(3*dm)));
					break;
				}                                              /* end else */
			}
			case -2:
			{
				if(dn==3)
				{
					return(-1.0);                            /* error code */
				}
				df1=(dn-1)*(dn-1);          /* random-walk frequency noise */
				df2=(dn-2)/((dm)*(dn-3)*(dn-3));               /* alpha=-2 */
				df3=(3*dm)*(1-dn)+(4*dm*dm);
				df=df2*(df1+df3);
				break;
			}
			default:                                           /* NG alpha */
			{
				df=-1.0;                                     /* error code */
			}
	}                                                        /* end switch */

	return(float)df;                     /* chi-squared degrees of freedom */
}                                                     /* end CalcDegFree() */

/***************************************************************************/
/*                                                                         */
/* 		                       TotvarEDF()                                 */
/*                                                                         */
/* 	Function to calc the estimated # of chi-squared degrees of freedom     */
/*  for the total variance (TOTVAR) of a power law noise process.          */
/*                                                                         */
/*  Params:			int   nAlpha = alpha                                   */
/*                  float fRatio = T/ç                                     */
/*                                                                         */
/*	Return:         float edf = estimated degrees of freedom               */
/*								or -1 if error                             */
/*                                                                         */
/*	Notes:			(1) This function applies only to -2 >= alpha >= 0     */
/*						(W FM, F FM and RW FM noise)                       */
/*					(2) This function applies only for T/ç  >= 2.0         */
/*					(3) T is defined as ç0 * # phase data points           */
/*                                                                         */
/* 	Reference:		D.A. Howe and C.A. Greenhall, "Total Variance: A       */
/* 					Progress Report on a New Frequency Stability           */
/* 					Characterization", Proc. 29th PTTI Meeting,            */
/* 					Dec. 2, 1997, (to be published).                       */
/*                                                                         */
/*  Revision record:                                                       */
/*      12/06/97    Created                                                */
/*		12/09/97	Added to FrequenC Library                              */
/*		01/01/98	Changes for MS VC++ compatibility & warnings           */
/*                                                                         */
/* (c) Copyright 1997-8	 Hamilton Technical Services  All Rights Reserved  */
/*                                                                         */
/***************************************************************************/

float WINAPI TotvarEDF(int nAlpha, float fRatio)
{
	// Local variables
	float fEDF;

	// Verify parameters
	if(fRatio<2.0)
	{
		return -1.0; // Error
	}

	if(nAlpha<-2 || nAlpha>0)
	{
		return -1.0; // Error
	}

	switch(nAlpha)
	{
		case 0:		// W FM noise
		{
			fEDF=(float)(1.500*fRatio);
		}
		break;

		case -1:	// F FM noise
		{
			fEDF=(float)(1.16832*fRatio - 0.222);
		}
		break;

		case -2:	// RW FM noise
		{
			fEDF=(float)(0.92715*fRatio -0.358);
		}
		break;
	}

	return fEDF;
}

/*****************************************************************************/
/*                                                                           */
/*                             HadTotvarEDF()                                */
/*                                                                           */
/*  Function to determine the estimated # of chi-squared degrees of freedom  */
/*  for the Hadamard total variance (HTOT) of a power law noise process.     */
/*                                                                           */
/*  Params:         int   nAlpha = alpha                                     */
/*                  float fRatio = T/tau = N/m                               */
/*                                 where: N = # phase data points            */
/*                                        m = avg factor = tau/tau0          */
/*                                                                           */
/*  Return:         float edf = estimated degrees of freedom                 */
/*                              or -1 if error                               */
/*                                                                           */
/*  The edf is modeled by the expression (T/tau)/[b0+b1(tau/T)], where:      */
/*                                                                           */
/*  Alpha    Noise     b0      b1       Values per D. Howe, et.al            */
/*   0        W FM    0.559   1.004     12/24/00 draft of PTTI'00 paper      */
/*  -1        F FM    0.868   1.140                                          */
/*  -2       RW FM    0.938   1.696                                          */
/*  -3       FW FM    0.974   2.554                                          */
/*  -4       RR FM    1.276   3.149                                          */
/*                                                                           */
/*  Notes:          (1) This function applies only to -4 >= alpha >= 0       */
/*                  (2) This function applies only for m > 16                */
/*                  (3) T is defined as tau0 * # phase data points           */
/*                                                                           */
/*  References:     (1) E-Mail from D.A. Howe/NIST, 10/27/00                 */
/*                  (2) Final draft of PTTI 2000 HTOT paper:                 */
/*                      Final3-PTTI-00.pdf rcv'd 12/24/00                    */
/*                                                                           */
/*  Revision record:                                                         */
/*      10/28/00    Created                                                  */
/*      12/24/00    Revised per Ref. 2                                       */
/*                                                                           */
/* (c) Copyright 2000  Hamilton Technical Services   All Rights Reserved     */
/*                                                                           */
/*****************************************************************************/

float WINAPI HadTotvarEDF(int nAlpha, float fRatio)
{
    // Local variables
    float fEDF;
    double b0;
    double b1;

	// Verify parameters
    if(nAlpha<-4 || nAlpha>0)
    {
        return -1.0; // Error
    }

    switch(nAlpha)
    {
        case 0: // W FM noise
        {
			b0=0.559;
            b1=1.004;
        }
        break;

        case -1: // F FM noise
        {
            b0=0.868;
            b1=1.140;
        }
        break;

        case -2: // RW FM noise
        {
            b0=0.938;
            b1=1.696;
        }
        break;

        case -3: // FW FM noise
        {
            b0=0.974;
            b1=2.554;
        }
        break;

        case -4: // RR FM noise
        {
            b0=1.276;
            b1=3.149;
        }
        break;
    }

    fEDF=(float)((fRatio)/(b0+(b1/fRatio)));
    return fEDF;
}

/*****************************************************************************/
/*                                                                           */
/* 		                        Theo1EDF()                                   */
/*                                                                           */
/* 	Function to determine the estimated # of chi-squared degrees of freedom  */
/*  for the Thêo1 variance of a power law noise process.                     */
/*                                                                           */
/*  Params:			int nAlpha   = Power law noise exponent                  */
/*					int nNum     = # phase analysis points                   */
/*                  float fRatio = Tau/Tau0 ratio = m = averaging factor     */
/*                                                                           */
/*	Return:         float EDF    = estimated degrees of freedom              */
/*								   or -1 if error                            */
/*                                                                           */
/*	Notes:			(1) This function applies to -2 >= alpha >= 2            */
/*						(W PM, F PM, W FM, F FM and RW FM noise)             */
/*					(2) fRatio is entered as AF=m, and is scaled to          */
/*						fR=0.75m in this function for use in the             */
/*						Thêo1 edf formuale                                   */
/*                                                                           */
/* 	References:		(1) D.A. Howe e-mail of 05/01/03                         */
/*					(2) D.A. Howe e-mail of 05/16/03 with "Thêo1 Summary"    */
/*					(3) D.A. Howe e-mail of 05/18/03 clarifications          */
/*					(4) D.A. Howe e-mail of 05/29/03 with F PM edf formula   */
/*						and fR=0.75m instead of previous fR=0.75(m-1)        */
/*                                                                           */
/*  Revision record:                                                         */
/*      05/01/03    Created                                                  */
/*		05/16/03	Updated per latest D. Howe inputs                        */
/*		05/18/03	Ditto                                                    */
/*		05/29/03	Ditto again                                              */
/*                                                                           */
/* (c) Copyright 2003	 Hamilton Technical Services   All Rights Reserved   */
/*                                                                           */
/*****************************************************************************/

float WINAPI Theo1EDF(int nAlpha, int nNum, float fRatio)
{
	// Local variables
	float fEDF;
	float fR;

	// Verify parameters
	if(nAlpha<-2 || nAlpha>2)
	{
		return -1.0; // Error
	}

	// Calc Thêo1 ratio 0.75m
	fR=(float)(0.75*fRatio);

	switch(nAlpha)
	{
		case 2:		// W PM noise
		{
			fEDF=(float)(((0.86*(nNum+1)*(nNum-((4.0/3.0)*fR)))/
				(nNum-fR))*(fR/(fR+1.14)));
		}
		break;

		case 1:		// F PM noise
		{
			fEDF=(float)((((4.798*nNum*nNum)-(6.374*nNum*fR)+(12.387*fR))/
				(sqrt(fR+36.6)*(nNum-fR)))*(fR/(fR+0.3)));
		}
		break;

		case 0:		// W FM noise
		{
			fEDF=(float)((((4.1*nNum+0.8)/fR)-((3.1*nNum+6.5)/nNum))*
				(pow(fR,1.5)/(pow(fR,1.5)+5.2)));
		}
		break;

		case -1:	// F FM noise
		{
			fEDF=(float)(((2*nNum*nNum-1.3*nNum*fR-3.5*fR)/
				(nNum*fR))*((fR*fR*fR)/
				(fR*fR*fR+2.3)));
		}
		break;

		case -2:	// RW FM noise
		{
			fEDF=(float)(((4.4*nNum-2)/(2.9*fR))*(((4.4*nNum-1)*
				(4.4*nNum-1)-8.6*fR*(4.4*nNum-1)+
				11.4*fR*fR)/((4.4*nNum-3)*(4.4*nNum-3))));
		}
		break;
	}

	return fEDF;
}

/***************************************************************************/
/*                                                                         */
/*                               CombinedEDF()                             */
/*                                                                         */
/*  Function to calc the estimated # of chi-squared degrees of freedom     */
/*  for the normal and overlapping, modified and unmodified, Allan and     */
/*  Hadamard variances of a power law noise process.                       */
/*                                                                         */
/*      Parameters:    int N = # phase data points at tau0                 */
/*                     int m = averaging factor = tau/tau0                 */
/*                     int a = alpha (-4 to +2)                            */
/*					   int d = order of phase difference                   */
/*                             1 = first difference (not used)             */
/*							   2 = Allan variance                          */
/*							   3 = Hadamard variance                       */
/*					   int F = filter factor                               */
/*							   1 = modified                                */
/*							   m = unmodified                              */
/*					   int S = stride factor                               */
/*							   1 = long (tau)                              */
/*							   m = short (tau0)                            */
/*					   int v = calc version type (0=simple, 1=full)        */
/*                                                                         */
/*      Return:        float = edf (or -1 if error)                        */
/*                                                                         */
/*      Reference:     C. Greenhall, "Normalized Uncertainty of            */
/*					   Difference Variances", June 4, 2003                 */
/*                     (private communication via e-mail)                  */
/*                                                                         */
/*		Notes:		   1. Applies for alpha between -4 and 2 (RR FM to     */
/*						  W PM) for long or short stride (tau or tau0),    */
/*						  modified or unmodified Allan or Hadamard         */
/*                        variances.  Does not apply to Total variances    */
/*                        or Thêo1.                                        */
/*                     2. Two algorithms are used, a simplified version    */
/*						  using a truncated summation, and a full version  */
/*                        using 4 cases and 3 tables.                      */
/*                                                                         */
/*      Revision Record:                                                   */
/*          06/08/03   Created                                             */
/*			06/17/03   Running                                             */
/*			06/23/03   Added full version                                  */
/*			06/24/03   Full version running                                */
/*          06/26/03   Looking pretty good                                 */
/*			06/27/03   Added calc version type parameter                   */
/*					   Added calc parameter recording                      */
/*			07/03/03   Made corrections per C. Greenhall e-mail of 7/1/03  */
/*			07/16/03   Made corrections per C/ Greenhall e-mail of 7/9/03  */
/*			                                                               */
/* (c) Copyright 2003 Hamilton Technical Services  All Rights Reserved     */
/*                                                                         */
/***************************************************************************/

float WINAPI CombinedEDF(int N, int m, int a, int d, int F, int S, int v)
{
	// Local macros
	#define RECORD 1	// Control recording of calculation: 0=No, 1=Yes
	
	// Notes re function inputs:
	// N is Greenhall Nx
	// m is Greenhall m
	// a is Greenhall alpha
	// d is Greenhall d
	// F is Greenhall F
	// S is Greenhall S
	
	// Supporting function prototypes
	// Notes re functions:
	// Function names are capitalized
	// Functions Sw(), Sx(), and Sz() are same as Greehall names
	// Function Sum is Greenhall BasicSum()
	double Sw(double t, int a);
	double Sx(double t, double F, int a);
	double Sz(double t, double F, int a, int d);
	double Sum(int J, double M, double S, double F, int a, int d);
	// NR binomial coefficient function protptypes
	double bico(int n, int k);
	double factln(int n);
	double gammln(double xx);

	// Local variables
	char szBuffer[BUFFER_SIZE+1];   // Text buffer
	int L;				// Same as Greenhall L
	double M;			// Same as Greenhall M
	int J;				// Summation limit
	int Jmax=100;		// Max J to use summation
	int k;				// Summation index
	int K;				// ceil(r)
	double r;			// Same as Greenhall r=M/S
	double s=(double)S; // double version of S
	double a0;			// Same as Greenhall a0 (from Table 2)
	double a1;			// Same as Greenhall a1 (from table 2)
	double b0;			// Same as Greenhall b0 (from Table 2)
	double b1;			// Same as Greenhall b1 (from table 2)
	double edf;			// Same as Greenhall edf

	// Lookup tables
	
	// Table 1: Coefficients for modified variances
	// a0 & a1 as a function of a (row) & d (column)
	//    d=1           d=2           d=3
	// a0      a1     a0     a1     a0     a1
	double T1[7][6]={
	{0.667, 0.333, 0.778, 0.500, 0.880, 0.667},  // a=+2
	{0.840, 0.345, 0.997, 0.616, 1.141, 0.843},  // a=+1
	{1.079, 0.368, 1.033, 0.607, 1.184, 0.848},  // a= 0
	{0.000, 0.000, 1.048, 0.534, 1.180, 0.816},  // a=-1
	{0.000, 0.000, 1.302, 0.535, 1.175, 0.777},  // a=-2
	{0.000, 0.000, 0.000, 0.000, 1.194, 0.703},  // a=-3
	{0.000, 0.000, 0.000, 0.000, 1.489, 0.702}}; // a=-4
	// Access as T1[2-a][i], where a=alpha
	// and i=2*(d-1)+j and j=a#
	
	// Table 2: Coefficients for unmodified variances
	// a0 & a1 as a function of a (row) & d (column)
	//    d=1           d=2           d=3
	// a0      a1     a0     a1     a0     a1
	double T2[7][6]={
	{1.500, 0.500, 1.944, 1.000, 2.310, 1.500},  // a=+2
	{78.60, 25.20, 790.0, 410.0, 9950., 6520.},  // a=+1
	{0.667, 0.167, 0.667, 0.333, 0.778, 0.500},  // a= 0
	{0.000, 0.000, 0.852, 0.375, 0.997, 0.617},  // a=-1
	{0.000, 0.000, 1.079, 0.368, 1.033, 0.607},  // a=-2
	{0.000, 0.000, 0.000, 0.000, 1.053, 0.553},  // a=-3
	{0.000, 0.000, 0.000, 0.000, 1.302, 0.535}}; // a=-4
	// Access as T2[2-a][i], where a=alpha
	// and i=2*(d-1)+j and j=a#

	// Table 3: Coefficients for log denominator,
	// unmodified variances, F FM
	// b0 & b1 as a function of b (row) & d (column)
	//    d=1           d=2           d=3
	// b0      b1     b0     b1     b0     b1
	double T3[6]=
	{6.000, 4.000, 15.23, 12.00, 47.80, 40.00};
	// Access as T3[2*(d-1)+j], where j=b#

	// Check arguments
	// d must be 1, 2 or 3
	if( (d<1) || (d>3))
	{
		return -1; // Error
	}

	// F must be either 1 or m
	if( (F!=1) && (F!=m) )
	{
		return -1; // Error
	}

	// m must be >=1
	if(m<1)
	{
		return -1; // Error
	}

	// a must be between -4 and +2
	if( (a<-4) || (a>2) )
	{
		return -1; // Error
	}
	
	// S must be either 1 or m
	if( (S!=1) && (S!=m) )
	{
		return -1; // Error
	}

	// Check for legal alpha: a+2d>1
	if(a+2*d<=1)
	{
		return -1; // Error - Illegal alpha
	}
	
	if(RECORD)
	{
		// Record N
		sprintf(szBuffer, "\r\n    N=%d\r\n", N);
		strcat(szParams, szBuffer);

		// Record m
		sprintf(szBuffer, "    m=%d\r\n", m);
		strcat(szParams, szBuffer);

		// Record a
		sprintf(szBuffer, "    a=%d\r\n", a);
		strcat(szParams, szBuffer);

		// Record d
		sprintf(szBuffer, "    d=%d\r\n", d);
		strcat(szParams, szBuffer);

		// Record F
		sprintf(szBuffer, "    F=%d\r\n", F);
		strcat(szParams, szBuffer);

		// Record S
		sprintf(szBuffer, "    S=%d\r\n", S);
		strcat(szParams, szBuffer);

		// Record Jmax
		sprintf(szBuffer, "    Jmax=%d\r\n", Jmax);
		strcat(szParams, szBuffer);
	}

	// Find # summands
	L=(m/F)+m*d; // L is always an int (F is either 1 or m)

	if(RECORD)
	{
		// Record L
		sprintf(szBuffer, "    L=(m/F)+m*d=%d\r\n", L);
		strcat(szParams, szBuffer);
	}

	// Calc M - not necessarily an int
	M=(double)(S); 
	M*=(double)(N-L);
	M/=(double)m;
	M=floor(M);
	M+=1.0;

	if(RECORD)
	{
		// Record M
		sprintf(szBuffer, "    M=1+floor(S*(N-L)/m)=%e\r\n", M);
		ShortenExp(szBuffer);
		strcat(szParams, szBuffer);
	}

	J=min((int)M, (d+1)*S);

	if(RECORD)
	{
		// Record J
		sprintf(szBuffer, "    J=min((int)M,(d+1)*S)=%d\r\n", J);
		strcat(szParams, szBuffer);
	}

	// Check # data points
	if(N<L)
	{
		return -1; // Error - N<L
	}

	if(v==0) // v=0 indicates simple version
	{
		if(RECORD)
		{
			// Record Sz()
			sprintf(szBuffer, "    Sz(0,F,a,d)=%e\r\n", Sz(0, F, a, d));
			ShortenExp(szBuffer);
			strcat(szParams, szBuffer);

			// Record Sum()
			sprintf(szBuffer, "    Sum(J,M,S,F,a,d)=%e\r\n",
				Sum(J, M, S, F, a, d));
			ShortenExp(szBuffer);
			strcat(szParams, szBuffer);
		}

		// Calc edf by simplified version
		edf=(Sz(0, F, a, d)*Sz(0, F, a, d)*M/Sum(J, M, S, F, a, d));
		return (float)edf;
	}

	// Calc edf by full version
	r=M/s;

	if(RECORD)
	{
		// Record r
		sprintf(szBuffer, "    r=M/S=%e\r\n", r);
		ShortenExp(szBuffer);
		strcat(szParams, szBuffer);
	}

	// Sort by case
	if(F==1) // Case 1. Modified variances, all alpha
			 // Note: Apparently this is also the code used by unmodified
			 //	variances when F=m=1
	{
		if(RECORD)
		{
			// Record case #
			strcat(szParams, "    Case 1\r\n");
		}

		if(J<=Jmax)
		{
			if(RECORD)
			{
				// Record test
				strcat(szParams, "    J<=Jmax\r\n");

				// Record Sz()
				sprintf(szBuffer, "    Sz(0,1,a,d)=%e\r\n", Sz(0, 1, a, d));
				ShortenExp(szBuffer);
				strcat(szParams, szBuffer);

				// Record Sum()
				sprintf(szBuffer, "    Sum(J,M,S,1,a,d)=%e\r\n",
					Sum(J, M, S, 1, a, d));
				ShortenExp(szBuffer);
				strcat(szParams, szBuffer);
			}

			// Calc edf
			edf=(Sz(0, 1, a, d)*Sz(0, 1, a, d)*M/
				Sum(J, M, S, 1, a, d));
		}
		else
		{
			if(RECORD)
			{
				// Record test
				strcat(szParams, "    J>Jmax\r\n");
			}

			if(r>=d+1)
			{
				if(RECORD)
				{
					// Record test
					strcat(szParams, "    r>=d+1\r\n");
				}

				// Get a0 & a1 from Table 1
				a0=T1[2-a][2*(d-1)+0];
				a1=T1[2-a][2*(d-1)+1];

				if(RECORD)
				{
					// Record a0 & a1
					strcat(szParams, "    a0,a1 from Table 1\r\n");
					sprintf(szBuffer, "    a0=%e, a1=%e\r\n", a0, a1);
					ShortenExp(szBuffer);
					strcat(szParams, szBuffer);
				}

				// Calc edf
				edf=((1.0/r)*(a0-(a1/r)));
				edf=(1.0/edf);
			}
			else 
			{
				if(RECORD)
				{
					// Record test
					strcat(szParams, "    r<d+1\r\n");

					// Record Sz()
					sprintf(szBuffer, "    Sz(0,1,a,d)=%e\r\n", Sz(0, 1, a, d));
					ShortenExp(szBuffer);
					strcat(szParams, szBuffer);

					// Record Sum()
					sprintf(szBuffer,
						"    Sum(Jmax,Jmax,(double)Jmax/r,1,a,d)=%e\r\n",
						Sum(Jmax, Jmax, (double)Jmax/r, 1, a, d));
					ShortenExp(szBuffer);
					strcat(szParams, szBuffer);
				}

				// Calc edf
				edf=(Sz(0, 1, a, d)*Sz(0, 1, a, d)*Jmax/
					Sum(Jmax, Jmax, (double)Jmax/r, 1, a, d));
			}
		}
		return (float)edf;  // EDF for Case 1
	}
	else // Unmodified variances: F=m
	{
		if(a<=0) // Case 2. W FM to RR FM
		{
			if(RECORD)
			{
				// Record case #
				strcat(szParams, "    Case 2\r\n");
			}

			if(J<=Jmax)
			{
				if(RECORD)
				{
					// Record test
					strcat(szParams, "    J<=Jmax\r\n");
				}

				if(m*(d+1)<=Jmax) // m'=m;
				{
					if(RECORD)
					{
						// Record test
						strcat(szParams, "    m*(d+1)<=Jmax\r\n");
						strcat(szParams, "    m'=m\r\n");

						// Record Sz()
						sprintf(szBuffer, "    Sz(0,m,a,d)=%e\r\n",
							Sz(0, m, a, d));
						ShortenExp(szBuffer);
						strcat(szParams, szBuffer);

						// Record Sum()
						sprintf(szBuffer, "    Sum(J,M,S,m,a,d)=%e\r\n",
							Sum(J, M, S, m, a, d));
						ShortenExp(szBuffer);
						strcat(szParams, szBuffer);
					}

					// Calc edf
					edf=(Sz(0, m, a, d)*Sz(0, m, a, d)*M/
						Sum(J, M, S, m, a, d));
				}
				else // m'=infinity, use F=m=0 as flag
				{
					if(RECORD)
					{
						// Record test
						strcat(szParams, "    m*(d+1)>Jmax\r\n");
						strcat(szParams, "    m'=infinity (F=0)\r\n");

						// Record Sz()
						sprintf(szBuffer, "    Sz(0,0,a,d)=%e\r\n",
							Sz(0, 0, a, d));
						ShortenExp(szBuffer);
						strcat(szParams, szBuffer);

						// Record Sum()
						sprintf(szBuffer, "    Sum(J,M,S,0,a,d)=%e\r\n",
							Sum(J, M, S, 0, a, d));
						ShortenExp(szBuffer);
						strcat(szParams, szBuffer);
					}

					edf=(Sz(0, 0, a, d)*Sz(0, 0, a, d)*M/
						Sum(J, M, S, 0, a, d));
				}
			}
			else // J>Jmax
			{
				if(RECORD)
				{
					// Record test
					strcat(szParams, "    J>Jmax\r\n");
				}

				if(r>=d+1)
				{
					if(RECORD)
					{
						// Record test
						strcat(szParams, "    r>=d+1\r\n");
					}

					// Get a0 & a1 from Table 2
					a0=T2[2-a][2*(d-1)+0];
					a1=T2[2-a][2*(d-1)+1];
					
					if(RECORD)
					{
						// Record a0 & a1
						strcat(szParams, "    a0,a1 from Table 2\r\n");
						sprintf(szBuffer, "    a0=%e, a1=%e\r\n", a0, a1);
						ShortenExp(szBuffer);
						strcat(szParams, szBuffer);
					}

					// Calc edf
					edf=((1.0/r)*(a0-(a1/r)));
					edf=(1.0/edf);
				}
				else // r<d+1
				{
					if(RECORD)
					{
						// Record test
						strcat(szParams, "    r<d+1\r\n");

						// Record Sz()
						sprintf(szBuffer, "    Sz(0,0,a,d)=%e\r\n",
							Sz(0, 0, a, d));
						ShortenExp(szBuffer);
						strcat(szParams, szBuffer);

						// Record Sum()
						sprintf(szBuffer,
							"    Sum(Jmax,Jmax,Jmax/r,0,a,d)=%e\r\n",
							Sum(Jmax, Jmax, (double)Jmax/r, 0, a, d));
						ShortenExp(szBuffer);
						strcat(szParams, szBuffer);
					}
					
					// Calc edf
					edf=(Sz(0, 0, a, d)*Sz(0, 0, a, d)*Jmax/
						Sum(Jmax, Jmax, (double)Jmax/r, 0, a, d));
				}
			}
			return (float)edf; // EDF for Case 2
		}
		else if(a==1) // Case 3. F PM
		{
			if(RECORD)
			{
				// Record case #
				strcat(szParams, "    Case 3\r\n");
			}

			if(J<=Jmax)
			{
				if(RECORD)
				{
					// Record test
					strcat(szParams, "    J<=Jmax\r\n");

					// Record Sz()
					sprintf(szBuffer, "    Sz(0,m,1,d)=%e\r\n",
						Sz(0, m, 1, d));
					ShortenExp(szBuffer);
					strcat(szParams, szBuffer);

					// Record Sum()
					sprintf(szBuffer, "    Sum(J,M,S,m,1,d)=%e\r\n",
						Sum(J, M, S, m, 1, d));
					ShortenExp(szBuffer);
					strcat(szParams, szBuffer);
				}

				// Calc edf
				// Note: m must be <1e6 to avoid roundoff error
				edf=(Sz(0, m, 1, d)*Sz(0, m, 1, d)*M/
					Sum(J, M, S, m, 1, d));
			}
			else
			{
				if(RECORD)
				{
					// Record test
					strcat(szParams, "    J>Jmax\r\n");
				}

				if(r>=d+1)
				{
					if(RECORD)
					{
						// Record test
						strcat(szParams, "    r>=d+1\r\n");
					}

					// Get a0 & a1 from Table 2
					a0=T2[2-a][2*(d-1)+0];
					a1=T2[2-a][2*(d-1)+1];

					// Get b0 & b1 from Table 3
					b0=T3[2*(d-1)+0];
					b1=T3[2*(d-1)+1];

					if(RECORD)
					{
						// Record a0 & a1
						strcat(szParams, "    a0,a1 from Table 2\r\n");
						sprintf(szBuffer, "    a0=%e, a1=%e\r\n", a0, a1);
						ShortenExp(szBuffer);
						strcat(szParams, szBuffer);

						// Record b0 & b1
						strcat(szParams, "    b0,b1 from Table 3\r\n");
						sprintf(szBuffer, "    b0=%e, b1=%e\r\n", b0, b1);
						ShortenExp(szBuffer);
						strcat(szParams, szBuffer);
					}

					// Calc edf
					edf=(((b0+b1*Log(m))*(b0+b1*Log(m))*r)/
						(a0-(a1/r)));
				}
				else
				{
					if(RECORD)
					{
						// Record test
						strcat(szParams, "    r<d+1\r\n");
					}

					// Get b0 & b1 from Table 3
					b0=T3[2*(d-1)+0];
					b1=T3[2*(d-1)+1];

					if(RECORD)
					{
						// Record b0 & b1
						strcat(szParams, "    b0,b1 from Table 3\r\n");
						sprintf(szBuffer, "    b0=%e, b1=%e\r\n", b0, b1);
						ShortenExp(szBuffer);
						strcat(szParams, szBuffer);

						// Record Sum()
						sprintf(szBuffer,
							"    Sum(Jmax,Jmax,Jmax/r,Jmax/r,1,d)=%e\r\n",
							Sum(Jmax, Jmax, (double)Jmax/r, (double)Jmax/r, 1, d));
						ShortenExp(szBuffer);
						strcat(szParams, szBuffer);
					}

					// Calc edf
					edf=(((b0+b1*Log(m))*(b0+b1*Log(m))*Jmax)/
						Sum(Jmax, Jmax, (double)Jmax/r, (double)Jmax/r, 1, d));
				}
			}
			return (float)edf; // EDF for Case 3
		}
		else // Case 4. W PM a=2
		{
			K=(int)ceil(r);

			if(RECORD)
			{
				// Record case #
				strcat(szParams, "    Case 4\r\n");

				// Record K
				sprintf(szBuffer, "    K=%d\r\n", K);
				strcat(szParams, szBuffer);
			}

			if(K<=d)
			{
				if(RECORD)
				{
					// Record test
					strcat(szParams, "    K<=d\r\n");
				}
				
				// Use a0 and a1 as working variables
				a0=bico(2*d, d);

				if(RECORD)
				{
					// Record a0
					strcat(szParams, "    a0,a1=binominial coeffs\r\n");
					sprintf(szBuffer, "    a0=%e, a1 inside loop\r\n", a0);
					ShortenExp(szBuffer);
					strcat(szParams, szBuffer);
				}
				
				// Calc sum
				edf=0.0;
				for(k=1; k<=K-1; k++)
				{
					a1=bico(2*d, d-k);
					// TRACE_EXP(a1);

					edf+=((1.0-((double)k/r))*a1*a1);
				}

				// Complete edf calc
				edf*=(2.0/(a0*a0));
				edf+=1.0;
				edf=M/edf;
			}
			else
			{
				if(RECORD)
				{
					// Record test
					strcat(szParams, "    K>d\r\n");
				}

				// Get a0 & a1 for a=2 from Table 2
				a0=T2[0][2*(d-1)+0];
				a1=T2[0][2*(d-1)+1];
				
				if(RECORD)
				{
					// Record a0 & a1
					strcat(szParams, "    a0,a1 from Table 2\r\n");
					sprintf(szBuffer, "    a0=%e, a1=%e\r\n", a0, a1);
					ShortenExp(szBuffer);
					strcat(szParams, szBuffer);
				}
				
				// Calc edf
				edf=(a0-(a1/r));
				edf=M/edf;
			}
			return (float)edf; // EDF for Case 4
		}
	}
}

// Supporting functions
double Sw(double t, int a)
{
	// Local variables
	double b;
	double sw;

	// Calc abs value
	b=fabs(t);
	
	// Calc sw
	switch(a)
	{
		case 2:
		{
			sw=-b;
		}
		break;

		case 1:
		{
			sw=t*t*Log(b);
		}
		break;

		case 0:
		{
			sw=b*b*b;
		}
		break;

		case -1:
		{
			sw=-t*t*t*t*Log(b);
		}
		break;

		case -2:
		{
			sw=-b*b*b*b*b;
		}
		break;

		case -3:
		{
			sw=t*t*t*t*t*t*Log(b);
		}
		break;

		case -4:
		{
			sw=b*b*b*b*b*b*b;
		}
		break;
	}
	return sw;
}

double Sx(double t, double F, int a)
{
	// Local variables
	double sx;

	// Calc sx
	if(F>0)
	{
		sx=F*F*(2*Sw(t, a) - Sw(t-1.0/F, a) - Sw(t+1.0/F, a));
	}
	else // F=0 is flag for F=infinity
	{
		// Note: This case applies only for a<=0
		sx=Sw(t, a+2);
	}
	return sx;
}

double Sz(double t, double F, int a, int d)
{
	// Local variables
	double sz;

	// Calc sz
	switch(d)
	{
		case 1:
		{
			sz=2*Sx(t, F, a) - Sx(t-1, F, a) - Sx(t+1, F, a);
		}
		break;
	
		case 2:
		{
			sz=6*Sx(t, F, a) - 4*Sx(t-1, F, a) - 4*Sx(t+1, F, a) +
				Sx(t-2,  F, a) + Sx(t+2, F, a);
		}
		break;

		case 3:
		{
			sz=20*Sx(t, F, a) - 15*Sx(t-1, F, a) - 15*Sx(t+1, F, a) +
				6*Sx(t-2, F, a) + 6*Sx(t+2, F, a) - Sx(t-3, F, a) -
				Sx(t+3, F, a);
		}
		break;
	}
	return sz;
}

double Sum(int J, double M, double S, double F, int a, int d)
{
	// Local variables
	int j;
	double sum;
	double z;

	// Initialize sum
	sum=Sz(0, F, a, d);
	sum*=sum;
	
	// Calc sum
	for(j=1; j<=J-1; j++)
	{
		z=Sz(j/S, F, a, d);
		z*=z;
		z*=(1-j/M);
		sum+=2*z;
	}

	// Add last term
	z=Sz(J/S, F, a, d);
	z*=z;
	z*=(1-J/M);
	sum+=z;

	return sum;
}

// NR functions to find binomial coefficients
double bico(int n, int k)
{
	return floor(0.5+exp(factln(n)-factln(k)-factln(n-k)));
}

double factln(int n)
{
	static double a[101];

	if (n < 0) return -1.0; // Error - negative factorial
	if (n <= 1) return 0.0;
	if (n <= 100) return a[n] ? a[n] : (a[n]=gammln(n+1.0));
	else return gammln(n+1.0);
}

double gammln(double xx)
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}

/*****************************************************************************/
