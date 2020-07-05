/*****************************************************************************/
/*                                                                           */
/*   edf2.c                                                                  */
/*                                                                           */
/*   Dialog box procedure for edf program                                    */
/*                                                                           */
/*   (c) 2000  W. Riley  Hamilton Technical Services  All Rights Reserved    */
/*                                                                           */
/*****************************************************************************/

// Includes
#include "edf.h"
#include "resource.h"

/*****************************************************************************/
/*                                                                           */
/*                   EDF Dialog Procedure                                    */
/*                                                                           */
/*      Revision record:                                                     */
/*          05/20/00    Created from HEDFDlgProc()                           */
/*			05/21/00	Revised # and AF entry tests                         */
/*			05/01/03	Added Thêo1                                          */
/*			06/14/03	Added controls for new Greenhall edf method          */
/*			06/27/03	Added View button                                    */
/*			06/28/03	Limited noise types to allowable values              */
/*						Added calc hourglass                                 */
/*			07/07/03	Added Hadamard total EDF                             */
/*                                                                           */
/*****************************************************************************/

BOOL CALLBACK EDFDlgProc(HWND hDlg, UINT msg, WPARAM wParam, LPARAM lParam)
{
    // Local macros
	#define OLD 0					// Classic edf formulae
	#define NEW 1					// Greenhall 2003 unified edf formulae
	#define SIMPLE 0				// Simple EDF calc version
	#define FULL 1					// Full EDF calc version
	#define TAU0 0					// Short (tau0) stride
	#define TAU 1					// Long (tau) stride 

	// Local variables
	char szBuffer[BUFFER_SIZE+1];   // Text buffer
	int nDummy;						// To avoid complier warning
	static nType=0;					// Variance type
	static int nNum=100;			// # phase data points
	static int nAF=1;				// Averaging factor
	static int nBeta;				// Beta of power law noise type
	static int nEntry;				// User entry to be validated
	static int nPrev;				// Previous value to be restored
	static int nMax=0x7FFFFFFF;		// Max allowable # data points (2^31)-1
	static int nMaxAF=1;			// Max allowable averaging factor
	static int nMinNum=1;			// Min allowable # data points
	static int nStride=1;			// S=AF (short) or 1 (long)
	static BOOL bErrorTrap=FALSE;	// Flag to handle entry errors
	static BOOL bType=OLD;			// EDF type (old or new)
	static BOOL bVersion=FULL;		// EDF calc version (full or simple)
	static BOOL bStride=TAU0;		// Stride type (long or short)
	static BOOL bCopy=FALSE;		// Flag that results copied to clipboard
	static float fEDF=0.0;			// EDF

	// Clipboard variables
	HGLOBAL hClipboard;				// Handle to clipboard text buffer
	LPSTR	szClipboard;			// Ptr to clipboard text string

    // Dialog box message switch statement
    switch(msg)
    {
        case WM_INITDIALOG:
        {
			// Set the application system menu icon
			SetClassLong(hDlg, GCL_HICON, (LONG)hIcon);

			// Variance types are as follows:
			// ALLAN_VARIANCE           0 (Overlapping)
			// MODALLAN_VARIANCE        1 (and Time)
			// HADAMARD_VARIANCE        2 (Overlapping)
			// TOTAL_VARIANCE           3 (and Time Total)
			// MODTOTAL_VARIANCE		4
			// HADTOTAL_VARIANCE		5
			// THEO1_VARIANCE			6

			// Load variance strings into the combo box
			// Note: Index is variance type
			nDummy=ComboBox_InsertString(GETHCTL(IDC_EDF_VARIANCE), 0,
					"Allan");
			nDummy=ComboBox_InsertString(GETHCTL(IDC_EDF_VARIANCE), 1,
					"Modified");
			nDummy=ComboBox_InsertString(GETHCTL(IDC_EDF_VARIANCE), 2,
					"Hadamard");
			nDummy=ComboBox_InsertString(GETHCTL(IDC_EDF_VARIANCE), 3,
					"Total");
			nDummy=ComboBox_InsertString(GETHCTL(IDC_EDF_VARIANCE), 4,
					"Mod Total");
			nDummy=ComboBox_InsertString(GETHCTL(IDC_EDF_VARIANCE), 5,
					"Had Total");
			nDummy=ComboBox_InsertString(GETHCTL(IDC_EDF_VARIANCE), 6,
					"Thêo1");

			// Select Allan variance
			nDummy=ComboBox_SetCurSel(GETHCTL(IDC_EDF_VARIANCE), 0);

			// Load noise strings into the combo box
			// Note: Read noise type as nBeta=index-4
			// FM FM and RR FM apply only to Hadamard variances
			// and are added to list if one of those variances is selected
			nDummy=ComboBox_InsertString(GETHCTL(IDC_EDF_NOISE), 0,
					"RW FM");
			nDummy=ComboBox_InsertString(GETHCTL(IDC_EDF_NOISE), 1,
					"FL FM");
			nDummy=ComboBox_InsertString(GETHCTL(IDC_EDF_NOISE), 2,
					"WH FM");
			nDummy=ComboBox_InsertString(GETHCTL(IDC_EDF_NOISE), 3,
					"FL PM");
			nDummy=ComboBox_InsertString(GETHCTL(IDC_EDF_NOISE), 4,
					"WH PM");

			// Select white FM noise
			nDummy=ComboBox_SetCurSel(GETHCTL(IDC_EDF_NOISE), 2);

			// Show the # phase data points value
			sprintf(szBuffer, "%d", nNum);
			SetWindowText(GETHCTL(IDC_EDF_NUM), szBuffer);

			// Show the averaging factor value
			sprintf(szBuffer, "%d", nAF);
			SetWindowText(GETHCTL(IDC_EDF_AF), szBuffer);

			// Show the edf type
			SendMessage(GETHCTL(IDC_EDF_OLDTYPE), BM_SETCHECK, !bType, 0);
			SendMessage(GETHCTL(IDC_EDF_FULLTYPE), BM_SETCHECK, bVersion && bType,
				0);
			SendMessage(GETHCTL(IDC_EDF_SIMPLETYPE), BM_SETCHECK, !bVersion && bType,
				0);

			// Show the stride type
			SendMessage(GETHCTL(IDC_EDF_SHORTSTRIDE), BM_SETCHECK, !bStride, 0);
			SendMessage(GETHCTL(IDC_EDF_LONGSTRIDE), BM_SETCHECK, bStride, 0);

			// Disable long stride selection if old edf type
			if(!bType)
			{
				EnableWindow(GETHCTL(IDC_EDF_LONGSTRIDE), FALSE);
			}

			// Disable the Copy button until calc is done
			EnableWindow(GETHCTL(IDC_EDF_COPY), FALSE);

			// Disable the View button until copy is done
			EnableWindow(GETHCTL(IDC_EDF_VIEW), FALSE);

			return TRUE;
        }

        case WM_COMMAND:
        {
            switch(GET_WM_COMMAND_ID(wParam, lParam))
            {
				// Respond to About button
				// Note: Could make this a real Help button
				case IDC_EDF_HELP:
                {
					// Create about dialog box
					// CreateDialog(hInst, "AboutBox", NULL, AboutDlgProc);

					// Invoke help
                    WinHelp(hDlg, "EDF.HLP", HELP_CONTENTS, 0L);
					
					return 0;
                }

                // Respond to change in variance type
                case IDC_EDF_VARIANCE:
                {
                    // Get new variance type
					nEntry=ComboBox_GetCurSel(GETHCTL(IDC_EDF_VARIANCE));

					// Has variance changed?
					if(nEntry!=nType)
					{
						// Disable the copy button until there are new results
						EnableWindow(GETHCTL(IDC_EDF_COPY), FALSE);

						// Disable the view button until new results copied
						EnableWindow(GETHCTL(IDC_EDF_VIEW), FALSE);

						// Clear all the results boxes
						SetWindowText(GETHCTL(IDC_EDF_EDF), "");

						// Clear the EDF value
						fEDF=0.0;

						// Set the averaging factor to 1
						nAF=1;

						// Set AF to 10 for Thêo1
						if(nEntry==THEO1_VARIANCE)
						{
							nAF=10;
						}

						// Show the averaging factor value
						sprintf(szBuffer, "%d", nAF);
						SetWindowText(GETHCTL(IDC_EDF_AF), szBuffer);

						// Save the new variance type
						nType=nEntry;

						// Clear noise type list
						nDummy=ComboBox_ResetContent(GETHCTL(IDC_EDF_NOISE));

						if(nType==HADAMARD_VARIANCE || nType==HADTOTAL_VARIANCE)
						{
							// Load noise strings into the combo box
							// Note: Read noise type as nBeta=index-6
							nDummy=ComboBox_InsertString(GETHCTL(IDC_EDF_NOISE),
								0, "RR FM");
							nDummy=ComboBox_InsertString(GETHCTL(IDC_EDF_NOISE),
								1, "FW FM");
							nDummy=ComboBox_InsertString(GETHCTL(IDC_EDF_NOISE),
								2, "RW FM");
							nDummy=ComboBox_InsertString(GETHCTL(IDC_EDF_NOISE),
								3, "FL FM");
							nDummy=ComboBox_InsertString(GETHCTL(IDC_EDF_NOISE),
								4, "WH FM");
							nDummy=ComboBox_InsertString(GETHCTL(IDC_EDF_NOISE),
								5, "FL PM");
							nDummy=ComboBox_InsertString(GETHCTL(IDC_EDF_NOISE),
								6, "WH PM");

							// Select white FM noise
							nDummy=ComboBox_SetCurSel(GETHCTL(IDC_EDF_NOISE), 4);
						}
						else // Non-Hadamard variances
						{
							// Load noise strings into the combo box
							// Note: Read noise type as nBeta=index-4
							nDummy=ComboBox_InsertString(GETHCTL(IDC_EDF_NOISE),
								0, "RW FM");
							nDummy=ComboBox_InsertString(GETHCTL(IDC_EDF_NOISE),
								1, "FL FM");
							nDummy=ComboBox_InsertString(GETHCTL(IDC_EDF_NOISE),
								2, "WH FM");
							nDummy=ComboBox_InsertString(GETHCTL(IDC_EDF_NOISE),
								3, "FL PM");
							nDummy=ComboBox_InsertString(GETHCTL(IDC_EDF_NOISE),
								4, "WH PM");

							// Select white FM noise
							nDummy=ComboBox_SetCurSel(GETHCTL(IDC_EDF_NOISE), 2);
						}
					}
					return 0;
				}

                // Respond to change in noise type
                case IDC_EDF_NOISE:
                {
                    // Get new beta
					if(nType==HADAMARD_VARIANCE || nType==HADAMARD_VARIANCE)
					{
						nEntry=ComboBox_GetCurSel(GETHCTL(IDC_EDF_NOISE))-6;
					}
					else // Non-Hadamard variances
					{
						nEntry=ComboBox_GetCurSel(GETHCTL(IDC_EDF_NOISE))-4;
					}

					// Has beta changed?
					if(nEntry!=nBeta)
					{
						// Disable the copy button until there are new results
						EnableWindow(GETHCTL(IDC_EDF_COPY), FALSE);

						// Disable the view button until new results copied
						EnableWindow(GETHCTL(IDC_EDF_VIEW), FALSE);

						// Clear all the results boxes
						SetWindowText(GETHCTL(IDC_EDF_EDF), "");

						// Clear the EDF value
						fEDF=0.0;

						// Save the new beta
						nBeta=nEntry;
					}

					return 0;
				}
 				
				// Respond to change in # data points
                case IDC_EDF_NUM:
                {
                    if(GET_WM_COMMAND_CMD(wParam, lParam)==EN_CHANGE)
                    {
						// Disable the copy button until there are new results
						EnableWindow(GETHCTL(IDC_EDF_COPY), FALSE);

						// Disable the view button until new results copied
						EnableWindow(GETHCTL(IDC_EDF_VIEW), FALSE);

                        // Clear all the results boxes
                        SetWindowText(GETHCTL(IDC_EDF_EDF), "");

						// Clear the EDF value
						fEDF=0.0;
					}

                    // Respond to edit control when it loses the input focus
					// unless it loses it to the cancel button
                    if((GET_WM_COMMAND_CMD(wParam, lParam)==EN_KILLFOCUS) &&
                        (!bErrorTrap) && (GetFocus()!=GETHCTL(IDC_EDF_CLOSE)))
					{
						// Set error flag to avoid multiple error messages
                        bErrorTrap=TRUE;

                        // Save previous value
                        nPrev=nNum;

                        // Get new # data points
						GetWindowText(GETHCTL(IDC_EDF_NUM), szBuffer,
							BUFFER_SIZE);

						// Bound size of num
						if(strlen(szBuffer)<=9)
						{
							nEntry=abs(atoi(szBuffer));
						}
						else
						{
							nEntry=999999999;
						}

						// Check entry - must be between 2*nAF and nMax for
						// AVAR nd TOTVAR, and between 3*nAF and nMax for
						// MVAR, HVAR and MTOT

						// Find nMinNum
						if(nType==ALLAN_VARIANCE || nType==TOTAL_VARIANCE)
						{
							nMinNum=2*nAF;
						}
						else if(nType==THEO1_VARIANCE)
						{
							nMinNum=nAF;
						}
						else
						{
							nMinNum=3*nAF;
						}

                        if(ValidateInteger(hDlg, nEntry, nMinNum, nMax)==TRUE)
                        {
                            // New value OK
                            nNum=nEntry;
                            sprintf(szBuffer, "%d", nNum);
							SetWindowText(GETHCTL(IDC_EDF_NUM), szBuffer);
						}
						else
						{
							// Invalid entry - restore previous value
							sprintf(szBuffer, "%d", nPrev);
							SetWindowText(GETHCTL(IDC_EDF_NUM), szBuffer);
							nNum=nPrev;
							SetFocus(GETHCTL(IDC_EDF_NUM));
						}
					}
					else
					{
						// Reset error flag
						bErrorTrap=FALSE;
					}
					return 0;
				}
				
				// Respond to change in averaging factor
                case IDC_EDF_AF:
                {
                    if(GET_WM_COMMAND_CMD(wParam, lParam)==EN_CHANGE)
                    {
						// Disable the copy button until there are new results
						EnableWindow(GETHCTL(IDC_EDF_COPY), FALSE);

						// Disable the view button until new results copied
						EnableWindow(GETHCTL(IDC_EDF_VIEW), FALSE);

                       // Clear the results box
                        SetWindowText(GETHCTL(IDC_EDF_EDF), "");

						// Clear the EDF value
						fEDF=0.0;
					}

                    // Respond to edit control when it loses the input focus
					// unless it loses it to the cancel button
                    if((GET_WM_COMMAND_CMD(wParam, lParam)==EN_KILLFOCUS) &&
                        (!bErrorTrap) && (GetFocus()!=GETHCTL(IDC_EDF_CLOSE)))
					{
						// Set error flag to avoid multiple error messages
                        bErrorTrap=TRUE;

                        // Save previous value
                        nPrev=nAF;

                        // Get new averaging factor
						GetWindowText(GETHCTL(IDC_EDF_AF), szBuffer,
							BUFFER_SIZE);
						nEntry=abs(atoi(szBuffer));

						// Check entry - must be between 1 and nNum/2
						// for AVAR and TOTVAR, and between 1 and nNum/3
						// for MVAR, HVAR & MTOT

						// Find nMaxAF
						if(nType==ALLAN_VARIANCE || nType==TOTAL_VARIANCE)
						{
							nMaxAF=nNum/2;
						}
						else if(nType==THEO1_VARIANCE)
						{
							nMaxAF=nNum;
						}
						else
						{
							nMaxAF=nNum/3;
						}

                        if(ValidateInteger(hDlg, nEntry, 1, nMaxAF)==TRUE)
                        {
                            // New value OK
                            nAF=nEntry;
                            sprintf(szBuffer, "%d", nAF);
							SetWindowText(GETHCTL(IDC_EDF_AF),	szBuffer);
						}
						else
						{
							// Invalid entry - restore previous value
							sprintf(szBuffer, "%d", nPrev);
							SetWindowText(GETHCTL(IDC_EDF_AF),	szBuffer);
							nAF=nPrev;
							SetFocus(GETHCTL(IDC_EDF_AF));
						}
					}
					else
					{
						// Reset error flag
						bErrorTrap=FALSE;
					}
					return 0;
				}
				
                case IDC_EDF_OLDTYPE:
                {
					// Set type to old
					bType=OLD;
					SendMessage(GETHCTL(IDC_EDF_OLDTYPE), BM_SETCHECK,
						TRUE, 0);
					SendMessage(GETHCTL(IDC_EDF_FULLTYPE), BM_SETCHECK,
						FALSE, 0);
					SendMessage(GETHCTL(IDC_EDF_SIMPLETYPE), BM_SETCHECK,
						FALSE, 0);

					// If bType is old, then bStride must be short
					SendMessage(hDlg, WM_COMMAND,
							GET_WM_COMMAND_MPS(IDC_EDF_SHORTSTRIDE, 0, 0));

					// Disable long stride selection
					EnableWindow(GETHCTL(IDC_EDF_LONGSTRIDE), FALSE);

					// Disable the copy button until there are new results
					EnableWindow(GETHCTL(IDC_EDF_COPY), FALSE);

					// Disable the view button until new results copied
					EnableWindow(GETHCTL(IDC_EDF_VIEW), FALSE);

                    // Clear the results box
                    SetWindowText(GETHCTL(IDC_EDF_EDF), "");

					return 0;
				}

                case IDC_EDF_FULLTYPE:
                {
					// Set type to new
					bType=NEW;

					// Set version to full
					bVersion=FULL;

					// Set radio buttons
					SendMessage(GETHCTL(IDC_EDF_OLDTYPE), BM_SETCHECK,
						!bType, 0);
					SendMessage(GETHCTL(IDC_EDF_FULLTYPE), BM_SETCHECK,
						bVersion && bType, 0);
					SendMessage(GETHCTL(IDC_EDF_SIMPLETYPE), BM_SETCHECK,
						!bVersion && bType, 0);

					// Enable long stride selection
					EnableWindow(GETHCTL(IDC_EDF_LONGSTRIDE), TRUE);

					// Disable the copy button until there are new results
					EnableWindow(GETHCTL(IDC_EDF_COPY), FALSE);

					// Disable the view button until new results copied
					EnableWindow(GETHCTL(IDC_EDF_VIEW), FALSE);

					// Clear the results box
                    SetWindowText(GETHCTL(IDC_EDF_EDF), "");

					return 0;
				}

                case IDC_EDF_SIMPLETYPE:
                {
					// Set type to new
					bType=NEW;

					// Set version to simple
					bVersion=SIMPLE;

					SendMessage(GETHCTL(IDC_EDF_OLDTYPE), BM_SETCHECK,
						!bType, 0);
					SendMessage(GETHCTL(IDC_EDF_FULLTYPE), BM_SETCHECK,
						bVersion && bType, 0);
					SendMessage(GETHCTL(IDC_EDF_SIMPLETYPE), BM_SETCHECK,
						!bVersion && bType, 0);

					// Enable long stride selection
					EnableWindow(GETHCTL(IDC_EDF_LONGSTRIDE), TRUE);

					// Disable the copy button until there are new results
					EnableWindow(GETHCTL(IDC_EDF_COPY), FALSE);

					// Disable the view button until new results copied
					EnableWindow(GETHCTL(IDC_EDF_VIEW), FALSE);

                    // Clear the results box
                    SetWindowText(GETHCTL(IDC_EDF_EDF), "");

					return 0;
				}
                case IDC_EDF_LONGSTRIDE:
                {
					// Set stride to long
					bStride=TAU;
					SendMessage(GETHCTL(IDC_EDF_SHORTSTRIDE), BM_SETCHECK,
						!bStride, 0);
					SendMessage(GETHCTL(IDC_EDF_LONGSTRIDE), BM_SETCHECK,
						bStride, 0);

					// Disable the copy button until there are new results
					EnableWindow(GETHCTL(IDC_EDF_COPY), FALSE);

					// Disable the view button until new results copied
					EnableWindow(GETHCTL(IDC_EDF_VIEW), FALSE);

                   // Clear the results box
                    SetWindowText(GETHCTL(IDC_EDF_EDF), "");

					return 0;
				}

                case IDC_EDF_SHORTSTRIDE:
                {
					// Set stride to short
					bStride=TAU0;
					SendMessage(GETHCTL(IDC_EDF_SHORTSTRIDE), BM_SETCHECK,
						!bStride, 0);
					SendMessage(GETHCTL(IDC_EDF_LONGSTRIDE), BM_SETCHECK,
						bStride, 0);

					// Disable the copy button until there are new results
					EnableWindow(GETHCTL(IDC_EDF_COPY), FALSE);

					// Disable the view button until new results copied
					EnableWindow(GETHCTL(IDC_EDF_VIEW), FALSE);

                   // Clear the results box
                    SetWindowText(GETHCTL(IDC_EDF_EDF), "");

					return 0;
				}

				// Respond to Calc button
				case IDC_EDF_CALC:
                {
					// TRACE_INT(bType);
					// TRACE_INT(bStride);
					
					// Show the wait cursor
					SetCursor(hWait);

					// Clear the EDF value
					fEDF=0.0;

					// Clear the calc record
					strcpy(szParams, "");
					
                    // Get variance type
					nType=ComboBox_GetCurSel(GETHCTL(IDC_EDF_VARIANCE));

                    // Get new beta
					if(nType==HADAMARD_VARIANCE || nType==HADTOTAL_VARIANCE)
					{
						nBeta=ComboBox_GetCurSel(GETHCTL(IDC_EDF_NOISE))-6;
					}
					else // Non-Hadamard variances
					{
						nBeta=ComboBox_GetCurSel(GETHCTL(IDC_EDF_NOISE))-4;
					}

					// Old edf type
					if(bType==OLD)
					{
						// Parse calculation based on variance type
						switch(nType)
						{
							case ALLAN_VARIANCE:
							{
								// Calculate Allan edf
								fEDF=CalcDegFree(nBeta+2, nNum, nAF);
							}
							break;

							case MODALLAN_VARIANCE:
							{
								// Calculate Mod Allan edf
								fEDF=HadamardEDF(nNum+1, nAF, nBeta-2);
							}
							break;

							case HADAMARD_VARIANCE:
							{
								// Calculate Hadamard edf
								fEDF=HadamardEDF(nNum, nAF, nBeta);
							}
							break;

							case TOTAL_VARIANCE:
							{
								// Calculate Total edf
								fEDF=TotvarEDF(nBeta+2, (float)nNum/(float)nAF);
							}
							break;

							case MODTOTAL_VARIANCE:
							{
								// Calculate Mod Total edf
								fEDF=ModTotvarEDF(nBeta+2, (float)nNum/(float)nAF);
							}
							break;

							case HADTOTAL_VARIANCE:
							{
								// Calculate Hadamard Total edf
								fEDF=HadTotvarEDF(nBeta+2, (float)nNum/(float)nAF);
							}
							break;

							case THEO1_VARIANCE:
							{
								// Calculate Thêo1 edf
								fEDF=Theo1EDF(nBeta+2, nNum, (float)nAF);
							}
							break;
						}
					}

					// New edf type
					if(bType==NEW)
					{
						if(bStride==TAU)
						{
							nStride=1; // Long stride
						}
						else
						{
							nStride=nAF; // Short stride
						}

						// Parse calculation based on variance type
						switch(nType)
						{
							case ALLAN_VARIANCE:
							{
								// Calculate Allan edf
								fEDF=CombinedEDF(nNum, nAF, nBeta+2, 2, nAF,
									nStride, bVersion);
							}
							break;

							case MODALLAN_VARIANCE:
							{
								// Calculate Mod Allan edf
								fEDF=CombinedEDF(nNum, nAF, nBeta+2, 2, 1,
									nStride, bVersion);
							}
							break;

							case HADAMARD_VARIANCE:
							{
								// Calculate Hadamard edf
								fEDF=CombinedEDF(nNum, nAF, nBeta+2, 3, nAF,
									nStride, bVersion);
							}
							break;

							case TOTAL_VARIANCE:
							{
								// Calculate Total edf
								fEDF=TotvarEDF(nBeta+2, (float)nNum/(float)nAF);
							}
							break;

							case MODTOTAL_VARIANCE:
							{
								// Calculate Mod Total edf
								fEDF=ModTotvarEDF(nBeta+2, (float)nNum/(float)nAF);
							}
							break;

							case THEO1_VARIANCE:
							{
								// Calculate Thêo1 edf
								fEDF=Theo1EDF(nBeta+2, nNum, (float)nAF);
							}
							break;
						}
					}

					// Display  edf
					if(fEDF>0)
                    {
						sprintf(szBuffer, "%e", fEDF);
						ShortenExp(szBuffer);
						SetWindowText(GETHCTL(IDC_EDF_EDF), szBuffer);
					}
					else
                    {
						SetWindowText(GETHCTL(IDC_EDF_EDF), "Error");
					}
					
					// Enable Copy button
					EnableWindow(GETHCTL(IDC_EDF_COPY), TRUE);

					// Show the normal cursor
					SetCursor(hArrow);

					return 0;
                }

				// Respond to Copy button
                case IDC_EDF_COPY:
                {
					// Set up clipboard usage
					hClipboard=GlobalAlloc(GMEM_DDESHARE,
						(DWORD)(10*80));
					szClipboard=GlobalLock(hClipboard);

					// Heading
					strcpy(szClipboard, "EDF Calculation");
					strcat(szClipboard, "\r\n\r\n");

					// Calc type
					strcat(szClipboard, "    Calculation Type: ");
					if(bType==OLD)
					{
						strcpy(szBuffer, "Old");
					}
					else
					{
						if(bVersion==FULL)
						{
							strcpy(szBuffer, "Full");
						}
						else
						{
							strcpy(szBuffer, "Simple");
						}
					}
					strcat(szClipboard, szBuffer);
					strcat(szClipboard, "\r\n");

					// Variance type
					strcat(szClipboard, "    Variance Type: ");
					GetWindowText(GETHCTL(IDC_EDF_VARIANCE), szBuffer,
						BUFFER_SIZE);
					strcat(szClipboard, szBuffer);
					strcat(szClipboard, "\r\n");

					// Noise type
					strcat(szClipboard, "    Noise Type: ");
					GetWindowText(GETHCTL(IDC_EDF_NOISE), szBuffer,
						BUFFER_SIZE);
					strcat(szClipboard, szBuffer);
					strcat(szClipboard, "\r\n");

					// # data points
					sprintf(szBuffer, "    # Data Points=%d", nNum);
					strcat(szClipboard, szBuffer);
					strcat(szClipboard, "\r\n");

					// Avg factor
					sprintf(szBuffer, "    Avg Factor=%d", nAF);
					strcat(szClipboard, szBuffer);
					strcat(szClipboard, "\r\n");

					// Stride
					if(bStride==TAU0)
					{
						strcpy(szBuffer, "    Stride=Short");
					}
					else
					{
						strcpy(szBuffer, "    Stride=Long");
					}
					strcat(szClipboard, szBuffer);
					strcat(szClipboard, "\r\n");

					// Calc params
					strcat(szClipboard, szParams);
					strcat(szClipboard, "\r\n");
					
					// edf
					sprintf(szBuffer, "    EDF=%e", fEDF);
					ShortenExp(szBuffer);
					strcat(szClipboard, szBuffer);
					strcat(szClipboard, "\r\n");

					// Write info to clipboard
					GlobalUnlock(hClipboard);
						OpenClipboard(NULL);
						EmptyClipboard();
						SetClipboardData(CF_TEXT, hClipboard);
						CloseClipboard();

					// Enable View button
					EnableWindow(GETHCTL(IDC_EDF_VIEW), TRUE);

					return 0;
                }

				// View clipboard
				case IDC_EDF_VIEW:
				{
					// WinExec("CLIPBRD.EXE", SW_SHOWNORMAL);
					// WinExec("PTClpVue.exe", SW_SHOWNORMAL);
					WinExec("CLIPVIEW.EXE", SW_SHOWNORMAL);
					return 0;
                }

				// Respond to CR in an edit control
				// Note: No button is a DEFPUSHBUTTON
				case IDOK:
				{
					// Does avg factor edit control have focus?
					if(GetFocus()==GETHCTL(IDC_EDF_AF))
					{
						// Send EN_KILLFOCUS message to avg factor edit control
						SendMessage(hDlg, WM_COMMAND,
							GET_WM_COMMAND_MPS(IDC_EDF_AF, 0,
							EN_KILLFOCUS));
					}

					// Does # data points edit control have focus?
					if(GetFocus()==GETHCTL(IDC_EDF_NUM))
					{
						// Send EN_KILLFOCUS msg to # data points edit control
						SendMessage(hDlg, WM_COMMAND,
							GET_WM_COMMAND_MPS(IDC_EDF_NUM, 0, EN_KILLFOCUS));
					}
				}
				break;

				// Respond to Close button, etc.
                case IDCANCEL:
                case IDC_EDF_CLOSE:
                {
					// Close program
					PostQuitMessage(0);
					return 0;
                }
                return TRUE;
            }
        }
        break;
    }
    return FALSE;
}

/*****************************************************************************/
/*                                                                           */
/*                   About Dialog Procedure                                  */
/*                                                                           */
/*      Revision record:                                                     */
/*          05/09/99    Created                                              */
/*                                                                           */
/*****************************************************************************/

BOOL CALLBACK AboutDlgProc(HWND hDlg, UINT message, WPARAM wParam,
	LPARAM lParam)
{
	// Local variables
	#if(0) // Not used
	char szBuffer[BUFFER_SIZE+1];   // Text buffer
	#endif

	switch(message)
	{
		case WM_INITDIALOG:
		{
			// Load and display any special strings here
			return TRUE;
		}

		case WM_COMMAND:
		{
			switch(GET_WM_COMMAND_ID(wParam, lParam))
            {
                case IDOK:
                case IDCANCEL:
                {
                    EndDialog (hDlg, 0);
                }
                return TRUE;
			}
        }
        break;
     }
     return FALSE;
}

/*****************************************************************************/
/*                                                                           */
/*							ShortenExp()                                     */
/*                                                                           */
/*		Function to shorten all 3-digit exponents in a string                */
/*      Use with Microsoft Visual C++ 5.0 for legacy display layouts         */
/*      Changes all numeric value formats of passed numeric value string     */
/*		from 3 digits to 2 digits                                            */
/*                                                                           */
/*      Parameters:                                                          */
/*      	szString		Numeric value string                             */
/*                                                                           */
/*      Return:                                                              */
/*			None	(void)                                                   */
/*                                                                           */
/*		Revision record:                                                     */
/*			01/16/98	Created                                              */
/*                                                                           */
/*****************************************************************************/
void WINAPI ShortenExp( char *szString)
{
	// Local variables
    int i;		// Main index
	int j;		// Secondary index
	int len;	// String length

	// Find string length
	len=strlen(szString);
	
	// Scan for next 'e' in the string followed by a + or - and then a zero
	// Looking for pattern e+0dd or e-0dd
	// Scan stops 4 chars before end of string
	for(i=0; i<len-4; i++)
	{
		if( (szString[i]=='e') && (szString[i+1]=='+' || szString[i+1]=='-') &&
			(szString[i+2]=='0') )
		{
			// Original format is e+0dd or e-0dd
			// Want to eliminate the leading exp digit (always zero)
			// Modified format will be e+dd or e-dd
			// Shift up the string chars (and EOS) after the unwanted 0
			for(j=i+2; j<len; j++)
			{
				szString[j]=szString[j+1];
			}
			len--;
		}
	}
}

/*****************************************************************************/
/*                                                                           */
/*							ValidateInteger()                                */
/*                                                                           */
/*		Function to validate integer input versus allowable limits.          */
/*                                                                           */
/*      Parameters:                                                          */
/*      	HWND	hwnd		Handle to window calling function            */
/*          int		nVal		Integer value to be checked                  */
/*          int		nMin 		Minimum allowable value                      */
/*			int		nMax  	    Maximum allowable value                      */
/*                                                                           */
/*      Return:                                                              */
/*			BOOL	status		TRUE=Value OK, FALSE=Value NG                */
/*                                                                           */
/*      Note 1:	Calling function S/B able to restore variable to its         */
/*              previous or default value before the bad entry was made.     */
/*      Note 2:	Call in response to EN_KILLFOCUS message when using with     */
/*              an edit control.                                             */
/*      Note 3:	Don't call in response to pressing the Cancel button.        */
/*                                                                           */
/*****************************************************************************/
BOOL WINAPI ValidateInteger(HWND hwnd, int nVal, int nMin,
	int nMax)
{
	// Local variables
	char szBuffer[BUFFER_SIZE+1];			// Text buffer

	if(nVal<nMin || nVal>nMax)
	{
		MessageBeep(0);
		sprintf(szBuffer, "Value %d must be between %d and %d",
			nVal, nMin, nMax);
		MessageBox(hwnd, szBuffer, "Invalid Entry", MB_ICONEXCLAMATION);
		return FALSE;
	}
	else
	{
		return TRUE;
	}
}

/*****************************************************************************/
