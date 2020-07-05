/*****************************************************************************/
/*                                                                           */
/*	 edf.c                                                                   */
/*                                                                           */
/*   Simple Win32 program for testing edf calculations.                      */
/*	 Calculates the following edfs:                                          */
/*                                                                           */
/*	 Variance				Listbox Label	EDF Function	Alpha Range      */
/*	 Overlapping Allan		Allan			CalcDegFree()   -2 to 2          */
/*	 Modified Allan			Modified        EDF()           -2 to 2          */
/*   Overlapping Hadamard   Hadamard        HadamardEDF()   -4 to 0          */
/*   Total					Total           TotvarEDF()     -2 to 0          */
/*   Modified Total         ModTotal        ModTotvarEDF()  -2 to 2          */
/*	 Thêo1					Thêo1			Theo1EDF()		-2 to 0          */
/*                                                                           */
/*	 These edf calculations do not include the following Stable32 features:  */
/*   Stable32 uses AVAR edf for HVAR where HadamardEDF() doesn't apply.      */
/*   It uses AVAR edf+2 for TOTVAR where TotvarEDF() doesn't apply.          */
/*	 Also, MTOT edf switches from ModTotvarEDF() to EDF() for m<=8.          */
/*                                                                           */
/*   Changed 05/22/00 per inputs from C. Greenhall/JPL:                      */
/*   Replace call to EDF() with call to HadamardEDF() for MVAR               */
/*   Let mvar_edf(N, m, beta) = Hvar_edf(N+1, m, beta-2)                     */
/*   where N = # phase samples, m = averaging factor, beta = alpha -2        */
/*   Range of validity: 3m <= N+1, beta = 0 thru -4 (alpha +2 to -2)         */
/*   m <= (N+1)/3 covers the entire span of allowable averaging factors      */
/*   and -2 <= alpha <= 2 covers all the usual power-law noise types         */
/*                                                                           */
/*   Uses dialog box as main window.                                         */
/*   No menu, toolbuttons or other windows.                                  */
/*                                                                           */
/*   Major revision record:                                                  */
/*      05/20/00  1.0 Created from hedf.c                                    */
/*		05/21/00  1.1 Revised MTOT coefficients and entry tests              */
/*		05/22/00  1.2 Reset AF=1 when variance type changed                  */
/*					  Changed MVAR edf calc from EDF() to HadamardEDF()      */
/*		05/01/03  1.3 Added Thêo1                                            */
/*		06/08/03  1.4 Added Greenhall combined edf algorithm                 */
/*                                                                           */
/*   (c) 2000-3  W. Riley Hamilton Technical Services  All Rights Reserved   */
/*                                                                           */
/*****************************************************************************/

// Includes
#include "edf.h"
#include "resource.h"

// Global variables
// Note: All string allocations have extra character for EOS
HWND hwndMain;								// Main window handle
HINSTANCE hInst;							// Application instance
HICON hIcon;								// Application icon
HCURSOR hWait;								// Handle to hourglass cursor
HCURSOR hArrow;								// Handle to normal cursor
char szTrace[BUFFER_SIZE+1];				// Trace buffer
char szParams[BIG_BUF_SIZE+1];				// Calc parameters

// WinMain() function
int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance,
	PSTR szCmdLine, int nCmdShow)
{
	// Local variables
    MSG         msg;						// Main message loop message

	// Load cursors
	hWait=LoadCursor(NULL, IDC_WAIT);
	hArrow=LoadCursor(NULL, IDC_ARROW);

	// Get handle to application icon
	hIcon=LoadIcon(hInstance, "EDFIcon");

	// Set global instance handle
	hInst=hInstance;

	// Create dialog box to be used as main window
	hwndMain=CreateDialog(hInstance, "EDFBox", NULL, EDFDlgProc);
	
	// Display main window
    ShowWindow(hwndMain, nCmdShow);
    UpdateWindow(hwndMain);

	// Main message loop
    while(GetMessage(&msg, NULL, 0, 0))
    {
		// Is this a dialog box message?
		if(IsDialogMessage(hwndMain, &msg))
		{
			// If so, skip normal message handling
			continue;
		}
		
		// Otherwise process message as normal
		TranslateMessage(&msg);
        DispatchMessage(&msg);
    }
    return msg.wParam;
}

/*****************************************************************************/
