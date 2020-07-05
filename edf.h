/*****************************************************************************/
/*                                                                           */
/*                             edf.h                                         */
/*                                                                           */
/*      Simple header file for edf program                                   */
/*                                                                           */
/*      Revision record:                                                     */
/*          05/20/00    Created from hedf.h                                  */
/*                                                                           */
/*      (c) 2000  W. Riley Hamilton Technical Services All Rights Reserved   */
/*                                                                           */
/*****************************************************************************/

// Include Windows header files
#include <windows.h>
#include <windowsx.h>
#include <commdlg.h>
#include <commctrl.h>

// Include standard C headers
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <commdlg.h>
#include <direct.h>
#include <math.h>
#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <dos.h>
#include <limits.h>

// BOOL macro
#define BOOL int

// Debugging header and macros
// Note: Define NDEBUG to disable asserts
//#define NDEBUG
#include <assert.h>

// Assert macro Ref: DDJ Feb 98 p.10
// Does nothing if x is TRUE or if NDEBUG is defined
// If x is FALSE, execution immediately stops and debugger opens at ASSERT line
#ifdef NDEBUG
    #define ASSERT(x)
#else
    #define ASSERT(x) if(x);else{*((void**)0)=0;}
#endif

// Macro to simplify getting handle to dialog box control
// Note: 1st arg of dialog box proc MUST be named hDlg
#define GETHCTL(id) GetDlgItem(hDlg, (id))

// Trace, debug & display macros
// Require global variables char szTrace[BUFFER_SIZE+1} and HWND hwndMain

// Trace int
// Usage  TRACE_INT(nTrace); where nTrace is int variable to be displayed
#define TRACE_INT(var) sprintf(szTrace, #var "=%d", var); \
    MessageBox(hwndMain, szTrace, "Trace Int", MB_OK)

// Debug int
// Usage  DEBUG_INT(nTrace); where nTrace is int variable to be displayed
#define DEBUG_INT(var) sprintf(szTrace, "\n" #var "=%d", var); \
    OutputDebugString(szTrace)

// Display int
// Usage  DISPLAY_INT(nTrace); where nTrace is int variable to be displayed
#define DISPLAY_INT(var) sprintf(szTrace, "int " #var "=%d", var); \
    DebugMessageBox(hwndMain, __FILE__, __LINE__, szTrace)

// Trace unsigned int
// Usage  TRACE_UINT(uTrace); where uTrace is unsigned int variable to disp
#define TRACE_UINT(var) sprintf(szTrace, #var "=%u", var); \
    MessageBox(hwndMain, szTrace, "Trace Unsigned Int", MB_OK)

// Debug unsigned int
// Usage  DEBUG_UINT(uTrace); where uTrace is uint variable to be displayed
#define DEBUG_UINT(var) sprintf(szTrace, "\n" #var "=%u", var); \
    OutputDebugString(szTrace)

// Display unsigned int
// Usage  DISPLAY_UINT(uTrace); where uTrace is uint variable to be displayed
#define DISPLAY_UINT(var) sprintf(szTrace, "uint " #var "=%u", var); \
    DebugMessageBox(hwndMain, __FILE__, __LINE__, szTrace)

// Trace long
// Usage  TRACE_LONG(lnTrace); where lnTrace is long int to be displayed
#define TRACE_LONG(var) sprintf(szTrace,#var "=%ld",var); \
    MessageBox(hwndMain, szTrace, "Trace Long", MB_OK)

// Debug long
// Usage  DEBUG_LONG(lnTrace); where lnTrace is long variable to be displayed
#define DEBUG_LONG(var) sprintf(szTrace, "\n" #var "=%ld", var); \
    OutputDebugString(szTrace)

// Display long
// Usage  DISPLAY_LONG(lnTrace); where lnTrace is long variable to be displayed
#define DISPLAY_LONG(var) sprintf(szTrace, "long " #var "=%ld", var); \
    DebugMessageBox(hwndMain, __FILE__, __LINE__, szTrace)

// Trace real (float or double)
// Usage  TRACE_EXP(fTrace); where fTrace is float or double to be displayed
#define TRACE_EXP(var) sprintf(szTrace,#var "=%e",var); \
    MessageBox(hwndMain, szTrace, "Trace Exp", MB_OK)

// Debug real
// Usage  DEBUG_EXP(fTrace); where fTrace is real variable to be displayed
#define DEBUG_EXP(var) sprintf(szTrace, "\n" #var "=%e", var); \
    OutputDebugString(szTrace)

// Display real
// Usage  DISPLAY_EXP(fTrace); where fTrace is real variable to be displayed
#define DISPLAY_EXP(var) sprintf(szTrace, "real " #var "=%e", var); \
    DebugMessageBox(hwndMain, __FILE__, __LINE__, szTrace)

// Trace address (pointer)
// Usage  TRACE_ADR(pTrace); where pTrace is pointer to be displayed
#define TRACE_ADR(var) sprintf(szTrace, #var "=%p", var); \
    MessageBox(hwndMain, szTrace, "Trace Adr", MB_OK)

// Debug pointer
// Usage  DEBUG_ADR(pTrace); where pTrace is pointer variable to be displayed
#define DEBUG_ADR(var) sprintf(szTrace, "\n" #var "=%p", var); \
    OutputDebugString(szTrace)

// Display int
// Usage  DISPLAY_ADR(pTrace); where pTrace is pointer variable to be displayed
#define DISPLAY_ADR(var) sprintf(szTrace, "pointer " #var "=%p", var); \
    DebugMessageBox(hwndMain, __FILE__, __LINE__, szTrace)

// Trace string
// Usage  TRACE_STR(szTrace); where szTrace is string to be displayed
#define TRACE_STR(var) MessageBox(hwndMain, var, "Trace Str", MB_OK)

// Debug string
// Usage  DEBUG_STR(szTrace); where szTrace is string variable to be displayed
#define DEBUG_STR(var) OutputDebugString(var)

// Display string
// Usage  DISPLAY_STR(szTrace); where szTrace is string var to be displayed
#define DISPLAY_STR(var) DebugMessageBox(hwndMain, __FILE__, __LINE__, var)

// Error string
// Usage  ERROR_STR(szTrace); where szTrace is error message to be displayed
#define ERROR_STR(var) MessageBox(hwndMain, var, "Simple Error", \
	MB_ICONEXCLAMATION)

// Data type
#define F_TYPE double

// Buffer sizes
#define BUFFER_SIZE          128
#define FILENAME_SIZE        256
#define BIG_BUF_SIZE        1024

// Global variables
extern HWND hwndMain;						// Main window handle
extern HINSTANCE hInst;						// Application instance
extern HICON hIcon;							// Application icon
extern HCURSOR hWait;						// Handle to hourglass cursor
extern HCURSOR hArrow;						// Handle to normal cursor
extern char szTrace[BUFFER_SIZE+1];			// Buffer for trace macro
extern char szParams[BIG_BUF_SIZE+1];	    // Calc parameters

// Variance types
#define ALLAN_VARIANCE          0
#define MODALLAN_VARIANCE       1
#define HADAMARD_VARIANCE       2
#define TOTAL_VARIANCE          3
#define MODTOTAL_VARIANCE		4
#define HADTOTAL_VARIANCE		5
#define THEO1_VARIANCE		    6

// Function prototypes
// edf.c - Main function
int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance,
	PSTR szCmdLine, int nCmdShow);

// edf2.c - Dialog box procedures
BOOL CALLBACK EDFDlgProc(HWND hDlg, UINT msg, WPARAM wParam, LPARAM lParam);
BOOL CALLBACK AboutDlgProc(HWND hDlg, UINT message, WPARAM wParam,
	LPARAM lParam);
void WINAPI ShortenExp( char *szString);
BOOL WINAPI ValidateInteger(HWND hwnd, int nVal, int nMin,
	int nMax);

// edf4.c - Calculation functions
float WINAPI HadamardEDF(int N, int m, int b);
double WINAPI Rx(double t, int b);
double WINAPI Log(double x);
float WINAPI ModTotvarEDF(int nAlpha, float fRatio);
float WINAPI EDF(int nBeta, int nAvgFactor, int nNum);
float WINAPI CalcDegFree(int alpha, int n, int m);
float WINAPI TotvarEDF(int nAlpha, float fRatio);
float WINAPI Theo1EDF(int nAlpha, int nNum,	float fRatio);
float WINAPI HadTotvarEDF(int nAlpha, float fRatio);
float WINAPI CombinedEDF(int N, int m, int a, int d, int F, int S, int v);
