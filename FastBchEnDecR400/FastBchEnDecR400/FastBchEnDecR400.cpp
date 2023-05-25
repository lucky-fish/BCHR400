//
// FastBchEnDecR400.cpp  ----- SEE REVISION HISTORY BELOW -----
//
// COPYRIGHT STATEMENT: Copyright (C) Neal Glover 1999,2000,2008,2010
// (boglover@msn.com 303-466-5434).
// LICENSE GRANT: A free non-exclusive, perpetual, irrevocable source
// and object code license, with the following restrictions and warranty
// disclaimer and without support or maintenance, is granted for the
// use, modification, reproduction, distribution and sub-licensing of
// this code.
// RESTRICTIONS: All parts of this statement including copyright statement,
// license grant, and warranty disclaimer must remain part of the source code
// or parts thereof that you use.  All source code comments associated with
// the parts of the code that you use must remain in the code.  The license
// grant does not include the right to sell this code or part thereof as a separate
// product.
// WARRANTY DISCLAIMER: The author of this code does not warrant
// that it is free of defects or that it is useful or fit for any purpose.
// The author does not warrant that this code will not violate
// intellectual property rights including patent rights or copyright rights
// of third parties.  The author assumes no risk for your use
// of this code.  If you use this code or parts thereof, it is at your own risk.
// By using this code, you agree to hold the author of this code free and
// harmless of any loss, liability, damage or costs that may arise from
// your use of this code.  This code is intended to be used for tutorial
// purposes only.  If you use it for any other purpose it is at your own risk.
// This code requires more testing.  You should not use it without first
// performing exhaustive testing on it.  If you are not qualified to make
// use of this code or parts thereof or to perform the required testing,
// seek the services of a qualified professional.
// -----Neal Glover  boglover[at]msn.com  303-466-5434-----
// -----June 6, 2010-----
//
// REVISION HISTORY
// --------------------------------------------
// Revision 100
// March 7, 2010 - Internal level
// --------------------------------------------
// Revision 200
// Abandoned
// --------------------------------------------
// Revision 300
// June 6, 2010  - First beta revision.
//               - Comment changes only.
// --------------------------------------------
// Revision 400
// Sept. 9, 2010  - Added the BTA root finder. Can now
//                - select between BTA and a Chien Search.
//                - Other small changes.
// Oct. 10, 2011  - Changed MAXCORR to 64 and retested
// --------------------------------------------
//
// NOTES:
// The level of flexibility in this code was motivated by the flexibiltiy
// in the Robert Morelos-Zaragoza bch decoder on the "ECC Page" web site.
//
// This code does not provide all the answers for a binary BCH code design.
// For example it does not make the code non-cyclic for the detection of
// incorrect synchronization framing and other errors.  If you are not aware
// of such considerations, you may need to consult a professional.
//
// This code will simulate most binary bch codes of parctical interest, but
// there are limitations.  If you enter a parameter that is rejected you can
// check the defines in the source code to see if the limitation might be
// overcome by changing a define.
//
// Global variables were used extensively in the test code.  The use of global
// varaibles was avoided in core BCH functions that might be used as a guide
// for designing functions of a firmware product.  There are a few global
// variables in bchDecode to assist testing.  They are marked "for testing only"
// and may be deleted in the function you write to integrate into a firmware
// product.
//
// It is recommended that you run lint or a lint like program to help find subtle
// errors in your code, especially array boundary problems.  Another way to
// help find such problems is to make a copy of your code and convert it to "C#" just
// for extensive testing. Converter programs exist to convert "C" to "C#".
//
// It is possible to add an (x-1) factor to the code generator poly of a
// binary bch code to increase the detection (not correction) capability
// by one.  Support for this is not included in this program.
//
// The inner most Chien search loop of this code is unrolled.  It is possible
// that speed can be improved further by unrolling some of the other loops.
// To unroll a loop, replicate the contents of the loop to eliminate loop
// overhead.  If for a particular loop the number of times the body
// is executed changes dynamically, just branch down into the
// replications to achieve the equivalent of looping the body the
// correct number of times.  For an example of loop unrolling see the Chien
// search function.  My original knowledge of unrolling loops came from
// the Earl Cohen Ph.D. thesis (Berkeley).  See also the Glover-Dudley
// patent 4,839,896.
//
// There are a number of places in the code where a mod operation is avoided
// because it is slow compared to the alternatives.
//
// There are still many places in the code where it can be speeded up some.
// Some of these places are identified in code comments.
//
// In this code the codeword is always stored in byte form even though the
// code is a binary BCH code and regardless of the size of the finite field.
// Finite field values such as syndromes and error location numbers are stored
// in "int" data types.  This is done to make the program flexible to handle
// a wide range of finite field sizes.  For a particular implementation
// you could try reducing the size of the values to save memory.  It is
// possibe that this would affect speed but I am not sure in which direction.
//
// This code supports two root finders, Chien and BTA.  You will need to
// test both to see which one best fits your requirements. For the parameter
// range that I have tested the BTA algorithm is much faster than the Chien
// search.  I have tested data block size 1024 bytes and GF(2^14) and
// 6 to 64 errors. I have not performed timing tests ourside that range.
// It should be noted that the BTA algorithm is far more complicated than the
// Chien search.
//
// ANY FEEDBACK WILL BE GREATLY APPRECIATED! ALL INPUT WILL BE ACCEPTED
// ESPECIALLY FEEDBACK POINTING OUT ERRORS IN THE CODE OR SUGGESTIONS
// FOR CODE IMPROVEMENT.
//
// I am considering writing a new book on Reed-Solomon, binary BCH, and
// LDPC codes.  If I do, it will cover the techniques of this program
// and many other practical aspects of error correcting code design.  It will
// also include extensions to many of the fundamental concepts from my earlier
// book.  I have been prepareing the material off and on for several years.
//
// -----Neal Glover  boglover[at]msn.com  303-466-5434-----
//
//
//#include "stdafx.h" // Not much in here.  Trying to be compiler independent
#include <stdio.h>  // Needed for printf and scanf
#include <tchar.h>  // Not needed right now
#include <math.h>   // Not needed right now
#include <time.h>	// Needed for time functions
#include <stdlib.h> // Needed for the rand and srand functions
//
// An instance of this structure is used to return 2 items from bchEval
struct statAndFCnt { // Status and failing pass number
	int stat;
	int FCnt;
};
//
// Defines for the functions that deal with codewords on disk
#define MAXFILESIZE  (5000000) // Maximum file size for reading codewords from disk
#define MAXLOOPALLCWSCNT (1000000) // For CWs from disk, max times to loop all CWs
// Definitions for evaluation code
#define MAXERRSTOSIM (200)     // Determines memory size for errors to simulate
#define ZERO			(0)			// Zero
// ################ DEFINITIONS AFFECTING STORAGE SPACE ################
// ***** IF YOU CHANGE MAXMPARM, YOU MUST CHANGE MAXFFSIZE AS WELL
#define MAXMPARM		(16)	// Max mParm - m of GF(2^m), Pgm not "designed" for m>16
// #define MINMPARM     (6)		// Min mParm, Pgm not "designed" to handle m<6
#define MINMPARM     (5)		// Min mParm, Pgm not "designed" to handle m<6
#define MAXCORR			(64)	// Max tParm, Pgm not "tested" for Max tParm> 64
//  Pgm not tested for MAXFFSIZE>65536
// ***** MAXFFSIZE MUST BE SET TO 2^MAXMPARM *****
#define MAXFFSIZE		(65536)	// Max finite field size - Sets storage requirement only
// Storage space but based on other definitions
#define MAXNUMSYN		(2*MAXCORR)	// Max num syndromes - Sets storage requirement only
#define MAXCODEWDBYTES	(MAXFFSIZE/8+1) // Max number of bytes in a codeword
#define MAXREDUNWDS (((MAXCORR*MAXMPARM)/8+1)/4+1)  //Max # redundancy words
// Definitions for clarity
#define BYTESTATES		(256)   // Number of states of a byte
// Init error definitions
//#define QUADBUILDERR  (2)		// Err building the quad table for special solutions
// Definition of error flag bits (uncorrectable errors)
#define DIVZRODIV	(0x0001)		// (1) Divide by 0 error in ffDiv
#define DIVZROINV	(0x0002)		// (2) Divide by 0 error in ffInv
#define ROOTSNEQLN	(0x0004)		// (4) Root count error in chienSearch
#define CGPFATAL    (0x0008)		// (8) Fatal error in cgp generation
//#define UNASSIGNED2 (0x0010)		// (16) Unasigned
//#define UNASSIGNED3 (0x0020)		// (32) Unasigned
#define CORROUTSIDE (0x0040)		// (64) Attempt to correct beyond last Codeword bit
#define BERMASERR	(0x0080)		// (128) ELP would exceed space allotted
#define TBLNOTINIT  (0x0100)		// (256) Tbls not init error in bchDecode
#define	E1QUADRATIC	(0x0200)		// (512) Divide by 0 err in ffQuadratic
#define E2QUADRATIC	(0x0400)		// (1024) Zero root of y^2+y+c in ffQuadratic
#define	E1CUBIC		(0x0800)		// (2048) Divide by 0 error in ffCubic
#define E2CUBIC		(0x1000)		// (4096) Zero root of y^2+y+c in ffCubic
#define E1QUARTIC	(0x2000)		// (8192) Sigma1=0 and Sigma3=0 in ffQuartic
#define E1CROOT	    (0x4000)		// (16384) "Not a cube" in ffCubeRoot
#define BTAFACTORTBLERR   (0x8000)	// BTA - factor table error
#define BTANUMROOTSERR    (0x010000)// BTA - # roots error
#define BTAREPEATEDROOTS  (0x020000)// BTA - repeated roots error
#define GCDZROONENTRYERR  (0x040000)// GCD - zero on entry error
#define MODZROONENTRYERR  (0x080000)// MOD - zero on entry error
#define QUOZROONENTRYERR  (0x100000)// QUOTIENT - zero on entry error
#define LOGALPHAIGTHMPARM (0x200000)// BTA - LOGALPHAi GTH mParm
#define LOGALOGBUILDERR   (0x400000)// Err building the log or alog table
//
// Definition of the status bits returned by eccDecode
#define CORR		(1)				// Correctable status
#define UNCORR		(2)				// Uncorrectable status
#define ERRFREE		(0)				// Error free status
//
// Definition of bchEval status bits - These added to status returned by bchDecode
#define UNCORRNOTEXPD (0x0010)		// (16) if dcdStatus==UNCORR && statusExpd<UNCORR
#define EXPDERRFREE   (0x0020)		// (32) if dcdStatus>ERRFREE && statusExpd==ERRFREE
#define EXPDERR		  (0x0040)		// (64) if dcdStatus==ERRFREE && statusExpd>ERRFREE
#define COMPAREERR    (0X0080)      // (128)Compare error
//
static int gblLogZVal,gblRootFindOption;
static int gblSigmaOrig[MAXCORR+1],gblSyndromes[MAXNUMSYN];
static int gblFFPoly,gblFFSize,gblLnOrig;
static int gblAlogTbl[2*MAXFFSIZE],gblLogTbl[MAXFFSIZE];
static int gblAppliedErrLocs[MAXERRSTOSIM],gblAppliedErrVals[MAXERRSTOSIM];
static int gblCodeword[MAXCODEWDBYTES], gblCodewordSav[MAXCODEWDBYTES];
static int gblRemainBytes[(MAXCORR*MAXMPARM)/8+1],gblLoc[MAXCORR];
static int gblKParm,gblMParm, gblNParm, gblMParmOdd,  gblTParm;
static int gblNumCodewordBytes,gblNumErrsApplied;
static int gblNumRedunBits, gblNumRedunBytes, gblNumDataBytes;
static int gblNumDataBits,gblNumRedunWords;
static int gblCgpBitArray[MAXCORR*MAXMPARM+1], gblCgpDegree;
static int gblMisCorrCnt,gblRawLoc[MAXERRSTOSIM];
static int gblBerMasUCECntr,gblRootFindUCECntr,gblFixErrorsUCECntr;
static int gblTraceTestVal,gblQuadCompTbl[MAXMPARM];
static unsigned int gblCgpFdbkWords[((MAXCORR*MAXMPARM)/8+1)/4+1];
static unsigned int gblEncodeTbl[BYTESTATES][MAXREDUNWDS];
static unsigned int gblRandomNum;
//
// Prototypes - If the functions are rearranged, more protypes will be required
static int ffInv(int opa,int *pErrFlg);
static int ffMult(int opa, int opb);
static int ffDiv(int opa, int opb,int *pErrFlg);
//
//****************************************************************
static void randomSetSeed(unsigned int mySeed)
{
	//****************************************************************
	//	Function: randomSetSeed
	//
	//	Function to initialize the random number generators.
	//  See comments for getRandom.
	//****************************************************************
	gblRandomNum=mySeed;
	srand(mySeed);
}

static unsigned int getRandom()
{
	//****************************************************************
	//	Function: getRandom
	//
	//	Function to get a random number for the simulator.  I use 2
	//  random number generators and add the 2 random numbers.  To
	//  generate one of the random numbers I use a degree 32 primitive
	//  GF(2) polynomial in order to make the length of the random
	//  number sequence very long before repeating starts.  The second
	//  random number generator is the "C" function rand().
	//
	//  fdbkConstant is the feedback pattern for the minimum polynomial
	//  of alpha 11 on page 492 of the Peterson and Weldon book
	//****************************************************************
	static const unsigned int fdbkConstant=0x76B553;

	if ((gblRandomNum & 0x80000000) > 0){
		gblRandomNum <<= 1;
		gblRandomNum ^= fdbkConstant;
	}
	else
		gblRandomNum <<= 1;
	// Below - divide by 2 because the add might affect number distribution in the
	// low digit.  This strategy probably ok for the small random numbers we want
	// in this program.  Non repeating is important and it should be achieved
	// with this strategy.
	//return (((gblRandomNum+(unsigned int)rand()) & 0x7fffffffU)/2);
	return (0x11111111);
}
static void pickFieldGenPoly()
{
	//****************************************************************
	//	Function: pickFieldGenPoly
	//
	//	If you choose not to enter your own GF(2) primitive polynomial,
	//  This function will use m of GF(2^m) to pick the first entry
	//  of the required degree from Peterson and Weldon's tables.
	//
	// This function was motivated by a similar function in the
	// Robert Morelos-Zaragoza bch decoder on the "ECC Page" web site.
	//****************************************************************
	// Source for field generator polynomials - Peterson and Weldon book.
#if 0
	int fieldPolyTbl[15] = {
		67, 137, 285, 529, 1033, 2053, 4179, 8219, 17475,
		32771, 69643, 131081, 262273, 524327, 1048585
	};
	gblFFPoly = fieldPolyTbl[gblMParm-6]; // m = 6 table entry is at location 0
#else
	int fieldPolyTbl[] = {
		37, 67, 137, 285, 529, 1033, 2053, 4179, 8219, 17475,
		32771, 69643, 131081, 262273, 524327, 1048585
	};
	gblFFPoly = fieldPolyTbl[gblMParm-5]; // m = 5 table entry is at location 0
#endif
}

static void buildLogAlogTbls()
//****************************************************************
//	Function: buildLogAlogTbls
//
//	This function is called during initialization to build finite field
//  log and alog tables for the encoder and decoder and for simulator
//  functions.
//****************************************************************
{
	int kx;
	unsigned int shiftReg,fdbkCon;
	int ffSizeDivTwo;

	// Construct the finite field log and alog tables
	shiftReg=1;
	fdbkCon=(unsigned int)(gblFFPoly-gblFFSize);
	ffSizeDivTwo=gblFFSize/2;
	// 9-1-10 Doubled size of alog tbl for speed
	for (kx=0;kx<2*gblFFSize;kx++){
		gblAlogTbl[kx] = (int)shiftReg;
		if (kx<gblNParm){
			gblLogTbl[shiftReg] = kx;
		}
		if (shiftReg>=(unsigned int)ffSizeDivTwo){
			shiftReg=((shiftReg<<1) & (unsigned int)gblNParm)^fdbkCon;
		}
		else {
			shiftReg<<=1;
		}
	}
	// 9-1-10 Changed value for log of zero
	gblLogTbl[0]=gblLogZVal; // This is the value for log of zero
	// 9/2010 Changed next line for double size alog table
	gblAlogTbl[gblLogZVal] = 0;
}
static int chkLogAlogTbls()
//****************************************************************
//	Function: chkLogAlogTbls
//
//	This function tests the finite field log and alog tables to
//  determine if they are correct.  If they are not, an error
//  code is returned.
//****************************************************************
{
	int kx;

	// Check log and alog tables
	for (kx=0;kx < gblNParm;kx++){
		if (kx != gblLogTbl[gblAlogTbl[kx]]){
			return(LOGALOGBUILDERR); // Flag error
		}
	}
	// Next line changed 9-5-10 for double size alog table
	if (gblLogTbl[0]!=gblLogZVal || gblAlogTbl[gblLogZVal] != 0){
		return(LOGALOGBUILDERR); // Flag error
	}
	return(ZERO);
}


static int ffSquareRoot(int opa)
{
	//****************************************************************
	//	Function: ffSquareRoot
	//
	//	This function returns the square root of its operand.
	//	If the finite field log of opa is even it returns -
	//	alog(log(opa)/2).  But if the finite field log of opa is
	//	odd it returns - alog((log(opa)+(finite field size -1))/2).
	//****************************************************************
	int logtmp;

	if (opa==0){
		return(0);
	}
	else
	{
		logtmp=gblLogTbl[opa];
		if (logtmp%2==1){
			logtmp=logtmp+gblNParm;
		}
		return(gblAlogTbl[logtmp/2]);
	}
}

static int ffCubeRoot(int opa, int *pErrFlg)
{
	//****************************************************************
	//	Function: ffCubeRoot
	//
	//	This function returns the cube root of its operand.
	//	If the finite field log of opa is divisible by three, it
	//	returns alog(log(opa)/3).  But if the finite field log of
	//	opa is "not" divisible by three, there are no cube roots
	//	of the operand in the field and the error flag is set.
	//	The returned value in this case is not important.  This
	//  function is used only with "m" even (m of GF(2^m))
	//****************************************************************
	int logtmp;

	if (opa==0){
		return(0);
	}
	else
	{
		logtmp=gblLogTbl[opa];
		if (logtmp%3 != 0) {
			*pErrFlg|=E1CROOT; // Or error into location pointed to
		}
		return(gblAlogTbl[logtmp/3]);
	}
}
static void genTraceTestVal()
{
	//****************************************************************
	// Function: genTraceTestVal
	//
	// Function to compute a trace test value (gblTraceTestVal).
	// This function is called at initialization
	// time only.  This value is used in fast root
	// finding.
	//
	// The trace test value is computed by computing
	// the trace of the first "m" elements of the field.
	// Let k be a number between 0 and gblNParm-1 then bit k
	// of the gblTraceTestVal will be "1" if the trace of alpha^k
	// is one and "0" otherwise.
	//****************************************************************
	int shifter,x,jx,kx,sum;

	gblTraceTestVal=0;
	shifter=1;
	for (jx=0;jx<gblMParm;jx++){ // Compute trace for alpha^jx
		x=gblAlogTbl[jx]; // x is element for which trace will be computed
		sum=0;
		for (kx=0;kx<gblMParm;kx++){ // This is the trace computing loop
			sum^=x;
			x=ffMult(x,x);
		}
		if (sum==1){
			gblTraceTestVal=gblTraceTestVal^shifter; // Set appropriate bit if trace is "1"
		}
		shifter*=2;
	}
}

static void genQuadCompTbl()
{
	//****************************************************************
	//	Function: genQuadCompTbl
	//
	//	Function to find a list of components of y to be used
	//	in finding solutions to Y^2+y+c using a table of
	//  linear components. This function is called at initialization
	//  time only.  This solution uses a small table with gblMParm
	//  entries.  A large table would be faster but would
	//  require 2^gblMParm entries.  It is possible to have an inbetween
	//  solution that uses two moderate size tables and is faster than
	//  this small table solution but not as fast as the one large table
	//  solution.  The reference for this solution is Dr. Berlekamp's 1968
	//  "Algebraic Coding Theory" book pages 243 and 244.
	//  The first loop of this function will put in the search table
	//  each "c" of y^2+y = c that has a single "1" bit or the xor of
	//  such a pattern with the value of a fixed element with
	//  trace = "1". The second loop searches for each element in the
	//  table within a set of solutions for
	//  y^2+y=c and if it finds one then it stores the associated
	//  value of y in the gblQuadCompTbl.
	//****************************************************************
	int c,y,kx,shifter,searchTbl[MAXMPARM],firstTraceOne;

	firstTraceOne=0;
	// This loop will put in the search table each "c" of y^2+y = c that has a
	// single "1" bit or the xor of such a pattern with a fixed pattern of an
	// element with trace = "1"
	shifter=1;
	for (kx=0;kx<gblMParm;kx++){
		searchTbl[kx]=0;
		if ((gblTraceTestVal&shifter)>0){
			if (firstTraceOne>0){
				searchTbl[kx]=shifter^firstTraceOne;
			}
			else{
				firstTraceOne=shifter;
			}
		}
		else{
			searchTbl[kx]=shifter;
		}
		shifter*=2;
	}
	for (kx=0;kx<gblMParm;kx++){
		gblQuadCompTbl[kx]=0; // Clear table
	}
	// This loop will search for each "c" of y^2+y = c that has a single bit
	// or is the xor of such a pattern with a fixed pattern of an
	// element with trace = "1"
	for (y=0;y<gblFFSize;y+=2){//Incr of 2 may not work with another basis
		c=ffMult(y,y)^y;
		for (kx=0;kx<gblMParm;kx++){
			if (c==searchTbl[kx]){
				gblQuadCompTbl[kx]=y;
			}
		}
	}
}

static int ffQuadFun(int c)
{
	//****************************************************************
	// Function: ffQuadFun
	//
	// This function is used during decode operations of
	// error correction.  When a solution is needed for the equation
	// y^2+y+c=0 the function first ANDs the trace
	// test value with "c" and then tests parity
	// of the resulting bits.  If parity is odd (trace = 1)
	// an error is posted and the function exits.  If parity is
	// even then the function continues on and XORs a set of values to
	// get y1, one solution of y^2+y=c.  The second solution is the
	// XOR of the first solution with a single one bit in the low
	// order position.
	//
	// This function can be made faster with a single large table.  This
	// table would be the same size as the log and alog tables.  So
	// for large "m" (m of GF(2^m)) the table will be large.  Using
	// lots of memory may not be good for a firmware application but
	// might be ok for a PC application.
	//****************************************************************
	int kx,y1,shifter,parityWord,parityBit;

	// "c" is "c" of y^2+y=c.  Given "c" we want one solution (y1) of y^2+y=c
	parityWord=gblTraceTestVal & c;  // Get bits that parity is to be checked on
	parityBit=0;
	y1=0;
	shifter=1;
	for (kx=0;kx<gblMParm;kx++){
		if ((c & shifter)>0){
			y1^=gblQuadCompTbl[kx];// XOR values from the quadratic component table
		}
		if ((parityWord & shifter)>0){
			parityBit ^= 1; // Accumulate parity to get trace of "c"
		}
		shifter*=2;
	}
	if (parityBit > 0){ // If trace of "c" is "1"
		return (0);
	}
	else{
		return(y1^1);// "1" Not for all basis.  y1 is one solution of y^2+y=c.
	}
}

static void bchInit()
{
	//****************************************************************
	//	Function: bchInit
	//
	//	Function to do initialization computations.
	//****************************************************************
	buildLogAlogTbls();
	genTraceTestVal();
	genQuadCompTbl();
}

static void linearElp(const int sigmaN[],int Loc[])
{
	//****************************************************************
	//	Function:  linearElp
	//
	//	Function to find the root of a degree one error locator
	//  polynomial (ELP)
	//****************************************************************
	Loc[0]=sigmaN[1];
}

static int quadraticElp(const int sigmaN[],int Loc[])
{
	//****************************************************************
	//	Function:  quadraticElp
	//
	//	Function to find the roots of a degree two error locator
	//  polynomial (ELP).  References for this algorithm is the 1978
	//  Flagg patent 4,099,160 and the paper "High Speed Interleaved
	//  Reed-Solomon Error Detection and Correction System" by
	//  Shirish Deodhar and E.J. Weldon, Journal of Photo-Optical
	//  Instrumentation Engineers (SPIE), 1983.  See also the Deodhar
	//  patent 4,567,594 and the Glover-Dudley patents 4,839,896 and
	//  5,280,488.
	//****************************************************************
	int c,y1,y2,errFlg;

	errFlg=0;
	if (sigmaN[1]==0){
		errFlg|=E1QUADRATIC;
	}
	c=ffDiv(sigmaN[2],ffMult(sigmaN[1],sigmaN[1]),&errFlg);
	y1=ffQuadFun(c);
	if (y1==0){
		errFlg|=E2QUADRATIC;
	}
	// In next line, if ever use basis other than poly then "1" may have
	// to change to alog[0]
	y2=y1^1;
	Loc[0]=ffMult(sigmaN[1],y1);
	Loc[1]=ffMult(sigmaN[1],y2);
	return (errFlg);  // Return error flag
}

static int cubicElp(int nParm, const int sigmaN[],int Loc[],const int alogTbl[])
{
	//****************************************************************
	//	Function:	cubicElp
	//
	//	Function to find the roots of a degree three error locator
	//  polynomial (ELP).  This code follows the algorithm of page 355
	//  of the 1991 Glover-Dudley book "Practical Error Correction Design
	//  for Engineers" REVISED SECOND EDITION.  My original reference for
	//  that book was the 1978 Flagg patent 4,099,160.  See also
	//  the paper "High Speed Interleaved Reed-Solomon Error
	//  Detection and Correction System" by Shirish Deodhar
	//  and E.J. Weldon, Journal of Photo-Optical
	//  Instrumentation Engineers (SPIE), 1983.  See also the Deodhar
	//  patent 4,567,594 and the Glover-Dudley patents 4,839,896 and
	//  5,280,488.  See also the paper by Chien, Cunningham, and
	//  Oldham titled "Hybrid Methods for Finding Roots of a Polynomial
	//  With Application to BCH Decoding" published in Transactions on
	//  Information Theory", March 1969 pages 329-335.
	//
	//  Note to Neal. ###### Next level - list here the conditions that
	//  result in an uncorrectable error and where each is detected.
	//***************************************************************
	int v1,u1,t1,t2,t3,n,d,n3,d2,c,errFlg;

	errFlg=0;
	// n for numerator, d for denominator
	n=sigmaN[2]^ffMult(sigmaN[1],sigmaN[1]);
	d=sigmaN[3]^ffMult(sigmaN[1],sigmaN[2]);
	// Note to Neal.  The error check on the next line is
	// redundant. This error would also get caught in "ffDiv" function.
	if (d==0){
		errFlg|=E1CUBIC; // Divide by "0" error
	}
	n3=ffMult(n,ffMult(n,n)); // Numerator cubed
	d2=ffMult(d,d);           // Denominator squared
	c=ffDiv(n3,d2,&errFlg);   // Finite field divide
	// Note to Neal.  ######## I think I put in the next decision
	// during debug in 1999.  The code could be extensively tested
	// without this decision to see if it can be left out.
	if (n==0){
		u1=d;
	}
	else
	{
		// The quad function is equiv to fetching from large table
		v1=ffQuadFun(c);
		if (v1==0){
			errFlg|=E2CUBIC;
		}
		u1=ffMult(v1,d);
	}
	// Roots of transformed cubic
	t1=ffCubeRoot(u1,&errFlg);
	t2=ffMult(t1,alogTbl[nParm/3]); // nParm/3 is 85 for gf(2^8)
	t3=t1^t2; // Equivalent to t2= line with nParm replaced by 2*nParm
	// Roots of original cubic
	Loc[0]=sigmaN[1]^t1^ffDiv(n,t1,&errFlg);
	Loc[1]=sigmaN[1]^t2^ffDiv(n,t2,&errFlg);
	Loc[2]=sigmaN[1]^t3^ffDiv(n,t3,&errFlg);
	return (errFlg);  // Return error flag
}

static int quarticElp(int nParm,int sigmaN[],int Loc[],const int alogTbl[])
{
	//****************************************************************
	//	Function:	quarticElp
	//
	//	Function to find the roots of a degree four error locator
	//  polynomial (ELP).  My code follows the quartic algorithm
	//  of the paper "High Speed Interleaved Reed-Solomon
	//  Error Detection and Correction System" by
	//  Shirish Deodhar and E.J. Weldon, Journal of Photo-Optical
	//  Instrumentation Engineers (SPIE), 1983.  See also the 1978
	//  Flagg patent 4,099,160, the Deodhar patent 4,567,594,
	//  and the Glover-Dudley patents 4,839,896 and
	//  5,280,488.  See also the paper by Chien, Cunningham, and
	//  Oldham titled "Hybrid Methods for Finding Roots of a Polynomial
	//  With Application to BCH Decoding" published in Transactions on
	//  Information Theory", March 1969 pages 329-335.
	//
	//  I am currently consolidating my notes on how the methods
	//  of fast root finding for cubic and quartic finite field
	//  polynomials were developed.
	//
	//  Note to Neal. ###### Next level - list here the conditions that
	//  result in an uncorrectable error and where each is detected.
	//****************************************************************
	int n,errFlg,Ln;
	//
	int sigbk[MAXCORR+1],b2,b3,b4,b4n,b4d;
	int qq,ss,tt,tmp;

	errFlg=0;
	Ln=4; // ELP degree is 4 for quartic
	for (n=0;n<=Ln;n++){
		sigbk[n]=sigmaN[n];
	}
	//	Test for special case - poly already in correct form
	if (sigbk[1]==0){
		// Note to Neal.  Re-test extensively without the next decision.
		// Uncorrectable errors caught by this check may be caught by
		// other checks.  Re-testing will require comparing miscorrection
		// rate with and without the decision.
		if (sigbk[3]==0){
			errFlg|=E1QUARTIC;
		}
		b2=sigbk[2];
		b3=sigbk[3];
		b4=sigbk[4];
	}
	else
	{
		// ---------- Step b of the Deodhar-Weldon paper ----------
		b4n=ffMult(sigbk[1],sigbk[1]);
		b4d=ffMult(sigbk[3],sigbk[3])
			^ffMult(sigbk[1],ffMult(sigbk[2],sigbk[3]))
			^ffMult(sigbk[4],b4n);
		b4=ffDiv(b4n,b4d,&errFlg);
		b3=ffMult(sigbk[1],b4);
		b2=ffMult(b4,
			ffSquareRoot(ffMult(sigbk[1],sigbk[3]))^sigbk[2]);
	}
	// ---------- Step c of the Deodhar-Weldon paper ----------
	//	Set up a cubic and find its 3 roots.
	sigmaN[0]=1;
	sigmaN[1]=0;
	sigmaN[2]=b2;
	sigmaN[3]=b3;
	errFlg|=cubicElp(nParm,sigmaN,Loc,alogTbl);
	qq=Loc[1];
	// ---------- Step d of the Deodhar-Weldon paper ----------
	//	Set up a quadratic and find its two roots
	sigmaN[0]=1;
	sigmaN[1]=ffDiv(b3,qq,&errFlg);
	sigmaN[2]=b4;
	errFlg|=quadraticElp(sigmaN,Loc);
	ss=Loc[0];
	tt=Loc[1];
	//	First 2 roots of quartic
	sigmaN[0]=1;
	sigmaN[1]=qq;
	sigmaN[2]=ss;
	errFlg|=quadraticElp(sigmaN,Loc);
	Loc[2]=Loc[0];
	Loc[3]=Loc[1];
	//	Last 2 roots of quartic
	sigmaN[0]=1;
	sigmaN[1]=qq;
	sigmaN[2]=tt;
	errFlg|=quadraticElp(sigmaN,Loc);
	//	Skip inverse substitution if special case
	if (sigbk[1]!=0){
		// ---------- Step e of the Deodhar-Weldon paper ----------
		//	Do inverse substitution
		tmp=ffSquareRoot(ffDiv(sigbk[3],sigbk[1],&errFlg));
		for (n=0;n<4;n++){
			Loc[n]=ffInv(Loc[n],&errFlg)^tmp;
		}
	}
	return (errFlg);  // Return error flag
}

static int genCodeGenPoly()
{
	//****************************************************************
	//	Function: genCodeGenPoly
	//
	//	Function to compute a code generator polynomial for a binary
	//  BCH code.  The action of this function is equivalent to finding
	//  the polynomial that is the least common multiple of a set of
	//  minimum polynomials.
	//
	//  This function was motivated by a similar function in the
	//  Robert Morelos-Zaragoza bch decoder on the "ECC Page" web site.
	//
	//  Generally, finding the code generator polynomial requires two
	//  steps  1) compute minimum polynomials and 2) find the LCM of
	//  the minimum polynomials to determine the code generator polynomial.
	//  But it is also possible to combine these steps.  Finding a
	//  minimum polynomial is accomplished by finding the product of a
	//  set of roots and finding the code generator polynomial is
	//  accomplished by finding the product of a set of minimum
	//  polynomials.  So the two steps can be combined to find the code
	//  generator polynomial in one step by finding the product of a
	//  larger set of roots.  This program uses this combined approach.
	//
	//  It might be possible though to pull a little time out of this function
	//  by using the more general two step approach.  By careful coding it may
	//  be possible to take advantage of the fact that once the minimum
	//  polynomials have been computed, computing the code generator polynomial
	//  from them involves only GF(2) finite field operations.
	//
	//  This function is used only during development, so it would not be part
	//  of an implementation in a product employing one fixed code.
	//
	//  NOTE: The code generator polynomial is stored in a array,
	//  one bit per int and Low order in address 0.
	//
	//  Note to Neal: On next update check value of root before using it
	//  as an index to store "1" in the flg array.  This will reduce
	//  significantly the size of the flg array and will reduce the time
	//  to initialize this array.
	//****************************************************************
	int flg[MAXFFSIZE], tmp[MAXCORR*MAXMPARM+1];
	int	kx, root, rootBase, wk, errFlg; // root and rootBase are in log form

	for (kx=0;kx<gblFFSize;kx++){
		flg[kx]=0;// Index to flg can have values root,2*root,4*root,8*root...
	}
	for (kx=0;kx<=gblMParm*gblTParm;kx++){
		gblCgpBitArray[kx]=0;// Initialize
		tmp[kx]=0;			 // Initialize
	}
	gblCgpDegree=0; // The degree of the code generator poly is initialized to "0"
	gblCgpBitArray[0]=1; // Now the initial code generator poly is "1" (degree "0")
	errFlg=0;
	for (rootBase=1;rootBase<=2*gblTParm-1;rootBase += 2){ // alpha 1,3,5,7 etc.
		if (flg[rootBase] == 0){ // If this root not already processed
			root = rootBase;
			for(;;){ // Infinite loop - Exit is by "break"
				// In loop - root will take values like 1,2,4,8... 3,6,12,24...etc
				if (gblCgpDegree+1>gblMParm*gblTParm){
					errFlg=CGPFATAL;  // Fatal problem of some type
					break;
				}
				// Start multiply this factor times current gblCgpBitArray
				// Intermediate poly products will have some values greater than one
				// The coefficients of the final product will have values 0 or 1 only.
				// Multiply code gen poly by (X - alpha^root)
				for (kx=gblCgpDegree+1;kx>=1;kx--){
					// Shift poly 1 place (low to high)
					tmp[kx]=gblCgpBitArray[kx-1];
				}
				tmp[0]=0;
				gblCgpDegree++;
				for (kx=gblCgpDegree;kx>=0;kx--){
					if (gblCgpBitArray[kx]>0){
						wk=ffMult(gblCgpBitArray[kx],gblAlogTbl[root]);
						gblCgpBitArray[kx]=tmp[kx]^wk;
					}
					else {
						gblCgpBitArray[kx]=tmp[kx]; // Move the shifted version of cgp
					}
				}
				// End multiply
				flg[root] = 1;
				root *= 2; // root is in finite field log form
				if (root>=gblNParm){			//
					root -= gblNParm; // These 3 lines do a fast mod op
				}							//
				if (root == rootBase){
					break;
				}
			} // end for(;;)
		} // end if
	} // end for
	return(errFlg);
} // end function
static void cvtCgpBitToCgpWord()
{
	//****************************************************************
	//	Function: cvtCgpBitToCgpWord
	//
	//	Function to convert, the code generator polynomial (cgp) stored one
	//  bit per int and organized low order to address 0, into a set of
	//  words that contain the shift register feedback pattern for the code
	//  generator polynomial.  The feedback words are stored high order to
	//  address 0.  SO THE ORGANIZATION CHANGES BETWEEN THE INPUT AND OUTPUT
	//  OF THIS FUNCTION.
	//****************************************************************
	int nnn, kkk, wordAddr;
	unsigned int bitMask;

	// Convert gblCgpBitArray to cgpWordArray
	// That is, 1 bit per word to 32 bits per word
	for (nnn=0; nnn < gblNumRedunWords;nnn++){
		gblCgpFdbkWords[nnn] = 0; // Clear fdbk array
	}
	wordAddr = 0;
	bitMask =0x80000000;
	// Place the code generator poly in the ints
	for (kkk=gblNumRedunBits-1;kkk>=0;kkk--){
		if (gblCgpBitArray[kkk]>0){
			gblCgpFdbkWords[wordAddr] ^= bitMask;
		}
		bitMask >>= 1; // bitMask must be unsigned
		if (bitMask ==0){
			bitMask =0x80000000;
			wordAddr++;
		}
	}
}

static void genEncodeTbls()
{
	//****************************************************************
	//	Function: genEncodeTbls
	//
	//  Function to generate encode tables.  The encode tables allow
	//  8 bits of data to be processed at a time.  A Left shifting shift
	//  reg is implemented in multiple 32-bit words.  A multiple word
	//  shift register value is stored in the table for each possible
	//  value of a byte.  A byte value is placed in the left most
	//  (high order position of lowest address word).  Then the shift
	//  register is shifted 8 times with feedback.  Then the multiple
	//  word shift register value is stored in the table at the addr
	//  specified by the byte.
	//
	//  The feedback words are highest order in lowest address and the
	//  resulting encode table is organized the same way.
	//****************************************************************
	unsigned int iii,fdbk,fdbkSav,SR[MAXREDUNWDS];
	int jjj,nnn;

	// Gen Encode Table
	for (iii = 0;iii<BYTESTATES;iii++){ // Encoding is 8 bits parallel
		for (jjj=0; jjj < gblNumRedunWords;jjj++){ // Clear shift register
			SR[jjj] = 0;
		}
		SR[0] ^= (iii << 24);
		for (jjj = 0;jjj<=7;jjj++){// 7 is the # of bits in a byte -1
			fdbk = 0;
			for (nnn=gblNumRedunWords-1; nnn >=0;nnn--){
				fdbkSav = fdbk;
				if (SR[nnn] & 0x80000000){
					fdbk = 1;
				}
				else {
					fdbk = 0;
				}
				SR[nnn] <<= 1;
				SR[nnn] ^= fdbkSav;
			}
			if (fdbk == 1){
				for (nnn=0; nnn < gblNumRedunWords;nnn++){
					SR[nnn] ^= gblCgpFdbkWords[nnn];
				}
			}
		}
		for (nnn=0; nnn < gblNumRedunWords;nnn++){
			gblEncodeTbl[iii][nnn] = SR[nnn];  // Move SR to encode table
		}
	}
}

static void clearWriteCW()
{
	//****************************************************************
	//	Function: clearWriteCW
	//
	//	Function to clear write data.
	//****************************************************************
	int kx;

	// Clears gblCodeword including pad bits
	for (kx=0;kx<gblNumCodewordBytes;kx++){
		gblCodeword[kx] = 0; // Clear gblCodeword
	}
}

static void genWriteData()
{
	//****************************************************************
	//	Function: genWriteData
	//
	//	Function to generate random data byte values for encoding.  No
	//  Matter what size of finite field we are using, the encoder
	//  operates on bytes in byte parallel.
	//****************************************************************
	int kx;

	// Clear gblCodeword including pad bits
	for (kx=0;kx<gblNumCodewordBytes;kx++){
		gblCodeword[kx] = 0; // Clear gblCodeword
	}
	// Put random user data in gblCodeword - there are no pad bits in user data
	for (kx=0;kx<gblNumDataBytes;kx++){
		gblCodeword[kx] = (getRandom() % BYTESTATES);// Reduce random # to byte size
	}
}

static void bchEncode(const unsigned int encodeTbl[BYTESTATES][MAXREDUNWDS],int codeword[],
					  int numRedunWords,int numRedunBytes,int numDataBytes)
{
	//****************************************************************
	//	Function: bchEncode
	//
	//  This is the encoder.  It performs its function byte parallel using
	//  an encode table that is built during initialization.  It processes
	//  one byte at a time regardless of the size of the finite field.
	//
	//  The parallel approach shifts a software shift register implementing the
	//  code generator polynomial once per codeword byte even though the binary
	//  BCH code generator polynomial is over GF(2).  The shift occurs once per
	//  byte regardless of the size of the finite field.  The value of the
	//  codeword byte being processed addresses a large table to get the next
	//  value of the software shift register.  See page 243 (k-bit serial techniques)
	//  of "Practical Error Correction Design for Engineers" (revised second
	//  edition 1991) by Neal Glover and Trent Dudley.  Page 347 of this book
	//  describes tables for an early implementation of k-bit serial for
	//  a binary code (a computer generated code).  The k-bit serial
	//  method described there would be the same as for a binary
	//  BCH code.  I learned early on that it was possible to implement
	//  a parallel shift register for a binary code.  I may have learned this
	//  from a paper by Hsiao and Sih titled "Serial-to-Parrallel Transformation
	//  of Linear-Feedback Shift-Register Circuits" which appeared in
	//  IEEE. Trans. on Elec. Comp., 738-740 (Dec. 1964).
	//****************************************************************
	unsigned int SR[MAXREDUNWDS], tmp, fdbk, fdbkSav;
	// +5 So that we can temporarily keep remainder bytes in whole words
	int redunByteArray[(MAXCORR*MAXMPARM)/8+5];
	int jx ,kx, writeCWAddr;

	for (kx=0; kx < numRedunWords;kx++){ // Clear encode shift register
		SR[kx] = 0;
	}
	for (writeCWAddr = 0; writeCWAddr < numDataBytes; writeCWAddr++) {// index to write data buffer
		fdbk = 0;
		for (jx=numRedunWords-1; jx >=0 ;jx--){
			// IT IS POSSIBLE TO UNROLL THIS LOOP USING A SWITCH STATEMENT.
			// WILL HAVE TO PROCESS FROM HIGH ADDRESS BACK TOWARDS ZERO
			fdbkSav = fdbk;
			// 32 bits per int - but process 8 bits at a time (8 bits unrelated to "m")
			fdbk = (SR[jx] >> 24) & 0x000000ff;
			SR[jx] = (SR[jx] << 8) ^ fdbkSav;  // 8 # bits in parallel - unrelated to "m"
		}
		fdbk ^= (unsigned int)codeword[writeCWAddr];
		// IT IS POSSIBLE TO UNROLL THIS LOOP USING A SWITCH STATEMENT
		for (jx=0; jx < numRedunWords;jx++){
			SR[jx] ^= encodeTbl[fdbk][jx];
		}
	}
	// Copy redundancy bytes from shift register (SR) word array
	for (kx=0;kx<numRedunWords;kx++){
		tmp=SR[kx]; // Fetch a 32 bit word
		redunByteArray[4*kx] = tmp>>24;
		redunByteArray[4*kx+1] = (tmp>>16) & 0x000000ff;
		redunByteArray[4*kx+2] = (tmp>>8) & 0x000000ff;
		redunByteArray[4*kx+3] = tmp & 0x000000ff;
	}
	// Copy redundancy to codeword.  Don't copy pad bytes of SR word array
	for (kx=0;kx<numRedunBytes;kx++){
		codeword[numDataBytes+kx]=redunByteArray[kx];
	}
}

static void savCodeword()
{
	//****************************************************************
	//	Function: savCodeword
	//
	//	Function to save the original codeword resulting from encoding,
	//  so that it can be compared later with the corrected codeword.
	//****************************************************************
	int kx;

	for (kx=0;kx<gblNumCodewordBytes;kx++){
		gblCodewordSav[kx]=gblCodeword[kx];
	}
}

static int applyErrors(int lowNumErrs,int hiNumErrs)
{
	//****************************************************************
	//	Function: applyErrors
	//
	//	Function to apply errors to the codeword before decoding.  The
	//  minimum and maximum number of errors to apply are established
	//  by input from the user.  Errors are being put in a codeword of
	//  bytes regardless of finite field size.
	//****************************************************************
	int jx,kx,bitLoc,byteLoc,byteValue;
	int byteBitNum,uniqueFlg;

	gblNumErrsApplied = ((int)getRandom()%((hiNumErrs-lowNumErrs)+1))+lowNumErrs;
	for (kx=0;kx<gblNumErrsApplied;kx++){
		// Highest bit err gblLoc will be gblNumDataBits+gblNumRedunBits-1,
		// which is the last gblCodeword bit location.  NOTE: Last
		// gblCodeword bit not necessarialy on a byte boundary.  But
		// this code works on bytes so there are possibly zero fill
		// bits in the last gblCodeword byte.  If we put errs in the
		// fill positions they could appear to be in the front
		// part of the gblCodeword and miscorrection could result.
		do {
			// Working with bits because of pad bits
			bitLoc = (int)getRandom()%(gblNumDataBits+gblNumRedunBits);
			uniqueFlg=1;
			for (jx=0;jx<kx;jx++){ // If kx==0 this loop will not do anything
				if (bitLoc==gblRawLoc[jx]){
					uniqueFlg=0;
					break;
				}
			}
		}while ( uniqueFlg ==0);
		gblRawLoc[kx]=bitLoc;
		byteLoc = bitLoc/8;
		byteBitNum = 7 - (bitLoc % 8);
		byteValue=1;
		for (jx=byteBitNum;jx>=1;jx--){
			byteValue *= 2;
		}
		gblCodeword[byteLoc] ^= byteValue;
		gblAppliedErrLocs[kx]=byteLoc;
		gblAppliedErrVals[kx]=byteValue;
	}
	return (gblNumErrsApplied);
}

static int computeRemainder(const int codeword[],int numRedunWords,int numRedunBytes,
							int numDataBytes,
							const unsigned int encodeTbl[BYTESTATES][MAXREDUNWDS],
							int remainBytes[])
{
	//****************************************************************
	//	Function: computeRemainder
	//
	//  This function computes a remainder.  This remainder is the remainder
	//  from dividing the received codeword by the full code generator polynomial
	//  (not by factors of it).  This is the key to speeding up syndrome computation.
	//  This remainder is computed using a parallel approach similar to that used
	//  in encoding.  Once we have a full remainder we can compute syndromes from it.
	//  This is faster than computing syndromes directly from the codeword because
	//  the remainder is much shorter than the codeword.
	//
	//  This also gives us a fast way to determine if the syndromes would all be zero.
	//  If the remainder is all zeros then there is no need to compute syndromes from
	//  it as they would all be zeros also.  If the remainder is zero then the
	//  codeword is assumed to be error free and no further decoding steps are
	//  necessary.
	//
	//  The parallel approach shifts a software shift register implementing the
	//  code generator polynomial once per codeword byte even though the binary
	//  BCH code generator polynomial is over GF(2).  The shift occurs once per
	//  byte regardless of the size of the finite field.  The value of the
	//  codeword byte being processed addresses a large table to get the next
	//  value of the software shift register.  See page 243 (k-bit serial techniques)
	//  of "Practical Error Correction Design for Engineers" (revised second
	//  edition 1991) by Neal Glover and Trent Dudley.  Page 347 of this book
	//  describes tables for an early implementation of k-bit serial for
	//  a binary code (a computer generated code).  The k-bit serial
	//  method described there would be the same as for a binary
	//  BCH code.  I learned early on that it was possible to implement
	//  a parallel shift register for a binary code.  I may have learned this
	//  from a paper by Hsiao and Sih titled "Serial-to-Parrallel Transformation
	//  of Linear-Feedback Shift-Register Circuits" which appeared in
	//  IEEE. Trans. on Elec. Comp., 738-740 (Dec. 1964).
	//****************************************************************
	unsigned int fdbk,fdbkSav,SR[MAXREDUNWDS];
	int nnn,readCWAddr,remainderDetdErr;

	for (nnn=0; nnn < numRedunWords;nnn++){ // Clear encode shift register
		SR[nnn] = 0;
	}
	// SHIFTS WITH FEEDBACK
	for (readCWAddr = 0; readCWAddr < numDataBytes; readCWAddr++) {// index to write data buffer
		fdbk = 0;
		for (nnn=numRedunWords-1; nnn >=0 ;nnn--){
			// IT IS POSSIBLE TO UNROLL THIS LOOP USING A SWITCH STATEMENT.
			// WILL HAVE TO PROCESS FROM HIGH ADDRESS BACK TOWARDS ZERO
			fdbkSav = fdbk;
			fdbk = (SR[nnn] >> 24); // 32 - # bits in word (unrelated to "m")
			SR[nnn] = (SR[nnn] << 8) ^ fdbkSav;  // 8 # bits in parallel - unrelated to "m"
		}
		fdbk ^= (unsigned int)codeword[readCWAddr];
		for (nnn=0; nnn < numRedunWords;nnn++){
			// IT IS POSSIBLE TO UNROLL THIS LOOP
			SR[nnn] ^= encodeTbl[fdbk][nnn];
		}
	}
	// SHIFTS WITHOUT FEEDBACK
	// Line below - This flag will be set later if the remainder is non zero.
	// Non-zero means either corr or uncorr err.  We will know which after decoding.
	remainderDetdErr=0;
	// index to codeword buffer
	for (readCWAddr = numDataBytes; readCWAddr < numDataBytes+numRedunBytes; readCWAddr++) {
		fdbk = 0;
		for (nnn=numRedunWords-1; nnn >=0 ;nnn--){
			fdbkSav = fdbk;
			fdbk = (SR[nnn] >> 24); // 32 - # bits in parallel (unrelated to "m")
			SR[nnn] = (SR[nnn] << 8) ^ fdbkSav;  // 8 # bits in parallel - unrelated to "m"
		}
		fdbk ^= (unsigned int)codeword[readCWAddr];;
		remainBytes[readCWAddr-numDataBytes]=(int)fdbk;
		if (fdbk!=0){
			remainderDetdErr = 1;
		}
	}
	return(remainderDetdErr);
}

static void computeSyndromes(int syndromes[],int numRedunBytes,const int remainBytes[],
							 const int alogTbl[],const int logTbl[],int nParm,int tParm)
{
	//****************************************************************
	//	Function: computeSyndromes
	//
	//  9-10-10 Made few changes for speed
	//
	//  Function to compute syndromes.  Syndromes are computed from a
	//  remainder not from the codeword.  See the comments
	//  for the ComputeRemainder function.
	//
	//  Note that even syndromes are computed from odd syndromes.
	//  I have know that this is possible for a binary
	//  BCH code for about 30 years.  I think that I learned
	//  that this is possible from one of the early ECC books.
	//  Possibly from an early book by Peterson, Berlekamp, or
	//  Shu Lin.  I just found an early reference to this on
	//  page 192 of "Algebraic Coding Theory" by Dr. Elwyn
	//  Berlekamp.
	//
	//  An early reference for computing syndromes from
	//  a remainder is page 160 of the Glover-Dudley 1991 book
	//  "Practical Error Correction Design for Engineers"
	//  REVISED SECOND EDITION.
	//****************************************************************
	int mask, data, evenSNum;
	int jjj,iii,kkk,x,numSyndromes,accumVal,bumpVal;

	// In a real implementation of one fixed code, numSyndromes would be a constant
	numSyndromes=2*tParm;
	// Clear syndromes
	for (kkk=0;kkk<numSyndromes;kkk++){
		syndromes[kkk] = 0;
	}
	for (iii=0;iii < numRedunBytes;iii++){
		mask=0x01;
		for (jjj=0;jjj<8;jjj++){
			data = remainBytes[iii] & mask;
			if (data>0){
				// 9-10-10 Changed for speed
				accumVal=(((numRedunBytes-1)-iii)*8+jjj); // Initialize
				bumpVal=2*accumVal; // Initialize
				// Note increment by 2
				for (kkk=0;kkk<numSyndromes;kkk+=2){
					//  "+1" is for syndrome offset
					//	CAN UNROLL THIS LOOP FOR A LITTLE MORE SPEED
					//  A SWITCH STATEMENT WILL BE REQUIRED.
					//  CAN GET RID OF MOD (FOR SPEED) BY CARRYING A SUM THAT STARTS
					//  OUT AS THE MAXIMUM VALUE OF ((kkk+1)*preComputeVal) AND
					//  SUBTRACTING PRECOMPUTEVAL EACH TIME AND IF VALUE GOES
					//  NEGATIVE ADD IN NPARM.  THE MAX VALUE ABOVE WILL VARY
					//  BASED ON NUMBER OF REDUN BYTES FOR SELECTED PARMS
					//  9-10-10 Changed next few lines for speed
					syndromes[kkk] ^= alogTbl[accumVal];
					accumVal+=bumpVal; // Add 2*preComputeVal
					if (accumVal>=nParm){
						accumVal-=nParm;
						if (accumVal>=nParm){
							accumVal-=nParm;
						}
					}
				}
			}
			mask *=2;
		}
	}
	// Compute even Syndromes from the odd Syndromes
	for (kkk=0;kkk<numSyndromes;kkk+=2){
		x=syndromes[kkk];
		evenSNum=(2*(kkk+1));
		while (evenSNum<=2*tParm){
			// Square "x"
			if (x>0){
				x=alogTbl[2*logTbl[x]];
			}
			syndromes[evenSNum-1]=x;
			evenSNum*=2;
		}
	}
}

static int berMas(int tParm,int sigmaN[],const int syndromes[],int *pErrFlg,
				  const int nParm,const int alogTbl[],const int logTbl[])
//***************************************************************
//	Function: BerMas
//
//  10-3-08 One line was changed and one line was added to change the
//  the code over from supporting a Reed-Solomon code to supporting a
//  binary BCH code.
//
//	This function computes the coefficient's of the ELP using the Berlekamp-
//	Massey algorithm. The reference for this algorithm is
//	"Shift-Register Synthesis and BCH Decoding", IEEE Transactions
//	on Information Theory, IT-15, 122-127, 1969.  This paper can also be
//  found in Blake's 1973 book, "Algebraic Coding Theory"
//  (a collection of papers).  When I wrote the 1999 version of this
//  function, I evaluated several variations of the Berlekamp-Massey
//  algorithm. The variations vary mostly in the initialization of "nn"
//  (0 or 1) and in the algebra used for the decisions.
//
//	I wrote this code originally in 1999 for a Reed-Solomon code.
//  I added the error evaluator polynomial at that time.  My reference
//  for adding the error evaluator polynomial was the thesis
//  by Doug Whiting.  I removed the error evaluator polynomial 10-3-08
//  for a binary BCH, since we do not compute error values for a binary
//  BCH code.  Variable names used here are related to the variable names
//  used in the original Massey article as shown below.
//
//   This Code   Massey Paper
//   ---------   ------------
//	 sigmaTmp    T(D)  D of T(D) is from signal processing and is equiv to using T(x)
//   SigmaN      C(D)
//   sigmaK      B(D)
//   n           N
//   nminusk     x
//   dk          b
//   dn          d
//   Ln          L
//
//  For implementations where the time spent in this function is a
//  dominating factor in performance, a technique due to C.L. Chen can
//  be used for a speed up. If you are interested see "High-Speed
//  Decoding of BCH Codes" by C.L. Chen - IEEE Info Theory VOL. IT-27,
//  NO. 2, March 1981.  Before implementing this technique you need
//  to understand its impact on your miscorrection rate.
//***************************************************************
{
	int logTmpQ,sigmaK[MAXCORR+1];
	int sigmaTmp[MAXCORR+1];
	int dn,dk,Ln;
	int nn,j,lk,nminusk;

	sigmaN[0]=1;sigmaK[0]=1;
	for (nn=1;nn<=tParm;nn++){
		sigmaN[nn]=0;
		sigmaK[nn]=0;
	}
	nminusk=1;dk=1;lk=0;Ln=0;
	// Next cmd line - 2*tParm is number of syndromes
	// Next cmd line - *****incr was changed from 1 (RS) to 2 (bin BCH)*****
	for (nn=0;nn<2*tParm;nn+=2){
		dn=0;
		for (j=0;j<=Ln;j++){
			dn^=ffMult(sigmaN[j],syndromes[nn-j]);
		}
		if (dn==0){
			nminusk++;
		}
		else {
			if (2*Ln>nn){
				if (nn-Ln>(tParm-1)){
					*pErrFlg|=BERMASERR;
					break;
				}
				// Next 2 "if" blks chgd to not use ffMult and ffDiv funs 9-9-2010
				if (dk==0){
					return (DIVZRODIV); // Divide by zero error
				}
				if (dn>0){
					logTmpQ=logTbl[dn]-logTbl[dk];
					if (logTmpQ<0){
						logTmpQ+=nParm;
					}
					for (j=0;j<=lk;j++) {
						if (sigmaK[j]>0){
							sigmaN[nminusk+j]^=alogTbl[logTbl[sigmaK[j]]+logTmpQ];
						}
					}
				}
				nminusk++;
			}
			else {
				if (nn-Ln>(tParm-1)) {
					*pErrFlg|=BERMASERR;break;
				}
				for (j=0;j<=Ln;j++){
					sigmaTmp[j]=sigmaN[j];
				}
				// Next 2 "if" blks chgd to not use ffMult and ffDiv funs 9-9-2010
				if (dk==0){
					return (DIVZRODIV); // Divide by zero error
				}
				if (dn>0){
					logTmpQ=logTbl[dn]-logTbl[dk];
					if (logTmpQ<0){
						logTmpQ+=nParm;
					}
					for (j=0;j<=lk;j++) {
						if (sigmaK[j]>0){
							sigmaN[nminusk+j]^=alogTbl[logTbl[sigmaK[j]]+logTmpQ];
						}					}
				}
				for (j=0;j<=Ln;j++){
					sigmaK[j]=sigmaTmp[j];
				}
				lk=Ln;
				Ln=nn+1-Ln;
				dk=dn;
				nminusk=1;
			}
		}
		nminusk++;// *****This line added for binary BCH*****
	} // end for
	if (sigmaN[Ln]==0){
		*pErrFlg|=BERMASERR;
	}
	return (Ln);
}

static int chienSearch(int sigmaN[],int Loc[],const int alogTbl[], const int logTbl[],
					   const int LnOrig,const int numCodewordBytes,const int nParm,
					   const int mParmOdd)
{
	//****************************************************************
	//  Function:	ChienSearch
	//
	//  This is an enhanced version of the Chien search algorithm.  As
	//  coded here the function finds roots until the degree (Ln) of
	//  the Error Locator Polynomial (ELP) has been reduced to four,
	//  if "m" is even, or to two, if "m" is odd.  The quartic function
	//  finds the last four roots if "m" is even or the quadratic function
	//  finds the last two roots if "m" is odd.  This function is not
	//  called if the degree of the original ELP is 4 or less("m" even) or 2
	//  or less ("m" odd).  To optimize speed two different techniques are
	//  used to advance the ELP for its next evaluation.
	//
	//  For this code a fast technique is used when there are no zero
	//  coefficients in the ELP.  A slower technique is used when there are
	//  zero coefficients in the ELP. The probability that there are zero
	//  coefficients in the ELP is relatively low.  Also for speed the ELP
	//  is divided down each time a root is found.  Also for speed the inner
	//  most loop of the faster technique is unrolled.  Again for speed, mod
	//  operations are replaced with equivalent but faster operations.
	//  My original knowledge of unrolling loops came from the Earl
	//  Cohen Ph.D. thesis (Berkeley).  See also the Glover-Dudley patent
	//  4,839,896.
	//
	//  The fastest Chien search technique that I am aware of is to
	//  keep the coefficient values in alog form and to use a table
	//  in memory to implement each constant multiplier of the Chien
	//  search.  With this technique combined with loop unrolling
	//  I think you would see a speed up by about a factor of two
	//  for root finding.  The memory space requirement will be very
	//  large for large finite fields.  The large memory requirement
	//  may not be a problem for decoders that run on a PC.
	//
	//  There are software root finding techniques that beat the
	//  speed of even very well written Chien search software.
	//  Some of these techniques are for particular cases such
	//  as for quintic and sextic polynomials.  But others are for
	//  the general case.  Contact the author for more information.
	//
	//  The original reference for the most basic Chien search is
	//  "Cyclic Decoding Procedures for Bose-Chaudhuri-
	//  Hocquenghem Codes", IEEE Transactions on Information Theory,
	//  vol. IT-10, pp 357-363, Oct. 1964.
	//
	//****************************************************************
	int nn,jj,kx,coeffContainsAZero;
	int accum,reg,tmp,errFlg,Ln;
	//
	errFlg=0;
	Ln=LnOrig;
	// Check for a zero coeff before converting to log domain and if find one, set flag
	coeffContainsAZero = 0;
	for (kx=1;kx<=Ln;kx++){
		if (sigmaN[kx]==0){
			coeffContainsAZero = 1;
			break;
		}
	}
	// Convert error locator poly to log domain for Chien Search
	for (nn=1;nn<=Ln;nn++){
		sigmaN[nn] = logTbl[sigmaN[nn]];
	}
	for (nn=0;nn<numCodewordBytes*8;nn++){
		accum = 0;
		// Do simple Chien Search if zero coeffs or if Ln>8.  YOU CANNOT
		// CHANGE THE # 8 IN "Ln>8" WITHOUT ADDING MORE ENTRIES FOR THE
		// UNROLLED LOOP
		if (coeffContainsAZero == 1 || Ln>8){
			for (jj=1;jj<=Ln;jj++){ // One step of Simple Chien Search in this loop
				if (sigmaN[jj] != gblLogZVal){ // Test for log of zero
					accum ^= alogTbl[sigmaN[jj]];
					sigmaN[jj] -= jj;
					if (sigmaN[jj] < 0){ // Compare & subtract is faster than mod
						sigmaN[jj] += nParm;
					}
				}
			}
		}
		else{
			// Unrolled loop
			switch (Ln)
			{
			case 8:accum ^= alogTbl[sigmaN[8]];  // accum is XOR sum of all alogs
				sigmaN[8] -= 8;
				if (sigmaN[8] < 0){
					sigmaN[8] += nParm;
				}
				//lint -fallthrough
			case 7:accum ^= alogTbl[sigmaN[7]];  // accum is XOR sum of all alogs
				sigmaN[7] -= 7;
				if (sigmaN[7] < 0){
					sigmaN[7] += nParm;
				}
				//lint -fallthrough
			case 6:accum ^= alogTbl[sigmaN[6]];  // accum is XOR sum of all alogs
				sigmaN[6] -= 6;
				if (sigmaN[6] < 0){
					sigmaN[6] += nParm;
				}
				//lint -fallthrough
			case 5:accum ^= alogTbl[sigmaN[5]];  // accum is XOR sum of all alogs
				sigmaN[5] -= 5;
				if (sigmaN[5] < 0){
					sigmaN[5] += nParm;
				}
				//lint -fallthrough
			case 4:accum ^= alogTbl[sigmaN[4]];  // accum is XOR sum of all alogs
				sigmaN[4] -= 4;
				if (sigmaN[4] < 0){
					sigmaN[4] += nParm;
				}
				//lint -fallthrough
			case 3:accum ^= alogTbl[sigmaN[3]];  // accum is XOR sum of all alogs
				sigmaN[3] -= 3;
				if (sigmaN[3] < 0){
					sigmaN[3] += nParm;
				}
				//lint -fallthrough
			case 2:accum ^= alogTbl[sigmaN[2]];  // accum is XOR sum of all alogs
				sigmaN[2] -= 2;
				if (sigmaN[2] < 0){
					sigmaN[2] += nParm;
				}
				//lint -fallthrough
			case 1:accum ^= alogTbl[sigmaN[1]];  // accum is XOR sum of all alogs
				sigmaN[1] -= 1;
				if (sigmaN[1] < 0){
					sigmaN[1] += nParm;
				}
				//lint -fallthrough
			default: break;
			}
		}
		if (accum==1){
			Loc[Ln-1]=alogTbl[nn];
			// Convert back to alog domain so we can divide down
			for (kx=1;kx<=Ln;kx++){
				sigmaN[kx] = alogTbl[sigmaN[kx]];
			}
			// Divide down the ELP to eliminate the root just found
			reg=0;
			for (kx=Ln;kx>=0;kx--){
				tmp=ffMult(reg,alogTbl[1]);// The number "1"
				reg=sigmaN[kx]^tmp;
				sigmaN[kx]=tmp;
			}
			Ln--; // Ln must be decremented right here - do not move
			// Check if any coefficients are zero and set flag appropriately
			coeffContainsAZero = 0;
			for (kx=1;kx<=Ln;kx++){
				if (sigmaN[kx]==0){
					coeffContainsAZero = 1;
					break;
				}
			}
			// If degree reduced, special cases will take it from here
			if ((Ln==4 && mParmOdd==0) || (Ln==2 && mParmOdd==1))  // 4 and 2
			{
				// We are still in alog domain so,
				// position the ELP back to its starting point for special cases
				for (kx=1;kx<=Ln;kx++){
					sigmaN[kx]=ffMult(sigmaN[kx],alogTbl[((nn+1)*kx)%nParm]);
				}
				break;
			}
			// Convert back to log domain so we can continue root search
			for (kx=1;kx<=Ln;kx++){
				sigmaN[kx] = logTbl[sigmaN[kx]];
			}
		}
	}
	// If degree of ELP has not been reduced properly
	if ((Ln!=4 && mParmOdd ==0) || (Ln!=2 && mParmOdd ==1)){// 4 and 2
		errFlg|=ROOTSNEQLN;
	}
	return (errFlg);
}

static int rootFindChien(int sigmaN[],int Loc[],const int alogTbl[], const int logTbl[],
						 const int LnOrig,const int numCodewordBytes,const int nParm,
						 const int mParmOdd)
{
	//****************************************************************
	//	Function: rootFindChien
	//
	//  NOTE: There are seperate rootFind functions for Chien and BTA.
	//
	//  This function is the highest level root finding function.  It calls
	//  lower level root finding functions.  There are special functions
	//  called for error locator polynomial (ELP) degrees 1 to 4 ("m" even)
	//  or 1 to 2 ("m" odd).  If the degree of the ELP is greater than
	//  4 ("m" even) or 2 ("m" odd) the Chien search function is called.
	//****************************************************************
	int errFlg;

	// Note that special cases for quartic and cubic used only for m even
	// while the special cases for linear and quadratic used for m odd or even
	if (LnOrig==1){
		linearElp(sigmaN,Loc); // No status returned - no errs detected in function
		errFlg=0;
	}
	else if (LnOrig==2){
		errFlg=quadraticElp(sigmaN,Loc);
	}
	else if (LnOrig==3 && mParmOdd==0){
		errFlg=cubicElp(nParm,sigmaN,Loc,alogTbl);
	}
	else if (LnOrig==4 && mParmOdd==0){
		errFlg=quarticElp(nParm,sigmaN,Loc,alogTbl);
	}
	else
	{
		errFlg=chienSearch(sigmaN,Loc,alogTbl,logTbl,
			LnOrig,numCodewordBytes,nParm,mParmOdd);
		if (errFlg==0 && mParmOdd==0){
			// m even - Chien will have divided down to quartic
			errFlg|=quarticElp(nParm,sigmaN,Loc,alogTbl);
		}
		else if (errFlg==0 && mParmOdd==1){
			// m odd - Chien will have divided down to quadratic
			errFlg|=quadraticElp(sigmaN,Loc);
		}
	}
	return (errFlg);
}
static int ffPFastQuotient(const int mIn[],int mDegIn,const int nIn[],
						   int nDegIn,int quotient[],
						   const int alogTbl[],const int logTbl[],
						   const int nParm)
{
	//****************************************************************
	//	Function: ffPFastQuotient
	//
	//  Function to compute a polynomial quotient in GF(2^m).
	//
	//  mDeg and nDeg are degrees.
	//
	//****************************************************************
	//
	int k,logQDigit,errFlg,mDeg,nDeg;
	int m[MAXCORR+1],n[MAXCORR+1];

	errFlg=0; // Clear error flag

	// Copy the input degrees and arrays so not to modify them
	mDeg=mDegIn;
	nDeg=nDegIn;
	for (k=0;k<=mDeg;k++){
		m[k]=mIn[k];
	}
	for (k=0;k<=nDeg;k++){
		n[k]=nIn[k];
	}
	if ((mDeg==0 && m[0]==0) || (nDeg==0 && n[0]==0)) {
		// '-----zero on entry to ffPFastQuotient-----'
		return (QUOZROONENTRYERR); // Return error
	}
	while (n[nDeg]==0 && nDeg>0) {
		nDeg=nDeg-1;
	}
	while (m[mDeg]==0 && mDeg>0) {
		mDeg=mDeg-1;
	}
	for (k=0;k<=mDeg-nDeg;k++) {  // For q array, # coeffs one gth degree
		quotient[k]=0; // This init loop necessary - do not remove
	}
	while (1) {
		if (mDeg<nDeg) {
			// Quotient will be returned
			break;
		}else{
			// Divide high coeffs to get tmp
			if (n[nDeg]==0){
				return (DIVZRODIV); // Divide by zero error
			}
			if (m[mDeg]>0){
				logQDigit=logTbl[m[mDeg]]-logTbl[n[nDeg]];
				if (logQDigit<0){
					logQDigit+=nParm;
				}
				quotient[mDeg-nDeg]=alogTbl[logQDigit];
				for (k=0;k<=mDeg;k++) {
					if ((k>=mDeg-nDeg) && (n[k-(mDeg-nDeg)]>0)) {
						m[k]=alogTbl[logTbl[n[k-(mDeg-nDeg)]]+logQDigit]^m[k];
					}
				}
			}
			while (m[mDeg]==0 && mDeg>0) {
				mDeg=mDeg-1;
			}
			if (mDeg==0 && m[0]==0) {
				// Quotient will be returned
				break;
			}
		}
	}
	return (0);
}

static int ffPFastMod(const int mIn[],int mDegIn,const int nIn[],
					  int nDegIn,int remainder[],int *pRDeg,
					  const int alogTbl[],const int logTbl[],
					  const int nParm)
{
	//****************************************************************
	//	Function: ffPFastMod
	//
	//  Function to compute a polynomial remainder in GF(2^m).
	//  Computes the poly m modulo the poly n over GF(2^m).
	//
	//  mDeg and nDeg are degrees.
	//
	//****************************************************************
	//
	int k,logQDigit,errFlg,mDeg,nDeg;
	// "+2" because input poly of deg MAXCORR is multiplied by x^2
	int m[MAXCORR+2];
	int n[MAXCORR+1];

	errFlg=0; // Clear error flag
	// Copy the input degrees and arrays so not to modify them
	mDeg=mDegIn;
	nDeg=nDegIn;
	for (k=0;k<=mDeg;k++){
		m[k]=mIn[k];
	}
	for (k=0;k<=nDeg;k++){
		n[k]=nIn[k];
	}
	if ((mDeg==0 && m[0]==0) || (nDeg==0 && n[0]==0)) {
		// '-----zero on entry to ffPFastMod-----'
		return (MODZROONENTRYERR); // Return error
	}
	while (n[nDeg]==0 && nDeg>0) {
		nDeg=nDeg-1;
	}
	while (m[mDeg]==0 && mDeg>0) {
		mDeg=mDeg-1;
	}
	while (1) {
		if (mDeg<nDeg) {
			for (k=0;k<=nDeg-1;k++) {
				if (k<=mDeg) {
					remainder[k]=m[k];
				}else{
					remainder[k]=0;
				}
			}
			*pRDeg=mDeg;
			break;
		}else{
			// Divide high coeffs to get tmp
			if (n[nDeg]==0){
				return (DIVZRODIV); // Divide by zero error
			}
			if (m[mDeg]>0){
				logQDigit=logTbl[m[mDeg]]-logTbl[n[nDeg]];
				if (logQDigit<0){
					logQDigit+=nParm;
				}
				for (k=0;k<=mDeg;k++) {
					if ((k>=mDeg-nDeg) && (n[k-(mDeg-nDeg)]>0)) {
						n[k-(mDeg-nDeg)]=
							alogTbl[logTbl[n[k-(mDeg-nDeg)]]+logQDigit];
						m[k]=(n[k-(mDeg-nDeg)]^m[k]); // add vectors
					}
				}
			}
			while (m[mDeg]==0 && mDeg>0) {
				mDeg=mDeg-1;
			}
			if (mDeg==0 && m[0]==0) {
				for (k=0;k<=nDeg-1;k++) {
					remainder[k]=0;
				}
				*pRDeg=0;
				break;
			}
		}
	}
	return (0);
}

static int ffPFastGcd(const int mIn[],int mDegIn,const int nIn[],
					  int nDegIn,int gcdresult[],int *pDegOut,
					  const int alogTbl[],const int logTbl[],const int nParm)
{
	//****************************************************************
	//	Function: ffPFastGcd
	//
	//  Function to compute a polynomial greatest common divisor
	//  in GF(2^m).  m and n are the GF(2^m) polynomials for which
	//  the gcd is required.
	//
	//  mDeg and nDeg are degrees.
	//
	//  This function implements Euclids algorithm.  MY reference is
	//  page 149 of the Berlekamp book.
	//
	//  E. Berlekamp (1968), Algebraic Coding Theory, McGraw-Hill.
	//
	//****************************************************************
	//
	int k,logQDigit,errFlg,mDeg,nDeg;
	int m[MAXCORR+1],n[MAXCORR+1];

	// Copy the input degrees and arrays so not to modify them
	mDeg=mDegIn;
	nDeg=nDegIn;
	errFlg=0;
	for (k=0;k<=mDeg;k++){
		m[k]=mIn[k];
	}
	for (k=0;k<=nDeg;k++){
		n[k]=nIn[k];
	}
	if ((mDeg==0 && m[0]==0) || (nDeg==0 && n[0]==0)) {
		// '-----zero on entry to ffPFastGcd-----'
		return (GCDZROONENTRYERR); // Return error
	}
	while (n[nDeg]==0 && nDeg>0) {
		nDeg=nDeg-1;
	}
	while (m[mDeg]==0 && mDeg>0) {
		mDeg=mDeg-1;
	}
	while (1) {
		if (nDeg<=mDeg) {
			// Divide high coeffs to get tmp
			if (n[nDeg]==0){
				return (DIVZRODIV); // Divide by zero error
			}
			if (m[mDeg]>0){
				logQDigit=logTbl[m[mDeg]]-logTbl[n[nDeg]];
				if (logQDigit<0){
					logQDigit+=nParm;
				}
				for (k=0;k<=mDeg;k++) {
					if ((k>=mDeg-nDeg) && (n[k-(mDeg-nDeg)]>0)) {
						n[k-(mDeg-nDeg)]=
							alogTbl[logTbl[n[k-(mDeg-nDeg)]]+logQDigit];
						m[k]=(n[k-(mDeg-nDeg)]^m[k]); // add vectors
					}
				}
			}
			while (m[mDeg]==0 && mDeg>0) {
				mDeg=mDeg-1;
			}
			if (mDeg==0 && m[0]==0) {
				for (k=0;k<=nDeg;k++) {
					gcdresult[k]=n[k];
				}
				*pDegOut=nDeg;
				break;
			}
		}else{
			// Divide high coeffs to get tmp
			if (m[mDeg]==0){
				return (DIVZRODIV); // Divide by zero error
			}
			if (n[nDeg]>0){
				logQDigit=logTbl[n[nDeg]]-logTbl[m[mDeg]];
				if (logQDigit<0){
					logQDigit+=nParm;
				}
				for (k=0;k<=nDeg;k++) {
					if ((k>=nDeg-mDeg) && (m[k-(nDeg-mDeg)])){
						m[k-(nDeg-mDeg)]=alogTbl[logTbl[m[k-(nDeg-mDeg)]]+logQDigit];
						n[k]=(m[k-(nDeg-mDeg)]^n[k]); // add vectors
					}
				}
			}
			while (n[nDeg]==0 && nDeg>0) {
				nDeg=nDeg-1;
			}
			if (nDeg==0 && n[0]==0) {
				for (k=0;k<=mDeg;k++) {
					gcdresult[k]=m[k];
				}
				*pDegOut=mDeg;
				break;
			}
		}
	}
	return (0);
}

static int BTA(const int sigmaN[],int Loc[],const int alogTbl[],
			   const int logTbl[],const int LnOrig,const int nParm,
			   const int mParmOdd,const int mParm,const int LogZVal,
			   const int ffSize)
{
	//****************************************************************
	//	Function: BTA
	//
	//  Function that implements the Berlekamp trace method for
	//	finding roots of finite field polynomials.  The finite
	//  field supported by this function is GF(2^m).
	//
	//  My references for the BTA algorithm are a paper [4] by
	//  Berlekamp and three papers [1,2,3] that discuss a BTA derivative
	//  algorithm called BTZ.  In this program I implemented much of the
	//  math from BTZ, but instead of implementing the special root
	//  finding methods of low degree polynomials by Zinoviev, I
	//  implemented alternative methods.  References [1-4] can be found
	//  on the internet, but there may be a fee for some of them.
	//
	//  1) V. Herbert. Efficient root finding of polynomials over fields
	//  of characteristic 2, WEWoRC 2009, INRIA Paris-Rocquencourt.
	//  2) V. Herbert. Efficient root finding of polynomials over fields
	//  of characteristic 2, INRIA Paris-Rocquencourt.
	//  3) B. Biswas, V. Herbert. Efficient root finding of polynomials
	//  over fields of characteristic 2, CRI INRIA Paris-Rocquencourt.
	//  4) E. Berlekamp (1970), Factoring polynomials over large finite
	//  fields, Mathematics of Computation, v. 24, 1970, pp. 713-735.
	//
	//  ---Overview of the Berlekamp Trace Algorithm (BTA)---. The BTA
	//  as defined by Berlekamp, splits the polynomial to be factored
	//  into two factors using a polynomial greatest common divisor
	//  (gcd) function.  Each resulting factor is also split into
	//  two factors and so on until there exist only degree one factors.
	//  To split a polynomial the greatest common divisor function is
	//  performed on the polynomial and a trace polynomial that has as
	//  its roots about half the elements of the finite field employed.
	//
	//  It is faster to stop splitting when the degree of a factor falls
	//  below a threshold and instead to find the roots of such factors
	//  by even faster methods for low degree polynomials. In this
	//  function I stop splitting at degree four for even “m” (“m” of
	//  GF(2^m)) and at degree two for odd “m”.  I am currently using
	//  special root finding algorithms for linear, quadratic, cubic
	//  and quartic polynomials.  Perhaps more time could be squeezed
	//  out of root finding by stopping the splitting process at degree
	//  six or less by using special root finding algorithms for quintic
	//  and sextic polynomials as well.
	//
	//  This function uses several techniques to achieve speed.  For
	//  example it does not call functions to do finite field multiplies
	//  or divides, they are done inline.  And at one point a loop is
	//  unrolled (straight line coded).  There are still some places in
	//  the function where speed can be increased a bit.
	//****************************************************************
	//
	int errFlg,degCF,tmpLogQ,tmpN,tmpLogD;
	int tmpVk,factorTblNxtEntryIdx,logALPHAi;
	int skipFactorCurrIdxInc,tmp[MAXCORR+1];
	int rootsFoundIdx,twoToKx1Pwr,specialCaseFlg;
	int jx,kx,kx0,kx1,kx2,kx3,kx4,degA,degB,junk;
	int TiCoeff,tmpDeg,tmpForSq;
	int factorTblCurrPosIdx;
	int tmpPoly[MAXCORR+1];
	int p[MAXCORR+1];
	int accumResidue[MAXCORR],v[MAXCORR];
	int MDblShiftRowHasAZero[MAXCORR];
	int MResiduesRowHasAZero[MAXMPARM];
	int degTbl[MAXCORR],alphaTbl[MAXCORR],alphaFlgs[MAXMPARM];
	int factorA[MAXCORR+1],factorB[MAXCORR+1],currFactor[MAXCORR+1];
	int workTiModP[MAXCORR],roots[MAXCORR];
	// Note +2.  Need 2 extra spaces cause poly multiplied by x^2
	int tmpV[MAXCORR+2],tmpVBack[MAXCORR+2];
	int MResidues[MAXMPARM][MAXCORR];
	int MDblShift[MAXCORR][MAXCORR];
	int TiModP[MAXCORR+1][MAXCORR+1];
	int factorTbl[MAXCORR][MAXCORR+1];

	errFlg=0; // Clear error flag
	//
	// Flip "sigmaN" & put in "p" to make input format compatible
	// with this function
	for (kx=0;kx<=LnOrig;kx++){
		tmpPoly[kx]=sigmaN[kx];
	}
	for (kx=0;kx<=LnOrig;kx++){
		p[kx]=tmpPoly[LnOrig-kx];
	}
	// ========== CONSTRUCT THE MDblShift MATRIX ===================
	// Using this matrix does the same thing as shifting twice a finite
	// field polynomial shift register implementing the poly "p".  Using
	// the matrix allows the number of multiplies to be cut in about half
	// compared to actual shifting of a software shift register twice.

	// Clear "row has a zero" flgs
	for (kx=0;kx<=LnOrig-1;kx++) {
		MDblShiftRowHasAZero[kx]=0;
	}
	// Create the matrix to be used for double shifts
	for (kx=0;kx<=LnOrig-1;kx++) {
		if (2*kx<LnOrig) {
			// Initialize row to "0"s, then stuff a "1"
			for (jx=0;jx<=LnOrig-1;jx++) {
				MDblShift[kx][jx]=LogZVal; // RHS is log of zero
			}
			MDblShift[kx][2*kx]=0; // RHS is log of one
			MDblShiftRowHasAZero[kx]=1; // Each of these rows has at least one "0"
		}else{
			// Next few lines - Start multiply residue poly by x^2
			tmpV[0]=0;
			tmpV[1]=0;
			// Finish multiply residue poly by x^2
			for (jx=0;jx<=LnOrig-1;jx++) {
				if (MDblShift[kx-1][jx]==LogZVal) {
					tmpV[jx+2]=0;
				}else{
					tmpV[jx+2]=alogTbl[MDblShift[kx-1][jx]];
				}
			}
			// --------------
			// Deg of residue is (LnOrig-1), "+2" is for the inserted zeros
			tmpDeg=LnOrig+1; // This is (LnOrig-1)+2
			errFlg |= ffPFastMod(tmpV,tmpDeg,p,LnOrig,tmpVBack,&junk,
				alogTbl,logTbl,nParm);
			if (errFlg>0){
				return (errFlg);
			}
			for (jx=0;jx<=LnOrig-1;jx++){
				tmpV[jx]=tmpVBack[jx];
			}
			for (jx=0;jx<=LnOrig-1;jx++) {
				if (tmpV[jx]==0) {
					MDblShiftRowHasAZero[kx]=1; // Set "has a zero" flg
				}
				MDblShift[kx][jx]=logTbl[tmpV[jx]];
			}
		}
	}
	// ========== BEGIN USING THE MDblShift MATRIX TO COMPUTE RESIDUES ====================
	// Residues means residues of Trace(Alpha*z) mod p.
	// Clear temp vector for first TiModP which is computed as residues
	// are computed. This is for first TiModP only, all other TiModP
	// vectors are computed in the factor loop.
	for (kx=0;kx<=LnOrig-1;kx++) {
		workTiModP[kx]=0;
	}
	// Clear "row has a zero" flgs
	for (kx=0;kx<=mParm-1;kx++) {
		MResiduesRowHasAZero[kx]=0;
	}
	// Initialize the "v" array
	for (kx=0;kx<=LnOrig-1;kx++) {
		v[kx]=0;
	}
	v[1]=1; // stuff a "1" in the x^1 position
	for (kx1=0;kx1<=mParm-1;kx1++) {  // Loop thru all the residues to be computed
		if (v[0]==0  ) {// This decision for speed
			MResiduesRowHasAZero[kx1]=1; // Set "has a zero" flg
			// Loop # one - faster, no need to check all entries for 0
			for (jx=0;jx<=LnOrig-1;jx++) {  // Add an entry to residue matrix
				workTiModP[jx]=(workTiModP[jx]^v[jx]); // Developing first TiModP
				// RHS is in log form
				MResidues[kx1][jx]=logTbl[v[jx]];
			}
		}else{
			// Loop # two - slower, must check all entries for 0
			for (jx=0;jx<=LnOrig-1;jx++) {  // Add an entry to residue matrix
				if (v[jx]==0) {
					MResiduesRowHasAZero[kx1]=1; // Set "has a zero" flg
				}
				workTiModP[jx]=(workTiModP[jx]^v[jx]); // Developing first TiModP
				// RHS is in log form
				MResidues[kx1][jx]=logTbl[v[jx]];
			}
		}
		// NEXT LINE - THE +1 NEEDED IN BOTH MATLAB AND "C"
		if (kx1+1<mParm ) {// This "if" just skips to end of loop on last loop pass
			for (jx=0;jx<=LnOrig-1;jx++) {
				if (v[jx]!=0 ) {// Test for zero
					// Compute square of v().  Input and output are in antilog form.
					// You can get a further speedup here by using a table to do the
					// squaring.  Use the element to address the table and read the
					// square of the element from the table.
					tmpForSq=logTbl[v[jx]]; // Fetch log
					v[jx]=alogTbl[tmpForSq+tmpForSq]; // Add logs and fetch alog
				}
			}
			// Initialize accumResidue to all zeros
			for (jx=0;jx<=LnOrig-1;jx++) {
				accumResidue[jx]=0;
			}
			for (kx2=0;kx2<=LnOrig-1;kx2++) {
				if (2*kx2<LnOrig) {
					accumResidue[2*kx2]=v[kx2];
				}else{
					if (v[kx2]>0  ) {// If this operand is "0" then block of code is skipped
						// Must take log of RHS and use it for the two loops
						tmpVk=logTbl[v[kx2]];
						if (MDblShiftRowHasAZero[kx2]==1) {
							// In this loop "need" to test for zero
							for (kx3=0;kx3<=LnOrig-1;kx3++) {
								if (MDblShift[kx2][kx3]!=LogZVal) {
									accumResidue[kx3]=
										(accumResidue[kx3]^alogTbl[tmpVk+MDblShift[kx2][kx3]]);
								}
							}
						}else{
							// In this loop do "not" need to test for zero
							// WARNING - Cannot change # in next line
							// without changing the number of lines of
							// replication in the unrolled loop which
							// follows the else
							if (LnOrig>14){
								for (kx3=0;kx3<=LnOrig-1;kx3++) {
									accumResidue[kx3]^=alogTbl[tmpVk+MDblShift[kx2][kx3]];
								}
							}else{
								switch (LnOrig-1)
								{
								case 13:accumResidue[13]^=alogTbl[tmpVk+MDblShift[kx2][13]];
									//lint -fallthrough
								case 12:accumResidue[12]^=alogTbl[tmpVk+MDblShift[kx2][12]];
									//lint -fallthrough
								case 11:accumResidue[11]^=alogTbl[tmpVk+MDblShift[kx2][11]];
									//lint -fallthrough
								case 10:accumResidue[10]^=alogTbl[tmpVk+MDblShift[kx2][10]];
									//lint -fallthrough
								case 9:accumResidue[9]^=alogTbl[tmpVk+MDblShift[kx2][9]];
									//lint -fallthrough
								case 8:accumResidue[8]^=alogTbl[tmpVk+MDblShift[kx2][8]];
									//lint -fallthrough
								case 7:accumResidue[7]^=alogTbl[tmpVk+MDblShift[kx2][7]];
									//lint -fallthrough
								case 6:accumResidue[6]^=alogTbl[tmpVk+MDblShift[kx2][6]];
									//lint -fallthrough
								case 5:accumResidue[5]^=alogTbl[tmpVk+MDblShift[kx2][5]];
									//lint -fallthrough
								case 4:accumResidue[4]^=alogTbl[tmpVk+MDblShift[kx2][4]];
									//lint -fallthrough
								case 3:accumResidue[3]^=alogTbl[tmpVk+MDblShift[kx2][3]];
									//lint -fallthrough
								case 2:accumResidue[2]^=alogTbl[tmpVk+MDblShift[kx2][2]];
									//lint -fallthrough
								case 1:accumResidue[1]^=alogTbl[tmpVk+MDblShift[kx2][1]];
									//lint -fallthrough
								case 0:accumResidue[0]^=alogTbl[tmpVk+MDblShift[kx2][0]];
									//lint -fallthrough
								default: break;
								}
							}
						}
					}
				}
			} // end kx2
			// This loop does "v=accumResidue"
			for (kx4=0;kx4<=LnOrig-1;kx4++) {
				v[kx4]=accumResidue[kx4];
			}
		}
	} // end of computing all residues
	// Place very first TiModP vector in the table
	// All other TiModP vectors will be computed in factor loop
	for (kx=0;kx<=LnOrig-1;kx++) {
		TiModP[0][kx]=workTiModP[kx];
	}
	// ========== ONCE PER CASE INITIALIZATION =========================
	// Initialize some arrays
	for (jx=0;jx<=LnOrig-1;jx++) {
		alphaTbl[jx]=0;
		degTbl[jx]=0;
	}
	// Initialize flags that tell if a TiModP has already been computed
	// The first entryTiModP vector was computed as residues were computed
	alphaFlgs[0]=1;
	for (jx=1;jx<=mParm-1;jx++) {
		alphaFlgs[jx]=0;
	}
	degTbl[0]=LnOrig;
	// Initialize factorTbl
	factorTblCurrPosIdx=0; // Point to 1st entry of factor table
	for (kx=0;kx<=LnOrig;kx++) {
		// Initialize first row of factorTbl to input poly
		factorTbl[0][kx]=p[kx];
	}
	factorTblNxtEntryIdx=1; // Init
	rootsFoundIdx=0;  // Init
	// ========== BEGIN PROCESSING FACTORS =============================
	while (factorTblCurrPosIdx<factorTblNxtEntryIdx) {
		// ---------- Once per factor initialization ---------------------
		// Get, from the table, a factor to split
		for (kx=0;kx<=degTbl[factorTblCurrPosIdx];kx++) { // n+1 coeffs in degree n poly
			currFactor[kx]=factorTbl[factorTblCurrPosIdx][kx];
		}
		logALPHAi=alphaTbl[factorTblCurrPosIdx];
		if (logALPHAi>=mParm) {
			// This is an uncorrectable error - BREAK
			return (LOGALPHAIGTHMPARM);
		}
		// ========== COMPUTE TiModP - JUST IN TIME, IF NEEDED ===============
		// Note that we already computed the first entry of the TiModP table.  We
		// did that as we were computing residues
		if (alphaFlgs[logALPHAi]==0 ) {// If have not computed this TiModP yet
			// Get TiModP by equiv of multiply vector times matrix
			// Initialize array
			for (kx0=0;kx0<=LnOrig-1;kx0++) {
				workTiModP[kx0]=0;
			}
			// TiCoeff must be in log form
			TiCoeff=logALPHAi;
			twoToKx1Pwr=1;
			for (kx1=0;kx1<=mParm-1;kx1++) {
				if (twoToKx1Pwr<LnOrig) {
					// Need to convert TiCoeff to alog form (will never be zero)
					workTiModP[twoToKx1Pwr]=alogTbl[TiCoeff];
				}else{
					if (MResiduesRowHasAZero[kx1]==1) {
						// In this loop "need" to test MResidues(,) for zero
						// TiCoeff will never be zero - no need to test
						for (kx2=0;kx2<=LnOrig-1;kx2++) {
							// The symbol "TiModP" stands for (Trace(ALPHAi*z)) mod p
							// TiCoeff and MResidues must be logs
							// Add logs and fetch alog
							if (MResidues[kx1][kx2]!=LogZVal) {
								workTiModP[kx2]=
									(workTiModP[kx2]^alogTbl[TiCoeff+MResidues[kx1][kx2]]);
							}
						}
					}else{
						// In this loop do "not" need to test MResidues(,) for zero
						// TiCoeff will never be zero - no need to test
						for (kx2=0;kx2<=LnOrig-1;kx2++) {
							// The symbol "TiModP" stands for (Trace(ALPHAi*z)) mod p
							// TiCoeff and MResidues must be logs
							// Add logs and fetch alog
							workTiModP[kx2]=
								(workTiModP[kx2]^alogTbl[TiCoeff+MResidues[kx1][kx2]]);
						}
					}
				}
				for (kx2=0;kx2<=LnOrig-1;kx2++) {
					TiModP[logALPHAi][kx2]=workTiModP[kx2];
				}
				// This code adds logs modulo (ffSize-1)
				TiCoeff=TiCoeff+TiCoeff; // Add logs to square
				if (TiCoeff>=(ffSize-1)) {
					TiCoeff=TiCoeff-(ffSize-1);
				}
				twoToKx1Pwr=2*twoToKx1Pwr;
			} // end for kx1
			alphaFlgs[logALPHAi]=1;
		}
		for (kx=0;kx<=LnOrig-1;kx++) {
			workTiModP[kx]=TiModP[logALPHAi][kx];
		}
		// Detect spceial case - Detect TiModP equal all zeros
		specialCaseFlg=1;
		for (jx=0;jx<=LnOrig-1;jx++) {
			if (workTiModP[jx]>0) {
				specialCaseFlg=0;
				break;
			}
		}
		// Handle special case of residues summing to all zeros
		// ========== EITHER SPLIT CURRENT FACTOR OR PASS IT ON ==========================
		if (specialCaseFlg==1) {
			// ========== CANNOT SPLIT BECAUSE OF SPECIAL CASE - PASS ON CURRENT FACTOR
			degA=degTbl[factorTblCurrPosIdx];
			degB=0;
			for (kx=0;kx<=degA;kx++) { // n+1 coefficients in a degree n polynomial
				factorA[kx]=currFactor[kx];
			}
			factorB[0]=1;
		}else{
			// ========== TRY SPLITTING THE FACTOR ========================================
			// Use gcd to split currFactor into two factors (factorA and factorB).
			// NOTE: Since TiModP is a residue its degree in the below
			// call is one less than the degree of the input poly to this function
			errFlg |= ffPFastGcd(currFactor,degTbl[factorTblCurrPosIdx],
				workTiModP,LnOrig-1,factorA,&degA,alogTbl,logTbl,nParm);
			if (errFlg>0){
				return (errFlg);
			}
			degCF=degTbl[factorTblCurrPosIdx];
			errFlg |= ffPFastQuotient(currFactor,degCF,factorA,degA,
				factorB,alogTbl,logTbl,nParm);
			if (errFlg>0){
				return (errFlg);
			}
			degB=degTbl[factorTblCurrPosIdx]-degA;
		}
		// ========== PUT factorA IN THE FACTOR TABLE IF DEGREE >(2,4) ====================
		skipFactorCurrIdxInc=0;
		if ((mParmOdd==1 && degA>2) || (mParmOdd==0 && degA>4)) {
			// Put factorA in curr pos of factorTbl
			for (kx=0;kx<=degA;kx++) { // n+1 coeffs for degree n poly
				factorTbl[factorTblCurrPosIdx][kx]=factorA[kx];
			}
			alphaTbl[factorTblCurrPosIdx]=logALPHAi+1;
			degTbl[factorTblCurrPosIdx]=degA;
			skipFactorCurrIdxInc=1;
		}else{
			// Flip L-R and div thru by highest coeff. -------------------
			// Note that loop starts at "1". This because we handle
			// hi coeff (lowest in array) seperately.
			// Does not use ffDiv for speed.
			if (factorA[degA]==0){
				return (DIVZRODIV); // Divide by zero error
			}
			tmpLogD=logTbl[factorA[degA]];
			for (kx=1;kx<=degA;kx++) {  // n+1 coefficients in a degree n polynomial
				tmpN=factorA[degA-kx];
				if ((tmpN)!=0){
					tmpLogQ=logTbl[tmpN]-tmpLogD;
					if (tmpLogQ<0){
						tmpLogQ=tmpLogQ+nParm;
					}
					tmp[kx]=alogTbl[tmpLogQ];
				}else{
					tmp[kx]=0;
				}
			}
			tmp[0]=1; // Force highest coefficient to "1"
			// End flip and divide through -------------------------------
			// Call quartic, cubic, quadratic, or linear to find roots
			if (degA==4) {
				errFlg|=quarticElp(nParm,tmp,roots,alogTbl); // 4 errors
			};
			if (degA==3) {
				errFlg|=cubicElp(nParm,tmp,roots,alogTbl); // 3 errors
			};
			if (degA==2) {
				errFlg|=quadraticElp(tmp,roots); // 2 errors
			};
			if (degA==1) {
				linearElp(tmp,roots); // 1 error
			};
			if (errFlg>0){
				return (errFlg);
			}
			// Put roots in Loc array
			for (kx=0;kx<=degA-1;kx++) {
				Loc[rootsFoundIdx]=roots[kx];
				rootsFoundIdx=rootsFoundIdx+1;
			}
		}
		// ========== PUT factorB IN THE FACTOR TABLE IF DEGREE >(2,4) ====================
		if ((mParmOdd==1 && degB>2) || (mParmOdd==0 && degB>4)) {
			if (skipFactorCurrIdxInc==0 ) {// If curr position of factor tbl not already taken for next loop
				// Put factorB in curr pos of factorTbl
				for (kx=0;kx<=degB;kx++) {  // n+1 coeffs for degree n poly
					factorTbl[factorTblCurrPosIdx][kx]=factorB[kx];
				}
				alphaTbl[factorTblCurrPosIdx]=logALPHAi+1;
				degTbl[factorTblCurrPosIdx]=degB;
				skipFactorCurrIdxInc=1;
			}else{
				// Put factorB in next position of factorTbl
				for (kx=0;kx<=degB;kx++) {  // n+1 coeffs for degree n poly
					factorTbl[factorTblNxtEntryIdx][kx]=factorB[kx];
				}
				alphaTbl[factorTblNxtEntryIdx]=logALPHAi+1;
				degTbl[factorTblNxtEntryIdx]=degB;
				factorTblNxtEntryIdx=factorTblNxtEntryIdx+1;
			}
		}else{
			if (factorB[degB]==0){
				return (DIVZRODIV); // Divide by zero error
			}
			tmpLogD=logTbl[factorB[degB]];
			for (kx=1;kx<=degB;kx++) {  // n+1 coefficients in a degree n polynomial
				tmpN=factorB[degB-kx];
				if ((tmpN)!=0){
					tmpLogQ=logTbl[tmpN]-tmpLogD;
					if (tmpLogQ<0){
						tmpLogQ=tmpLogQ+nParm;
					}
					tmp[kx]=alogTbl[tmpLogQ];
				}else{
					tmp[kx]=0;
				}
			}
			tmp[0]=1; // Force lowest coefficient to "1"
			// Call quartic, cubic, quadratic, or linear to find roots
			if (degB==4) {
				errFlg|=quarticElp(nParm,tmp,roots,alogTbl); // 4 errors
			};
			if (degB==3) {
				errFlg|=cubicElp(nParm,tmp,roots,alogTbl); // 3 errors
			};
			if (degB==2) {
				errFlg|=quadraticElp(tmp,roots); // 2 errors
			};
			if (degB==1) {
				linearElp(tmp,roots); // 1 error
			};
			if (errFlg>0){
				return (errFlg);
			}
			// Put roots in Loc array
			for (kx=0;kx<=degB-1;kx++) {
				Loc[rootsFoundIdx]=roots[kx];
				rootsFoundIdx=rootsFoundIdx+1;
			}
		}
		// This "if" decides if use curr factor tbl pos again or advance
		if (skipFactorCurrIdxInc==0) {
			factorTblCurrPosIdx=factorTblCurrPosIdx+1;
		}
		if (factorTblCurrPosIdx>LnOrig) {
			// Uncorrectable Err - factorTblCurrPosIdx exceeds LnOrig
			return (BTAFACTORTBLERR);
		}
	} // ========== END WHILE - END OF FACTORING FOR CURRENT CASE =======================
	// Check if correct number of roots found
	if ((errFlg==0) && (rootsFoundIdx)!=LnOrig) {
		errFlg |= BTANUMROOTSERR;
	}
	if (errFlg==0){
		for (kx1=0;kx1<=LnOrig-1;kx1++) {
			for (kx2=kx1+1;kx2<=LnOrig-1;kx2++) {
				if (Loc[kx1]==Loc[kx2]) {
					return (BTAREPEATEDROOTS);
				}
			}
		}
	}
	return (errFlg);
}

static int rootFindBTA(int sigmaN[],int Loc[],const int alogTbl[], const int logTbl[],
					   const int LnOrig,const int nParm,const int mParmOdd,
					   const int mParm,const int LogZVal,const int ffsize)
{
	//****************************************************************
	//	Function: rootFindBTA
	//
	//  NOTE: There are seperate rootFind functions for Chien and BTA.
	//
	//  This function is the highest level root finding function.  It calls
	//  lower level root finding functions.  There are special functions
	//  called for error locator polynomial (ELP) degrees 1 to 4 ("m" even)
	//  or 1 to 2 ("m" odd).  If the degree of the ELP is greater than
	//  4 ("m" even) or 2 ("m" odd) the Chien search function is called.
	//****************************************************************
	int errFlg;

	// Note that special cases for quartic and cubic used only for m even
	// while the special cases for linear and quadratic used for m odd or even
	if (LnOrig==1){
		linearElp(sigmaN,Loc); // No status returned - no errs detected in function
		errFlg=0;
	}
	else if (LnOrig==2){
		errFlg=quadraticElp(sigmaN,Loc);
	}
	else if (LnOrig==3 && mParmOdd==0){
		errFlg=cubicElp(nParm,sigmaN,Loc,alogTbl);
	}
	else if (LnOrig==4 && mParmOdd==0){
		errFlg=quarticElp(nParm,sigmaN,Loc,alogTbl);
	}
	else
	{
		errFlg=BTA(sigmaN,Loc,alogTbl,logTbl,
			LnOrig,nParm,mParmOdd,mParm,LogZVal,ffsize);
	}
	return (errFlg);
}

static int fixErrors(int codeword[],const int Loc[],const int logTbl[],int Ln,int numCodewordBytes,
					 int numDataBits, int numRedunBits,int nParm)
{
	//****************************************************************
	//	Function: fixErrors
	//
	//	Function to do the actual correction of errors after error
	// locations have been found by the decode function.
	//****************************************************************
	int kx,jx,bitLoc,byteLoc,byteValue,byteBitNum,errFlg;

	errFlg=0;
	for (kx=0;kx<Ln;kx++){
		bitLoc = (((numCodewordBytes*8 - logTbl[Loc[kx]])-1)%nParm);
		// Bounds check fwd displacement because pad bits at end.
		// Note to Neal.
		if (bitLoc>=0 && bitLoc<numDataBits+numRedunBits){
			byteLoc = bitLoc/8;
			byteBitNum = 7 - (bitLoc % 8);
			byteValue=1;
			for (jx=byteBitNum;jx>=1;jx--){
				byteValue *=2;
			}
			codeword[byteLoc] ^= byteValue;
		}
		else{
			errFlg |= CORROUTSIDE;
		}
	}
	return (errFlg);
}

static int bchDecode(int Loc[],const int alogTbl[], const int logTbl[],int FFSize,
					 int tParm,int numCodewordBytes,int nParm,int mParmOdd,
					 int numDataBits,int numRedunBits,int codeword[],
					 int numRedunWords,int numRedunBytes,int numDataBytes,
					 const unsigned int encodeTbl[BYTESTATES][MAXREDUNWDS],
					 int remainBytes[],int syndromes[],int *pErrFlg,
					 const int LogZVal,const int mParm,const int ffsize)
{
	//****************************************************************
	//	Function: bchDecode
	//
	//	This function performs decoding by calling -
	//  - a function to compute a remainder
	//  - a function to compute syndromes from the remainder
	//  - a Berlekamp-Massey function to compute an error
	//    locator polynomial
	//  - a root finding function to find error locations
	//    (roots of the error locator polynomial)
	//
	//  There is some test code in this function.  If you plan to implement
	//  this function then you may want to remove the test code after
	//  testing is complete.  All such code is commented "For testing only"
	//
	//  Array addresses for syndrome symbols and remainder bytes are passed
	//  to this function but on entry to this function these arrays do not
	//  contain useful data.
	//****************************************************************
	int status,remainderDetdErr,Ln,kx;
	int sigmaN[MAXCORR+1];

	*pErrFlg=0;
	status=0;
	//	Tests three entries of the decode tables to determine
	//	if they have been initialized.
	for (kx=2;kx<=4;kx++){
		if (logTbl[alogTbl[FFSize-kx]] != FFSize-kx){
			*pErrFlg|=TBLNOTINIT;
			return(UNCORR); // Return tables not initialized status
		}
	}
	for (kx=0;kx<MAXCORR;kx++){
		// Changed to "LogZVal" 9-9-10
		Loc[kx]=LogZVal; // Set to log of zero
	}
	for(;;){ // Infinite loop - Exit is by "break"
		remainderDetdErr=computeRemainder(codeword,numRedunWords,numRedunBytes,
			numDataBytes,encodeTbl,remainBytes);
		// If remainderDetdErr not 0, CW is not err free - could be corr or uncorr
		if (remainderDetdErr!=0){
			// GET HERE IF REMAINDER INDICATES AN ERROR (NON-ZERO REMAINDER)
			status=CORR;
			computeSyndromes(syndromes,numRedunBytes,remainBytes,
				alogTbl,logTbl,nParm,tParm);
			//	Compute coeff's of ELP using Berlekamp/Massey
			Ln=berMas(tParm,sigmaN,syndromes,pErrFlg,nParm,alogTbl,logTbl);
			gblLnOrig=Ln; // Line for testing only ################################
			for (kx=0;kx<=Ln;kx++){
				gblSigmaOrig[kx]=sigmaN[kx]; // Loop for testing only #############
			}
			if (*pErrFlg!=0){
				gblBerMasUCECntr++; // Line for testing only ######################
				break;
			}
			//	Find the roots of the ELP
			if (gblRootFindOption==1){  // "1" - BTA, "0" - Chien
				*pErrFlg|=rootFindBTA(sigmaN,Loc,alogTbl,logTbl,Ln,nParm,
					mParmOdd,mParm,LogZVal,ffsize);
			}else{
				*pErrFlg|=rootFindChien(sigmaN,Loc,alogTbl,logTbl,
					Ln,numCodewordBytes,nParm,mParmOdd);
			}
			if (*pErrFlg!=0){
				gblRootFindUCECntr++; // Line for testing only ####################
				break;
			}
			//	Fix the errors in the data buffer
			*pErrFlg|=fixErrors(codeword,Loc,logTbl,Ln,numCodewordBytes,
				numDataBits, numRedunBits,nParm);
			if (*pErrFlg!=0){
				gblFixErrorsUCECntr++; // Line for testing only ###################
				break;
			}
		}
		else {
			for (kx=0;kx<2*tParm;kx++){
				syndromes[kx]=0; // Loop for testing only #########################
			}
		}
		break;
	}
	if (*pErrFlg!=0){
		status=UNCORR;
	}
	// Status 0, CORR, or UNCORR (for UNCORR, *pErrFlg further defines FOR TESTING)
	return(status);
}
static void printAppliedErrs()
{
	//****************************************************************
	//	Function: printAppliedErrs
	//****************************************************************
	int kx,bitLocFromEnd;

	printf("\nRaw bit Locs in applyErrors (log form)(from FRONT of gblCodeword)\n");
	for (kx=0;kx<gblLnOrig;kx++){
		printf("%d ",gblRawLoc[kx]);
	}
	printf("\nRaw bit Locs in applyErrors (log form)(from END of gblCodeword)\n");
	for (kx=0;kx<gblLnOrig;kx++){
		bitLocFromEnd=(gblNumCodewordBytes*8-1)-gblRawLoc[kx];
		printf("%d ",bitLocFromEnd);
	}
	printf("\nApplied errors - byte Locs (from FRONT of gblCodeword) and Values");
	for (kx=0;kx<gblNumErrsApplied;kx++){
		printf("\nLoc %d  Val %d ",gblAppliedErrLocs[kx],gblAppliedErrVals[kx]);
	}
}

static void printSigmaN()
{
	//****************************************************************
	//	Function: printSigmaN
	//****************************************************************
	int kx;

	printf("\n***** gblSigmaOrig (Alog format), gblLnOrig*****\n");
	for (kx=0;kx<=gblLnOrig;kx++){
		printf("%d ",gblSigmaOrig[kx]);
	}
	printf("   (%d)",gblLnOrig);
}

static int compareResults()
{
	//****************************************************************
	//	Function: compareResults
	//
	//	This function compares the saved codeword with the codeword
	//  after correction and returns the number of miscompares
	//****************************************************************
	int kx;
	int misCompareCnt;

	// compare the saved gblCodeword with the gblCodeword after correction
	// and return the number of miscompares
	misCompareCnt=0;
	for (kx=0;kx<gblNumCodewordBytes;kx++){
		if (gblCodeword[kx]!=gblCodewordSav[kx]){
			misCompareCnt++;
		}
	}
	return (misCompareCnt);
}

static int ffMult(int opa, int opb)
{
	//****************************************************************
	//	Function: ffMult
	//
	//	This is the finite field multiply function.  Logs of the
	//	operands are fetched from a table.  The logs are added and
	//	the result is used to fetch an alog.
	//****************************************************************
	int tmp;

	if (opa==0 || opb==0){
		return(0);
	}
	tmp =  gblLogTbl[opa]+gblLogTbl[opb];
	if (tmp >= gblNParm){
		tmp -= gblNParm;
	}
	return (gblAlogTbl[tmp]);
}

static int ffInv(int opa, int *pErrFlg)// Pointer to status
{
	//****************************************************************
	//	Function: ffInv
	//
	//	This is the finite field inversion (1/opa) function.  The
	//	log of the operand is fetched from a table.  The log is
	//	subtracted from field size minus 1 and the result is
	//	used to fetch an alog.
	//****************************************************************
	if (opa==0){
		*pErrFlg|=DIVZROINV; // Or into location pointed to
		return(0);
	}
	if (gblLogTbl[opa]==0){
		return(gblAlogTbl[0]);// gblAlogTbl[0] for general case - other types of finite fields
	}
	else {
		return(gblAlogTbl[gblNParm-gblLogTbl[opa]]);
	}
}

static int ffDiv(int opa, int opb, int *pErrFlg)
{
	//****************************************************************
	//	Function: ffDiv
	//
	//	This is the finite field division (opa/opb) function.  The
	//	logs of the operands are fetched.  The logs are then subtracted
	//	and the result is used to fetch an alog.
	//****************************************************************
	int tmp;

	if (opb==0) {
		*pErrFlg|=DIVZRODIV; // Or into location pointed to
		return(0);
	}
	if (opa==0){
		return(0);
	}
	tmp = gblLogTbl[opa]-gblLogTbl[opb];
	if (tmp<0){
		tmp += gblNParm;
	}
	return (gblAlogTbl[tmp]);
}

static void printLogAlogTbls()
{
	//****************************************************************
	//	Function: printLogAlogTbls
	//****************************************************************
	int kx;

	printf("\n***** gblAlogTbl ***** = \n");
	for (kx = 0; kx <gblFFSize; kx++) {
		printf("%d ", gblAlogTbl[kx]);
		if (kx && ((kx % 20) == 0)){
			printf("\n");
		}
	}
	printf("\n");
	printf("\n***** gblLogTbl ***** = \n");
	for (kx = 0; kx <gblFFSize; kx++) {
		printf("%d ", gblLogTbl[kx]);
		if (kx && ((kx % 20) == 0)){
			printf("\n");
		}
	}
}

static void printCgpBits()
{
	//****************************************************************
	//	Function: printCgpBits
	//
	//  Prints code generator polynomial bits
	//****************************************************************
	int kx;

	printf("\n");
	printf("\n***** gblCgpBitArray ***** = \n");
	for (kx = 0; kx <= gblCgpDegree; kx++) {
		printf("%d ", gblCgpBitArray[kx]);
		if ((kx % 20) == 19){
			printf("\n");
		}
	}
	printf("\n");
}

static void printCgpWords()
{
	//****************************************************************
	//	Function: printCgpWords
	//
	//  Prints code generator polynomial feedback constant
	//****************************************************************
	int kx;

	// Print the ints containing the code generator poly (32 bits per int)
	printf("\n***** cgpWordArray ***** = \n");
	for (kx = 0; kx < gblNumRedunWords; kx++){
		printf("%8x ", gblCgpFdbkWords[kx]);
	}
	printf("\n");
	//printf("\n");
}

static void printOrigCodeword()
{
	//****************************************************************
	//	Function: printOrigCodeword
	//
	//  Print the encoded codeword as it was before errors were applied
	//****************************************************************
	int kx;

	printf("\n\nCodeword after encoding = \n");
	for (kx=0;kx<gblNumDataBytes+gblNumRedunBytes;kx++){
		printf("%x-",gblCodewordSav[kx]);
		if (kx % 20 == 19){
			printf("\n");
		}
	}
	printf("\n");
}

static void printCodeword()
{
	//****************************************************************
	//	Function: printCodeword
	//****************************************************************
	int kx;

	printf("\n\nCodeword after correction = \n");
	for (kx=0;kx<gblNumDataBytes+gblNumRedunBytes;kx++){
		printf("%x-",gblCodeword[kx]);
		if (kx % 20 == 19){
			printf("\n");
		}
	}
	printf("\n");
}

static void printRemainBytes()
{
	//****************************************************************
	//	Function: printRemainBytes
	//
	//  Print the remainder bytes
	//****************************************************************
	int kx;

	printf("\nRemainder bytes = \n");
	for (kx=0;kx<gblNumRedunBytes;kx++){
		printf("%d-",gblRemainBytes[kx]);
		if (kx % 20 == 19){
			printf("\n");
		}
	}
	printf("\n");
}

static void printSyndromes()
{
	//****************************************************************
	//	Function: printSyndromes
	//****************************************************************
	int kx;

	printf("\nSyndromes = \n");
	for (kx=0;kx<2*gblTParm;kx++){
		printf("%d-",gblSyndromes[kx]);
		if (kx % 20 == 19){
			printf("\n");
		}
	}
}

static void printLocs()
{
	//****************************************************************
	//	Function: printLocs
	//
	//	Function to print error locations
	//****************************************************************
	int kx,locTmp;

	printf("\nRaw err locations (alog) (found by root search)(from end of CW) = \n");
	for (kx=0;kx<gblLnOrig;kx++){
		printf("%d-",gblLoc[kx]);
		if (kx % 20 == 19){
			printf("\n");
		}
	}
	printf("\nRaw err locations (log) (found by root search)(from end of CW) = \n");
	for (kx=0;kx<gblLnOrig;kx++){
		printf("%d-",gblLogTbl[gblLoc[kx]]);
		if (kx % 20 == 19){
			printf("\n");
		}
	}
	printf("\nAdjusted Err Byte Locations (found by root search)(from front of CW) = \n");
	for (kx=0;kx<gblLnOrig;kx++){
		locTmp = ((((gblNumCodewordBytes)*8
			- gblLogTbl[gblLoc[kx]])-1)%gblNParm)/8;
		printf("%d-",locTmp);
		if (kx % 20 == 19){
			printf("\n");
		}
	}
}

static void printMiscompares()
{
	//****************************************************************
	//	Function: printMiscompares
	//****************************************************************
	int kx;

	for (kx=0;kx<gblNumCodewordBytes;kx++){
		if (gblCodeword[kx]!=gblCodewordSav[kx]){
			printf("\nMiscompare - Loc %d  -  Rcvd %x  Expd %x",
				kx,gblCodeword[kx],gblCodewordSav[kx]);
		}
	}
}
static struct statAndFCnt bchEval(int CWsPerPass,
								  unsigned int mySeed, int randomDataFlg,
								  int doCompareFlg,int *pErrFlg,
								  int minErrsToSim,int maxErrsToSim)
{
	//***************************************************************
	//	Function: bchEval
	//
	//	Function to test the bch encoding and decoding functions
	//***************************************************************
	int dcdStatus,statusExpd,evalStatus;
	int numErrsSimed,misCompareCnt,CWCntr;
	struct statAndFCnt statusAndFCnt;

	//*********************BEGIN PASS LOOP**********************

	evalStatus = 0;
	gblMisCorrCnt = 0;
	misCompareCnt = 0;
	randomSetSeed(mySeed);
	for (CWCntr=0;CWCntr<CWsPerPass;CWCntr++){
		if (randomDataFlg==1){
			// Generate a random data record
			genWriteData();
			// Encode the random data record
			bchEncode(gblEncodeTbl,gblCodeword,gblNumRedunWords,
				gblNumRedunBytes,gblNumDataBytes);
		}
		else {
			// Clear the write codeword
			clearWriteCW();//The all 0's CW is a valid CW, no need encode this path
		}
		// Save the encoded data buffer
		savCodeword();
		// Set expected status
		statusExpd=0;
		// Go pick and apply random errors
		numErrsSimed=applyErrors(minErrsToSim,maxErrsToSim);
		if (numErrsSimed>gblTParm){
			statusExpd=UNCORR;
		}
		else if (numErrsSimed>=1 && numErrsSimed<=gblTParm){
			statusExpd=CORR;
		}
		else {
			statusExpd=ERRFREE;
		}
		dcdStatus=bchDecode(gblLoc,gblAlogTbl,gblLogTbl,gblFFSize,gblTParm,
			gblNumCodewordBytes,gblNParm,gblMParmOdd,
			gblNumDataBits,gblNumRedunBits,gblCodeword,
			gblNumRedunWords,gblNumRedunBytes,gblNumDataBytes,gblEncodeTbl,
			gblRemainBytes,gblSyndromes,pErrFlg,gblLogZVal,gblMParm,gblFFSize);   //  *****DECODE*****
		if (dcdStatus==UNCORR && statusExpd<UNCORR){
			//		Return error
			evalStatus=(UNCORRNOTEXPD+dcdStatus);
			break;
		}
		else if (dcdStatus>ERRFREE && statusExpd==ERRFREE){
			//		Return error
			evalStatus=(EXPDERRFREE+dcdStatus);
			break;
		}
		else if (dcdStatus==ERRFREE && statusExpd>ERRFREE){
			//		Return error
			evalStatus=(EXPDERR);// No need to include dcdStatus here because it is "0"
			break;
		}
		else if (dcdStatus<UNCORR && statusExpd==UNCORR){
			// This is a miscorrection, count it
			// We need to make sure miscorrection rate is not higher than expected
			gblMisCorrCnt++;
			// Note: No return in this path because we will keep looping
		}
		if (doCompareFlg==1){
			misCompareCnt=compareResults();
			if (misCompareCnt>0 && statusExpd<UNCORR){
				//		Return error - 7x
				evalStatus=(COMPAREERR+dcdStatus);
				break;
			}
		}
	}
	//**********************END LOOP**********************
	statusAndFCnt.stat=evalStatus;
	statusAndFCnt.FCnt=CWCntr;  // Failure count
	return(statusAndFCnt);
}

static void correctCWsFromDisk(int toDoCode,int loopAllCWsCnt)
{
	//****************************************************************
	//	Function: correctCWsFromDisk
	//
	//	This function is used to correct codewords from disk.
	//  And on option to write corrected codewords back to disk.
	//****************************************************************
	static unsigned char fileBuff[MAXFILESIZE];
	int numDiskCodewords;
	char inFileName[100];
	char outFileName[100];
	time_t timeStart,timeEnd; // "time_t" is a "typedef" defined in "time.h"
	//                This is in "time.h" -> "typedef long time_t;"
	FILE *infp,*outfp; // --type FILE--    File pointers
	int k1,k2,loops,junk,tmp,errFlg;
	int dcdStatus,errFreeCnt,correctableCnt,unCorrectableCnt;
	size_t count,readLength,writeLength,numElements;
	//
	errFlg=0;
	// These three initializations are to make "PC lint" happy
	errFreeCnt=0;
	correctableCnt=0;
	unCorrectableCnt=0;
	do{
		printf("\nEnter # CWs to read and correct from disk.");
		printf("\nMust be less than or eq # CWs on disk and # CWs * CW length");
		printf("\nin bytes must be less than or eq to %d.\n",MAXFILESIZE);
		(void)scanf_s("%d", &numDiskCodewords);
		tmp=numDiskCodewords*gblNumCodewordBytes;
	}while (numDiskCodewords<1 || tmp>MAXFILESIZE);
	do {
		printf("\nEnter file path and name for READING - Example - C://Folder/File.bin.\n");
		// Was unsuccessful in using scanf_s at this point
		scanf("%s", inFileName); // No "&" - already addr
		infp=fopen(inFileName,"rb");
		if (infp==0){
			printf("*****OPEN ERROR ON READ INPUT FILE*****\n");
		}
	} while (infp==0);
	printf("The program will read from - %s.\n",inFileName);
	if (fseek(infp,0,0)!=0)	{
		printf("\nFile seek error\n");
		fclose(infp);
		printf("\n-----ENTER ANY NUMBER TO EXIT-----\n");
		(void)scanf_s("%d",&junk);
		return;
	}
	readLength=(size_t)(numDiskCodewords*gblNumCodewordBytes);
	numElements=1;
	count=fread(fileBuff,readLength,numElements,infp);//No "&" - its already an addr
	if (count!=numElements){
		printf("\nFile read error");
		printf("\nMake sure you entered the correct number of codewords.");
		printf("\nAlso make sure you entered the correct code parameters.\n");
		fclose(infp);
		printf("\n-----ENTER ANY NUMBER TO EXIT-----\n");
		(void)scanf_s("%d",&junk);
		return;
	}
	fclose(infp);
	printf( "\nBUSY - Correcting codewords from disk.\n\n");
	(void)time( &timeStart ); // Get current time
	for (loops=1;loops<=loopAllCWsCnt;loops++){
		errFreeCnt=0;
		correctableCnt=0;
		unCorrectableCnt=0;
		for (k1=0;k1<=numDiskCodewords-1;k1++){
			for (k2=0;k2<gblNumCodewordBytes;k2++){
				// Copy CW from gblFileBuff to global CW array
				gblCodeword[k2]=fileBuff[k1*gblNumCodewordBytes+k2];
			}
			dcdStatus=bchDecode(gblLoc,gblAlogTbl,gblLogTbl,gblFFSize,gblTParm,
				gblNumCodewordBytes,gblNParm,gblMParmOdd,
				gblNumDataBits,gblNumRedunBits,gblCodeword,
				gblNumRedunWords,gblNumRedunBytes,gblNumDataBytes,gblEncodeTbl,
				gblRemainBytes,gblSyndromes,&errFlg,gblLogZVal,gblMParm,gblFFSize);
			if (toDoCode==2) { // If to write corrected CWs back to disk
				for (k2=0;k2<gblNumCodewordBytes;k2++){
					// Copy corrected CW array back to gblFileBuff
					fileBuff[k1*gblNumCodewordBytes+k2]=(unsigned char)gblCodeword[k2];
				}
			}
			if (loops==loopAllCWsCnt && dcdStatus==0){
				errFreeCnt++;
			}
			if (loops==loopAllCWsCnt && dcdStatus==1){
				correctableCnt++;
			}
			if (loops==loopAllCWsCnt && dcdStatus==2){
				unCorrectableCnt++;
			}
		}
	}
	(void)time( &timeEnd ); // Get current time
	printf("Elapsed Time in Seconds   - %d\n\n",(int)(timeEnd-timeStart));
	printf("-----Error counts for the last loop follow.-----\n");
	printf("Error Free Count = %d \n",errFreeCnt);
	printf("Correctable Error Count = %d \n",correctableCnt);
	printf("unCorrectable Error Count = %d \n",unCorrectableCnt);
	if (toDoCode==2){
		// ############ WRITE THE OUTPUT FILE #####################################
		do {
			printf("\nEnter file path and name for WRITING - Example - C://Folder/File.bin.");
			printf("\nThe file must be a new file - existing files will not be overwritten.\n");
			// Was unsuccessful in using scanf_s at this point
			(void)scanf("%s", outFileName);// No "&" - already addr of array
			outfp=fopen(outFileName,"rb"); // See if file for writing exists already
			if (outfp!=0){ // If the file for writing already exists
				fclose(outfp);
				printf("\n*****FILE FOR WRITING EXISTS ALREADY*****\n");
				outfp=0; // Do this so we will repeat the do-while loop
			}
			else {
				outfp=fopen(outFileName,"wb");
				if (outfp==0){
					printf("\n*****OPEN ERROR ON FILE FOR WRITING*****\n");
				}
			}
		} while (outfp==0);
		printf("\n**** YOU ARE ABOUT TO WRITE THE FILE - %s ****",outFileName);
		printf("\n**** IF YOU DO NOT WANT TO WRITE THIS FILE, TERMINATE THIS PROGRAM ****");
		printf("\n**** IF YOU WISH TO WRITE THE FILE - ENTER ANY NUMBER ****\n");
		(void)scanf_s("%d", &junk);
		writeLength=readLength;
		numElements=1;
		count=fwrite(fileBuff,writeLength,numElements,outfp);//No "&" - already addr
		if (count!=numElements){
			printf("\nFile write error");
			fclose(outfp);
			printf("\n-----ENTER ANY NUMBER TO EXIT-----\n");
			(void)scanf_s("%d",&junk);
			return;
		}
		fclose(outfp);
	}
	// ########################################################################
	printf("\n************ DONE - ENTER ANY NUMBER TO EXIT ***********\n");
	(void)scanf_s("%d",&junk);
	return;
}
static void wrtTestCWsToDisk(int randomDataFlg,int minErrsToSim, int maxErrsToSim)
{
	//****************************************************************
	//	Function: wrtTestCWsToDisk()
	//
	//	This function is used to write test codewords to disk.
	//****************************************************************
	static unsigned char fileBuff[MAXFILESIZE];
	unsigned int seed;
	int numDiskCodewords;
	char outFileName[100];
	FILE *outfp; // --type FILE--    Output file pointer
	int k1,k2,junk,tmp;
	size_t count,writeLength,numElements;
	time_t timeForSeed;
	do{
		printf("\nEnter # test CWs to generate, apply errors, & wrt to disk.");
		printf("\n(The # of CWs) * (CW length)in bytes must be less than");
		printf("\nor eq to %d.\n",MAXFILESIZE);
		(void)scanf_s("%d", &numDiskCodewords);
		tmp=numDiskCodewords*gblNumCodewordBytes;
	}while (numDiskCodewords<1 || tmp>MAXFILESIZE);
	// #########################################################################
	// THIS LOOP ENCODES ALL THE TEST CODEWORDS TO WRITE TO DISK
	(void)time(&timeForSeed); //timeForSeed is a type "time_t" which is type long
	seed=(unsigned int)(timeForSeed % 2147483647); // Constant is 2^31-1
	randomSetSeed(seed);
	for (k1=0;k1<=numDiskCodewords-1;k1++){
		if (randomDataFlg==1){
			// Generate a random data record
			genWriteData();
			// Encode the random data record
			bchEncode(gblEncodeTbl,gblCodeword,gblNumRedunWords,
				gblNumRedunBytes,gblNumDataBytes);
		}
		else {
			// Clear the write codeword
			clearWriteCW();//The all 0's CW is a valid CW, no need encode this path
		}
		// Go pick and apply random errors
		(void)applyErrors(minErrsToSim,maxErrsToSim);
		// We have a test codeword, now put it in the file buffer
		for (k2=0;k2<gblNumCodewordBytes;k2++){
			// Copy test CW array to the file buffer
			fileBuff[k1*gblNumCodewordBytes+k2]=(unsigned char)gblCodeword[k2];
		}
	}
	// Finished putting all the test CWs in the file buffer
	// ############ WRITE THE FILE BUFFER TO FILE #####################################
	do {
		printf("\nEnter file path and name for WRITING - Example - C://Folder/File.bin.");
		printf("\nThe file must be a new file - existing files will not be overwritten.\n");
		// Was unsuccessful in using scanf_s at this point
		(void)scanf("%s", outFileName);//No "&" - already addr
		outfp=fopen(outFileName,"rb"); // See if file for writing exists already
		if (outfp!=0){ // If the file for writing already exists
			fclose(outfp);
			printf("\n*****FILE FOR WRITING EXISTS ALREADY*****\n");
			outfp=0; // Do this so we will repeat the while loop
		}
		else {
			outfp=fopen(outFileName,"wb");
			if (outfp==0){
				printf("\n*****OPEN ERROR ON FILE FOR WRITING*****\n");
			}
		}
	} while (outfp==0);
	printf("\n**** YOU ARE ABOUT TO WRITE THE FILE - %s ****",outFileName);
	printf("\n**** IF YOU DO NOT WANT TO WRITE THIS FILE, TERMINATE THIS PROGRAM ****");
	printf("\n**** IF YOU WISH TO WRITE THE FILE - ENTER ANY NUMBER ****\n");
	(void)scanf_s("%d", &junk);
	writeLength=(size_t)(numDiskCodewords*gblNumCodewordBytes);
	numElements=1;
	count=fwrite(fileBuff,writeLength,numElements,outfp);//No "&" - already addr
	if (count!=numElements){
		printf("\nFile write error");
		fclose(outfp);
		printf("\n-----ENTER ANY NUMBER TO EXIT-----\n");
		(void)scanf_s("%d",&junk);
		return;
	}
	fclose(outfp);
	// ########################################################################
	printf("\n************ DONE - ENTER ANY NUMBER TO EXIT ***********\n");
	(void)scanf_s("%d",&junk);
	return;
}
int main()
{
	//****************************************************************
	//	Function:  main
	//
	//	This is the main function.  It obtains most of the options
	//  from the user.  It calls "correctCWsFromDisk" if the options
	//  indicate that codewords from disk should be corrected.
	//  Otherwise it calls "bchEval" to encode random data, add errors,
	//  and perform correction.  Timing data and statistics are printed.
	//
	//  In this function the little code snippets that get user inputs
	//  was motivated by similar code in the Robert Morelos-Zaragoza
	//  bch decoder on the "ECC Page" web site.
	//****************************************************************
	struct statAndFCnt statusAndFCnt;
	time_t timeForSeed,startTime,stopTime;
	int kx,junk,initStatus,evalStatus,toDoCode,failCWCnt;
	int accumMisCorrCnt,passesToDo,CWsPerPass,tmp,tmpMax;
	int randomDataFlg,doCompareFlg,printTblsFlg,passCntr;
	int loopAllCWsCnt,errFlg;
	int minErrsToSim,maxErrsToSim;

	unsigned int seed,userSeed;

	seed=0;
	userSeed=0;
	// Some of these initializations are to make compiler and lint happy
	errFlg=0;
	minErrsToSim=0;
	maxErrsToSim=0;
	passesToDo=0;
	CWsPerPass=0;
	doCompareFlg=0;
	randomDataFlg=0;
	printf("\n\nFast Encoder-Decoder for binary BCH codes using parallel");
	printf("\ntechniques for encode and syndrome generation.\n");
	printf("\nCopyright (C) Neal Glover 1999,2000,2008,2010,2011 (boglover@msn.com).\n");
	printf("\nSee the source code for the license grant and warranty disclaimer");
	printf("\nand for restrictions on the free license.");
	printf("\n-----------------------------------------------------------------");
	printf("\nCodewords are in byte format on disk.  Highest order byte is");
	printf("\nin location 0.");
	printf("\n-----------------------------------------------------------------");
	printf("\nNOTE:SOME INVALID ENTRIES WILL CAUSE FAST SCROLLING OF THE SCREEN.");
	printf("\nIF THIS HAPPENS YOU WILL NEED TO CLOSE THE PROGRAM AND START OVER.");
	printf("\n-----------------------------------------------------------------");
	printf("\n\nEnter any number to continue.\n");
	(void)scanf_s("%d",&junk);
	// Get the major function the program is to perform on this run
	do{
		// 0,11,22,33 Trying to avoid toDoCode being confused with any other entry
		printf("\nSELECT THE MAJOR FUNCTION THE PROGRAM IS TO PERFORM ON THIS RUN.");
		printf("\nEnter ---0--- to generate CWs, apply errors, perform correction,");
		printf("\nand report elapsed time.");
		printf("\nEnter ---11--- to read all CWs from disk to a buffer, loop");
		printf("\ncorrection of CWs in the buffer, and report elapsed time.");
		printf("\nEnter ---22--- to read CWs from disk, perform correction, and");
		printf("\nwrite the corrected CWs back to a new file on disk.");
		printf("\nEnter ---33--- to generate test CWs and write them to a new");
		printf("\nfile on disk.\n");
		(void)scanf_s("%d", &toDoCode);
	}while ((toDoCode/10>3 || toDoCode/10<0) || ((toDoCode/10)!=(toDoCode % 10)));
	toDoCode=toDoCode % 10;
	if (toDoCode==0){
		printf("\nAt the end of each pass the pgm will print pass info.  This info");
		printf("\nincludes 4 counts that are accumulated over all passes in a run.");
		printf("\nThese counts include uncorrectable error counts for berMas,");
		printf("\nrootFind, and fixErrors.  If you elect to apply more than");
		printf("\n't' errors, it is normal for these counts to be non-zero.");
		printf("\nAlso printed is a count of mis-corrections.  This count");
		printf("\nis very important and must be understood.  This count may also");
		printf("\nbe non-zero if you elect to apply more that 't' errors.");
		printf("\nCompute mis-correction rate for your code and make sure this");
		printf("\nrate is not significantly greater than the rate you compute.");
		printf("\nIf you don't understand, consult a professional.");
		printf("\nEnter any number to continue.\n");
		(void)scanf_s("%d",&junk);
	}
	if (toDoCode!=3){
		do{
			printf("\nEnter 0 to use the fast Chien search root finder.");
			printf("\nEnter 1 to use the Berlekamp trace (BTA) root finder.");
			printf("\nYou may need to test with both to determine which one best");
			printf("\nfits your requirement. As an example, the BTA algorithm");
			printf("\nis faster for a data block size of 1024 bytes and");
			printf("\nGF(2^14)and for between 6 and 64 errors occurring. I have");
			printf("\nnot performed timing tests ourside that range.\n");
			(void)scanf_s("%d", &gblRootFindOption);
		}while (gblRootFindOption < 0 || gblRootFindOption>1);
	}else{
		gblRootFindOption=1; // Should not matter, but set to something
	}
	do{
		printf("\nEnter m of GF(2^m), must be between %d & %d.\n",MINMPARM,MAXMPARM);
		(void)scanf_s("%d", &gblMParm);
	}while (gblMParm<MINMPARM || gblMParm > MAXMPARM);
	// Compute finite field size
	gblFFSize=1;
	for (kx = 1; kx <= gblMParm; kx++){
		gblFFSize *= 2;
	}
	if (gblFFSize > MAXFFSIZE){
		printf("\nffSize > MAXFFSIZE - change the m parameter you entered or ");
		printf("\nchange MAXMPARM and MAXFFSIZE in the source code and restart.\n");
		printf("\n************ ENTER ANY NUMBER TO EXIT ***********\n");
		(void)scanf_s("%d",&junk);
		return(0);
	}
	gblMParmOdd = gblMParm % 2;
	gblNParm=gblFFSize - 1;	// n = (2^m)-1
	// Ask if user wishes to specify the primitive polynomial for the finite field
	// Added next line 9-9-2010
	gblLogZVal=2*gblNParm;
	do{
		printf("\nEnter primitive polynomial in decimal (Example- enter 67 for 1000011).");
		printf("\nEnter 0 to have pgm select the primitive poly.\n");
		(void)scanf_s("%d",&gblFFPoly);
	}while ((gblFFPoly<=gblFFSize || gblFFPoly>=gblFFSize*2) && gblFFPoly!=0);
	if (gblFFPoly>0 && (gblFFPoly<=gblFFSize || gblFFPoly>=2*gblFFSize)){
		printf("\nThe primitive polynomial you entered is not a valid polynomial\n");
		printf("\nof degree m (you entered m above).  Restart and try again.\n");
		printf("\n************ ENTER ANY NUMBER TO EXIT ***********\n");
		(void)scanf_s("%d",&junk);
		return(0);
	}
	// Get correction capability
	do{
		printf("\nEnter 't', the max # of single bit errors the code will be");
		printf("\ndesigned to correct.  t*m must be < 2^m-1 and t must be less");
		printf("\nthan or equal to MAXCORR (a define in the source code).\n");
		(void)scanf_s("%d", &gblTParm);
	}while (gblTParm > MAXCORR || gblTParm*gblMParm>=gblFFSize);
	if (gblFFPoly==0){ // gblFFPoly will be greater than 0 if user entered a poly above
		pickFieldGenPoly(); // PICK FINITE FIELD GENERATOR POLYNOMIAL
	}
	printf("\ngblFFPoly=%d  gblFFSize=%d\n",gblFFPoly,gblFFSize);
	bchInit();  // GENERATE LOG AND ALOG TABLES
	initStatus = chkLogAlogTbls();
	if (initStatus>ZERO){
		printf("\n***** ERROR IN bchInit. *****  initStatus %d.",initStatus);
		printf("\nThe polynomial you entered several steps above may be");
		printf("\nNON-primitive or it may have been entered incorrectly.");
		printf("\nEXIT and try again.\n");
		printf("\n************ ENTER ANY NUMBER TO EXIT ***********\n");
		(void)scanf_s("%d",&junk);
		return(0);
	}
	tmp=genCodeGenPoly();// GENERATE CODE GENERATOR POLYNOMIAL (CGP)
	// Call or previous line will set gblCgpDegree
	if (tmp!=0){
		printf("\n***** Fatal error in cgp generation *****  error flag %d.",tmp);
		printf("\n************ ENTER ANY NUMBER TO EXIT ***********\n");
		(void)scanf_s("%d",&junk);
		return(0);
	}
	gblNumRedunBits = gblCgpDegree; //Do not try to make this # evenly divisible by 8
	gblNumRedunBytes = (gblCgpDegree)/8;
	if ((gblCgpDegree % 8) >0){
		gblNumRedunBytes++;
	}
	// Get number of data bytes
	do{
		tmp=(gblNParm-gblCgpDegree)/8;
		printf("\n\nEnter data length in bytes, must be <= %d. \n",tmp);
		(void)scanf_s("%d", &gblNumDataBytes);
	}while (gblNumDataBytes>tmp || gblNumDataBytes<0);
	gblNumCodewordBytes = gblNumDataBytes+gblNumRedunBytes;
	gblNumDataBits = gblNumDataBytes*8;
	gblKParm = gblNumDataBits;
	gblNumRedunWords = gblNumRedunBytes/4;
	if ((gblNumRedunBytes % 4) >0){
		gblNumRedunWords++;
	}
	if (toDoCode==0 || toDoCode==3){//If pgm to gen CWs, apply errs, correct, report time
		do{
			tmpMax=2*gblTParm;
			if (MAXERRSTOSIM<tmpMax){
				tmpMax=MAXERRSTOSIM; // tmpMax is min(MAXERRSTOSIM,2*gblTParm)
			}
			printf("\nEnter max # errors to apply.");
			printf("\nMust be less than or equal to %d.\n",tmpMax);
			(void)scanf_s("%d", &maxErrsToSim);
		}while (maxErrsToSim>tmpMax || maxErrsToSim<0);
		if(maxErrsToSim>gblTParm){
			printf("\nYou have elected to sim UNCORR as well as CORR errors.\n");
		}
		do{
			printf("\nEnter min # errors to apply.");
			printf("\nMust be less than or equal to %d.\n",maxErrsToSim);
			(void)scanf_s("%d", &minErrsToSim);
		}while (minErrsToSim>maxErrsToSim || minErrsToSim<0);
		if (toDoCode==0){
			do{
				printf("\nEnter 1 for random data.");
				printf("\nEnter 0 for all zeros data and to skip encoding.");
				printf("\nChoose 0 if you want to time decode only.\n");
				(void)scanf_s("%d", &randomDataFlg);
			}while (randomDataFlg!=0 && randomDataFlg!=1);
		}
		else{ // Do this if toDoCode==3
			do{
				printf("\nEnter 1 for random data.");
				printf("\nEnter 0 for all zeros data.\n");
				(void)scanf_s("%d", &randomDataFlg);
			}while (randomDataFlg!=0 && randomDataFlg!=1);
		}
	}
	if (toDoCode==0){
		if (randomDataFlg==1){
			printf("\nEnter 0 to use a random seed, or");
			printf("\nEnter the seed from a failure printout to");
			printf("\ntry to repeat a failure.  In this case");
			printf("\nseed must be >0 and <=4294967295.  This will");
			printf("\ncontrol the seed for the 1st pass only, after");
			printf("\nthat it will be random.\n");
			(void)scanf_s("%d", &userSeed);
		}
		do{
			printf("\nEnter # of codewords to generate per pass.");
			printf("\nThis # should be large for timing accuracy.");
			printf("\nRecommend first try 100000\n");
			(void)scanf_s("%d", &CWsPerPass);
		}while (CWsPerPass<1);
		do{
			printf("\nEnter # of passes to do.");
			printf("\nA report will be printed after each pass.");
			printf("\nTo verify code changes make this # large enough");
			printf("\nto prove the code for the parameters you entered.\n");
			(void)scanf_s("%d", &passesToDo);
		}while (passesToDo<1);
		do{
			printf("\nEnter 1 to compare corrected CW with original.");
			printf("\nEnter 0 to skip the compare for MORE ACCURATE timings.\n");
			(void)scanf_s("%d", &doCompareFlg);
		}while (doCompareFlg!=0 && doCompareFlg!=1);
	}
	loopAllCWsCnt=1; // For all cases except toDoCode==1
	if (toDoCode==1){//If to rd CWs from disk
		do{
			printf("\nEnter # times to loop all codewords.");
			printf("\nMust be less than or equal to %d.",MAXLOOPALLCWSCNT);
			printf("\nIf you have a large number of codewords enter a");
			printf("\nsmall number for loops.\n");
			(void)scanf_s("%d", &loopAllCWsCnt);
		}while (loopAllCWsCnt<1 || loopAllCWsCnt>MAXLOOPALLCWSCNT);
	}
	printf("\nnParm=%d,gblKParm=%d,gblCgpDegree=%d",gblNParm,gblKParm,gblCgpDegree);
	printf("\n\nnumRedunBits=%d,gblNumRedunBytes=%d,gblNumRedunWords=%d",
		gblNumRedunBits,gblNumRedunBytes,gblNumRedunWords);
	printf("\n\nnumDataBytes=%d,gblNumCodewordBytes=%d\n",gblNumDataBytes,gblNumCodewordBytes);
	cvtCgpBitToCgpWord();// CONVERT CODE GENERATOR POLY (CGP) FROM BIT TO WORD FORMAT
	genEncodeTbls();	 // GENERATE ENCODE TABLES
	printCgpBits();
	printCgpWords();
	accumMisCorrCnt = 0;
	gblBerMasUCECntr=0;
	gblRootFindUCECntr=0;
	gblFixErrorsUCECntr=0;
	if (toDoCode==1 || toDoCode==2){ // Rd and corr CWs from disk, on option write back
		correctCWsFromDisk(toDoCode,loopAllCWsCnt); // Load and correct codewords from disk
		return(0); // Done EXIT the program
	}
	if (toDoCode==3){ // Write test CWs to disk
		wrtTestCWsToDisk(randomDataFlg,minErrsToSim,maxErrsToSim); // Load and correct codewords from disk
		return(0); // Done EXIT the program
	}
	if (toDoCode==0){ // Gen CWs, apply errors, and correct (do timings)
		printf( "\nBUSY - Gen'ing and corr'ing CWs - 1st printout will occur shortly.\n");
		for (passCntr=1;passCntr<=passesToDo;passCntr++){
			(void)time(&timeForSeed); //timeForSeed is a type "time_t" which is type long - see above
			if (passCntr==1 && userSeed!=0){
				seed=userSeed; // User can control seed for first pass only
			}
			else {
				do {
					// Constant in next line is 2^31-1
					seed=(unsigned int)((timeForSeed % 2147483647)+127*passCntr);
				}while (seed==0); // The do-while loop is for a very rare case
			}
			(void)time(&startTime); // Start time
			// **** CALL EVAL FUNCTION *****
			statusAndFCnt=bchEval(CWsPerPass,seed,randomDataFlg,doCompareFlg,
				&errFlg,minErrsToSim,maxErrsToSim);
			evalStatus=statusAndFCnt.stat;
			failCWCnt=statusAndFCnt.FCnt;
			(void)time(&stopTime);  // End time
			accumMisCorrCnt += gblMisCorrCnt; // On return add miscorr count
			if (evalStatus>0){
				printf("\n*******************************************************");
				printf("\nERROR RETURNED FROM BCHEVAL.  EVALSTATUS(hex)= %x",evalStatus);
				printf("\n# CWs per pass = %d",CWsPerPass);
				printf("\nPass # = %d ",passCntr);
				printf("\nFailing codeword # = %d",failCWCnt);
				printf("\nRandom number seed = %u",seed);
				printf("\nquadCompTbl = ");
				for (kx=0;kx<gblMParm;kx++){
					printf("%d ",gblQuadCompTbl[kx]);
				}
				printf("\ngblTraceTestVal = %d ",gblTraceTestVal);
				printAppliedErrs();
				printRemainBytes();
				printSyndromes();
				printSigmaN();
				printLocs();
				printMiscompares();
				printf("\nerrFlg(hex) %x  gblMParmOdd %d",
					errFlg,gblMParmOdd);
				do{
					printf("\nEnter 1 to print codeword and log/alog tables.");
					printf("\nIf the finite field is very large, the print out");
					printf("\nwill be very large.");
					printf("\nEnter 0 to skip this printing.\n");
					(void)scanf_s("%d", &printTblsFlg);
				}while (printTblsFlg!=0 && printTblsFlg!=1);
				if (printTblsFlg==1){
					printLogAlogTbls();
					printOrigCodeword();
					printCodeword();
				}
				printf("\nEnter any number to continue testing.\n");
				(void)scanf_s("%d", &junk);
			}
			else
			{
				printOrigCodeword();
				printCodeword();
				printf("\nNO FAILURES OF THE SOFTWARE ON THIS PASS.");
				printf("\nFOR EACH CASE THE STATUS RECEIVED MATCHED THE EXPECTED STATUS.");
			}
			if (randomDataFlg==1){
				printf("\nTime for %d loops of encode and decode = %d seconds.",
					CWsPerPass,(int)(stopTime-startTime));
				printf("\n(Time includes the time to generate random data).");
			}
			else{
				printf("\nTime for %d loops of decode only = %d seconds.",
					CWsPerPass,(int)(stopTime-startTime));
				printf("\n(Time includes the time to clear the codeword).");
				printf("\n(It does not include encode time since data is all zeros).");
			}
			if (doCompareFlg==1) {
				printf("\n(Time includes time for a full codeword compare).");
			}
			else {
				printf("\n(Time does not include time for a codeword compare).");
			}
			printf("\nPass # %d  # CWs per pass %d  accumMisCorrCnt %d",
				passCntr,CWsPerPass,accumMisCorrCnt);
			printf("\ngblBerMasUCECntr %d gblRootFindUCECntr %d gblFixErrorsUCECntr %d",
				gblBerMasUCECntr, gblRootFindUCECntr, gblFixErrorsUCECntr);
			printf("\nRandom number seed = %u\n",seed);
			printf("\nYou are using the following parameters -- ");
			printf("\nm %d t %d max errs to sim %d min errs to sim %d # data bytes %d",
				gblMParm,gblTParm,maxErrsToSim,minErrsToSim,gblNumDataBytes);
			printf("\nCompare Flg %d Random Flg %d",doCompareFlg,randomDataFlg);
			if (gblRootFindOption==0){
				printf("\nYou are using the Chien Search root finder\n");
			}else{
				printf("\nYou are using the Berlekamp trace (BTA) root finder\n");
			}
		}
	}
	printf("\n************ DONE - PROGRAM FINISHED ************");
	printf("\n************ ENTER ANY NUMBER TO EXIT ***********\n");
	(void)scanf_s("%d", &junk);
	return(0);
}











