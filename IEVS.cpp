/* why I am editting this code (a note so I remember):

   I want to be able to quickly add additional methods I find and
   compare them to others, especially identifying the percentages when
   a particular system comes out "on top". I also want to understand
   what Warren means by "RandWinner". I could ask Him by e-mail but I
   believe He is retired and He owes Me nothing; therefore, I do not
   want to impose unnecessarily.
 */
#if defined(MSWINDOWS) && MSWINDOWS
  #include "stdafx.h"  /*needed if Microsoft Windows OS, unwanted if LINUX; will be more of that*/
#endif
#include <algorithm>
#include <array>
#include <climits>
#include <cmath>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <ctype.h>
#include <sys/types.h>
#if defined(MSWINDOWS) && MSWINDOWS
  #include <time.h>
  #include <sys/timeb.h>
  #include <process.h>
#else
  #include <sys/time.h>
  #include <unistd.h>
#endif
#include <stdint.h>
#include <string.h>
#include <typeinfo>

void ensure(bool good, int number);

/* #define NDEBUG    uncomment this line if want to turn off asserts for speed */

/**IEVS ***** Warren D. Smiths's 2007 infinitely extendible voting system comparison engine ****
Copyright (C): Warren D. Smith 2007.  Anybody can use this non-commercially
provided they acknowledge me, notify me (and make available for my unrestricted use)
re any & all code improvements created or bugs found.  If you want commercial use,
then we'll negotiate.

On Unix, Mac OSX, and Linux: Compilation:  gcc IEVS.c -Wall -lm -O6 -o IEVS
    for extra speed (but less safety) add the #define NDEBUG line at top.
On MSWINDOWS:   #define MSWINDOWS to true.  Should work, see notes below by Cary -
    thanks much to David Cary for the MSWINDOWS port!!  (Except that he introduced a bunch
    of bugs with careless typecasts, which hopefully I now have fixed...)
What election methods are available?:    fgrep EMETH IEVS.c
What utility-generators now available?:  fgrep UTGEN IEVS.c

This program should also work (with minor modifications) under microsoft windows OS.
David Cary successfully ported it and his notes on how he did it, are below.
As a side effect, Cary also noticed how to speed up Yee-picture-generation significantly
and fixed a bug in OutputCompressedBarray(). :)
If you do such porting, please report your experience, how you compiled it,
etc.   warren.wds AT gmail.com.

Software Architecture:
I. voting methods.
II. voting strategies.
III. ignorance generators
IV. utility generators

Anybody can add new voting methods, new strategies, or new utility generators.
The idea is to build a "Chinese menu" system which can investigate A*B*C*D kinds of scenarios
BUT the effort to write it is only A+B+C+D.

The information-flow-direction is  IV-->III-->II-->I-->winners-->regrets.

I am initially writing this with not very many of each.  It is a "skeleton"
system which can be fleshed out later by adding more of each of I,II,III,IV.
Please send your contributed routines to warren.wds AT gmail.com .
You should be able to make an implementation of a new method pretty easily by
imitating some already-implemented methods similar to yours.
BUGGY: Rouse??, Copeland??
As-yet unimplemented voting methods include:
ICA, Maxtree, range+runoff(based on range ballots with ties discarded),
RRV+runoff, Sarvo-methods, Banks set, Boehm, Dodgson.
Asset & Candidate-withdrawal-option-IRV (both don't really fit into our model).
To add a new voting method you need to:
add your new EMETH subroutine, and you need to change 3 more lines:
  #define NumMethods PutCorrect#Methodshere
  and the case lines in GimmeWinner() and PrintMethName().
(If you add a new "core" method also must change NumCoreMethods and
perhaps renumber the methods...)

Currently only rating-based, strict-ranking-based, and approval-based methods (& combinations)
are supported - equalities are not allowed in rankings.  ("3-slot" methods
also supported.)  To extend the code to allow
equalities, I suggest adding a boolean vector
with v[x]=true meaning x is really tied with, not below,
the candidate immediately preceding.  That is a later goal.

HonestyStrat() now supports an arbitrary mixture of honest and (one kind of) strategic voters.
More strategies can/will be added later.

3 utility generators now implemented (some parameterized), plus
another "reality-based" utility generator based on
Tideman-et-al's real-world election collection.

Future plans/ideas:  translate from C into D?
(Templates & automatic array bounds checking would really help.)
auto-searcher for property violations?
L1 utils should be skewed distn.
Allow a user-selected voting-method subset.
Targeted BR finder (new user interface to BR finder).
Other strat generators.
Use of 2-candidate elections is knd of silly.
Need to make better controls on summarizers so you can ask it to summarize only
SUBSETS of the data (such as, honest voters only, large #s of candidates only, etc)...
More realistic strategies for plur+runoff, MCA, IRNR, Benham2AppRunoff?
Build good hybrid voting methods?
Add feature to input votes from user and tally election.
Two-humped Gaussian distributions?
RandomWinner regret mean & variance should be computed more exactly.
DMC & the like - do the approvals and rankings correspond with strategy? Should they?
Summarizer: Also look at dominance relations, "approval".
***************************

Now notes by David Cary, 2007-02-14 re modifications to port IEVS 2.58 from LINUX to
Windows XP using Microsoft Visual Studio .Net 2003.

David Cary's Changes (not listing ones WDS did anyhow) include:
 1) Fix compile errors due to passing an int[] to a function that expects a uint[],
 2) Fix compile errors due to passing a uint[] to a function that expects an int[],
 3) Fix runtime error caused by opening the .bmp output file in text mode, to avoid
      LF bytes being written as CR LF combinations.
 4) Add an include for "stdafx.h".

 This program was successfully compiled and run as a C++ program after making the following
 adjustments to the default Visual Studio Project:
 1)  Increase the stack size to 4 MB (under Linker/System)
        (the default is 1 MB, but sizeof(edata) is about 2.2 MB and lives on the stack).

 All reported warnings were eliminated by:
 1) Updating the Visual Studio Project to not report 64-bit issues
      (under "C/C++" / General)
 2) Demoting warning messages #4018 & 4244 to warning level 4 by adding command line options
      /w44018 /w44244 under "C/C++" / Command Line.
      These messages warn about implicit int/uint conversions and conversions from __int64
      to unsigned int.
****************************/

#define until(x) while(!(x))
typedef double real;
typedef void (*driver_t)(void);

const real PI = 3.14159265358979323844;

const int MaxNumCands = 32;
const int MaxNumVoters = 2048;
const int MaxNumIssues = 8;
const int MaxScenarios = 99999;//17424
const int NumMethods = 68;  /*EMETH the # in this line needs to change if add new voting method*/
const int NumSlowMethods = 8;
const int NumFastMethods = (NumMethods-NumSlowMethods);
const int NumUtilGens = 11;  /*UTGEN the # in this line needs to change if add new util gen*/
const int NumCoreMethods = 12;
/* #define CWSPEEDUP to be 1 causes to cause it to be faster; but 0 or leaving it undefined is better for diagnosing bugs */

const uint HTMLMODE = 1U; /*output adding "HTML table" formatting commands*/
const uint TEXMODE = 2U;  /*output adding "TeX table" formatting commands*/
const uint NORMALIZEREGRETS = 4U;  /*Output "normalized" Regrets so RandomWinner=1, SociallyBest=0. */
const uint SORTMODE = 8U;          /*Output with voting methods in sorted order by increasing regrets.*/
const uint SHENTRUPVSR = 16U;  /*Output "Shentrup normalized" Regrets so RandomWinner=0%, SociallyBest=100%. */
const uint OMITERRORBARS = 32U;
const uint VBCONDMODE = 64U; /*vote-based condorcet winner agreement counts; versus true utility based CWs*/
const uint DOAGREETABLES = 128U;
const uint ALLMETHS = 256U;
const uint TOP10METHS = 512U;
uint BROutputMode=0;
struct oneCandidateToTheVoter;
struct oneVoter;

struct oneVotingMethod
{
	std::array<uint, NumMethods> agreementCountWithMethod;
	uint regCount;
	real meanRegret;
	real sRegret;
	uint trueCondorcetAgreementCount;
	uint CondorcetAgreementCount;
	int Winner;
};

typedef std::array<real,MaxNumCands> ArmytageData;
typedef std::array<int,MaxNumCands> ArmytageDefeatsData;
typedef std::array<real,MaxNumCands> ArmytageMarginData;
typedef std::array<int,MaxNumCands> DefeatsData;
typedef std::array<int64_t,MaxNumCands> MarginsData;
typedef std::array<int64_t,MaxNumCands> PairApprovalData;
typedef std::array<int,MaxNumCands> TrueDefeatsData;

/*	oneCandidate:	information about a particular Candidate
 */
struct oneCandidate
{
	DefeatsData DefeatsMatrix;
	TrueDefeatsData TrueDefeatsMatrix;
	ArmytageData ArmytageMatrix;
	ArmytageDefeatsData ArmytageDefeatsMatrix;
	ArmytageMarginData ArmytageMarginsMatrix;
	MarginsData margins;
	uint64_t Ibeat;
	PairApprovalData alsoApprovedWith;
	uint64_t antiPluralityVotes;
	uint approvals;
	int64_t electedCount;
	uint64_t drawCount;
	int64_t BordaVotes;
	uint64_t pluralityVotes;
	real rangeVote;
	int64_t lossCount;
	bool uncovered;
	bool IsASchwartzMember;
	bool IsASmithMember;
	real normalizedRatingSum;
	real utilitySum;
	int64_t voteCountForThisRound;
	bool eliminated;
};

typedef std::array<oneCandidate,MaxNumCands> CandidateSlate;

template< class T >
		int ArgMinArr(uint64_t N, const T Arr[], int RandPerm[]);
template< class T >
		int ArgMaxArr(uint64_t N, const T Arr[], int RandPerm[]);
template< class T >
		int ArgMaxArr(uint64_t N, const oneCandidateToTheVoter (&Candidates)[MaxNumCands], int RandPerm[]);
/*	reset(iteratableCollection):	resets 'iteratableCollection'
 *							with default values
 *	iteratableCollection:	the collection to reset
 */
template <typename T> void reset(T& iteratableCollection) {
	std::fill(iteratableCollection.begin(), iteratableCollection.end(), typename T::value_type());
}
int flipACoin(int choice1, int choice2);
template<class T>
		void PermShellSortDown( uint64_t N, int Perm[], const T Key[] );
template<class T>
		void PermShellSortDown(uint64_t N, int Perm[], const oneCandidateToTheVoter (&Candidates)[MaxNumCands]);
template<class T>
		void PermShellSortDown(uint64_t N, const CandidateSlate& Candidates, const T oneCandidate::*member);
void PrintBRPreamble(void);
void printName(const char *name, bool padding, int spaces);
void PrintSummaryOfNormalizedRegretData(uint scenarios);
void RandomTest(real &s, real &mn, real &mx, real &v, int (&ct) [10], real (*func1)(void), real (*func2)(void));
void RandomTestReport(const char *mean_str, const char *meansq_str, real s, real mn, real mx, real v, int (&ct)[10]);
void resetFavorites(oneVoter (&Voters)[MaxNumVoters], const int64_t& Candidate = 0);
template<class T>
		int Sign(T x);
template<class T>
		int SortedKey(uint64_t N, const int Arr[], const T Key[]);
template<class T>
		int SortedKey(uint64_t N, const int Arr[], const oneCandidateToTheVoter (&Candidates)[MaxNumCands]);
int SortedKey(const int Arr[], const std::array<oneVotingMethod, NumMethods>& methods);
void Test(const char *name, const char *direction, real (*func1)(void), real (*func2)(void), const char *mean_str, const char *meansq_str);
template< class T1 >
		T1 TwiceMedian(uint N, T1 A[] );
#define NullFunction ((real(*)(void))NULL)

/*	oneCandidateToTheVoter:	information about a particular Candidate
 *				from the perspective of a particular
 *				Voter
 */
struct oneCandidateToTheVoter
{
	real actualUtility;
	bool approve;
	bool approve2;
	real perceivedUtility;
	uint64_t ranking;
	real score;
};

/******************** GENERAL PURPOSE routines having nothing to do with voting: ******/

/******** Fns to deal with sets represented by bit-vectors in 1 machine word: ******/
bool SingletonSet(uint x){ return ((x&(x-1))==0); } /*assumes non-empty*/
bool StrictSuperset(uint64_t x, uint64_t y){ return ((x&y)==y && x!=y); }
bool EmptySet(uint x){ return (x==0); }
/****** convenient fns: *******/
/*	Square(x):	returns the square of 'x'
 *	x:		the number to square
 */
template <class T> T Square(const T &x) {
	return (T)pow(x, 2);
}
/*	PosInt(x):	returns 'x' if 'x' is positive and 0 otherwise
 *	x:	the value to examine
 */
uint64_t PosInt(int64_t x)
{
	if(x>0) {
		return x;
	} else {
		return 0;
	}
}
int MaxInt(int a, int b){ return (((a)>(b)) ? (a):(b)); }

/*	GCD(a, b):	returns the greatest common divisor (or factor) of 'a' and 'b'
 *	a:	one value to examine
 *	b:	the other value to examine
 */
uint GCD(uint a, uint b)
{ /*Euclid alg*/
	if(a>b) {
		a %= b;
		if(a==0) {
			return b;
		}
	}
	for(;;) {
		b %= a;
		if(b==0) {
			return a;
		}
		a %= b;
		if(a==0) {
			return b;
		}
	}
}

uint LCMfact[32];

void BuildLCMfact()
{
	uint j,x;
	printf("\nComputing sequence LCM(1,2,...N) for N=1..22:\n");
	LCMfact[0] = 1;
	for(j=1; j<23; j++) { /*LCMfact[23]=5354228880 won't fit in 32 bit word*/
		LCMfact[j] = LCMfact[j-1];
		x = LCMfact[j]%j;
		if(x) {
			LCMfact[j] *= j/GCD(x,j);
		}
	}
	printf("LCMfact[%d]=%u\n", j, LCMfact[j]);
}


/***************************************** other "Artin primes" are
3, 5, 11, 13, 19, 29, 37, 53, 59, 61, 67, 83, 101, 107, 131, 139, 149, 163, 173, 179, 181,
197, 211, 227, 269, 293, 317, 347, 349, 373, 379, 389, 419, 421, 443, 461, 467, 491, 509,
523, 541, 547, 557, 563, 587, 613, 619, 653, 659, 661, 677, 701, 709, 757, 773, 787, 797
i.e. these are primes which have 2 as a primitive root.
This routine finds the greatest Artin prime P with P<=x.  (Not intended to be fast.)
  *******************************************/
uint ARTINPRIME;
/*	FindArtinPrime(x): returns the largest Artin prime less than or equal to 'x'
 *	x:	the maximum possible value for searching for Artin primes
 */
uint FindArtinPrime(uint x)
{
	uint j,p,k;
	int mod2;
	p = x;
	mod2 = p%2;
	if(mod2==0) {
		p--;  /* make p be odd */
	}
	for(; p>2; p -= 2) {
		for((j=4), (k=3); j!=2; j=(j*2)%p, (k += 1)) {
			if(k >= p) {
				return(p);
			}
		}
	}
	printf("FindArtinPrime failed - terminating\n");
	exit(EXIT_FAILURE);
}

/*	PrintNSpaces(N):	prints 'N' space characters to stdout
 *	N:	the number of space characters to print
 */
void PrintNSpaces(int N)
{
	int printed;
	for(printed=0; printed<N; printed++) {
		putchar(' ');
	}
}
/****** convenient constants: *******/
#define BIGINT 0x7FFFFFFF
#define MAXUINT ((uint)-1)
#define MAXUINT64 ((uint64_t)-1)
/* defn works on 8,16,32, and 64-bit machines */


uint32_t BLC32x[60];  /* 32*60=1920 bits of state. Must be nonzero mod P. */
int BLC32NumLeft;

#define lohalf(x) ((uint32_t)(x))
#define A1 ((uint64_t)1284507170)
#define A2 ((uint64_t)847441413)
#define A3 ((uint64_t)650134147)
#define AllF 0xffffffffU
/********************************************************
Warren D. Smith 2001
**********************************************************
Linear congruential pseudo-random number generator mod P,
where P is the enormous prime (578 base-10 digits long;
60 words long if each word is 32 bits)
  P = [2^(48*32) - 1] * 2^(12*32) + 1.
This prime can yield PRNGs suitable for use on
computers with w-bit words, w=8,16,32,48,64,96,128.
The following code is intended for w=32.
The fact that 2^(60*32) = 2^(12*32) - 1 (mod P)
makes modular arithmetic mod P easy, and is the
reason this particular P was chosen.
The period of our generator is P-1.
***********************************************************
Although it is usually easy to detect the nonrandomness of
linear congruential generators because they generate d-tuples
lying on a lattice, in the present case the "grain size" of the
lattice is invisibly small (far less than a least significant bit),
for 32-bit words, if 1<=d<=180. If you want 64-bit words, then need
to concatenate two 32-bit words, and then grain size invisibly small
if 1<=d<=90. These bounds are conservative; in fact I suspect
one is ok, even for 64-bit words, even in up to 1000 dimensions.
***********************************************************
Some even more enormous primes which we have not used are:
[2^{59*32} - 1] * 2^{8 *32} + 1,
[2^{63*32} - 1] * 2^{24*32} + 1,
[2^{69*32} - 1] * 2^{14*32} + 1,
[2^{95*32} - 1] * 2^{67*32} + 1,
[2^{99*32} - 1] * 2^{35*32} + 1;
these would also be suitable for (8,16,32)-bit computers,
and the second of them would also be good for (48,96)-bit computers.
Unfortunately the factorization of P-1 is not known for the last 3
I've listed here, preventing you from being certain you've found a
primitive root mod that P. A still more enormous prime is
  [2^4253 - 1] * 2^4580 + 1    [2659 digits long!]
(where note 2^4253-1 is also prime so that factorizing P-1 is
trivial) but doing arithmetic mod this P is (although still fairly
easy) less pleasant because bit-shifting is required.
*************************************************************/
uint32_t BigLinCong32()
{
	uint32_t y[120];
	int i;
	uint64_t u;

	if(BLC32NumLeft==0) {
		/* Need to refill BLC32x[0..59] with 60 new random numbers: */

 /****************************************************************
 * If BLC32x[0..59] is the digits, LS..MS, of a number in base 2^w,
 * then the following code fragment puts A times that number
 * in y[0..119].  Here
 *  A = 1284507170 * 2^(w*3) + 847441413 * 2^(w*44) + 650134147 * 2^(w*59)
 * is a "good" primitive root mod P, if w=32.
 *****************************************************************/
		for(i=0; i<3; i++) {
			y[i] = 0;
		}
		u=0;
		for(/*i=3*/; i<44; i++) {
			u += A1 * BLC32x[i-3];
			y[i] = lohalf(u);
			u = u>>32;
		}
		for(/*i=44*/; i<59; i++) {
			u += A1 * BLC32x[i-3];
			u += A2 * BLC32x[i-44];
			y[i] = lohalf(u);
			u = u>>32;
		}
		for(/*i=59*/; i<60+3; i++) {
			u += A1 * BLC32x[i-3];
			u += A2 * BLC32x[i-44];
			u += A3 * BLC32x[i-59];
			y[i] = lohalf(u);
			u = u>>32;
		}
		for(/*i=60+3*/; i<60+44; i++) {
			u += A2 * BLC32x[i-44];
			u += A3 * BLC32x[i-59];
			y[i] = lohalf(u);
			u = u>>32;
		}
		for(/*i=60+44*/; i<60+59; i++) {
			u += A3 * BLC32x[i-59];
			y[i] = lohalf(u);
			u = u>>32;
		}
		/*i=60+59=119*/
		y[i] = lohalf(u);
 /*************************************************************
 * If y[0..119] is the digits, LS..MS, of a number in base 2^w,
 * then the following code fragment replaces that number with
 * its remainder mod P in y[0..59]  (conceivably the result will
 * be >P, but this does not matter; it will never be >=2^(w*60)).
 **************************************************************/
		u=1; /*borrow*/
		/* Step 1: y[0..72] = y[0..59] + y[60..119]shift12 - y[60..119]: */
		for(i=0; i<12; i++) {
			u += y[i];
			u += (uint64_t)(uint32_t)~y[60+i];
			y[i] = lohalf(u);
			u = u>>32;
		}
		for(/*i=12*/; i<60; i++) {
			u += y[i];
			u += y[48+i];
			u += (uint64_t)(uint32_t)~y[60+i];
			y[i] = lohalf(u);
			u = u>>32;
		}
		for(/*i=60*/; i<72; i++) {
			u += AllF;
			u += y[48+i];
			y[i] = lohalf(u);
			u = u>>32;
		}
		assert(u>0);
		y[72] = (uint32_t)(u-1); /*unborrow*/

		/*  Step 2: y[0..60] = y[0..59] + y[60..72]shift12  - y[60..72]: */
		u=1; /*borrow*/
		for(i=0; i<12; i++) {
			u += y[i];
			u += (uint64_t)(uint32_t)~y[60+i];
			y[i] = lohalf(u);
			u = u>>32;
		}
		/*i=12*/
		u += y[i] + y[48+i];
		u += (uint64_t)(uint32_t)~y[60+i];
		y[i] = lohalf(u);
		u = u>>32;
		i++;
		for(/*i=13*/; i<25; i++) {
			u += AllF;
			u += y[i];
			u += y[48+i];
			y[i] = lohalf(u);
			u = u>>32;
		}
		for(/*i=25*/; i<60; i++) {
			u += AllF;
			u += y[i];
			y[i] = lohalf(u);
			u = u>>32;
		}
		/*i=60*/
		assert(u>0);
		y[i] = (uint32_t)(u-1); /*unborrow*/

	/*It is rare that any iterations of this loop are needed:*/
		while(y[60]!=0) {
			printf("rare loop\n");
			/*Step 3+:  y[0..60] = y[0..59] + y[60]shift12 - y[60]:*/
			u=1; /*borrow*/
			u += y[0];
			u += (uint64_t)(uint32_t)~y[60];
			y[0] = lohalf(u);
			u = u>>32;
			for(i=1; i<12; i++) {
				u += AllF;
				u += y[i];
				y[i] = lohalf(u);
				u = u>>32;
			}
			/*i=12*/
			u += AllF;
			u += y[i];
			u += y[60];
			y[i] = lohalf(u);
			u = u>>32;
			i++;
			for(/*i=13*/; i<60; i++) {
				u += AllF;
				u += y[i];
				y[i] = lohalf(u);
				u = u>>32;
			}
			/*i=60*/
			assert(u>0);
			y[i] = (uint32_t)(u-1); /*unborrow*/
		}

		/* Copy y[0..59] into BLC32x[0..59]: */
		for(i=0; i<60; i++) {
			BLC32x[i] = y[i];
		}
		/*printf("%u\n", BLC32x[44]);*/
		BLC32NumLeft=60;
	}
	/* (Else) We have random numbers left, so return one: */
	BLC32NumLeft--;
	return BLC32x[BLC32NumLeft];
}

/********************************
MAPLE script to check it works:

w := 32;
A := 1284507170 * 2^(w*3) + 847441413 * 2^(w*44) + 650134147 * 2^(w*59);
P := (2^(48*32) - 1) * 2^(12*32) + 1;
for i from 1 to 11 do
   print( floor((A &^ i mod P) / 2^(44*w)) mod (2^w) );
od;
ap := A &^ 11 mod P;
ap mod (2^w);
quit;
********************
MAPLE output:
847441413
4038410930
102374915
470100141
3896743552
243412576
1911259311
1640083353
4014446395
2679952782
4087228849
and
2475219032

C output:
847441413
4038410930
102374915
470100141
3896743552
243412576
oops
2990118053
2614294516
3539391572
1589778147
1758817216
2847725135 1364008005 3563722108 2020382641 1091616930
*************************/

real Rand01(){ /* returns random uniform in [0,1] */
	return ((BigLinCong32()+0.5)/(1.0 + MAXUINT) + BigLinCong32())/(1.0 + MAXUINT);
}

/*	InitRand(seed):	initializes the psuedo-random number generator
 *	seed:	the seed value for the generator
 */
void InitRand(uint seed)
{ /* initializes the randgen */
	int i;
	typeof(timeval::tv_sec) seed_sec=0;
	int processId=0;
	uint seed_usec=0;
#if defined(MSWINDOWS) && MSWINDOWS
	tm* locTime;
	_timeb currTime;
	time_t now;
#else
	struct timeval tv;
#endif
	if(seed==0) {
		printf("using time of day and PID to generate seed");
#if defined(MSWINDOWS) && MSWINDOWS
		now = time(NULL);
		locTime = localtime(&now);
		_ftime(&currTime);
		seed_usec = currTime.millitm*1001;
		seed_sec = locTime->tm_sec + 60*(locTime->tm_min + 60*locTime->tm_hour);
		processId = _getpid();
#else
		gettimeofday(&tv,0);
		seed_sec = tv.tv_sec;
		seed_usec = tv.tv_usec;
		processId = getpid();
#endif
		seed = (1000000*(uint)seed_sec) + (uint)seed_usec +
			(((uint)processId)<<20) + (((uint)processId)<<10);
		printf("=%u\n", seed);
	}
	for(i=0; i<60; i++) {
		BLC32x[i]=0;
	}
	BLC32x[0] = seed;
	BLC32NumLeft = 0;
	for(i=0; i<599; i++) {
		BigLinCong32();
	}
	printf("Random generator initialized with seed=%u:\n", seed);
	for(i=0; i<7; i++) {
		printf("%.6f ", Rand01());
	}
	printf("\n");
}

/*	TestRand01():	performs 100000 calls to 'randgen()' to test if 'randgen[0,1]'
 *			is behaving as expected
 */
void TestRand01()
{
	int i,y, ct[10];
	real x,s,mx,mn,v;
	s=0.0; v=0.0;
	mn = HUGE;
	mx = -HUGE;
	for(i=0; i<10; i++) {
		ct[i]=0;
	}
	printf("Performing 100000 randgen calls to test that randgen[0,1] behaving ok:\n");
	for(i=0; i<100000; i++) {
		x = Rand01();
		s += x;
		if(mx<x) {
			mx=x;
		}
		if(x<mn) {
			mn=x;
		}
		v += pow((x-0.5), 2);
		y = (int)(x*10.0);
		if((x>=0) && (y<10)) {
			ct[y]++;
		}
	}
	printf("mean=%g(should be 1/2)  min=%f  max=%g   variance=%g(should be 1/12=%g)\n",
			s/100000.0, mn, mx, v/100000.0, 1/12.0);
	printf("counts in 10 bins 0-0.1, 0.1-0.2, etc: ");
	for(i=0; i<10; i++) {
		printf(" %d", ct[i]);
	}
	printf("\n");
}

/*	RandBool():	returns a random bool value
 */
bool RandBool()
{ /* returns random boolean */
	if( Rand01() > 0.5 ) {
		return true;
	}
	return false;
}

real RandExpl(){ /* returns standard exponential (density=exp(-x) for x>0) random deviate */
	real x;
	do{
		x = Rand01();
	}while( x==0.0 );
	return -log(x);
}

void TestRandExpl()
{
	Test("expl randgen", "", RandExpl, NullFunction, "1", "2?");
}

real RandNormal(){ /* returns standard Normal (gaussian variance=1 mean=0) deviate */
	real w, x1;
	static real x2;
	static bool ready = false;
	if(ready){
		ready = false;
		return x2;
	}
	do{
		x1 = 2*Rand01() - 1.0;
		x2 = 2*Rand01() - 1.0;
		w = x1*x1 + x2*x2;
	}while ( w > 1.0 || w==0.0 );
	w = sqrt( (-2.0*log(w)) / w );
	x1 *= w;
	x2 *= w;  /* Now x1 and x2 are two indep normals (Box-Muller polar method) */
	ready = true;
	return x1;
}

void TestNormalRand()
{
	Test("normal randgen", "", RandNormal, NullFunction, "0", "1");
}

/* If a 2D normal [x & y coords of which are i.i.d. standard normals]
 * is known to lie on a ray thru center; this returns the distance along the ray.
 * Equivalently, sin(2*pi*T)*R and cos(2*pi*T)*R are iid standard normals if T is
 * random uniform in [0,1] and if R is the independent output of RandRadialNormal().  ****/
real RandRadialNormal(){
	real w;
	do{ w = Rand01(); }while(w==0.0);
	w = sqrt( -2.0*log(w) );
	return w;
}

void TestRadialNormalRand()
{
	Test("radial normal randgen", "", RandRadialNormal, NullFunction, "~1.25", "2");
}

void TestRadialNormalRand2()
{
	Test("normal randgen", " radially", RandNormal, RandNormal, "~1.25", "2");
}

/*	RandSkew():	returns the mean=0 skewed deviate
 */
real RandSkew()
{ /* returns mean=0 skewed deviate; the fact the mean is really 0, has been confirmed by
   * Monte Carlo to +-0.0001
   */
	const real RECIPRTPI = 0.564189583547756286948079451560772585844050;   /* 1/sqrt(pi) */
	real x,y;
	x = RandNormal();
	y = RandNormal();
	if(x<y) {
		x=y;
	}
	return( 1.21129 * (x-RECIPRTPI) );  /*1.21129 chosen so that variance is 1*/
}

void GenRandNormalArr(uint64_t N, real Arr[]){ /* returns Arr[0..N-1] of standard normal randoms */
	int i;
	for(i=0; i<N; i++){
		Arr[i] = RandNormal();
	}
}

/*	RandInt(N):	returns a random integer less than 'N' but at
 *			least 0
 *	N:		the number of possible integers from which to
 *			choose
 */
uint64_t RandInt(uint64_t N) {
	return (uint64_t)(Rand01()*N);
}

real wkc(int a, int b){
	int si;
	assert(b>0);
	si = (a%2) ? 1 : -1; /*si=pseudorandom sign*/
	return  si*cos( (2*a+1)*PI / (4*b) );
}

real wks(int a, int b){
	assert(b>0);
	return sin( (2*a+1)*PI / (4*b) );
}

/* This is a goofy probability density on N-vectors.  Its mean is the zero vector.
 * The variance of each component is 1.  Its support set is all of N-space.
 * It is NOT centrosymmetric.    It has 3 or fewer modes.
 * Centrosymmetric densities (e.g. normal) have the disadvantage that they cause every Condorcet
 * voting method to be identical (and cause Condorcet winner always to exist)
 * with distance-based utilities...  which is not realistic.  Hence this.
 **********/
void GenRandWackyArr(int N, real Arr[])
{ /* returns Arr[0..N-1] of skew randoms */
	const real RT85 = 1.264911064067351732799557417773;  /* sqrt(8/5) */
	const real RT25 = 0.632455532033675866399778708886;  /* sqrt(2/5) */
	int k, which;
	which = (int)RandInt(3);
	switch(which) {
	case(0) :
		for(k=0; k<N; k++) {
			Arr[k] = RandSkew();
		}
		break;
	case(1) :
		for(k=0; k<N; k++) {
			Arr[k] =  wkc(k, N) + RT85*wks(k, N) * RandNormal();
		}
		break;
	case(2) :
		for(k=0; k<N; k++) {
			Arr[k] = -wkc(k, N) + RT25*wks(k, N) * RandNormal();
		}
		break;
	default :
		printf("impossible case in GenRandWackyArr\n");
		exit(EXIT_FAILURE);
	}
}

void TestsOfRand(){
	TestRand01(); TestNormalRand(); TestRandExpl();
	TestRadialNormalRand();
	TestRadialNormalRand2();
}

/*	IsPerm(N, Perm)	returns true if each integer value from 0 to
 *			N-1 appears exactly once in 'Perm' and false
 *			otherwise
 *	N:	the number of unique non-negative integers, starting
 *		with 0, expected to appear in Perm
 *	Perm:	an array of integers
 */
bool IsPerm( uint64_t N, const uint Perm[] )
{ /* true if is a perm of [0..N-1] */
	int i;
	int ct[MaxNumCands];
	assert(N<MaxNumCands);
	for(i=0; i<(int)N; i++) {
		ct[i] = 0;
	}
	for(i=0; i<(int)N; i++) {
		ct[ Perm[i] ]++;
	}
	for(i=0; i<(int)N; i++) {
		if(ct[i]!=1) {
			return false;
		}
	}
	return true;
}

/*	IsPerm(N, Perm)	returns true if each integer value from 0 to
 *			N-1 appears exactly once in 'Perm' and false
 *			otherwise
 *	N:	the number of unique non-negative integers, starting
 *		with 0, expected to appear in Perm
 *	Perm:	an array of integers
 */
bool IsPerm( uint64_t N, const oneCandidateToTheVoter (&Candidates)[MaxNumCands] )
{ /* true if is a perm of [0..N-1] */
	int i;
	int ct[MaxNumCands];
	assert(N<MaxNumCands);
	for(i=0; i<(int)N; i++) {
		ct[i] = 0;
	}
	for(i=0; i<(int)N; i++) {
		ct[ Candidates[i].ranking ]++;
	}
	for(i=0; i<(int)N; i++) {
		if(ct[i]!=1) {
			return false;
		}
	}
	return true;
}

void MakeIdentityPerm( uint64_t N, uint Perm[] ){
	int i;
	for(i=0; i<(int)N; i++){ Perm[i] = i; }
}

/*	RandomlyPermute(N, RandPerm)	randomly permutes the each
 *					element, 'i' with a random
 *					element, 'j', with 0<=j<=i,
 *					for i = N-1, N-2, N-3, ...,
 *					2, 1
 *	N:		one more than the number of times to permute RandPerm
 *	RandPerm:	the array to permute
 */
void RandomlyPermute( uint64_t N, uint RandPerm[] ){ /* randomly permutes RandPerm[0..N-1] */
	uint64_t i;
	uint64_t j;
	for(i=N-1; i>0; i--){
		j = RandInt((uint)i);
		std::swap(RandPerm[i], RandPerm[j]);
	}
	assert(IsPerm(N,RandPerm));
}

/******* vector handling: **********/
/*	CopyArray(N, src, dest):	copies 'N' elements from 'src[]' into 'dest[]';
 *					only elements from index 0 to 'N-1' are copied;
 *					both 'src[]' and 'dest[]' are expected to have
 *					at least 'N' elements
 *	N:	number of elements to copy
 *	src:	the source array
 *	dest:	the destination array
 */
template <class T> void CopyArray(uint64_t N, const T src[], T dest[] ) {
	uint64_t i;
	for(i=0; i<N; i++) {
		dest[i] = src[i];
	}
}

/******* vector handling: **********/
/*	CopyArray(N, theCandidates, dest, member):	copies 'N'
 *							elements from
 *							'theCandidates[0..N-1].member'
 *							into 'dest[]';
 *							only elements
 *							from index 0
 *							to 'N-1' are
 *							copied; both
 *							'theCandidates'
 *							and 'dest[]'
 *							are expected
 *							to have at
 *							least 'N' elements
 *	N:		number of elements to copy
 *	theCandidates:	the slate of Candidates from which to copy
 *	dest:		the destination array
 *	member:		the member of Each Candidate to copy
 */
template <class T> void CopyArray(uint64_t N, const CandidateSlate& theCandidates, T dest[], T oneCandidate::*member ) {
	uint64_t i;
	for(i=0; i<N; i++) {
		dest[i] = theCandidates[i].*member;
	}
}

/*	FillArray(N, Arr, Filler):	assigns 'Filler' to 'N' elements of 'Arr[]';
 *					only elements from index 0 to 'N-1' are filled;
 *					Arr[]' is expected to have at least 'N' elements
 *	N:	number of elements to copy
 *	Arr:	the array to fill
 *	Filler:	the value to assign to the alements of 'Arr[]'
 */
template <class T> void FillArray(uint64_t N, T Arr[], T Filler) {
	uint64_t i;
	for(i=0; i<N; i++) {
		Arr[i] = Filler;
	}
}

/*	ZeroArray(N, Arr):	sets the first 'N' elements of 'Arr'
 *				to 0
 *	N:	number of elements to set to 0
 *	Arr:	array of elements to set
 */
template <class T>
void ZeroArray(uint64_t N, T Arr[]) {
	FillArray(N, Arr, T(0));
}

/* Assumes RandPerm[0..N-1] contains perm. Returns index of random max entry of Arr[0..N-1]. */
int ArgMaxUIntArr(uint64_t N, const uint Arr[], int RandPerm[] ){
	uint maxc;
	int i, r, winner;
	winner = -1;
	maxc = 0;
	RandomlyPermute( N, (uint*)RandPerm );
	for(i=0; i<(int)N; i++){
		r = RandPerm[i];
		if(Arr[r]>=maxc){
			maxc=Arr[r];
			winner=r;
		}
	}
	assert(winner>=0);
	return(winner);
}

/* Assumes RandPerm[0..N-1] contains random perm and MaxInd is index for Arr[] yielding max.
 * Returns index of second-max. */
int Arg2MaxUIntArr(uint64_t N, const uint Arr[], const int RandPerm[], int MaxInd ){
	uint maxc;
	int i, r, winner;
	winner = -1;
	maxc = 0;
	for(i=0; i<(int)N; i++){
		r = RandPerm[i];
		if(Arr[r]>=maxc && r!=MaxInd){
			maxc=Arr[r];
			winner=r;
		}
	}
	assert(winner>=0);
	return(winner);
}

uint RandCandPerm[MaxNumCands]; /* should initially contain 0..NumCands-1 */

/*	SecondMaximum(count, allCandidates, member, maximumIndex):	returns the index
 *									of the entry in
 *									'allCandidates' with
 *									the second highest
 *									'member'
 *	count:		the maximum number of entries in 'allCandidates' to examine
 *	allCandidates:	the slate of Candidates to examine
 *	member:		the data member to use for comparison
 *	maximumIndex:	the index into 'allCandidates' presumed to have the hightest value
 *			of 'member'
 */
template<class T> int SecondMaximum(uint64_t count, const CandidateSlate& allCandidates, const T oneCandidate::*member, int maximumIndex ) {
	T maxc;
	int i;
	int r;
	int winner;
	bool typeIsUnsigned = (typeid(T) == typeid(uint) || typeid(T) == typeid(uint64_t));
	winner = -1;
	if(typeid(T)==typeid(uint)) {
		maxc = (T) 0;
	} else {
		maxc = (T)(-HUGE);
	}
	for(i=0; i<(int)count; i++){
		r = RandCandPerm[i];
		if(((typeIsUnsigned && allCandidates[r].*member>=maxc) || allCandidates[r].*member>maxc) && r!=maximumIndex){
			maxc=allCandidates[r].*member;
			winner=r;
		}
	}
	assert(winner>=0);
	return(winner);
}

/*	ScaleVec(N, a, scalefac):	scales the first 'N' elements of 'a[]' by
 *					'scalefac'; 'a[]' is expected to have at least
 *					'N' elements
 *	N:		the number of elements to scale
 *	a:		the array of values to scale
 *	scalefac:	the scaling factor
 */
template <class T> void ScaleVec(uint64_t N, T a[], T scalefac) {
	int i;
	for(i=0; i<(int)N; i++) {
		a[i] *= scalefac;
	}
}

class range
{
private:
	int64_t last;
	int64_t iter;

public:
	range(int64_t end, int64_t start = 0):
	last(end),
	iter(start)
	{
		ensure(end>=start,30);
	}

	// Iterable functions
	const range& begin() const { return *this; }
	const range& end() const { return *this; }

	// Iterator functions
	bool operator!=(const range&) const { return iter < last; }
	void operator++() { ++iter; }
	int64_t operator*() const { return iter; }
};

/*	ScaleRegrets(methods, scalefac):	scales the 'sRegret'
 *						of each member of
 *						'methods' by 'scalefac'
 *	methods:	the voting methods with 'sRegret' values to
 *			scale
 *	scalefac:	the scaling factor
 */
void ScaleRegrets( std::array<oneVotingMethod, NumMethods>& methods, real scalefac) {
	for(oneVotingMethod& m : methods) {
		m.sRegret *= scalefac;
	}
}

/*Donald Knuth's "The Art of Computer Programming, Volume 2: Seminumerical Algorithms", section 4.2.2
describes how to compute mean and standard deviation using a recurrence relation, like this:
M(1) = x(1),   M(k) = M(k-1) + (x(k) - M(k-1)) / k
S(1) = 0,      S(k) = S(k-1) + (x(k) - M(k-1)) * (x(k) - M(k))
for 2 <= k <= n, then
sigma = sqrt(S(n) / (n - 1))
Attributes this method to B.P. Welford, Technometrics 4 (1962) 419-420.
*******/
void WelfordUpdateMeanSD(real NewDatum, int *Count, real *M, real *S){
	real OldMean;
	OldMean = *M;
	(*Count)++;
	*M += (NewDatum - OldMean)/(*Count);
	*S += (NewDatum - OldMean)*(NewDatum - *M);
	return;
}

void WelfordUpdateMeanSD(real newDatum, oneVotingMethod& theMethod) {
	real& meanRegret = theMethod.meanRegret;
	uint& regCount = theMethod.regCount;
	real& sRegret = theMethod.sRegret;
	real existingMean = meanRegret;

	regCount++;
	meanRegret += (newDatum - existingMean)/regCount;
	sRegret += (newDatum - existingMean)*(newDatum - meanRegret);
	return;
}

/*	DistanceSquared(N, a, b):	returns the sum of the square of the differences
 *					between corresponding elements of 'a[]' and
 *					'b[]'; both arrays are expected to have at least
 *					'N' elements
 *	N:	the number of elements to square and sum, in that order
 *	a:	the first array of values
 *	b:	the second array of values
 */
real DistanceSquared(uint N, const real a[], const real b[] )
{
	real d = 0.0;
	int i;
	for(i=0; i<(int)N; i++) {
		d += pow(( a[i]-b[i] ), 2);
	}
	return d;
}

/*	L1Distance(N, a, b):	returns the sum of the magnitude of the differences
 *				between corresponding elements of 'a[]' and 'b[]'; both
 *				arrays are expected to have at least 'N' elements
 *	N:	the number of magnitudes to sum
 *	a:	the first array of values
 *	b:	the second array of values
 */
real L1Distance(uint64_t N, const real a[], const real b[] )
{
	real d = 0.0;
	int i;
	for(i=0; i<(int)N; i++) {
		d += fabs( a[i]-b[i] );
	}
	return d;
}

real LpDistanceSquared(uint N, const real a[], const real b[], real Lp )
{
	real d = 0.0;
	int i;
	assert(Lp >= 1.0);
	if(Lp==1.0) {
		return pow((L1Distance(N,a,b)), 2);
	}
	if(Lp==2.0) {
		return DistanceSquared(N,a,b);
	}
	for(i=0; i<(int)N; i++) {
		d += pow( fabs( a[i]-b[i] ), Lp );
	}
	return pow( d, 2.0/Lp );
}

/*	SumRealArray(N, a):	returns the sum of the first 'N' elements of 'a[]';
 *				'a[]' is presumed to have at least 'N' elements
 *	N:	the number of elements to sum
 *	a:	the array of values to sum
 */
real  SumRealArray( uint64_t N, const real a[] )
{
	real s = 0.0;
	int i;
	for(i=0; i<(int)N; i++) {
		s += a[i];
	}
	return s;
}

/*	DotProd(N, a, b):	returns the dot product of the first 'N' elements of
 *				'a[]' and 'b[]'; both 'a[]' and 'b[]' are presumed to
 *				have at least 'N' elements
 *	N:	the number of elements to use for the dot product
 *	a:	the first array of values
 *	b:	the second array of values
 */
real DotProd(uint N, const real a[], const real b[] )
{
	real d = 0.0;
	int i;
	for(i=0; i<(int)N; i++) {
		d += a[i]*b[i];
	}
	return d;
}

const int ShellIncs[] = {1750, 701, 301, 132, 57, 23, 10, 4, 1, 0};
/* Marcin Ciura: Best Increments for the Average Case of Shellsort,
   13th International Symposium on Fundamentals of Computation Theory,
   Riga, Latvia, 22-24 August 2001;
   Springer Lecture Notes in Computer Science #2138, pp.106-117.
Here 1750 is unsure and how the sequence continues past 1750 is unknown.
 ***/

/* Rearranges Perm[0..N-1] so Key[Perm[0..N-1]] is in increasing order: */
void RealPermShellSortUp(uint N, int Perm[], const real Key[])
{
	int h,i,j,k;
	int x;
	int sortedTrend;
	for(k=0; (h=ShellIncs[k])>0; k++) {
		for(i=h; i<(int)N; i++) {
			x=Perm[i];
			for(j=i-h; j>=0 && Key[Perm[j]]>Key[x]; j -= h) {
				Perm[j+h]=Perm[j];
			}
			Perm[j+h]=x;
		}
	}
	sortedTrend = SortedKey<real>(N,Perm,Key);
	assert(((sortedTrend == 0) || (sortedTrend == 1)));

}

/*	RealPermShellSortUp(Perm, methods)	rearranges Perm[0..N-1],
 *						where 'N' is the number
 *						of entries in 'methods',
 *						so methods[Perm[0..N-1]].meanRegret
 *						is in increasing order
 *	Perm:					an array to permute
 *	methods:				a set of voting methods
 */
void RealPermShellSortUp(int Perm[], const std::array<oneVotingMethod, NumMethods>& methods)
{
	int h,k;
	int64_t j;
	int x;
	int sortedTrend;
	const range& theRange = methods.size();

	for(k=0; (h=ShellIncs[k])>0; k++) {
		for(const size_t& i : theRange) {
			x=Perm[i];
			for(j=i-h; j>=0 && methods[Perm[j]].meanRegret>methods[x].meanRegret; j -= h) {
				Perm[j+h]=Perm[j];
			}
			Perm[j+h]=x;
		}
	}
	sortedTrend = SortedKey(Perm, methods);
	assert(((sortedTrend == 0) || (sortedTrend == 1)));

}

/*************************** VOTING METHODS:  ********
all are subroutines with a common format - here is the format (which is all subsumed in
the convenient data structure "edata"):
input:
uint NumVoters;
uint64_t NumCands;

uint TopDownPrefs[NumVoters*NumCands];
Entry x*NumCands+y says the candidate who is the (y-1)th choice of voter x, x=0..NumVoters-1.
Candidates are numbered 0..NumCands-1.

uint CandRankings[NumVoters*NumCands];
Entry x*NumCands+y says the ranking (0=top) of the yth candidate (y=0..NumCands-1)
according to voter x, x=0..NumVoters-1.

real Score[NumVoters*NumCands];
Entry x*NumCands+y says the score (a floating point real from 0 to 1)
of the yth candidate (y=0..NumCands-1) according to voter x, x=0..NumVoters-1.

bool Approve[NumVoters*NumCands];
Entry x*NumCands+y says the Approval (a boolean)
of the yth candidate (y=0..NumCands-1) according to voter x, x=0..NumVoters-1.

bool Approve2[NumVoters*NumCands];
Entry x*NumCands+y says the second-level Approval (a boolean)
of the yth candidate (y=0..NumCands-1) according to voter x, x=0..NumVoters-1.
Used in MCA.  A higher level of approval.

real PerceivedUtility[NumVoters*NumCands];
Entry x*NumCands+y says the utility (a floating point real)
of the yth candidate (y=0..NumCands-1) according to voter x, x=0..NumVoters-1.
Note, this is "honest", i.e. undistorted by strategy,
although it IS distorted by ignorance.

real Utility[MaxNumCands*MaxNumVoters];
Entry x*NumCands+y says the utility (a floating point real)
of the yth candidate (y=0..NumCands-1) according to voter x, x=0..NumVoters-1.
Note, this is "honest", i.e. undistorted by strategy, AND undistorted by ignorance.

semi-input:  (used as input by some election methods, but computed from the preceding
genuine input by BuildDefeatsMatrix)
int DefeatsMatrix[MaxNumCands*MaxNumCands];
DefeatsMatrix[i*NumCands+j] is #voters who ranked i above j  (always nonnegative)

real Armytage[MaxNumCands*MaxNumCands];
Armytage[i*NumCands+j] is sum of rating[i]-rating[j] for voters who rated i above j  (always nonnegative)

int MarginsMatrix[MaxNumCands*MaxNumCands];
MarginsMatrix[i*NumCands+j] is margin of i over j  (either sign)

real MargArmy[MaxNumCands*MaxNumCands];
MargArmy[i*NumCands+j] is Armytage[i*NumCands+j]-Armytage[j*NumCands+i]


output (the returned value):
int Winner;       Number in {0,1,...,NumCands-1}  saying who won.
                  (Negative number indicates error condition.)
                  Note: tie breaking is random.  A global array RandPerm[0..NumCands-1] is
                  assumed available to help with that task.

CORE METHODS: Plurality, Borda, IRV, SociallyBest, SociallyWorst, RandomWinner, Approval, Range,
  SmithSet, SchwartzSet, Antiplurality
are "core" voting methods which must always be run (and first), otherwise
some other methods may not work.

side-effects:  can alter these global variables (list follows):
***************************/

 int PlurWinner;
 int AntiPlurWinner;
 int PSecond;
 int ApprovalWinner;
 int IRVwinner;
 int SmithWinner;
 int SchwartzWinner;
 int RangeWinner;
 int BordaWinner;
 int WorstWinner;
 int BestWinner;
 int RandomUncoveredMemb;
uint RangeGranul;
 int IRVTopLim = BIGINT;
 int CondorcetWinner;   /*negative if does not exist*/
 int TrueCW;   /*cond winner based on undistorted true utilities; negative if does not exist*/
 int CopeWinOnlyWinner;
 int SmithIRVwinner;
 int FavListNext[MaxNumVoters];
 int HeadFav[MaxNumCands];
bool CoverMatrix[MaxNumCands*MaxNumCands];

void InitCoreElState(){ /*can use these flags to tell if Plurality() etc have been run*/
	PlurWinner = -1;
	BordaWinner = -1;
	ApprovalWinner = -1;
	RangeWinner = -1;
	IRVwinner = -1;
	BestWinner = -1;
	WorstWinner = -1;
	AntiPlurWinner = -1;
	PSecond = -1;
	SmithWinner = -1;
	SchwartzWinner = -1;
	CopeWinOnlyWinner = -1;
	SmithIRVwinner = -1;
	RandomUncoveredMemb = -1;
	IRVTopLim = BIGINT;
}

/*	oneVoter:	information about a particular Voter
 */
struct oneVoter
{
	uint topDownPrefs[MaxNumCands];
	oneCandidateToTheVoter Candidates[MaxNumCands];
	int64_t favoriteCandidate;
};

template< class T >
		int Minimum(uint64_t N, const CandidateSlate& allCandidates, T oneCandidate::*member, const bool& permute = true, const bool& checkElimination = false);
template< class T >
		int Maximum(uint64_t N, const CandidateSlate& allCandidates, T oneCandidate::*member, const bool& permute = true, const bool& checkElimination = false);

typedef struct dum1 {
	uint NumVoters;
	uint64_t NumCands;
	oneVoter Voters[MaxNumVoters];
	CandidateSlate Candidates;
} edata;

int calculateForRunoff(const edata& E, int first, int second);

/*	PrintEdata(F, E):	prints election data to a file
 *	F:	a pointer to a FILE structure representing the file to which the data is
 *		to be written
 *	E:	the election data to write
 */
void PrintEdata(FILE *F, const edata& E)
{ /* prints out the edata */
	int v;
	uint j;
	const oneVoter (&allVoters)[MaxNumVoters] = E.Voters;
	const uint &numberOfVoters = E.NumVoters;
	const uint64_t& numberOfCandidates = E.NumCands;
	fprintf(F, "numberOfVoters=%d  numberOfCandidates=%lld\n", numberOfVoters, numberOfCandidates);
	for(v=0; v<(int)numberOfVoters; v++){
		const oneVoter& theVoter = allVoters[v];
		const oneCandidateToTheVoter (&allCandidates)[MaxNumCands] = theVoter.Candidates;
		fprintf(F, "Voter %2d:\n",v);
		fprintf(F, "Utility: ");
		for(j=0; j < numberOfCandidates; j++) {
			fprintf(F, "%6.3f", theVoter.Candidates[j].actualUtility);
		}
		fprintf(F, "\n");
		fprintf(F, "PercUti: ");
		for(j=0; j < numberOfCandidates; j++){  fprintf(F, "%6.3f", theVoter.Candidates[j].perceivedUtility);  }
		fprintf(F, "\n");
		fprintf(F, "RangeScore: ");
		for(j=0; j < numberOfCandidates; j++) {
			fprintf(F, "%6.3f", allCandidates[j].score);
		}
		fprintf(F, "\n");
		fprintf(F, "CandRank: ");
		for(j=0; j < numberOfCandidates; j++){  fprintf(F, "%2llu", allCandidates[j].ranking);  }
		fprintf(F, "\n");
		fprintf(F, "TopDown: ");
		for(j=0; j < numberOfCandidates; j++) {
			fprintf(F, "%2d", theVoter.topDownPrefs[j]);
		}
		fprintf(F, "\n");
		fprintf(F, "Approve: ");
		for(j=0; j < numberOfCandidates; j++) {
			fprintf(F, "%2d", allCandidates[j].approve ? 1:0);
		}
		fprintf(F, "\n");
		fprintf(F, "Approve2: ");
		for(j=0; j < numberOfCandidates; j++) {
			fprintf(F, "%2d", allCandidates[j].approve2 ? 1:0);
		}
		fprintf(F, "\n");
	}
	/*???more?*/
}

typedef int EMETH;  /* allows fgrep EMETH IEVS.c to find out what Election methods now available */

EMETH runoffForApprovalVoting(const edata& E);

/*	IBeatsWho(margins, numberOfCandidates):	returns a hash of
 *						Candidates beaten by
 *						the Candidates Whose
 *						margins are given by
 *						'margins'
 *	margins:			array of margins of the given
 *				Candidate with respect to other
 *				Candidates
 *	numberOfCandidates:	the number of Candidates in the
 *				current election
 */
uint64_t IBeatsWho(const MarginsData& margins, const uint64_t& numberOfCandidates)
{
	uint64_t beatsHash = 0;
	for(uint64_t j=0; j<numberOfCandidates; j++) {
		if( margins[j] > 0 ) {
			beatsHash |= (1U<<j);
		}
	}
	return beatsHash;
}

void BuildDefeatsMatrix(edata& E)
{ /* initializes  E->DefeatsMatrix[], E->MarginsMatrix[], RandCandPerm[], NauruWt[], Each Candidate's 'electedCount', Each Candidate's 'drawCount', CondorcetWinner, CopeWinOnlyWinner, TrueCW */
	int k,i,j;
	int64_t y;
	real t;
	bool CondWin, TrueCondWin;
	const oneVoter (&allVoters)[MaxNumVoters] = E.Voters;
	const uint &numberOfVoters = E.NumVoters;
	const uint64_t& numberOfCandidates = E.NumCands;
	CandidateSlate& allTheCandidates = E.Candidates;
	MakeIdentityPerm( numberOfCandidates, RandCandPerm );
	RandomlyPermute( numberOfCandidates, RandCandPerm );

	assert(numberOfCandidates <= MaxNumCands);
	for(oneCandidate& eachCandidate : allTheCandidates) {
		const real utilitySum = eachCandidate.utilitySum;
		const auto pluralityVotes = eachCandidate.pluralityVotes;
		eachCandidate = oneCandidate();
		eachCandidate.utilitySum = utilitySum;
		eachCandidate.pluralityVotes = pluralityVotes;
	}
	for(k=0; k<(int)numberOfVoters; k++) {
		const oneVoter& theVoter = allVoters[k];
		const oneCandidateToTheVoter (&allCandidatesToTheVoter)[MaxNumCands] = theVoter.Candidates;
		for(i=0; i<numberOfCandidates; i++) {
			const oneCandidateToTheVoter &firstCandidateToTheVoter = allCandidatesToTheVoter[i];
			oneCandidate& firstCandidate = allTheCandidates[i];
			for(j=0; j<i; j++) {
				const oneCandidateToTheVoter &secondCandidateToTheVoter = allCandidatesToTheVoter[j];
				oneCandidate& secondCandidate = allTheCandidates[j];
				if( secondCandidateToTheVoter.ranking>firstCandidateToTheVoter.ranking ) {
					firstCandidate.DefeatsMatrix[j]++;	/*i preferred above j*/
				}else{
					secondCandidate.DefeatsMatrix[i]++;	/*j preferred above i*/
				}
				t = firstCandidateToTheVoter.score - secondCandidateToTheVoter.score;
				if(t > 0.0) {
					firstCandidate.ArmytageMatrix[j] += t;
					firstCandidate.ArmytageDefeatsMatrix[j]++;
				}else{
					secondCandidate.ArmytageMatrix[i] -= t;
					secondCandidate.ArmytageDefeatsMatrix[i]++;
				}
				if(firstCandidateToTheVoter.actualUtility > secondCandidateToTheVoter.actualUtility) {
					firstCandidate.TrueDefeatsMatrix[j]++;
				}else{
					secondCandidate.TrueDefeatsMatrix[i]++;
				}
			}
		}
	}
	CondorcetWinner = -1;
	TrueCW = -1;
	for(i=0; i<numberOfCandidates; i++) {
		oneCandidate& firstCandidate = allTheCandidates[i];
		firstCandidate.electedCount = 0;
		firstCandidate.drawCount = 0;
		CondWin = true;
		TrueCondWin = true;
		for(j=0; j<numberOfCandidates; j++) {
			const int &iDefeatsJ = firstCandidate.DefeatsMatrix[j];
			oneCandidate& secondCandidate = allTheCandidates[j];
			const int &jDefeatsI = secondCandidate.DefeatsMatrix[i];
			int64_t& marginOfTheFirstCandidateComparedToTheSecond = firstCandidate.margins[j];
			assert( iDefeatsJ <= (int)numberOfVoters );
			assert( iDefeatsJ >= 0 );
			assert( iDefeatsJ + jDefeatsI <= (int)numberOfVoters );
			y = firstCandidate.ArmytageDefeatsMatrix[j];
			y -= secondCandidate.ArmytageDefeatsMatrix[i];
			if(y>0) {
				firstCandidate.ArmytageMarginsMatrix[j] = firstCandidate.ArmytageMatrix[j];
			} else {/*y<=0*/
				firstCandidate.ArmytageMarginsMatrix[j] = 0;
			}
			y = iDefeatsJ;
			y -= jDefeatsI;
			marginOfTheFirstCandidateComparedToTheSecond = y;
			assert(i!=j || marginOfTheFirstCandidateComparedToTheSecond == 0);
			if(y > 0) {
				firstCandidate.electedCount++;
			}
			if(y==0) {
				firstCandidate.drawCount++;
			}
			if(y<=0 && j!=i) {
				CondWin = false;
			} /* if beaten or tied, not a CondorcetWinner by this defn */
			y = firstCandidate.TrueDefeatsMatrix[j];
			y -= secondCandidate.TrueDefeatsMatrix[i];
			if(y<=0 && j!=i) { TrueCondWin = false; }
		}
		if( CondWin ) {
			CondorcetWinner = i;  /* i beats all opponents (unique Condorcet winner) */
		}
		if( TrueCondWin ) {
			TrueCW = i;  /* i beats all opponents (unique Condorcet winner) */
		}
	}
	/* find who-you-beat sets: */
	for(i=0; i<numberOfCandidates; i++) {
		oneCandidate& firstCandidate = allTheCandidates[i];
		uint64_t& beat = firstCandidate.Ibeat;
		const MarginsData& marginsOfI = firstCandidate.margins;
		beat = IBeatsWho(marginsOfI, numberOfCandidates);
	}
	CopeWinOnlyWinner = Maximum<int64_t>(numberOfCandidates, allTheCandidates, &oneCandidate::electedCount);
	for(i=0; i<(int)numberOfVoters; i++) {
		const oneCandidateToTheVoter (&allCandidatesToTheVoter)[MaxNumCands] = allVoters[i].Candidates;
		for(j=0; j<numberOfCandidates; j++) {
			const oneCandidateToTheVoter &firstCandidateToTheVoter = allCandidatesToTheVoter[j];
			oneCandidate& firstCandidate = allTheCandidates[j];
			for(k=j+1; k<numberOfCandidates; k++) {
				const oneCandidateToTheVoter &secondCandidateToTheVoter = allCandidatesToTheVoter[k];
				oneCandidate& secondCandidate = allTheCandidates[k];
				y = ((firstCandidateToTheVoter.approve && secondCandidateToTheVoter.approve) ? 1:0);
				firstCandidate.alsoApprovedWith[k] += y; /* count of voters who approve of both j and k */
				secondCandidate.alsoApprovedWith[j] += y;
			}
		}
	}
}

template <class T>
	void Zero(const uint64_t& number, CandidateSlate& allCandidates, T oneCandidate::*member);

EMETH SociallyBest(edata& E  /* greatest utility-sum winner */)
{ /* side effects: Each Candidate's 'utilitySum', BestWinner */
	int i,j;
	const uint64_t& numberOfCandidates = E.NumCands;
	const uint& numberOfVoters = E.NumVoters;
	const oneVoter (&allVoters)[MaxNumVoters] = E.Voters;
	CandidateSlate& allCandidates = E.Candidates;
	Zero(numberOfCandidates, allCandidates, &oneCandidate::utilitySum);
	for(i=0; i<(int)numberOfVoters; i++) {
		const oneVoter& theVoter = allVoters[i];
		const oneCandidateToTheVoter (&allCandidatesToTheVoter)[MaxNumCands] = theVoter.Candidates;
		for(j=0; j<numberOfCandidates; j++) {
			allCandidates[j].utilitySum += allCandidatesToTheVoter[j].actualUtility;
		}
	}
	BestWinner = Maximum(numberOfCandidates, allCandidates, &oneCandidate::utilitySum);
	for(j=0; j<numberOfCandidates; j++) {
		assert( allCandidates[BestWinner].utilitySum>= allCandidates[j].utilitySum );
	}
	return BestWinner;
}

/*	Determine(Winner, theMethod, E):	determines the
 *						Winner of an
 *						election
 *						conducted using
 *						a given voting
 *						method
 *	Winner:		the Winner to determine; receives the
 *			results, if needed
 *	theMethod:	the voting method to determine the
 *			Winner
 *	E:		the election data to use to determine
 *			the Winner
 */
template<class T> void Determine(int& Winner, T theMethod, edata& E)
{
	if(0 > Winner) {
		theMethod(E);
	}
}

/*	SociallyWorst(E):	returns the Candidate with the lowest utility sum
 *	E:	the election data used to determine the socially worst Candidate
 */
EMETH SociallyWorst(edata& E   /* least utility-sum winner */)
{ /* side effects: Each Candidate's 'utilitySum', WorstWinner */
	Determine(BestWinner, SociallyBest, E);
	WorstWinner = Minimum(E.NumCands, E.Candidates, &oneCandidate::utilitySum);
	return WorstWinner;
}

EMETH RandomWinner(const edata& E){ return (int)RandInt( E.NumCands ); }

/*	RandomBallot(E):	returns the honest top Candidate of a randomly selected
 *				Voter; such a technique is considered 'strategy-proof'
 *	E:	the election data used to determine thw Winner
 */
EMETH RandomBallot(const edata& E)
{ /*honest top choice of a random voter is elected. Strategyproof.*/
	int winner;
	const uint64_t& numberOfCandidates = E.NumCands;
	const uint64_t VoterIndex = RandInt(E.NumVoters);
	winner = ArgMaxArr<real>(numberOfCandidates, E.Voters[VoterIndex].Candidates, (int*)RandCandPerm);
	return winner;
}

/*	minimumActualUtility(Candidates, numberOfCandidates):	returns the minimum
 *								actual utility of All
 *								Candidates from the
 *								perspective of a
 *								particular Voter
 *	Candidates:		an array of Candidates from the perspective
 *				of a particular Voter
 *	numberOfCandidates:	the number of Candidates in the current
 *				election
 */
real minimumActualUtility(const oneCandidateToTheVoter (&Candidates)[MaxNumCands], const uint64_t& numberOfCandidates)
{
	real lowest = HUGE;

	for(uint64_t j=0; j<numberOfCandidates; j++) {
		real utility = Candidates[j].actualUtility;
		if(lowest > utility) {
			lowest = utility;
		}
	}
	return lowest;
}

EMETH Hay(const edata& E /*Strategyproof. Prob of election proportional to sum of sqrts of [normalized] utilities*/)
{
	real UtilityRootSum[MaxNumCands];
	int i,j,winner;
	real minu, sumrts, t;
	const uint64_t& numberOfCandidates = E.NumCands;
	const uint& numberOfVoters = E.NumVoters;
	const oneVoter (&allVoters)[MaxNumVoters] = E.Voters;
	ZeroArray( numberOfCandidates, UtilityRootSum );
	for(i=0; i<(int)numberOfVoters; i++) {
		const oneVoter& theVoter = allVoters[i];
		const oneCandidateToTheVoter (&allCandidatesToTheVoter)[MaxNumCands] = theVoter.Candidates;
		minu = minimumActualUtility(allCandidatesToTheVoter, numberOfCandidates);
		sumrts = 0.0;
		for(j=0; j<numberOfCandidates; j++) {
			sumrts += sqrt(allCandidatesToTheVoter[j].actualUtility - minu);
		}
		if(sumrts>0.0) {
			for(j=0; j<numberOfCandidates; j++) {
				UtilityRootSum[j] += sqrt(allCandidatesToTheVoter[j].actualUtility - minu)/sumrts;
			}
		}
	}
	sumrts = 0.0;
	for(j=0; j<numberOfCandidates; j++) {
		sumrts += UtilityRootSum[j];
	}
	t = Rand01() * sumrts;
	sumrts = 0.0;
	winner = -1;
	for(j=(int)numberOfCandidates-1; j>=0; j--) {
		sumrts += UtilityRootSum[j];
		if( t<sumrts ) { winner=j; break; }
	}
	assert(winner >= 0);
	return winner;
}

EMETH Plurality(edata& E   /* canddt with most top-rank votes wins */)
{ /* side effects: Each Candidate's plurality vote count, PlurWinner */
	int i;
	const uint64_t& numberOfCandidates = E.NumCands;
	const uint& numberOfVoters = E.NumVoters;
	const oneVoter (&allVoters)[MaxNumVoters] = E.Voters;
	CandidateSlate& allCandidates = E.Candidates;
	Zero(numberOfCandidates, allCandidates, &oneCandidate::pluralityVotes);
	for(i=0; i<(int)numberOfVoters; i++) {
		allCandidates[ allVoters[i].topDownPrefs[0] ].pluralityVotes++;
	}
	PlurWinner = Maximum(numberOfCandidates, allCandidates, &oneCandidate::pluralityVotes);
	return PlurWinner;
}

EMETH AntiPlurality(edata& E   /* canddt with fewest bottom-rank votes wins */)
{ /* side effects: each Candidates anti-plurality vote counts, AntiPlurWinner */
	int i;
	const uint64_t& numberOfCandidates = E.NumCands;
	const uint64_t& lastCandidateIndex = numberOfCandidates - 1;
	const uint& numberOfVoters = E.NumVoters;
	const oneVoter (&allVoters)[MaxNumVoters] = E.Voters;
	CandidateSlate& allCandidates = E.Candidates;
	Zero(numberOfCandidates, allCandidates, &oneCandidate::antiPluralityVotes);
	for(i=0; i<numberOfVoters; i++) {
		allCandidates[allVoters[i].topDownPrefs[lastCandidateIndex]].antiPluralityVotes++;
	}
	AntiPlurWinner = Minimum<uint64_t>(numberOfCandidates, allCandidates, &oneCandidate::antiPluralityVotes);
	return AntiPlurWinner;
}

/* Plurality needs to have already been run before running Dabagh.  */
EMETH Dabagh(const edata& E   /* canddt with greatest Dabagh = 2*#top-rank-votes + 1*#second-rank-votes, wins */)
{
	uint64_t DabaghVoteCount[MaxNumCands];
	int i, winner;
	const uint64_t& numberOfCandidates = E.NumCands;
	const uint& numberOfVoters = E.NumVoters;
	const oneVoter (&allVoters)[MaxNumVoters] = E.Voters;
	const CandidateSlate& allCandidates = E.Candidates;
	CopyArray(numberOfCandidates, allCandidates, DabaghVoteCount, &oneCandidate::pluralityVotes);
	ScaleVec( numberOfCandidates, DabaghVoteCount, (uint64_t)2 );
	for(i=0; i<(int)numberOfVoters; i++) { /*add 2nd-pref votes with weight=1:*/
		DabaghVoteCount[ allVoters[i].topDownPrefs[1] ]++;
	}
	winner = ArgMaxArr(numberOfCandidates, DabaghVoteCount, (int*)RandCandPerm);
	return winner;
}

/* Plurality needs to have already been run before running VtForAgainst.  */
EMETH VtForAgainst(const edata& E   /* canddt with greatest score = #votesFor - #votesAgainst,  wins */)
{
	int64_t VFAVoteCount[MaxNumCands];
	int i, winner;
	uint64_t last;
	const uint64_t& numberOfCandidates = E.NumCands;
	const uint& numberOfVoters = E.NumVoters;
	const oneVoter (&allVoters)[MaxNumVoters] = E.Voters;
	const CandidateSlate& allCandidates = E.Candidates;
	last = numberOfCandidates - 1;
	CopyArray(numberOfCandidates, allCandidates, (uint64_t*)VFAVoteCount, &oneCandidate::pluralityVotes);
	for(i=0; i<(int)numberOfVoters; i++) {
		VFAVoteCount[ allVoters[i].topDownPrefs[last] ]--;
	}
	winner = ArgMaxArr(numberOfCandidates, VFAVoteCount, (int*)RandCandPerm);
	return winner;
}

/*	Top2Runoff(E):	returns the index of the Winner of a runoff between the top 2
 *			Candidates; this implementation only compares perceived
 *			utilities; if the Candidates are perceived to have the same
 *			utility to the same Voter, this procedure acts as if said Voter
 *			did not vote; this procedure also presumes honest voting, which
 *			makes sense for this scenario; in the event of a tie between the
 *			top 2 Candidates, the Winner is decided with a 'coin flip'
 *	E:	the election data used to determine the Winner
 */
EMETH Top2Runoff(edata& E    /* Top2Runoff=top-2-runoff, 2nd round has fully-honest voting */)
{ /* side effects: PSecond */
	PSecond = -1;
	Determine(PlurWinner, Plurality, E);
	PSecond = SecondMaximum(E.NumCands, E.Candidates, &oneCandidate::pluralityVotes, PlurWinner);
	assert(PSecond>=0);
	return calculateForRunoff(E, PlurWinner, PSecond);
}

/*	VenzkeDisqPlur(E):	returns the Venzke disqualified plurality Winner; with
 *				this method, the Candidate with the most for votes is
 *				elected unless that Candidate receives more than half of
 *				the against votes, in which case said Candidate is not
 *				eligible to be elected, in which case the 2nd place
 *				Winner of the plurality vote wins
 *	E:	the election data used to determine the Winner
 */
EMETH VenzkeDisqPlur(edata& E   /* Plurality winner wins unless over 50% vote him bottom; then plurality 2nd-winner wins. */)
{ /* side effects: VFAVoteCount[] */
	Determine(PlurWinner, Plurality, E);
	Determine(AntiPlurWinner, AntiPlurality, E);
	if( (2*E.Candidates[PlurWinner].antiPluralityVotes) > E.NumVoters ) {
		Determine(PSecond, Top2Runoff, E);
		return PSecond;
	}
	return PlurWinner;
}

EMETH PlurIR(edata& E    /* PlurIR=plur+immediate runoff (give ranking as vote) */)
{
	int64_t i;
	Determine(CopeWinOnlyWinner, BuildDefeatsMatrix, E);
	Determine(PSecond, Top2Runoff, E);
	Determine(PlurWinner, Plurality, E);
	Determine(PSecond, Top2Runoff, E);
	assert(PSecond>=0);
	i = E.Candidates[PlurWinner].margins[PSecond];
	if(i>0) {
		return(PlurWinner);
	}
	if(i<0) {
		return(PSecond);
	}else if(RandBool()) {
		return(PSecond);
	} else {
		return(PlurWinner);
	}
}

EMETH Borda(edata& E  /* Borda: weighted positional with weights N-1, N-2, ..., 0  if N-canddts */)
{ /* side effects: Each Candidate's 'BordaVotes', BordaWinner */
	int i,j,t;
	const uint64_t& numberOfCandidates = E.NumCands;
	CandidateSlate& allCandidates = E.Candidates;
	Determine(CopeWinOnlyWinner, BuildDefeatsMatrix, E);
	for(i=0; i<numberOfCandidates; i++) {
		const MarginsData& marginsOfI = allCandidates[i].margins;
		t=0;
		for(j=0; j<numberOfCandidates; j++) {
			t += marginsOfI[j];
		}
		allCandidates[i].BordaVotes = t;
	}
	BordaWinner = Maximum(numberOfCandidates, allCandidates, &oneCandidate::BordaVotes);
	return BordaWinner;
}

/*	Black(E):	returns the Condorcet Winner if One exists and the Borda Winner
 *			otherwise
 *	E:	the election data used to determine the Winner
 */
EMETH Black(edata& E  /* Condorcet winner if exists, else use Borda */)
{
	Determine(CopeWinOnlyWinner, BuildDefeatsMatrix, E);
	if(CondorcetWinner >= 0) {
		return CondorcetWinner;
	}
	Determine(BordaWinner, Borda, E);
	return BordaWinner;
}

/*	RandomPair(E):	returns a pairwise honest-utility Winner from 2 randomly
 *			selected Candidates; if the 2 Candidates have the same
 *			honest-utility, the Winner is determined by a 'coin flip'
 *	E:	the election data used to determine the Winner
 */
EMETH RandomPair(const edata& E)
{ /*pairwise honest-util winner among 2 random candidates is elected*/
	int x,y;
	const uint64_t& numberOfCandidates = E.NumCands;
	const CandidateSlate& allCandidates = E.Candidates;
	x = (int)RandInt(numberOfCandidates);
	y = (int)RandInt(numberOfCandidates);
	if( allCandidates[x].utilitySum > allCandidates[y].utilitySum ) {
		return x;
	}
	if( allCandidates[x].utilitySum < allCandidates[y].utilitySum ) {
		return y;
	}
	return flipACoin(y, x);
}

/*	findLoser(count, Candidates, weights):	returns the index of
 *						the Candidate which is
 *						both non-eliminated
 *						and has the lowest
 *						corresponding entry in
 *						'weights'
 *	count:		the number of Candidates in the current
 *			election
 *	Candidates:	the set of Candidates in the current election
 *	weights:	the set of weights used to find the 'Loser'
 */
int findLoser(const uint64_t& count, const CandidateSlate& Candidates, const int64_t (&weights)[MaxNumCands])
{
	int Loser = -1;
	int64_t minc = BIGINT;
	int r;
	for(int i=(int)count-1; i>=0; i--) {
		r = RandCandPerm[i];
		if(!Candidates[r].eliminated && (weights[r]<minc)) {
			minc=weights[r];
			Loser=r;
		}
	}
	return Loser;
}

/*	subtract(voteCount, theCandidate, BordaLoser):	subtracts from 'voteCount'
 *							the margin of 'theCandidate'
 *							with respect to the 'BordaLoser'
 *							iff 'theCandidate' has not
 *							yet been eliminated from
 *							contention
 *	voteCount:	the value to adjust
 *	theCandidate:	the Candidate to consider
 *	BordaLoser:	the Borda count Loser
 */
void subtract(int64_t& voteCount, const oneCandidate& theCandidate, const int& BordaLoser)
{
	if(!theCandidate.eliminated) {
		voteCount -= theCandidate.margins[BordaLoser];
	}
}

/*	NansonBaldwin(E):	returns the index of the Baldwin Winner, not necessarily
 *				the Nanson Winner, and, even then, it looks as if the
 *				Baldwin Winner may be calculated incorrectly by this
 *				function; this function should be reviewed at some point
 *				for algorithmic correctness; also, a separate procedure
 *				should be created to determine the Nanson Winner; the
 *				procedure returns -1 if an error occurs
 *	E:	the election data used to determine the Baldwin Winner
 */
EMETH NansonBaldwin(edata& E  /* repeatedly eliminate Borda loser */)
{ /* side effects: Each Candidte's 'eliminated' member */
	int64_t NansonVoteCount[MaxNumCands];
	int i, BordaLoser, rnd;
	const uint64_t& numberOfCandidates = E.NumCands;
	CandidateSlate& allCandidates = E.Candidates;
#if defined(CWSPEEDUP) && CWSPEEDUP
	if(CondorcetWinner >= 0) return CondorcetWinner;
#endif
	Determine(BordaWinner, Borda, E);
	Zero(numberOfCandidates, allCandidates, &oneCandidate::eliminated);
	CopyArray(numberOfCandidates, allCandidates, NansonVoteCount, &oneCandidate::BordaVotes);
	RandomlyPermute( numberOfCandidates, RandCandPerm );
	for(rnd=1; rnd < (int)numberOfCandidates; rnd++) {
		BordaLoser = findLoser(numberOfCandidates, allCandidates, NansonVoteCount);
		assert(BordaLoser>=0);
		ensure(BordaLoser>=0, 7);
		allCandidates[BordaLoser].eliminated = true;
		for(i=0; i<numberOfCandidates; i++) {
			subtract(NansonVoteCount[i], allCandidates[i], BordaLoser);
		}
	} /* end of for(rnd) */
	for(i=0; i<numberOfCandidates; i++) { /* find non-eliminated candidate... */
		if(!allCandidates[i].eliminated) {
			return i; /*NansonBaldwin winner*/
		}
	}
	return(-1); /*error*/
}

/*	Rouse(E):	returns the Rouse Winner; Rouse is like Nanson-Baldwin but with
 *			an extra level of recursion: it successively pseudo-eliminates
 *			the candidate with the highest Borda score until one is left,
 *			then it genuinely-eliminates that one from the original list;
 *			this step is repeated until a single candidate is left; this
 *			implementation is both slow, O(N^4) steps to find winner for N
 *			candidates from pairwise matrix, and according to commentary,
 *			"buggy"
 *	E:	the election data to used to determine the Winner
 */
EMETH Rouse(edata& E  /*like Nanson-Baldwin but with an extra level of recursion BUGGY*/)
{
	bool Rmark[MaxNumCands];
	bool rRmark[MaxNumCands];
	int i,j,k,m,r,highestb,bordsum, maxb, winner;
	const uint64_t& numberOfCandidates = E.NumCands;
	const CandidateSlate& allCandidates = E.Candidates;
	Determine(CopeWinOnlyWinner, BuildDefeatsMatrix, E);
	FillArray(numberOfCandidates, rRmark, true); /* nobody eliminated initially */
	for(m= (int)numberOfCandidates; m>1; m--) { /* NumCands-1 elimination-rounds */
		for(k=1; k<m; k++) {  /* m-1 pseudo-elimination rounds */
			maxb = -BIGINT; highestb = -1;
			FillArray(numberOfCandidates, Rmark, true); /* nobody pseudo-eliminated initially */
			RandomlyPermute( numberOfCandidates, RandCandPerm );
			for(i=(int)numberOfCandidates-1; i>=0; i--) {
				r = RandCandPerm[i];
				if(rRmark[r] && Rmark[r]) {
					const MarginsData& marginsOfR = allCandidates[r].margins;
					bordsum = 0;
					for(j=0; j<numberOfCandidates; j++) {
						if(Rmark[j] && rRmark[j]) {
							bordsum += marginsOfR[j];
						}
					}
					if(maxb < bordsum) { maxb = bordsum;  highestb = r; }
				}
			}
			assert(highestb >= 0);
			ensure(highestb >= 0, 8);
			assert(rRmark[highestb]);
			assert(Rmark[highestb]);
			Rmark[highestb] = false; /* pseudo-eliminate borda-winner */
		}
		/* Find the non-psu-eliminated canddt i: */
		for(i=(int)numberOfCandidates-1; i>=0; i--) { if(Rmark[i] && rRmark[i]) { break; }}
		assert(i>=0);
		assert(i<(int)numberOfCandidates);
		ensure(i>=0, 9);
		rRmark[i] = false; /* (genuinely) eliminate it */
	}
	winner = -1;
	for(k=0; k<numberOfCandidates; k++) {
		if(rRmark[k]) {
			winner = k;
			break;
		}
	}
	assert((winner==0) || !rRmark[0]);
	/***???
	if(winner != CondorcetWinner && CondorcetWinner>=0) { *???*
		printf("Rouse aint Condorcet (RW=%d CW=%d):\n", winner, CondorcetWinner);
		PrintEdata(stdout, E);
	}***/
	return winner;
}

/*	IterCopeland(E):	iterates Copeland on tied-winner set from previous
 *				iteration until either no ties remain or no further
 *				reduction of ties occurs; either the iterated Copeland
 *				winner is returned or -1 is returned to indicate an
 *				error
 *	E:	the election data to use for iteration
 */
EMETH IterCopeland( const edata& E  /*iterate Copeland on tied-winner set from previous iter until fixpt*/)
{ /*Idea of this is due to Alex Small. */
	int summ[MaxNumCands];
	int i, r, j, z, maxc, winner, tiect, oldtiect;
	int mxs;
	const uint64_t& numberOfCandidates = E.NumCands;
	const CandidateSlate& allCandidates = E.Candidates;
	assert(numberOfCandidates >= 2);
#if defined(CWSPEEDUP) && CWSPEEDUP
	if(CondorcetWinner >= 0) return CondorcetWinner;
#endif
	winner = -1; oldtiect = -1;
	ZeroArray(numberOfCandidates, summ);
	RandomlyPermute( numberOfCandidates, RandCandPerm );
	mxs = 0;
	for(z=0; z<numberOfCandidates; z++) {
		tiect = 0;
		for(i=0;i<numberOfCandidates; i++) {
			if(summ[i] < mxs) {
				summ[i] = -BIGINT;
			}
		}
		maxc = -BIGINT;
		for(i=(int)numberOfCandidates-1; i>=0; i--) {
			r = RandCandPerm[i];
			if(summ[r] >= mxs) {
				const MarginsData& marginsOfR = allCandidates[r].margins;
				for(j=0; j<numberOfCandidates; j++) {
					if((j!=r) && (summ[j]>=mxs)) {
						summ[r] += 1 + Sign<int64_t>(marginsOfR[j]);
					}
				}
				if(summ[r] >= maxc) {
					if(summ[r]==maxc) { tiect++; }
					else{ maxc=summ[r]; winner=r;  tiect=1; }
				}
			}
		}
		assert(tiect>=1);
		if((tiect<=1) || (tiect==oldtiect)) {
			return(winner);
		}
		oldtiect = tiect;
	}
	return(-1);
}

EMETH Nauru(const edata& E  /* weighted positional with weights 1, 1/2, 1/3, 1/4,... */)
{
	uint NauruVoteCount[MaxNumCands];
	int i;
	int j;
	int winner;
	const uint64_t& numberOfCandidates = E.NumCands;
	uint64_t x = LCMfact[numberOfCandidates];
	uint64_t NauruWt[numberOfCandidates];
	const uint& numberOfVoters = E.NumVoters;
	ZeroArray( numberOfCandidates, (int*)NauruVoteCount );
	for(j=1; j<=numberOfCandidates; j++) {
		NauruWt[j-1] = x / j;
	}
	for(i=0; i<(int)numberOfVoters; i++) {
		const oneVoter& theVoter = E.Voters[i];
		for(j=0; j<numberOfCandidates; j++) {
			const uint64_t& theRank = theVoter.Candidates[j].ranking;
			ensure(theRank < numberOfCandidates, 50);
			NauruVoteCount[j] += NauruWt[theRank];
		}
	}
	winner = ArgMaxUIntArr( numberOfCandidates, NauruVoteCount, (int*)RandCandPerm );
	return(winner);
}

EMETH HeismanTrophy(const edata& E  /* Heisman: weighted positional with weights 3, 2, 1, 0,...,0 */)
{
	uint HeismanVoteCount[MaxNumCands];
	int i,j;
	int winner;
	const uint64_t& numberOfCandidates = E.NumCands;
	const uint& numberOfVoters = E.NumVoters;
	const oneVoter (&allVoters)[MaxNumVoters] = E.Voters;
	ZeroArray( numberOfCandidates, (int*)HeismanVoteCount );
	for(i=0; i<(int)numberOfVoters; i++) {
		const oneVoter& theVoter = allVoters[i];
		for(j=0; j<std::min<typeof(numberOfCandidates)>(numberOfCandidates, 3); j++) {
			HeismanVoteCount[ theVoter.topDownPrefs[j] ] += 3-j;
		}
	}
	winner = ArgMaxUIntArr( numberOfCandidates, HeismanVoteCount, (int*)RandCandPerm );
	return(winner);
}

EMETH BaseballMVP(const edata& E  /* weighted positional with weights 14,9,8,7,6,5,4,3,2,1,0,...,0 */)
{
	uint BaseballVoteCount[MaxNumCands];
	static const uint BaseballWt[] = {14,9,8,7,6,5,4,3,2,1};
	static const uint64_t numberOfBaseballWeightings = sizeof(BaseballWt)/sizeof(BaseballWt[0]);
	int i,j;
	int winner;
	const uint64_t& numberOfCandidates = E.NumCands;
	const uint64_t& maximumNumberOfLoops = std::min<typeof(numberOfCandidates)>( numberOfCandidates, numberOfBaseballWeightings);
	const uint& numberOfVoters = E.NumVoters;
	const oneVoter (&allVoters)[MaxNumVoters] = E.Voters;
	ZeroArray( numberOfCandidates, (int*)BaseballVoteCount );
	for(i=0; i<(int)numberOfVoters; i++) {
		const oneVoter& theVoter = allVoters[i];
		for(j=0; j<maximumNumberOfLoops; j++) {
			BaseballVoteCount[ theVoter.topDownPrefs[j] ] += BaseballWt[j];
		}
	}
	winner = ArgMaxUIntArr( numberOfCandidates, BaseballVoteCount, (int*)RandCandPerm );
	return(winner);
}

EMETH CondorcetLR(edata& E   /* candidate with least sum-of-pairwise-defeat-margins wins */)
{
	uint SumOfDefeatMargins[MaxNumCands]={0};
	int i,j,t,winner;
	const uint64_t& numberOfCandidates = E.NumCands;
	const CandidateSlate& allCandidates = E.Candidates;
	Determine(CopeWinOnlyWinner, BuildDefeatsMatrix, E);
#if defined(CWSPEEDUP) && CWSPEEDUP
	if(CondorcetWinner >= 0) return CondorcetWinner;
#endif
	for(i=0; i<numberOfCandidates; i++) {
		t = 0;
		for(j=0; j<numberOfCandidates; j++) {
			t += PosInt( allCandidates[j].margins[i] );
		}
		SumOfDefeatMargins[i] = t;
	}
	winner = ArgMinArr<int>(numberOfCandidates, (int*)SumOfDefeatMargins, (int*)RandCandPerm);
	return winner;
}

/*	Sinkhorn(E):	returns the index of the Candidate with the maximum Sinkhorn
 *			rating from an all-positive DefeatsMatrix+1; the Sinkhorn rating
 *			is described in "Sinkhorn ratings, and new strongly polynomial
 *			time algorithms for Sinkhorn balancing, Perron eigenvectors, and
 *			Markov chains" by Warren D. Smith, 21 Shore Oaks Drive, Stony
 *			Brook, New York 11790
 *	E:	the election data used to determine the Winner
 */
EMETH Sinkhorn(edata& E  /* candidate with max Sinkhorn rating (from all-positive DefeatsMatrix+1) */)
{
	real SinkRat[MaxNumCands]={0};
	real SinkRow[MaxNumCands], SinkCol[MaxNumCands];
	real SinkMat[MaxNumCands*MaxNumCands]={0};
	int j,k,winner;
	const uint64_t& numberOfCandidates = E.NumCands;
	const CandidateSlate& allCandidates = E.Candidates;
	real t,maxsum,minsum,sum,maxminRatio;
	Determine(CopeWinOnlyWinner, BuildDefeatsMatrix, E);
	FillArray( numberOfCandidates, SinkRow, 1.0 );
	FillArray( numberOfCandidates, SinkCol, 1.0 );
	do{
		for(k=0; k < (int)numberOfCandidates; k++) {
			const DefeatsData& defeats = allCandidates[k].DefeatsMatrix;
			for(j=0; j<numberOfCandidates; j++) {
				SinkMat[k*numberOfCandidates + j] =
					SinkRow[k]*SinkCol[j] * (defeats[j] + 1.0);
			}
		}
		maxsum = -HUGE; minsum = HUGE;
		for(k=0; k < (int)numberOfCandidates; k++) {
			sum = 0.0;
			for(j=(int)numberOfCandidates-1; j>=0; j--) {
				sum += SinkMat[k*numberOfCandidates + j];
			}
			minsum = std::min(minsum, sum);
			maxsum = std::max(maxsum, sum);
			ensure((sum!=0.0), 26);
			SinkRow[k] /= sum;
		}
		ensure((minsum != 0.0), 1);
		maxminRatio = maxsum/minsum;
		maxsum = -HUGE; minsum = HUGE;
		for(k=0; k < (int)numberOfCandidates; k++) {
			sum = 0.0;
			for(j=(int)numberOfCandidates-1; j>=0; j--) {
				sum += SinkMat[j*numberOfCandidates + k];
			}
			minsum = std::min(minsum, sum);
			maxsum = std::max(maxsum, sum);
			ensure((sum!=0.0), 27);
			SinkCol[k] /= sum;
		}
		ensure((minsum != 0.0), 2);
		t = maxsum/minsum;
		maxminRatio = std::max(maxminRatio, t);
	}until(maxminRatio < 1.000003);
	for(j=0; j<numberOfCandidates; j++) {
		SinkRat[j] = SinkCol[j]/SinkRow[j];
	}
	winner = ArgMaxArr<real>(numberOfCandidates, SinkRat, (int*)RandCandPerm);
	return winner;
}

EMETH KeenerEig(edata& E  /* winning canddt has max Frobenius eigenvector entry (of all-positive DefeatsMatrix+1) */)
{
	real EigVec[MaxNumCands], EigVec2[MaxNumCands];
	int j,k,winner;
	real t,sum,dist;
	const uint64_t& numberOfCandidates = E.NumCands;
	const CandidateSlate& allCandidates = E.Candidates;
	Determine(CopeWinOnlyWinner, BuildDefeatsMatrix, E);
	FillArray( numberOfCandidates, EigVec, 1.0 );
	FillArray( numberOfCandidates, EigVec2, 1.0 );
	do{
		for(k=0; k < (int)numberOfCandidates; k++) {
			const DefeatsData& defeats = allCandidates[k].DefeatsMatrix;
			t = 0;
			for(j=(int)numberOfCandidates-1; j>=0; j--) {
				t += (defeats[j] + 1.0) * EigVec[j];
			}
			EigVec2[k] = t;
		}
		sum = SumRealArray( numberOfCandidates, EigVec2 );
		ScaleVec( numberOfCandidates, EigVec2, 1.0/sum );
		dist = L1Distance(numberOfCandidates, EigVec, EigVec2);
		CopyArray( numberOfCandidates, EigVec2, EigVec );
	}until( dist < 0.00001 );
	winner = ArgMaxArr<real>(numberOfCandidates, EigVec, (int*)RandCandPerm);
	return winner;
}

EMETH SimpsonKramer(edata& E  /* candidate with mildest worst-defeat wins */)
{
	int64_t WorstDefeatMargin[MaxNumCands]={0};
	int i,r;
	int64_t x;
	int64_t t;
	int j,winner;
	const uint64_t& numberOfCandidates = E.NumCands;
	const CandidateSlate& allCandidates = E.Candidates;
	Determine(CopeWinOnlyWinner, BuildDefeatsMatrix, E);
#if defined(CWSPEEDUP) && CWSPEEDUP
	if(CondorcetWinner >= 0) return CondorcetWinner;
#endif
	for(i=(int)numberOfCandidates-1; i>=0; i--) {
		t = 0;
		RandomlyPermute( numberOfCandidates, RandCandPerm );
		for(j=(int)numberOfCandidates-1; j>=0; j--) {
			r = RandCandPerm[j];
			x = allCandidates[r].margins[i];
			t = std::max(t, x);
		}
		WorstDefeatMargin[i] = t;
	}
	winner = ArgMinArr<int64_t>(numberOfCandidates, WorstDefeatMargin, (int*)RandCandPerm);
	return winner;
}

/*	RaynaudElim(E):	repeatedly eliminates the Candidate suffering the
 *			worst-margin-defeat until only one Candidate remains; the index
 *			of the remaining Candidate is returned or -1 if an error occurs
 *	E:	the election data used to determine the Raynaud Winner
 */
EMETH RaynaudElim(edata& E  /* repeatedly eliminate canddt who suffered the worst-margin-defeat */)
{ /* side effects: Each Candidate's 'eliminated' member */
	int64_t RayDefeatMargin[MaxNumCands]={0};
	int RayBeater[MaxNumCands]={0};
	int i, j;
	int64_t x;
	int64_t t;
	int RayLoser, rnd;
	int64_t maxc;
	int r, beater;
	const uint64_t& numberOfCandidates = E.NumCands;
	CandidateSlate& allCandidates = E.Candidates;
	Determine(CopeWinOnlyWinner, BuildDefeatsMatrix, E);
#if defined(CWSPEEDUP) && CWSPEEDUP
	if(CondorcetWinner >= 0) return CondorcetWinner;
#endif
	for(i=(int)numberOfCandidates-1; i>=0; i--) {
		t = -BIGINT; beater = -1;
		RandomlyPermute( numberOfCandidates, RandCandPerm );
		for(j=(int)numberOfCandidates-1; j>=0; j--) {
			r = RandCandPerm[j];
			x = allCandidates[r].margins[i];
			if(x>t) { t=x; beater=r; }
		}
		assert(beater >= 0);
		RayDefeatMargin[i] = t; /*worst margin of defeat of i, nonpositive if undefeated */
		RayBeater[i] = beater; /*who administered that beating*/
	}
	Zero(numberOfCandidates, allCandidates, &oneCandidate::eliminated);
	for(rnd=1; rnd < (int)numberOfCandidates; rnd++) {
		RayLoser = -1;
		maxc = -BIGINT;
		RandomlyPermute( numberOfCandidates, RandCandPerm );
		for(i=(int)numberOfCandidates-1; i>=0; i--) {
			r = RandCandPerm[i];
			if(!allCandidates[r].eliminated && (RayDefeatMargin[r]>maxc)) {
				maxc=RayDefeatMargin[r];
				RayLoser=r;
			}
		}
		assert(RayLoser >= 0);
		ensure(RayLoser >= 0, 10);
		if( maxc <= 0 ) { return RayLoser; } /*"loser" is undefeated*/
		allCandidates[RayLoser].eliminated = true;
		for(i=(int)numberOfCandidates-1; i>=0; i--) {
			if(!allCandidates[i].eliminated && (RayBeater[i]==RayLoser)) {
				t = -BIGINT;
				beater = -1;
				RandomlyPermute( numberOfCandidates, RandCandPerm );
				for(j=(int)numberOfCandidates-1; j>=0; j--) {
					r = RandCandPerm[j];
					if(!allCandidates[r].eliminated) {
						x = allCandidates[r].margins[i];
						if(x>t) { t=x; beater=r; }
					}
				}
				assert(beater >= 0);
				RayDefeatMargin[i] = t;
				RayBeater[i] = beater;
			}
		}
	} /* end of for(rnd) */
	for(i=0; i<numberOfCandidates; i++) { /* find non-eliminated candidate... */
		if(!allCandidates[i].eliminated) {
			return i; /*Raynaud winner*/
		}
	}
	return(-1); /*error*/
}

/*	ArrowRaynaud(E):	repeatedly eliminates the Candidate with smallest
 *				largest-margin-of-victory; the index of the remaining
 *				Candidate, the Arrow-Raynaud Winner, is returned or -1
 *				if an error occurs; cf., "Solving the Discrete
 *				MultiCriterion Problem: An Overview of Lessons Learned
 *				from Social Choice Literature" by Giuseppe Munda,
 *				section 5, for a complete description of the algorithm;
 *				note: this code may not accurately apply the
 *				Arrow-Raynaud algorithm and should be examined for any
 *				needed alterations
 *	E:	the election data used to determine the Arrow-Raynaud Winner
 */
EMETH ArrowRaynaud(edata& E  /* repeatedly eliminate canddt with smallest {largest margin of victory, which is <=0 if never won} */)
{ /* side effects: Each Candidate's 'eliminated' member, ArrowRaynaud can eliminate a Condorcet Winner in round #1. */
	int64_t ARVictMargin[MaxNumCands]={0};
	int ARchump[MaxNumCands]={0};
	int i, j;
	int64_t x;
	int64_t t;
	int ARLoser, rnd;
	int r, chump;
	const uint64_t& numberOfCandidates = E.NumCands;
	CandidateSlate& allCandidates = E.Candidates;
	Determine(CopeWinOnlyWinner, BuildDefeatsMatrix, E);
	for(i=(int)numberOfCandidates-1; i>=0; i--) {
		const MarginsData& marginsOfI = allCandidates[i].margins;
		t = -BIGINT; chump = -1;
		RandomlyPermute( numberOfCandidates, RandCandPerm );
		for(j=(int)numberOfCandidates-1; j>=0; j--) {
			r = RandCandPerm[j];
			x = marginsOfI[r];
			if(x>t) { t=x; chump=r; }
		}
		assert(chump >= 0);
		ARVictMargin[i] = t; /*largest margin of victory of i, nonpositive if never won*/
		ARchump[i] = chump; /*who suffered that beating*/
	}
	Zero(numberOfCandidates, allCandidates, &oneCandidate::eliminated);
	for(rnd=1; rnd < (int)numberOfCandidates; rnd++) {
		RandomlyPermute( numberOfCandidates, RandCandPerm );
		ARLoser = findLoser(numberOfCandidates, allCandidates, ARVictMargin);
		assert(ARLoser >= 0);
		ensure(ARLoser >= 0, 11);
		allCandidates[ARLoser].eliminated = true;
		for(i=(int)numberOfCandidates-1; i>=0; i--) {
			const oneCandidate& CandidateI = allCandidates[i];
			if(!CandidateI.eliminated && (ARchump[i]==ARLoser)) {
				t = -BIGINT;
				chump = -1;
				RandomlyPermute( numberOfCandidates, RandCandPerm );
				for(j=(int)numberOfCandidates-1; j>=0; j--) {
					r = RandCandPerm[j];
					if(!allCandidates[r].eliminated) {
						x = CandidateI.margins[r];
						if(x>t) {
							t=x;
							chump=r;
						}
					}
				}
				assert(chump >= 0);
				ARVictMargin[i] = t;
				ARchump[i] = chump;
			}
		}
	} /* end of for(rnd) */
	for(i=0; i<numberOfCandidates; i++) { /* find non-eliminated candidate... */
		if(!allCandidates[i].eliminated) {
			return i; /*ArrowRaynaud winner*/
		}
	}
	return(-1); /*error*/
}

/*	beatPathWinnerExists(bpsArray, k, count):	examines
 *							'bpsArray'
 *							to see if
 *							a Schulze
 *							Beat Path
 *							Winner exists;
 *							if so, 'true'
 *							is returned;
 *							otherwise,
 *							'false'
 *							is
 *	bpsArray:	a beat path strength array
 *	k:		the index of a Candidate to consider
 *	count:		the number of Candidates in the current
 *			election
 */
template<class T> bool beatPathWinnerExists(const T (&bpsArray)[MaxNumCands*MaxNumCands], const int& k, const uint64_t& count)
{
	for(int j=0; j<count; j++) {
		if(k != j) {
			if( bpsArray[j*count +k] > bpsArray[k*count +j] ) {
				return false;
			}
		}
	}
	return true;
}

/*	SchulzeBeatpaths(E):	returns the index of the Schulze beatpaths Winner; this
 *				implementation is an O(N^3) algorithm, but it is known
 *				how to speed it up to O(N^2); the Schulze beatpath
 *				algorithm is described below
 *	E:	the election data used to determine the Schulze beatpaths Winner
 *
 *	===Algorith===
 *	(1) "votes": are orderings of the candidates of the form A>B>C=D>E (equalities
 *		are permitted). All candidates must be ordered and none omitted (if any
 *		are omitted, the system either refuses to accept your vote, or tells you
 *		it will assume all the unlisted candidates are ranked coequal last –
 *      	either way you have no way to express ignorance about any candidate).
 *	(2) For each pair of candidates (a,b) we compute the "winning vote count" for a
 *		over b, which is the number of voters whose vote said a>b, assuming this
 *		number exceeds the count for b>a. It would also be possible to base the
 *		method on "margins" (the number of votes saying A>B, minus the number
 *		saying B>A) but Schulze himself does not recommend that. The two are
 *		equaivalent if there are no candidate-equalities in any votes. The
 *		result is a "directed graph" with N vertices wheree winning-vote count
 *		we just described.
 *	(3) Now, in this directed graph, a "beatpath" is a directed path of edges,
 *		always walking in the direction of an arrow, which leads from some
 *		candidate L to some other W. The "strength" of this beatpath is the
 *		minimum value of the numerical labels on its arrows. Non-obviously, it
 *		is possible to compute all the beatpath strengths in polynomial time by
 *		a modification, invented by Schulze, of R.W.Floyd's dynamic programming
 *		algorithm for digraph all-shortest-paths.
 *	(4) If the strongest path from L to W, is stronger than, or at least as strong
 *		as, the strongest path from W to L, and if this is simultaneously true
 *		for every L, then W is a "Schulze winner." Schulze proved the theorem
 *		that such a W always exists (at least using "margins"; I am confused re
 *		the "winning-votes" enhancement).
 *	(5) Usually (i.e. in the absence of unlikely "ties"), W will be unique, but ties
 *		for winner are possible. In that case, Schulze advocates breaking the
 *		ties by a rather complicated procedure, below.
 *	(6) In this procedure, it perhaps is possible or desirable to further enhance it
 *		by permitting "partial-order-type" votes so voters can express
 *		ignorance/no-opinion. Some such votes would be illegal since they
 *		contain a "preference cycle" and hence make no logical sense. They would
 *		have to be detected and rejected. That would add considerable additional
 *		complexity to the algorithm.
 *	===Tie Breaking===
 *	If there is only one potential winner, then this potential winner is the unique
 *	winner. If there is more than one potential winner, then a Tie Breaking Ranking
 *	of the Candidates (TBRC) is calculated as follows:
 *	1. Pick a random ballot and use its rankings; consider ties as unsorted with
 *		regard to each other.
 *	2. Continue picking ballots randomly from those that have not yet been picked.
 *		When you find one that orders previously unsorted candidates, use the
 *		ballot to sort them. Do not change the order of the already sorted.
 *	3. If you go through all ballots, and some candidates are still not sorted,
 *		order them randomly.
 *	The winner is that potential winner who is ranked highest in this TBRC.
 *
 *	Note: However, if I read this implementation correctly, the TBRC approach is not
 *	yet implemented.
 */
EMETH SchulzeBeatpaths(edata& E  /* winner = X so BeatPathStrength over rivals Y exceeds strength from Y */)
{
	int64_t BeatPathStrength[MaxNumCands*MaxNumCands]={0};
	int i,j,k;
	int64_t minc;
	int winner;
	const uint64_t& numberOfCandidates = E.NumCands;
	const CandidateSlate& allCandidates = E.Candidates;
	Determine(CopeWinOnlyWinner, BuildDefeatsMatrix, E);
	for(i=0; i<numberOfCandidates; i++) {
		const MarginsData& marginsOfI = allCandidates[i].margins;
		for(j=0; j<numberOfCandidates; j++) {
			if(i != j) {
				BeatPathStrength[i*numberOfCandidates +j] = marginsOfI[j];
			}
		}
	}
	for(i=0; i<numberOfCandidates; i++) {
		for(j=0; j<numberOfCandidates; j++) {
			if(i != j) {
				for(k=0; k<(int)numberOfCandidates; k++) {
					if((k != j) && (k != i)) {
						minc = BeatPathStrength[j*numberOfCandidates+i];
						if( BeatPathStrength[i*numberOfCandidates +k] < minc ) {
							minc = BeatPathStrength[i*numberOfCandidates +k];
						}
						if( BeatPathStrength[j*numberOfCandidates +k] < minc ) {
							BeatPathStrength[j*numberOfCandidates +k] = minc;
						}
					}
				}
			}
		}
	}
	for(i=(int)numberOfCandidates-1; i>=0; i--) {
		bool haveAWinner;
		k = RandCandPerm[i];
		haveAWinner = beatPathWinnerExists(BeatPathStrength, k, numberOfCandidates);
		if( haveAWinner ) {
			winner = k;   return winner;
		}
	}
	return(-1);
}

/* Note about Smith and Schwartz sets:
BeatPathStrength[k*numberOfCandidates+j] > 0   for all k in the "Smith Set" and j outside it.
BeatPathStrength[k*numberOfCandidates+j] >= 0  for all k in the "Schwartz Set" and j outside it.
*******/

void beatDFS( const int& x, const int& diff, CandidateSlate& relationships, uint64_t N, bool oneCandidate::*member )
{
	int i;
	for(i=(int)N-1; i>=0; i--) {
		if(i!=x) {
			oneCandidate& theCandidate = relationships[i];
			if( theCandidate.margins[x]>=diff ) {
				bool& theMember = theCandidate.*member;
				if( !theMember ) {
					theMember = true;
					beatDFS( i, diff, relationships, N, member );
				}
			}
		}
	}
}

/*	SmithSet(E):	returns the index of a randomly selected Member of the Smith set
 *			or -1 if an error occurs
 *	E:	the election data used to determine the Smith set
 */
EMETH SmithSet(edata& E  /* Smith set = smallest nonempty set of canddts that pairwise-beat all nonmembers */)
{ /* side effects: Each Candidate's Smith member status */
	int i,r;
	const uint64_t& numberOfCandidates = E.NumCands;
	CandidateSlate& allCandidates = E.Candidates;
	Zero(numberOfCandidates, allCandidates, &oneCandidate::IsASmithMember);
	assert(CopeWinOnlyWinner>=0);
	assert(CopeWinOnlyWinner < (int)numberOfCandidates);
	allCandidates[CopeWinOnlyWinner].IsASmithMember = true;
	beatDFS(CopeWinOnlyWinner, 1, allCandidates, numberOfCandidates, &oneCandidate::IsASmithMember);
	RandomlyPermute( numberOfCandidates, RandCandPerm );
	for(i=(int)numberOfCandidates-1; i>=0; i--) {
		r = RandCandPerm[i];
		if(allCandidates[r].IsASmithMember) {
			return r; /*return random set member*/
		}
	}
	return(-1);
}

//	Function: SchwartzSet
//
//	Returns:
//		an index representing a randomly selected
//		Candidate which is also a Member of the Schwartz set
//		or -1 if an error occurs
//
//	Parameter:
//		E - the election data to use for determining the
//		    Schwartz set
EMETH SchwartzSet(edata& E  /* Schwartz set = smallest nonempty set of canddts undefeated by nonmembers */)
{ /* side effects: Each Candidate's Schwartz Member state */
	int i,r;
	const uint64_t& numberOfCandidates = E.NumCands;
	CandidateSlate& allCandidates = E.Candidates;
	Zero(numberOfCandidates, allCandidates, &oneCandidate::IsASchwartzMember);
	assert(CopeWinOnlyWinner>=0);
	assert(CopeWinOnlyWinner < (int)numberOfCandidates);
	allCandidates[CopeWinOnlyWinner].IsASchwartzMember = true;
	beatDFS(CopeWinOnlyWinner, 0, allCandidates, numberOfCandidates, &oneCandidate::IsASchwartzMember);
	RandomlyPermute( numberOfCandidates, RandCandPerm );
	for(i=(int)numberOfCandidates-1; i>=0; i--) {
		r = RandCandPerm[i];
		if(allCandidates[r].IsASchwartzMember) {
			return r; /*return random set member*/
		}
	}
	return(-1);
}

/* Uncovered set.  A "covers" B if the candidates A beats pairwise
are a superset of those B beats pairwise.  "Landau's theorem":
All Copeland winners are members of the uncovered set.
Now an order of magnitude faster assuming MaxNumCands<sizeof(uint)*4 because set inclusion can
be tested in 1 step using wordwide operations.
****/

//	Function: UncoveredSet
//
//	Returns:
//		an index representing a randomly selected
//		Candidate from the set of All "uncovered"
//		Candidates; a Candidate "covers" another if the
//		set of Candidates which the first Candidate
//		beats pairwise is a superset of the set which
//		the second Candidate beats pairwise; a Candidate
//		is "uncovered" if no Candidate covers said
//		Candidate; if no uncovered Candidate exists, -1
//		is returned
//
//	Parameter:
//		E - the election data to use for determining the
//		    uncovered Winner
EMETH UncoveredSet(edata& E)
{ /* side effects: Each Candidate's 'uncovered' state, CoverMatrix[] */
	int A,B,i,r;
	const uint64_t& numberOfCandidates = E.NumCands;
	CandidateSlate& allCandidates = E.Candidates;
	if( numberOfCandidates > 4*sizeof(numberOfCandidates) ) {
		printf("UncoveredSet: too many candidates %lld to use machine words(%d) to represent sets\n",
			numberOfCandidates,
			(int)(4*sizeof(numberOfCandidates)) );
		printf("You could rewrite the code to use uint128s to try allow up to 128 canddts\n");
		exit(EXIT_FAILURE);
	}
	Determine(CopeWinOnlyWinner, BuildDefeatsMatrix, E);
	Determine(SchwartzWinner, SchwartzSet, E);
	/*find cover relation:*/
	for(A=0; A < (int)numberOfCandidates; A++) {
		const uint64_t& Abeats = allCandidates[A].Ibeat;
		for(B=0; B < (int)numberOfCandidates; B++) {
			if(B!=A) {
				CoverMatrix[A*numberOfCandidates + B] = StrictSuperset(Abeats, allCandidates[B].Ibeat);
			}
		}
		allCandidates[A].uncovered = true;	/*initialization*/
	}
	/*find 'uncovered' Candidates:*/
	for(A=0; A < (int)numberOfCandidates; A++) {
		for(B=0; B < (int)numberOfCandidates; B++) {
			if(B!=A) {
				if( CoverMatrix[B*numberOfCandidates+A] ) {
					allCandidates[A].uncovered = false;
					break;
				}
			}
		}
	}
	/*select random uncovered winner:*/
	RandomlyPermute( numberOfCandidates, RandCandPerm );
	for(i=(int)numberOfCandidates-1; i>=0; i--) {
		const oneCandidate& CandidateI = allCandidates[i];
		const bool& IIsUncovered = CandidateI.uncovered;
		const bool& IIsASchwartz = CandidateI.IsASchwartzMember;
		assert( !IIsUncovered || IIsASchwartz );
		ensure( !IIsUncovered || IIsASchwartz, 55 );
		r = RandCandPerm[i];
		if(allCandidates[r].uncovered){
			RandomUncoveredMemb = r;
			return r;
		}
	}
	return(-1);
}

//	Function: Bucklin
//
//	Returns:
//		an index representing a randomly selected Candidate
//		which is also a Bucklin Winner or -1 if an error
//		occurs; Bucklin Winners are determine by first
//		counting the first choice votes of All Voters; if a
//		Candidate has a majority, that Candidate wins;
//		otherwise, the second choices are added to the first
//		choices; again, if a Candidate with a majority vote is
//		found, the Winner is the Candidate with the most votes
//		accumulated; lower rankings are added as needed; yes,
//		it is possible to find multiple Candidates with a
//		majority from the 2nd round on
//
//	Parameter:
//		E - the election data to use for determining the
//		    Bucklin Winner(s)
EMETH Bucklin(edata& E)
{ /* side effects: Each Candidates 'voteCountForThisRound' */
	int rnd,i,winner;
	int64_t best;
	const uint64_t& numberOfCandidates = E.NumCands;
	const uint& numberOfVoters = E.NumVoters;
	const oneVoter (&allVoters)[MaxNumVoters] = E.Voters;
	CandidateSlate& allCandidates = E.Candidates;
	winner = -1;
	Zero(numberOfCandidates, allCandidates, &oneCandidate::voteCountForThisRound);
	for(rnd=0; rnd<(int)numberOfCandidates; rnd++) {
		for(i=0; i<(int)numberOfVoters; i++) {
			const oneVoter& theVoter = allVoters[i];
			allCandidates[ theVoter.topDownPrefs[rnd] ].voteCountForThisRound++;
		}
		winner = Maximum(numberOfCandidates, allCandidates, &oneCandidate::voteCountForThisRound);
		best = allCandidates[winner].voteCountForThisRound;
		if((best*2) > numberOfVoters) {
			break;
		}
	}
	return winner;
}

/******
James Green-Armytage:  Definition of the cardinal-weighted pairwise comparison method
I. Voters rate the candidates, e.g. on a scale from 0 to 100.
Equal ratings are allowed.   There then is an implied ranking.
II. Tally
1. Determine the direction of the pairwise defeats by using the rankings for a standard pairwise
comparison tally.
2. Determine the strength of the pairwise defeats by finding the weighted magnitude as follows.
Suppose candidate A pairwise beats B, and we want to know the strength of the defeat.
For each voter who ranks A over B, and only for voters who rank A over B,
subtract their rating of B from their rating of A, to get the rating differential.
The sum of these individual winning rating differentials is the total weighted magnitude of the defeat.
(Note that, because voters who rank B over A do not contribute to the weighted magnitude of the A>B defeat,
it cannot be negative.)
3. Now that the direction of the pairwise defeats have been determined (in step 1) and
the strength of the defeats have been determined (in step 2), you can choose from a
variety of Condorcet completion methods to determine the winner.
I recommend the ranked pairs, beatpath, and river methods.
WDS: IEVS will use beatpaths.  What is good voter strategy in this method?
*********/

EMETH ArmytagePCSchulze(edata& E  /*Armytage pairwise comparison based on Schulze*/)
{
	real ArmyBPS[MaxNumCands*MaxNumCands]={0};
	int i,j,k,winner;
	real minc;
	const uint64_t& numberOfCandidates = E.NumCands;
	const CandidateSlate& allCandidates = E.Candidates;
	Determine(CopeWinOnlyWinner, BuildDefeatsMatrix, E);
	for(i=0; i<numberOfCandidates; i++) {
		const ArmytageMarginData& ArmytageMarginsOfI = allCandidates[i].ArmytageMarginsMatrix;
		for(j=0; j<numberOfCandidates; j++) {
			if(i != j) {
				ArmyBPS[i*numberOfCandidates +j] = ArmytageMarginsOfI[j];
			}
		}
	}
	for(i=0; i<numberOfCandidates; i++) {
		for(j=0; j<numberOfCandidates; j++) {
			if(i != j) {
				for(k=0; k<(int)numberOfCandidates; k++) {
					if(k != j && k != i) {
						minc = ArmyBPS[j*numberOfCandidates+i];
						if( ArmyBPS[i*numberOfCandidates +k] < minc ) {
							minc = ArmyBPS[i*numberOfCandidates +k];
						}
						if( ArmyBPS[j*numberOfCandidates +k] < minc ) {
							ArmyBPS[j*numberOfCandidates +k] = minc;
						}
					}
				}
			}
		}
	}
	for(i=(int)numberOfCandidates-1; i>=0; i--) {
		bool haveAWinner;
		k = RandCandPerm[i];
		haveAWinner = beatPathWinnerExists(ArmyBPS, k, numberOfCandidates);
		if(haveAWinner) {
			winner = k;   return winner;
		}
	}
	return(-1);
}

/*	Copeland(E):	returns the Copland Winner; Copeland's method or Copeland's
 *			pair-wise aggregation method is a Condorcet method in which
 *			Candidates are ordered by the number of pair-wise victories,
 *			minus the number of pair-wise defeats; it is not entirely clear
 *			whether this implementation takes into account pair-wise defeats
 *	E:	the election data used to determine the Winner
 */
EMETH Copeland(edata& E   /* canddt with largest number of pairwise-wins elected (tie counts as win/2) BUGGY */)
{
	int CopelandWinner;
	uint64_t CopeScore[MaxNumCands]={0};
	int i;
	const uint64_t& numberOfCandidates = E.NumCands;
	const CandidateSlate& allCandidates = E.Candidates;
	Determine(CopeWinOnlyWinner, BuildDefeatsMatrix, E);
#if defined(CWSPEEDUP) && CWSPEEDUP
	if(CondorcetWinner >= 0) {
		CopelandWinner = CondorcetWinner;
		return CopelandWinner;
	}
#endif
	for(i=0; i<numberOfCandidates; i++) {
		const oneCandidate& theCandidate = allCandidates[i];
		CopeScore[i] = (2*theCandidate.electedCount)+theCandidate.drawCount;
	}
	CopelandWinner = ArgMaxArr<uint64_t>(numberOfCandidates, CopeScore, (int*)RandCandPerm);
	/* Currently just break ties randomly, return random highest-scorer */
	return CopelandWinner;
}

typedef EMETH (methodFunction)(edata&);

/*	Determine(Winner, theMethod, E):	determines the
 *						Winner of an
 *						election
 *						conducted using
 *						a given voting
 *						method
 *	Winner:		the Winner to determine; receives the
 *			results, if needed
 *	theMethod:	the voting method to determine the
 *			Winner
 *	E:		the election data to use to determine
 *			the Winner
 */
void Determine(int& Winner, methodFunction theMethod, edata& E)
{
	if(0 > Winner) {
		theMethod(E);
	}
}

//	Function: calculateSimmonsVotesAgainst
//
//	Returns:
//		the sum of all plurality votes for Each Rival of
//		the given Candidate which beats said Candidate
//		pairwise; to help break ties in the
//		"SimmonsCond" system, an additional 1/3 of a
//		vote is added in certain cases
//	Parameters:
//		theCandidate       - the Candidate against which
//		                     to compare Rivals
//		allCandidates      - a slate of Candidates for
//		                     this election
//		numberOfCandidates - the number of Candidates in
//		                     the current election
int64_t calculateSimmonsVotesAgainst(const oneCandidate& theCandidate,
                                     const CandidateSlate& allCandidates,
                                     const uint64_t& numberOfCandidates)
{
	int64_t t=0;
	const MarginsData& margins = theCandidate.margins;
	for(uint64_t j=0; j<numberOfCandidates; j++) {
		if(margins[j]<0) {
			const auto& CandidateJ = allCandidates[j];
			const auto& pluralityVotes = CandidateJ.pluralityVotes;
			int64_t SchwartzFactor;
			int64_t SmithFactor;
			if(CandidateJ.IsASchwartzMember) {
				SchwartzFactor = 1;
			} else {
				SchwartzFactor = 0;
			}
			if(CandidateJ.IsASmithMember) {
				SmithFactor = 1;
			} else {
				SmithFactor = 0;
			}
			t += 3*pluralityVotes + SchwartzFactor + SmithFactor;
		}
	}
	return t;
}

//	Function: SimmonsCond
//
//	Returns:
//		the index of a randomly selected Candidate from
//		All Candidates with the least sum of top-rank-votes
//		for Rivals pairwise-beating said Candidate; note:
//		this implementation also adds 1/3 of a vote in certain
//		cases in order to help break ties
//	Parameters:
//		E	- the election data to use for determining
//			  the Simmons Condorect Winner
EMETH SimmonsCond(edata& E)
{
	int64_t SimmVotesAgainst[MaxNumCands]={0};
	int i;
	int winner;
	const uint64_t& numberOfCandidates = E.NumCands;
	const CandidateSlate& allCandidates = E.Candidates;
	Determine(CopeWinOnlyWinner, BuildDefeatsMatrix, E);
	Determine(PlurWinner, Plurality, E);
	Determine(SmithWinner, SmithSet, E);
	Determine(SchwartzWinner, SchwartzSet, E);
#if defined(CWSPEEDUP) && CWSPEEDUP
	if(CondorcetWinner>=0) return(CondorcetWinner);
#endif
	for(i=0; i<numberOfCandidates; i++) {
		SimmVotesAgainst[i] = calculateSimmonsVotesAgainst(allCandidates[i], allCandidates, numberOfCandidates);
	}
	winner = ArgMinArr(numberOfCandidates, SimmVotesAgainst, (int*)RandCandPerm);
	return winner;
}

//	Function: countLosses
//
//	Returns:
//		the number of Candidates against which a given
//		Candidate has a negative electoral margin
//
//	Parameters:
//		theCandidate       - the Candidate for which to
//		                     count losses
//		numberOfCandidates - the number of Candidates in
//		                     the current election
int64_t countLosses(const  oneCandidate& theCandidate, const uint64_t& numberOfCandidates)
{
	const MarginsData& margins = theCandidate.margins;
	int64_t t=0;
	for(uint64_t j=0; j<numberOfCandidates; j++) {
		if(margins[j] < 0) {
			t++;
		}
	}
	return t;
}

/*******
 * The following IRV (instant runoff voting) algorithm has somewhat long code, but
 * by using linked lists achieves fast (very sublinear) runtime.
 * This code can also be used to compute the winner in "Top3" (bastardized) IRV where
 * voters only allowed to indicate their top 3 choices in order..
 * For normal IRV (or SmithIRV) set IRVTopLim=BIGINT before run.
 * For Top3-IRV set IRVTopLim=3 (or any other integer N for TopN-IRV). ***/
//	Function: IRV
//
//	Returns:
//		the index of the instant runoff voting Winner or
//		-1 if an error occurs; if SmithIRVwinner<0, the
//		SmithIRVwinner is also determined as a side effect
//	Parameters:
//		E	- the election data to use for determining
//			  the instant runoff voting Winner
EMETH IRV(edata& E   /* instant runoff; repeatedly eliminate plurality loser */)
{ /* side effects: Each Candidate's 'eliminated' member, 'favoriteCandidate's of Each Voter, Each Candidate's 'voteCountForThisRound', FavListNext[], HeadFav[], Each Candidate's 'lossCount' member, SmithIRVwinner, IRVwinner  */
	int Iround,i,RdLoser,NextI,j;
	int64_t t;
	int x,stillthere,winner;
	const uint64_t& numberOfCandidates = E.NumCands;
	const uint& numberOfVoters = E.NumVoters;
	oneVoter (&allVoters)[MaxNumVoters] = E.Voters;
	CandidateSlate& allCandidates = E.Candidates;
	assert(numberOfCandidates <= MaxNumCands);
	if((SmithIRVwinner<0) && (IRVTopLim==BIGINT) && (CopeWinOnlyWinner<0)) {
		BuildDefeatsMatrix(E);
	}
	RandomlyPermute( numberOfCandidates, RandCandPerm );
	for(i=0; i<numberOfCandidates; i++) {
		oneCandidate& CandidateI = allCandidates[i];
		CandidateI.eliminated =false;
		HeadFav[i] = -1; /*HeadFav[i] will be the first voter whose current favorite is i*/
		if((SmithIRVwinner<0) && (IRVTopLim==BIGINT)) {
			CandidateI.lossCount = countLosses(CandidateI, numberOfCandidates);
		}
	} /*end for(i)*/
	Zero(numberOfCandidates, allCandidates, &oneCandidate::voteCountForThisRound);
	resetFavorites(allVoters);
	/* 'favoriteCandidate' is the rank of the 1st noneliminated canddt in voter i's topdownpref list (initially 0) */
	FillArray(numberOfVoters, FavListNext, -1);
	/* FavListNext is "next" indices in linked list of voters with common current favorite; -1 terminated. */
	if((SmithIRVwinner<0) && (IRVTopLim==BIGINT)) {
		SmithIRVwinner = CondorcetWinner;
	}
	/* compute vote totals for 1st round and set up forward-linked lists (-1 terminates each list): */
	for(i=0; i<numberOfVoters; i++) {
		const oneVoter& theVoter = allVoters[i];
		x = theVoter.topDownPrefs[theVoter.favoriteCandidate]; /* the initial favorite of voter i */
		assert(x >= 0);
		assert(x < (int)numberOfCandidates);
		allCandidates[x].voteCountForThisRound++;
		FavListNext[i] = HeadFav[x];
		HeadFav[x] = i;
	}
	RandomlyPermute( numberOfCandidates, RandCandPerm );
	for(Iround=1; Iround < (int)numberOfCandidates; Iround++) { /*perform IRV rounds*/
		RdLoser = Minimum(numberOfCandidates, allCandidates, &oneCandidate::voteCountForThisRound, false, true);
		assert(RdLoser>=0);
		assert(RdLoser < (int)numberOfCandidates);
		ensure(RdLoser>=0, 12);
		allCandidates[RdLoser].eliminated = true; /* eliminate RdLoser */
		if((IRVTopLim==BIGINT) && (SmithIRVwinner < 0)) {
			const MarginsData& marginsOfRdLoser = allCandidates[RdLoser].margins;
			for(j=(int)numberOfCandidates-1; j>=0; j--) {
				if(!allCandidates[j].eliminated) { /* update j's 'lossCount' member */
					int64_t& lossCount = allCandidates[j].lossCount;
					t = marginsOfRdLoser[j];
					if(t>0) {
						lossCount--;
					}
					if( lossCount <= 0 ) { SmithIRVwinner = j; break; }
				}
			}
		}
		for(i=HeadFav[RdLoser]; i>=0; i=NextI) {/*Go thru linked list of voters with favorite=RdLoser, adjust:*/
			oneVoter& theVoter = allVoters[i];
			int64_t& favorite = theVoter.favoriteCandidate;
			NextI =  FavListNext[i];
			ensure(favorite >= 0, 42);
			ensure(favorite < (int)numberOfCandidates, 43);
			ensure( theVoter.topDownPrefs[ favorite ] == RdLoser, 34 );
			do{
				favorite++;
				x = theVoter.topDownPrefs[favorite];
			} while( allCandidates[x].eliminated);
			/* x is new favorite of voter i (or ran out of favorites) */
			ensure( favorite < (int)numberOfCandidates, 44 );
			assert(x >= 0);
			assert(x < (int)numberOfCandidates);
			/* update favorite-list: */
			FavListNext[i] = HeadFav[x];
			HeadFav[x] = i;
			/* update vote count totals: */
			if(favorite < IRVTopLim) {
				allCandidates[x].voteCountForThisRound++;
			}
		} /*end for(i)*/
	}  /* end of for(Iround) */
	stillthere = 0;
	if(IRVTopLim >= (int)numberOfCandidates) {
		IRVwinner = -1;
	}
	winner = -1;
	for(i=0; i<numberOfCandidates; i++) { /* find non-eliminated candidate... */
		if(!allCandidates[i].eliminated) {
			winner=i;
			stillthere++;
		}
	}
	if(IRVTopLim >= (int)numberOfCandidates) {
		IRVwinner=winner;
	}
	assert(stillthere==1);
	ensure((stillthere==1), 3);
	return(winner);
}

//	Function: SmithIRV
//
//	Returns:
//		the index of a Candidate after eliminate the plurality
//		losing Candidate until an unbeaten Candidate exists;
//		the index returned is the index of the unbeaten
//		Candidate; note: this function must be called after
//		'IRV()'
//	Parameters:
//		E	- the election data to use for determining
//			  the instant runoff voting Winner
EMETH SmithIRV(edata& E)
{ /* must be run after IRV. */
	if(IRVwinner<0){
		SmithIRVwinner = -1;
		IRV(E);
	}
	return SmithIRVwinner;
}

//	Function: Top3IRV
//
//	Returns:
//		the index of the instant runoff voting Winner or
//		-1 if an error occurs; Voters indicate only Their
//		top 3 choices in order; if SmithIRVwinner<0, the
//		SmithIRVwinner is also determined as a side effect
//	Parameters:
//		E	- the election data to use for determining
//			  the instant runoff voting Winner
EMETH Top3IRV(edata& E)
{
	int w;
	IRVTopLim = 3;
	w = IRV(E);
	IRVTopLim = BIGINT;
	return w;
}

enum Direction {Up, Down};

//	Function: findNewFavorite
//
//	Returns:
//		the index of a given Voter's next favorite Candidate
//		which has not yet been eliminated by 'BTRIRV()'
//	Parameters:
//		favorite	- the index of the Voter's current favorite
//				  Candidate
//		preferences	- an array of Candidate indices ranked
//				  according to a given Voter's preferences;
//				  Candidates with indices located at
//				  lower indices into this array are more
//				  prefered than Candidates Whose indices
//				  are located at higher indices into
//				  this array
//		allCandidates	- the suite of All Candidates in this
//				  this election
//		CandidateCount	- the number of Candidates in this election
//		direction	- whether to bump 'favorite' up or down
uint findNewFavorite(int64_t& favorite,
		     const uint (&preferences)[MaxNumCands],
		     const CandidateSlate& allCandidates,
		     const uint64_t& CandidateCount,
		     const Direction& direction)
{
	uint indexOfNewFavorite;
	ensure(favorite >= 0, 36);
	ensure(favorite < (int)CandidateCount, 37);
	ensure((direction == Up) || (direction == Down), 51);
	int64_t delta;
	if(direction == Up) {
		delta = +1;
	} else {
		delta = -1;
	}
	do {
		favorite += delta;
		indexOfNewFavorite = preferences[favorite];
	} while( allCandidates[indexOfNewFavorite].eliminated );
	ensure( favorite < (int)CandidateCount, 38 );
	return indexOfNewFavorite;
}

//	Function: BTRIRV
//
//	Returns:
//		the index of the instant runoff voting Winner or
//		-1 if an error occurs; Candidates eliminated are
//		either the plurality Loser or the plurality "2nd Loser",
//		whichever loses pairwise to the Other
//	Parameters:
//		E	- the election data to use for determining
//			  the instant runoff voting Winner
EMETH BTRIRV(edata& E)
{ /* side effects: Each Candidate's 'eliminated' member, 'favoriteCandidate's, Each Candidate's 'voteCountForThisRound', FavListNext[], HeadFav[], */
	int Iround,i,RdLoser,RdLoser2,NextI;
	const uint64_t& numberOfCandidates = E.NumCands;
	const uint& numberOfVoters = E.NumVoters;
	oneVoter (&allVoters)[MaxNumVoters] = E.Voters;
	CandidateSlate& allCandidates = E.Candidates;
	uint64_t x;
	assert(numberOfCandidates <= MaxNumCands);
#if defined(CWSPEEDUP) && CWSPEEDUP
	if(CondorcetWinner>=0) return(CondorcetWinner);
#endif
	for(i=0; i<numberOfCandidates; i++) {
		allCandidates[i].eliminated = false;
		HeadFav[i] = -1;
	}
	Zero(numberOfCandidates, allCandidates, &oneCandidate::voteCountForThisRound);
	resetFavorites(allVoters);
	/* compute vote totals for 1st round and set up forward-linked lists (-1 terminates each list): */
	for(i=0; i<numberOfVoters; i++) {
		const oneVoter& theVoter = allVoters[i];
		const int64_t& favorite = theVoter.favoriteCandidate;
		ensure(favorite >= 0, 45);
		ensure(favorite < (int)numberOfCandidates, 46);
		x = theVoter.topDownPrefs[favorite]; /* the favorite of voter i */
		assert(x >= 0);
		assert(x < (int)numberOfCandidates);
		allCandidates[x].voteCountForThisRound++;
		FavListNext[i] = HeadFav[x];
		HeadFav[x] = i;
	}
	RandomlyPermute( numberOfCandidates, RandCandPerm );
	for(Iround=1; Iround<(int)numberOfCandidates; Iround++) {
		RdLoser = Minimum(numberOfCandidates, allCandidates, &oneCandidate::voteCountForThisRound, false, true);
		assert(RdLoser>=0);
		bool eliminationState = allCandidates[RdLoser].eliminated;
		allCandidates[RdLoser].eliminated = true;
		RdLoser2 = Minimum(numberOfCandidates, allCandidates, &oneCandidate::voteCountForThisRound, false, true);
		allCandidates[RdLoser].eliminated = eliminationState;
		assert(RdLoser2>=0);
		if( allCandidates[RdLoser].margins[RdLoser2] > 0 ) {
			RdLoser = RdLoser2;
		}
		ensure(RdLoser>=0, 13);
		allCandidates[RdLoser].eliminated = true; /* eliminate RdLoser */
		for(i=HeadFav[RdLoser]; i>=0; i=NextI){ /* Go thru list of voters with favorite=RdLoser, adjust: */
			oneVoter& theVoter = allVoters[i];
			int64_t& favorite = theVoter.favoriteCandidate;
			const uint (&preferences)[MaxNumCands] = theVoter.topDownPrefs;
			ensure( preferences[favorite] == RdLoser, 35 );
			x = findNewFavorite(favorite, preferences, allCandidates, numberOfCandidates, Up);
			NextI =	FavListNext[i];
			/* update favorite-list: */
			FavListNext[i] = HeadFav[x];
			HeadFav[x] = i;
			/* update vote count totals: */
			assert(x < (int)numberOfCandidates);
			allCandidates[x].voteCountForThisRound++;
		}
	}
	for(i=0; i<numberOfCandidates; i++){ /* find the non-eliminated candidate... */
		if(!allCandidates[i].eliminated){
			return i; /*IRV winner*/
		}
	}
	return(-1); /*error*/
}

//	Function: Coombs
//
//	Returns:
//		the index of the Coombs Winner or
//		-1 if an error occurs; Each Voter rank-orders All
//		of the Candidates on Their ballot; if at any time
//		1 Candidate is ranked first (among non-eliminated
//		Candidates) by an absolute majority of the Voters,
//		that Candidate wins; otherwise, the Candidate ranked
//		last (again among non-eliminated Candidates) by
//		the largest number of (or a plurality of) Voters
//		is eliminated; conversely, under instant runoff
//		voting, the Candidate ranked first (among non-eliminated
//		Candidates) by the fewest Voters is eliminated
//	Parameters:
//		E	- the election data to use for determining
//			  the instant runoff voting Winner
EMETH Coombs(edata& E)
{ /*side effects: Each Candidate's 'eliminated' member, 'favoriteCandidate's, Each Candidates 'voteCountForThisRound', FavListNext[], HeadFav[] */
	int Iround,i,RdLoser,NextI,x;
	const uint64_t& numberOfCandidates = E.NumCands;
	const uint& numberOfVoters = E.NumVoters;
	oneVoter (&allVoters)[MaxNumVoters] = E.Voters;
	CandidateSlate& allCandidates = E.Candidates;
	assert(numberOfCandidates <= MaxNumCands);
	for(i=0; i<numberOfCandidates; i++){
		allCandidates[i].eliminated = false;
		HeadFav[i] = -1; /*HeadFav[i] will be the first voter whose current most-hated is i*/
	}
	Zero(numberOfCandidates, allCandidates, &oneCandidate::voteCountForThisRound);
	assert(numberOfCandidates >= 2);
	assert(numberOfCandidates <= MaxNumVoters);
	resetFavorites(allVoters,numberOfCandidates-1);
	/* 'favoriteCandidate' is the rank of the last noneliminated canddt in voter i's topdownpref list
	 * (initially NumCands-1) */
	FillArray(numberOfVoters, FavListNext, -1);
	/* FavListNext is "next" indices in linked list of voters with common current favorite; -1 terminated. */
	/* compute vote totals for 1st round and set up forward-linked lists (-1 terminates each list): */
	for(i=0; i<numberOfVoters; i++) {
		const oneVoter& theVoter = allVoters[i];
		const int64_t& favorite = theVoter.favoriteCandidate;
		ensure(favorite >= 0, 47);
		ensure(favorite < (int)numberOfCandidates, 48);
		x = theVoter.topDownPrefs[favorite]; /* the initial most-hated of voter i */
		assert(x >= 0);
		assert(x < (int)numberOfCandidates);
		allCandidates[x].voteCountForThisRound++; /*antiplurality voting*/
		FavListNext[i] = HeadFav[x];
		HeadFav[x] = i;
	}
	RandomlyPermute( numberOfCandidates, RandCandPerm );
	for(Iround=1; Iround < (int)numberOfCandidates; Iround++) {
		RdLoser = Maximum(numberOfCandidates, allCandidates, &oneCandidate::voteCountForThisRound, false, true);
		assert(RdLoser>=0);
		ensure(RdLoser>=0, 14);
		allCandidates[RdLoser].eliminated = true; /* eliminate RdLoser */
		for(i=HeadFav[RdLoser]; i>=0; i=NextI){/*Go thru linked list of voters with favorite=RdLoser, adjust:*/
			oneVoter& theVoter = allVoters[i];
			const uint (&preferences)[MaxNumCands] = theVoter.topDownPrefs;
			int64_t& favorite = theVoter.favoriteCandidate;
			ensure( preferences[favorite] == RdLoser, 49 );
			x = findNewFavorite(favorite, preferences, allCandidates, numberOfCandidates, Down);
			NextI =	FavListNext[i];
			/* update favorite-list: */
			FavListNext[i] = HeadFav[x];
			HeadFav[x] = i;
			/* update vote count totals: */
			assert(x>=0);
			assert(x < (int)numberOfCandidates);
			allCandidates[x].voteCountForThisRound++;
		}
	}	/* end of for(Iround) */
	for(i=0; i<numberOfCandidates; i++) { /* find the non-eliminated candidate... */
		if(!allCandidates[i].eliminated) {
			return i; /*Coombs winner*/
		}
	}
	return(-1); /*error*/
}


//	Function: addTheApprovalOfTheVoter
//
//	increments the approval count of Each Candidate a given Voter
//	approves
//
//	Parameters:
//		allCandidatesToTheVoter - a collection of Candidates
//		                          as perceived by 1 particular
//		                          Voter
//		allCandidates           - a collection of Candidates
//		                          in general for the current
//		                          election
//		numberOfCandidates      - the number of Candidates in
//		                          the current election
void addTheApprovalOfTheVoter(const oneCandidateToTheVoter (&allCandidatesToTheVoter)[MaxNumCands],
                              CandidateSlate& allCandidates,
                              const uint64_t& numberOfCandidates)
{
	for(int j=0; j<numberOfCandidates; j++) {
		if(allCandidatesToTheVoter[j].approve) {
			allCandidates[j].approvals++;
		}
	}
}

//	Function: Approval
//
//	Returns:
//		the index corresponding to the Winner according to the
//		approval method; the Candidate with the most approvals wins
//	Parameters:
//		E	- the election data used to determine the Winner
EMETH Approval(edata& E)
{ /* side effects: Each Candidate's 'approval' member, ApprovalWinner */
	int i;
	CandidateSlate& allCandidates = E.Candidates;
	const oneVoter (&allVoters)[MaxNumVoters] = E.Voters;
	const uint64_t& numberOfCandidates = E.NumCands;
	const uint& numberOfVoters = E.NumVoters;
	Zero(numberOfCandidates, allCandidates, &oneCandidate::approvals);
	for(i=0; i<(int)numberOfVoters; i++) {
		const oneCandidateToTheVoter (&allCandidatesToTheVoter)[MaxNumCands] = allVoters[i].Candidates;
		addTheApprovalOfTheVoter(allCandidatesToTheVoter, allCandidates, numberOfCandidates);
	}
	RandomlyPermute( numberOfCandidates, RandCandPerm );
	ApprovalWinner = Maximum(numberOfCandidates, allCandidates, &oneCandidate::approvals);
	return(ApprovalWinner);
}

//	Function: App2Runoff
//
//	Returns:
//		the index of the Winner from a run-off between the top
//		2 Candidates as determined by approval voting; the 1st
//		round uses the approval method; the 2nd round has fully-honest
//		voting
//	Parameters:
//		E	- the election data used to determine the Winner
EMETH App2Runoff(edata& E)
{ /* side effects: ASecond */
	EMETH winner;
	Determine(ApprovalWinner, Approval, E);
	winner = runoffForApprovalVoting(E);
	return winner;
}

EMETH HeitzigDFC(edata& E)
{ /*winner of honest runoff between ApprovalWinner and RandomBallot winner*/
	int i, Rwnr;
	uint pwct=0, wct=0;
	const uint64_t& numberOfCandidates = E.NumCands;
	const uint& numberOfVoters = E.NumVoters;
	const uint64_t VoterIndex = RandInt(E.NumVoters);
	const oneVoter (&allVoters)[MaxNumVoters] = E.Voters;
	Rwnr = ArgMaxArr<real>(numberOfCandidates, E.Voters[VoterIndex].Candidates, (int*)RandCandPerm);
	Determine(ApprovalWinner, Approval, E);
	for(i=0; i<(int)numberOfVoters; i++) {
		const oneVoter& theVoter = allVoters[i];
		const oneCandidateToTheVoter (&allCandidatesToTheVoter)[MaxNumCands] = theVoter.Candidates;
		const real perceivedUtilityOfTheApprovalWinner = allCandidatesToTheVoter[ApprovalWinner].perceivedUtility;
		const real perceivedUtilityOfTheRandomBallotWinner = allCandidatesToTheVoter[Rwnr].perceivedUtility;
		if( perceivedUtilityOfTheApprovalWinner > perceivedUtilityOfTheRandomBallotWinner ) {
			pwct++;
		}else if( perceivedUtilityOfTheApprovalWinner < perceivedUtilityOfTheRandomBallotWinner ) {
			wct++;
		} else {
			/* tie, do nothing */
		}
	}
	if( pwct > wct ) {
		return(ApprovalWinner);
	}else if (pwct == wct) {
		if (RandBool()) {
			return(ApprovalWinner);
		}
	} else {
		/* do nothing */
	}
	return(Rwnr);
}

EMETH HeitzigLFC(edata& E)
{ /*random canddt who is not "strongly beat" wins; y strongly beats x if Approval(y)>Approval(x) and less than Approval(x)/2 voters prefer x>y.*/
	/*Side effects: Each Candidate's 'eliminated' member */
	int i,j;
	const uint64_t& numberOfCandidates = E.NumCands;
	CandidateSlate& allCandidates = E.Candidates;
	Zero(numberOfCandidates, allCandidates, &oneCandidate::eliminated);
	for(j=0; j<numberOfCandidates; j++) {
		const oneCandidate& CandidateJ = allCandidates[j];
		for(i=0; i<numberOfCandidates; i++) {
			oneCandidate& CandidateI = allCandidates[i];
			if(CandidateJ.approvals > CandidateI.approvals) {
				if(2*CandidateI.DefeatsMatrix[j] < CandidateI.approvals) {
					/*Candidate i is "strongly beat"*/
					allCandidates[i].eliminated = true;
				}
			}
		}
	}
	RandomlyPermute( numberOfCandidates, RandCandPerm );
	for(i=(int)numberOfCandidates-1; i>=0; i--) { /* find random non-eliminated candidate... */
		j = RandCandPerm[i];
		if(!allCandidates[j].eliminated) {
			return j; /*winner*/
		}
	}
	return(-1); /*error*/
}


/*
reset All Candidates' elimination status
until there is 1 Candidate remaining:
	reset All Candidate's normalized rating sum (NRS)
	for Each Voter:
		normalize each vote vector
		add the normalized vote for Each Candidate to the Candidate's NRS
	randomly permute the Candidates
	by random permutation, find the Candidate with the lowest NRS
	eliminate that Candidate
by random permutation, find the 1 uneliminated Candidate
return that Candidate if found
else return -1
*/

typedef std::array<real, MaxNumCands> voteVector;
typedef voteVector (*normalizationFunction)(const oneVoter&, const CandidateSlate&, const uint64_t&);
real IRNRPOWER=2.0;
/*	IRNR(E, mode):	returns the instant-runoff normalized-ratings Winner or -1 if an
 *			error occurs
 *	E:		the election data used to determine the Winner
 *	normalizer:	the normalization function to use
 */
EMETH IRNR(edata& E, normalizationFunction normalizer /*Brian Olson's voting method described above*/)
{ /* side effects: Each Candidate's 'eliminated' member, SumNormedRatings[]*/
	uint64_t rd;
	int i;
	int j;
	int loser;
	const oneVoter (&allVoters)[MaxNumVoters] = E.Voters;
	const uint64_t& numberOfCandidates = E.NumCands;
	const uint& numberOfVoters = E.NumVoters;
	CandidateSlate& allCandidates = E.Candidates;
	extern void addToNormalizedRatingSum(const voteVector&, CandidateSlate&, const uint64_t&);
	Zero(numberOfCandidates, allCandidates, &oneCandidate::eliminated);
	for(rd=numberOfCandidates; rd>1; rd--) {
		Zero(numberOfCandidates, allCandidates, &oneCandidate::normalizedRatingSum);
		for(i=0; i<(int)numberOfVoters; i++) {
			voteVector normalizedVoteVector = normalizer(allVoters[i], allCandidates, numberOfCandidates);
			addToNormalizedRatingSum(normalizedVoteVector, allCandidates, numberOfCandidates);
		}
		loser = Minimum(numberOfCandidates, allCandidates, &oneCandidate::normalizedRatingSum, true, true);
		assert(loser>=0);
		ensure(loser>=0, 15);
		allCandidates[loser].eliminated = true;
	}
	for(i=(int)numberOfCandidates-1; i>=0; i--) { /* find random non-eliminated candidate... */
		j = RandCandPerm[i];
		if(!allCandidates[j].eliminated) {
			return j; /*winner*/
		}
	}
	return(-1); /*error*/
}

EMETH MCA(edata& E  /*canddt with most-2approvals wins if gets >50%, else regular approval-winner wins*/)
{
	uint MCAVoteCount[MaxNumCands];
	int i;
	int j;
	int winner;
	const oneVoter (&allVoters)[MaxNumVoters] = E.Voters;
	const uint64_t& numberOfCandidates = E.NumCands;
	const uint& numberOfVoters = E.NumVoters;
	Determine(ApprovalWinner, Approval, E);
	ZeroArray( numberOfCandidates, (int*)MCAVoteCount );
	for(i=0; i<(int)numberOfVoters; i++) {
		const oneCandidateToTheVoter (&allCandidates)[MaxNumCands] = allVoters[i].Candidates;
		for(j=0; j<numberOfCandidates; j++) {
			if(allCandidates[j].approve2) {
				MCAVoteCount[j] += 1;
			}
		}
	}
	winner = ArgMaxUIntArr( numberOfCandidates, MCAVoteCount, (int*)RandCandPerm );
	if(2*MCAVoteCount[winner] > numberOfVoters) {
		return(winner);
	}
	return(ApprovalWinner);
}

/***Chris Benham:
This is my suggested rule for picking the two second-round qualifiers
(for use in a top2 runoff) using Approval for the first round:
  "The two finalists are the Approval winner A; and of those candidates
X who would have more approval than A if ballots that make no
approval distinction between A and X were altered to exclusively
approve X, the second qualifier is the X who is most approved on
ballots that don't approve A."
[There might not be a "second qualifier"; then we just make the approval winner win.]
   This is the "instant" version of a (delayed) 2-round system I
suggested a while back.
  This is a big improvement over simply promoting the two most approved
candidates because (a)it kills the strategy of richer factions trying
to capture both runoff spots by running a pair of clones, and (b) it
makes the method much less vulnerable to "turkey raising" strategy
(because voters in the first round can't have their votes do both of
promote one of their sincerely approved candidates and also
a "turkey").

Commentary by Me:
	Benham2AppRunoff(E, alwaysRunoff):	see above
 *
 *	E:		the election data
 *	alwaysRunoff:	whether to always do a run-off election or not
***/
EMETH Benham2AppRunoff(edata& E, bool alwaysRunoff)
{
	int i,j,r;
	int64_t maxc;
	int64_t y;
	const uint64_t& numberOfCandidates = E.NumCands;
	const CandidateSlate& allCandidates = E.Candidates;
	Determine(ApprovalWinner, Approval, E);
	RandomlyPermute( numberOfCandidates, RandCandPerm );
	maxc = -BIGINT;
	j = -1;
	for(i=0; i<(int)numberOfCandidates; i++) {
		r = RandCandPerm[i];
		y = allCandidates[ApprovalWinner].alsoApprovedWith[r];
		if( (allCandidates[r].approvals + y) > allCandidates[ApprovalWinner].approvals ) {
			if( (allCandidates[r].approvals - y) > maxc ) {
				maxc = allCandidates[r].approvals - y;
				j = r;
			}
		}
	}
	if(j<0) {
		if(alwaysRunoff) {
			EMETH winner;
			winner = runoffForApprovalVoting(E);
			return winner;
		} else {
			return(ApprovalWinner);
		}
	}
	/* now for honest runoff between ApprovalWinner and j */
	return calculateForRunoff(E, ApprovalWinner, j);
}

/*	CondorcetApproval(E):	returns the index of the Condorcet Winner, if any, and
 *				the Approval Winner otherwise
 *	E:	the election data used to determine the Winner
 */
EMETH CondorcetApproval(edata& E  /*Condorcet winner if exists, else use Approval*/)
{
	Determine(ApprovalWinner, Approval, E);
	Determine(CopeWinOnlyWinner, BuildDefeatsMatrix, E);
	if(CondorcetWinner >= 0) {
		return CondorcetWinner;
	}
	return ApprovalWinner;
}

EMETH Range(edata& E    /* canddt with highest average Score wins */)
{ /* side effects:   Each Candidate's range vote, RangeWinner  */
	int i;
	int j;
	const oneVoter (&allVoters)[MaxNumVoters] = E.Voters;
	const uint64_t& numberOfCandidates = E.NumCands;
	const uint& numberOfVoters = E.NumVoters;
	CandidateSlate& allCandidates = E.Candidates;
	Zero(numberOfCandidates, allCandidates, &oneCandidate::rangeVote);
	for(i=0; i<numberOfCandidates; i++) {
		allCandidates[i].rangeVote = 0;
	}
	for(i=0; i<(int)numberOfVoters; i++) {
		const oneCandidateToTheVoter (&allCandidatesToTheVoter)[MaxNumCands] = allVoters[i].Candidates;
		for(j=0; j<numberOfCandidates; j++) {
			allCandidates[j].rangeVote += allCandidatesToTheVoter[j].score;
		}
	}
	RangeWinner = Maximum(numberOfCandidates, allCandidates, &oneCandidate::rangeVote);
	return(RangeWinner);
}

EMETH RangeN(const edata& E /*highest average rounded Score [rded to integer in 0..RangeGranul-1] wins*/)
{ /* uses global integer RangeGranul  */
	uint RangeNVoteCount[MaxNumCands];
	int i;
	int j;
	int winner;
	const oneVoter (&allVoters)[MaxNumVoters] = E.Voters;
	const uint64_t& numberOfCandidates = E.NumCands;
	const uint& numberOfVoters = E.NumVoters;
	assert(RangeGranul>=2);
	assert(RangeGranul<=10000000);
	ZeroArray( numberOfCandidates, (int*)RangeNVoteCount );
	for(i=0; i<(int)numberOfVoters; i++) {
		const oneCandidateToTheVoter (&allCandidates)[MaxNumCands] = allVoters[i].Candidates;
		for(j=0; j<numberOfCandidates; j++) {
			RangeNVoteCount[j] += (uint)( (allCandidates[j].score)*(RangeGranul-0.0000000001) );
		}
	}
	winner = ArgMaxUIntArr( numberOfCandidates, RangeNVoteCount, (int*)RandCandPerm );
	return(winner);
}

/*	Range2Runoff(E):	returns the index of the Winner from a run-off between
 *				the top 2 Candidates as determined by range voting
 *	E:	the election data used to determine the Winner
 */
EMETH Range2Runoff(edata& E    /*top-2-runoff, 1stRd=range, 2nd round has fully-honest voting*/)
{
	int RSecond;
	const uint64_t& numberOfCandidates = E.NumCands;
	const CandidateSlate& allCandidates = E.Candidates;
	Determine(RangeWinner, Range, E);
	RandomlyPermute( numberOfCandidates, RandCandPerm );
	RSecond = SecondMaximum(numberOfCandidates, allCandidates, &oneCandidate::rangeVote, RangeWinner);
	assert(RSecond>=0);
	return calculateForRunoff(E, RangeWinner, RSecond);
}

EMETH ContinCumul(const edata& E    /* Renormalize scores so sum(over canddts)=1; then canddt with highest average Score wins */)
{
	real CCumVoteCount[MaxNumCands];
	int i;
	int j;
	int winner;
	real sum;
	const oneVoter (&allVoters)[MaxNumVoters] = E.Voters;
	const uint64_t& numberOfCandidates = E.NumCands;
	const uint& numberOfVoters = E.NumVoters;
	ZeroArray( numberOfCandidates, CCumVoteCount );
	for(i=0; i<(int)numberOfVoters; i++) {
		const oneCandidateToTheVoter (&allCandidates)[MaxNumCands] = allVoters[i].Candidates;
		sum = 0.0;
		for(j=0; j<numberOfCandidates; j++) {
			sum += allCandidates[j].score;
		}
		if(sum > 0.0) {
			sum = 1.0/sum;
		} else {
			sum=0.0;
		}
		for(j=0; j<numberOfCandidates; j++) {
			CCumVoteCount[j] += sum * allCandidates[j].score;
		}
	}
	winner = ArgMaxArr<real>(numberOfCandidates, CCumVoteCount, (int*)RandCandPerm);
	return(winner);
}

EMETH TopMedianRating(const edata& E    /* canddt with highest median Score wins */)
{
	real MedianRating[MaxNumCands]={0};
	real CScoreVec[MaxNumVoters];
	int i;
	int j;
	int winner;
	const oneVoter (&allVoters)[MaxNumVoters] = E.Voters;
	const uint64_t& numberOfCandidates = E.NumCands;
	const uint& numberOfVoters = E.NumVoters;
	for(j=0; j<numberOfCandidates; j++) {
		for(i=0; i<numberOfVoters; i++) {
			const oneCandidateToTheVoter &theCandidate = allVoters[i].Candidates[j];
			CScoreVec[i] = theCandidate.score;
		}
		MedianRating[j] = TwiceMedian<real>(numberOfVoters, CScoreVec);
	}
	winner = ArgMaxArr<real>(numberOfCandidates, MedianRating, (int*)RandCandPerm);
	return(winner);
}

EMETH LoMedianRank(const edata& E    /* canddt with best median ranking wins */)
{
	int64_t MedianRank[MaxNumCands]={0};
	int64_t CRankVec[MaxNumVoters];
	int i,j;
	int winner;
	const uint64_t& numberOfCandidates = E.NumCands;
	const uint& numberOfVoters = E.NumVoters;
	for(j=0; j<numberOfCandidates; j++) {
		for(i=0; i<numberOfVoters; i++) {
			CRankVec[i] = E.Voters[i].Candidates[j].ranking;
		}
		MedianRank[j] = TwiceMedian<int64_t>(numberOfVoters, CRankVec);
		assert( MedianRank[j] >= 0 );
		assert( MedianRank[j] <= 2*((int)numberOfCandidates - 1) );
	}
	winner = ArgMinArr<int64_t>(numberOfCandidates, MedianRank, (int*)RandCandPerm);
	return(winner);
}

/* Tideman ranked pairs with Honest voters.
 * Tideman = Condorcet variant in which you
 * pick the A>B comparison with the largest margin and "lock it in".
 * Then you  pick the next largest one available
 * ("available" means: not already used and not creating a cycle),
 * and continue on.  This creates an ordering of the candidates. The
 * topmost in the ordering wins.  It is a bit tricky to spot
 * the cycles as we go...
 * (This code based on email from Blake Cretney bcretney@postmark.net
 * and runs in order N^4 time with N candidates, i.e slow.  It is possible to speed it
 * up to Otilde(N^2) steps with the use of fast data structures...):
 ******************************************/
//	Function: TidemanRankedPairs
//
//	Returns:
//		an index of a Candidate randomly selected from
//		all Winners according to Tideman ranked pair method
//		or -1 if an error occurs
//
//	Parameters:
//		E	- the election data to use for determining
//			  the Winner
EMETH TidemanRankedPairs(const edata& E  /*lock in comparisons with largest margins not yielding cycle*/)
{  /*side effects: Tpath[] is used as a changeable copy of MarginsMatrix.*/
	CandidateSlate Tpath = E.Candidates;
	int i,r,j,winner;
	const uint64_t& numberOfCandidates = E.NumCands;
#if defined(CWSPEEDUP) && CWSPEEDUP
	if(CondorcetWinner>=0) return(CondorcetWinner);
#endif
	RandomlyPermute( numberOfCandidates, RandCandPerm );
	for(i=0; i < (int)numberOfCandidates; i++) { Tpath[i].margins[i]=BIGINT; }
	/* Whenever a victory
	 * is locked in, the appropriate Tpath[] cell is set to BIGINT.
	 * pi,pj are used with randperm to give precedence to victories higher
	 * in the random permutation (when tie-breaking).
	 * This loop finds the next pair (i,j) to lock in: ********/
	for(;;) {
		int64_t maxp;
		int oi,oj,pi,pj;
		maxp = -BIGINT; j = -1;
		for(pi=0; pi < (int)numberOfCandidates; pi++) {
			oi=RandCandPerm[pi];
			for(pj=pi+1; pj < (int)numberOfCandidates; pj++) {
				oj=RandCandPerm[pj];
				if((Tpath[oi].margins[oj]!=BIGINT)
					&& (Tpath[oj].margins[oi]!=BIGINT)) {/*not locked-out*/
						if(Tpath[oi].margins[oj]>maxp) { maxp=Tpath[oi].margins[oj]; i=oi; j=oj; }
						if(Tpath[oj].margins[oi]>maxp) { maxp=Tpath[oj].margins[oi]; i=oj; j=oi; }
				}
			}
		}
		if(maxp == -BIGINT) {
			break;
		}
		assert(j>=0);
		/********* print the pair and its margin:
			   printf("(%d %d) %d\n",i,j,maxp);
		***********************/
		/*lock in the pair and clobber future no-good pairs:*/
		for(oi=0; oi < (int)numberOfCandidates; oi++) {
			oneCandidate& relationshipsOfOI = Tpath[oi];
			for(oj=0; oj < (int)numberOfCandidates; oj++) {
				if((relationshipsOfOI.margins[i]==BIGINT) && (Tpath[j].margins[oj]==BIGINT)) {
					relationshipsOfOI.margins[oj]=BIGINT;
				}
			}
		}
	}
	/* The above code assumes that pairels has been set properly.
	 * Tpath[*numberOfCandidates +] ends up with the winning row having all cells
	 * set to BIGINT.  In fact, a complete ranking is given,
	 * where Tpath[i*numberOfCandidates +j]==BIGINT means that i is
	 * ranked over j (where i!=j).   So to find the winner: ****/
	winner = -99;
	for(i=0; i < (int)numberOfCandidates; i++) {
		r = RandCandPerm[i];
		for(j=0; j < (int)numberOfCandidates; j++) {
			if(Tpath[r].margins[j] != BIGINT) {
				break;
			}
		}
		if(j >= (int)numberOfCandidates) { winner=r; break; }
	}
	assert(winner >= 0);
	return(winner);
}

//	Function: updateTheRoot
//
//	updates all instances of 'priorRoot' in the given Heitzig River root array to 'newRoot'
//
//	Parameters:
//		theRoots           - a collection of Heitzig River
//		                     root values
//		priorRoot          - the value in 'theRoots' to
//		                     replace
//		newRoot            - the value with which to replace
//		                     'priorRoot'
//		numberOfCandidates - the number of Candidates in
//		                     the current election
void updateTheRoot(int (&theRoots)[MaxNumCands],
                   const int priorRoot,
                   const int& newRoot,
                   const uint64_t& numberOfCandidates)
{
	for(int i=0; i < (int)numberOfCandidates; i++) {
		int& ithRoot = theRoots[i];
		if(ithRoot == priorRoot) {
			ithRoot = newRoot;
		}
	}
}

//	Function: HeitzigRiver
//
//	Returns:
//		the index of a Candidate randomly selected from
//		all Winners according to Jobst Heitzig's River
//		method
//
//	Parameters:
//		E	- the election data to use for determining
//			  the Winner
EMETH HeitzigRiver(const edata& E /*http://lists.electorama.com/pipermail/election-methods-electorama.com/2004-October/013971.html*/)
{
	int Hpotpar[MaxNumCands]={0};
	int Hpar[MaxNumCands]={0};
	int Hroot[MaxNumCands]={0};
	int r;
	uint64_t z;
	int i,j,k,pp,newroot;
	int64_t maxc;
	const uint64_t& numberOfCandidates = E.NumCands;
	const CandidateSlate& allCandidates = E.Candidates;
#if defined(CWSPEEDUP) && CWSPEEDUP
	if(CondorcetWinner>=0) return(CondorcetWinner);
#endif
	/* find potential (and actual) parents: */
	for(i=0; i < (int)numberOfCandidates; i++) {
		maxc = -BIGINT; pp = -1;
		RandomlyPermute( numberOfCandidates, RandCandPerm );
		for(j=0; j < (int)numberOfCandidates; j++) {
			r = RandCandPerm[j];
			if(r!=i && allCandidates[r].margins[i]>maxc) {
				maxc = allCandidates[r].margins[i];
				pp = r;
			}
		}
		assert(pp>=0);
		Hpotpar[i] = pp;
		Hpar[i] = i;
		Hroot[i] = i;
	}
	for(z = numberOfCandidates; ; ) { /* loop that adds arcs: */
		/* scan tree roots to find best new arc to add: */
		maxc = -BIGINT; k = -1;
		RandomlyPermute( numberOfCandidates, RandCandPerm );
		for(i=0; i < (int)numberOfCandidates; i++) {
			r = RandCandPerm[i];
			if(Hpar[r]==r) { /* tree root */
				pp = Hpotpar[r];
				if(maxc < allCandidates[pp].margins[r]) {
					maxc = allCandidates[pp].margins[r];
					k = r;
				}
			}
		}/*end for(i)*/
		assert(k>=0);
		ensure(k>=0, 18);
		/* add it: */
		pp = Hpotpar[k];
		Hpar[k] = pp;
		/* update roots: */
		newroot = Hroot[pp];
		z--;
		if(z<=1) {
			break;
		}
		updateTheRoot(Hroot, Hroot[k], newroot, numberOfCandidates);
		/* update potential parent of newroot: */
		maxc = -BIGINT; k = -1;
		for(j=0; j < (int)numberOfCandidates; j++) {
			r = RandCandPerm[j];
			if(Hroot[r]!=newroot && allCandidates[r].margins[newroot]>maxc) {
				maxc = allCandidates[r].margins[newroot];
				k = r;
			}
		}
		assert(k>=0);
		Hpotpar[newroot] = k;
	}
	return newroot;
}

//	Function: updateLossCount
//
//	decrements the loss count of Each Candidate from 0 thru 'i-1'
//	if Candidate 'i' has a margin over said Candidate
//
//	Parameters:
//		allCandidates - the collection of Candidates in this
//		                election
//		margins       - the margins of a particular Candidate
//		                with respect to All Other Candidates
//		i             - the index of the Candidate to Whom 'margins'
//		                belongs
void updateLossCount(CandidateSlate& allCandidates,
                     const MarginsData& margins,
                     const int& i)
{
	for(int j=0; j<i; j++) {
		if(margins[j]>0) {
			allCandidates[j].lossCount--;
		}
	}
}

//	Function: DMC
//
//	Returns:
//		the index of the Definite Majority Choice (a.k.a.
//		Ranked Approval Voting) Winner or -1 if an error
//		occurs; the method is summarized as "While no
//		undefeated candidates exist, eliminate the least-approved
//		candidate"
//
//	Parameters:
//		E - the election data used to determine the Winner
EMETH DMC(edata& E  /* eliminate least-approved candidate until unbeaten winner exists */)
{ /* side effects: Each Candidate's 'lossCount' member */
	int i;
	const uint64_t& numberOfCandidates = E.NumCands;
	CandidateSlate& allCandidates = E.Candidates;
	Determine(CopeWinOnlyWinner, BuildDefeatsMatrix, E);
#if defined(CWSPEEDUP) && CWSPEEDUP
	if(CondorcetWinner>=0) return(CondorcetWinner);
#endif
	for(i=0; i<numberOfCandidates; i++) {
		oneCandidate& CandidateI = allCandidates[i];
		CandidateI.lossCount = countLosses(CandidateI, numberOfCandidates);
	}
	RandomlyPermute( numberOfCandidates, RandCandPerm );
	PermShellSortDown(numberOfCandidates, allCandidates, &oneCandidate::approvals);
	for(i=(int)numberOfCandidates-1; i>0; i--) {
		const oneCandidate& CandidateI = allCandidates[i];
		const MarginsData& marginsOfI = CandidateI.margins;
		if( CandidateI.lossCount <= 0 ){  return(i); /*winner*/ }
		updateLossCount(allCandidates, marginsOfI, i);
	}
	return(i);
}

void BSbeatDFS( const int& x, const int& diff, bool Set[], const bool (&OK)[MaxNumCands], const CandidateSlate& relationships, uint64_t N )
{
	int i;
	for(i=0; i<N; i++) {
		if(OK[i] && i!=x) {
			if( relationships[i].margins[x]>=diff ) {
				if( !Set[i] ) {
					Set[i] = true;
					BSbeatDFS( i, diff, Set, OK, relationships, N );
				}
			}
		}
	}
}

/*	BramsSanverPrAV(E):	Returns the index of the Winner according to Preference
 *				Approval Voting, as described in the paper "Voting
 *				Systems That Combine Approval and Preference" by Steven
 *				J. Brams and M. Remzi Sanver (March 2006) or -1 if an
 *				error occurs; the Winner under PAV is determined by two
 *				rules, the second comprising two cases:
 *					1. If no candidate, or exactly one candidate,
 *						receives a majority of approval votes,
 *						the PrAV winner is the AV Winner; that
 *						is, the Candidate who receives the most
 *						approval votes.
 *					2. If two or more Candidates receive a majority
 *						of approval votes, then
 *						(i) If One of these Candidates is
 *							preferred by a majority to every
 *							other majority approved
 *							Candidate, then He or She is the
 *							PrAV Winner, even if not the AV
 *							or Condorcet Winner.
 *						(ii) If there is not One
 *							majority-preferred Candidate
 *							because of a cycle among the
 *							majority-approved Candidates,
 *							then the AV Winner among them is
 *							the PrAV Winner, even if not the
 *							AV or Condorcet Winner.
 *				Note: Should see if the Brams/Sanver "Fallback Voting"
 *				method from the same paper is implemented here and, if
 *				not, do so.
 *	E:	the election data used to determine the Winner
 */
EMETH BramsSanverPrAV(edata& E  /*SJ Brams & MR Sanver: Voting Systems That Combine Approval and Preference,2006*/)
{
	bool MajApproved[MaxNumCands]={false};
	bool BSSmithMembs[MaxNumCands]={false};
	int i,j,winner,ctm,CopeWinr,r;
	uint t,maxt;
	const uint64_t& numberOfCandidates = E.NumCands;
	const uint& numberOfVoters = E.NumVoters;
	const CandidateSlate& allCandidates = E.Candidates;
	Determine(CopeWinOnlyWinner, BuildDefeatsMatrix, E);
	Determine(ApprovalWinner, Approval, E);
	ctm=0;
	for(i=0; i<numberOfCandidates; i++) {
		if((2*allCandidates[i].approvals) > numberOfVoters) {
			MajApproved[i] = true;
			ctm++;
		} else {
			MajApproved[i] = false;
		}
	}
	if(ctm<=1) {
		/*if exactly 0 or 1 canddt majority-approved, ApprovalWinner wins*/
		return ApprovalWinner;
	}
	for(i=0; i<numberOfCandidates; i++) {
		if(MajApproved[i]) {
			bool haveAWinner = true;
			const MarginsData& marginsOfI = allCandidates[i].margins;
			for(j=0; j<numberOfCandidates; j++) {
				if((j!=i) && MajApproved[j]) {
					if( marginsOfI[j] <= 0 ) {
						haveAWinner = false;
						break;
					}
				}
			}
			if(haveAWinner) {
				winner = i; /*beats-all winner among >=2 majority-approved canddts wins*/
				return winner;
			}
		}
	}
	/*now Brams&Sanver want the most-approved member of the
	 *top-cycle among the majority-approved canddts, to win: */
	maxt = 0;
	CopeWinr = -1;
	for(i=(int)numberOfCandidates-1; i>=0; i--) {
		if(MajApproved[i]) {
			const MarginsData& marginsOfI = allCandidates[i].margins;
			t = 0;
			for(j=0; j<numberOfCandidates; j++) {
				if((j!=i) && MajApproved[j]) {
					if( marginsOfI[j] > 0 ) { t++; }
				}
			}
			if(t >= maxt) { maxt=t; CopeWinr=i; }
		}
	}
	assert(CopeWinr >= 0);
	ensure(CopeWinr >= 0, 19);
	assert(MajApproved[CopeWinr]);
	FillArray(numberOfCandidates, BSSmithMembs, false);
	BSSmithMembs[CopeWinr] = true;
	BSbeatDFS(CopeWinr, 1, BSSmithMembs, MajApproved, E.Candidates, numberOfCandidates);
	assert(BSSmithMembs[CopeWinr]);
	RandomlyPermute( numberOfCandidates, RandCandPerm );
	winner = -1;
	maxt = 0;
	for(i=(int)numberOfCandidates-1; i>=0; i--) {
		r = RandCandPerm[i];
		if(BSSmithMembs[r] && (allCandidates[r].approvals>maxt)) {
			maxt=allCandidates[r].approvals;
			winner=r;
		}
	}
	return winner;
}

/*	MDDA(E):	returns the majority defeat disqualification approval Winner or
 *			-1 if an error occurs; the procedure for this voting method is
 *			as follows:
 *				(1) The Voter submits a ranking of the Candidates. The
 *					Candidates explicitly ranked are considered
 *					approved by that Voter.
 *				(2) A Candidate is 'dominated' if more than half of the
 *					Voters rank some other Candidate strictly above
 *					Him/Her.
 *				(3) All dominated Candidates are eliminated, unless this
 *					would eliminate all the Candidates.
 *				(4) Of remaining Candidates, the One approved by the
 *					most Voters is elected.
 *	E:	the election data used to determine the Winner
 */
EMETH MDDA(edata& E  /* approval-count winner among canddts not majority-defeated (or all, if all maj-defeated) */)
{
	bool MDdisquald[MaxNumCands]={false};
	int i,j,r,dqct,thresh,maxc,winner;
	const uint64_t& numberOfCandidates = E.NumCands;
	const CandidateSlate& allCandidates = E.Candidates;
	/*if(CWSPEEDUP && CondorcetWinner >=0 ) return(CondorcetWinner); valid??*/
	Determine(CopeWinOnlyWinner, BuildDefeatsMatrix, E);
	Determine(ApprovalWinner, Approval, E);
	dqct=0;
	thresh = (E.NumVoters)/2;
	for(i=0; i<numberOfCandidates; i++) {
		MDdisquald[i] = false;
		for(j=0; j<numberOfCandidates; j++) {
			if(allCandidates[j].DefeatsMatrix[i] > thresh) {
				MDdisquald[i] = true;
				dqct++;
				break;
			}
		}
	}
	if( dqct >= (int)numberOfCandidates ) { /*If all disqualified, none are. */
		return(ApprovalWinner);
	}
	winner = -1;
	maxc = 0;
	RandomlyPermute( numberOfCandidates, RandCandPerm );
	for(i=(int)numberOfCandidates-1; i>=0; i--) {
		r = RandCandPerm[i];
		if(!MDdisquald[r] && (allCandidates[r].approvals >= maxc)) {
			maxc=allCandidates[r].approvals;
			winner=r;
		}
	}
	return(winner);
}

/***Forest Simmons 27 Feb 2007:
UncAAO stands for Uncovered, Approval, Approval Opposition.
Here's how it works:
For each candidate X, if X is uncovered, then let f(X)=X,
else let f(X) be the candidate against which X has the least approval
opposition, among those candidates that cover X.
["Approval opposition" of X against Y is the number of ballots on which
X but not Y is approved.]
Start with the approval winner A and apply the function f repeatedly
until the output equals the input.  This "fixed point" of f is the
method winner.
  This method requires a tally of both pairwise approval and pairwise
ordinal information, but both are efficiently summable in N by N
matrices, where N is the number of candidates.
  This method (UncAAO) is monotone, clone free, and always picks from the
uncovered set, which is a subset of the Smith set.
  Zero info strategy is sincere.
Even perfect info incentives for burial and betrayal are practically  nil.
   As near as I can tell, the only bad thing about the method is the
"tyranny of the majority" problem shared by most, if not all, deterministic methods.
*****/
/*	UncAAO(E):	returns the index of the "uncovered, approval, approval
 *			opposition" method Winner, described above
 *	E:	the election data used to determine the Winner
 */
EMETH UncAAO(edata& E)
{
	int UncAAOF[MaxNumCands]={0};
	int i,j,ff,r,winner;
	const uint64_t& numberOfCandidates = E.NumCands;
	const CandidateSlate& allCandidates = E.Candidates;
	Determine(ApprovalWinner, Approval, E);
	Determine(RandomUncoveredMemb, UncoveredSet, E);
	RandomlyPermute( numberOfCandidates, RandCandPerm );
	for(i=(int)numberOfCandidates -1; i>=0; i--) {
		if(allCandidates[i].uncovered) {
			UncAAOF[i] = i;
		}else{
			int64_t MnAO = BIGINT;
			ff = -1;
			for(j=(int)numberOfCandidates -1; j>=0; j--) {
				r = RandCandPerm[j];
				if( CoverMatrix[r*numberOfCandidates+i]  ) {
					const int64_t& AppOpp = allCandidates[r].approvals - allCandidates[r].alsoApprovedWith[i];
					if(AppOpp < MnAO) { MnAO = AppOpp; ff = r; }
				}
			}
			assert(ff >= 0);
			UncAAOF[i] = ff;
		}
	}
	winner = ApprovalWinner;
	do{ winner =  UncAAOF[winner]; }until( winner == UncAAOF[winner] );
	return(winner);
}


/*  Woodall-DAC:
1. Rank-order ballots.
2. A voter "acquiesces"' to a set of candidates if he does not rank any
candidate outside the set higher than any inside the set.
(Every voter acquiesces to the full candidate-set.)
3. Sort all possible sets from most acquiescing voters to fewest.
Going down the list, disqualify every candidate not found in each set (i.e. take
set intersections) unless that
would disqualify all remaining candidates (i.e. would result in the empty set).
Continue until only one candidate
is not disqualified; he is the winner.
*****************************************/
/*	WoodallDAC(E):	returns an index corresponding to the Winner according to
 *			Woodall's "descending acquiescing coalitions" method or -1 if an
 *			error occurs
 *	E:	the election data used to determin the Winner
 *
 *	According to "Monotonicity and Single-Seat Election Rules" by Douglas R. Woodall, the
 *	Descending Acquiescing Coalitions (DAC) approach works like so:
 *		1. Voters rank each Candidate in decreasing order of preference.
 *		2. The number of times each possible combination of rank appears amongst all
 *		   votes is counted.
 *		3. Define the "set of eligible Candidates" to be all Candidates.
 *		4. Determine all possible subsets of Candidates without regards to ordering.
 *		5. Define an "acquiescing coalition" of any subset of Candidates as "all Votes
 *		   which do not indicated a preference for any candidate not in that set to any
 *		   candidate in that set"; by definition, all Votes are within the acquiescing
 *		   coalition with respect to all Candidates.
 *		6. The number of occurences of each possible acquiescing coalition is
 *		   determined.
 *		7. The acquiescing coalitions are ordered by decreasing coalition size; by
 *		   definition, the coalition acquiescing to all Candidates is the largest.
 *		8. Begining with the largest coalition, each next largest coalition is used to
 *		   eliminated Candidates from the set of eligible Candidates until only 1
 *		   Candidate remains in the set; if the subset of Candidates associated with a
 *		   given coalition contains no Candidates in the set of eligible Candidates,
 *		   that coalition is ignored; otherwise, all Candidates not in both the set of
 *		   eligible Candidates and the given subset are removed from the set of eligible
 *		   Candidates.
 *
 *	For illustration purposes, consider the following election involing 4 Candidates and 30
 *	Voters, casting the following rankings:
 *
 *			Rankings	Votes
 *			--------	-----
 *			  adcb		  5
 *			  bcad		  5
 *			  cabd		  8
 *			  dabc		  4
 *			  dbca		  8
 *
 *	The set of eligible Candidates are {a, b, c, d}.
 *	The possible subsets are:
 *
 *		{a}
 *		{b}
 *		{c}
 *		{d}
 *		{a, b}
 *		{a, c}
 *		{a, d}
 *		{b, c}
 *		{b, d}
 *		{c, d}
 *		{a, b, c}
 *		{a, b, d}
 *		{a, c, d}
 *		{b, c, d}
 *		{a, b, c, d}
 *
 *	The acquiescing coalition for each subset is:
 *
 *		   Subset	Coalition
 *		------------	---------
 *		{a}		5
 *		{b}		5
 *		{c}		8
 *		{d}		12
 *		{a, b}		0
 *		{a, c}		8
 *		{a, d}		9
 *		{b, c}		5
 *		{b, d}		8
 *		{c, d}		0
 *		{a, b, c}	13
 *		{a, b, d}	4
 *		{a, c, d}	5
 *		{b, c, d}	8
 *		{a, b, c, d}	30
 *
 *	Ordering the acquiescing coalitions by descending size gives:
 *
 *		   Subset	Coalition
 *		------------	---------
 *		{a, b, c, d}	30
 *		{a, b, c}	13
 *		{d}		12
 *		{a, d}		9
 *		{c}		8
 *		{a, c}		8
 *		{b, d}		8
 *		{b, c, d}	8
 *		{a}		5
 *		{b}		5
 *		{b, c}		5
 *		{a, c, d}	5
 *		{a, b, d}	4
 *		{a, b}		0
 *		{c, d}		0
 *
 *	(Woodall's paper does not specify what to do in the event coalitions are the same size.
 *	This may or may not be an issue.)
 *
 *	Beginning with the set of eligible Candidates, {a, b, c, d}, and the largest acquiescing
 *	coalitions, We compare for Candidates to eliminated. The subset of Candidates
 *	corresponding to the first coalition consists of all Candidates, leaving the set of
 *	eligible Candidates unchanged.
 *
 *	The next largest coalition corresponds to subset {a, b, c}. Since Candidate 'd' is not
 *	part of this subset, that Candidate is eliminated from the set of eligible Candidates,
 *	leaving {a, b, c}.
 *
 *	The next largest coalition corresponds to subset {d}. Since no eligible Candidates are
 *	in this subset, it is ignored.
 *
 *	The next largest coalition corresponds to subset {a, d}. Of the 3 remaining eligible
 *	Candidates, 'a', 'b', and 'c', only 'a' appears in the subset, resulting in the
 *	elimination of 'b' and 'c' and electing 'a'.
 */
EMETH WoodallDAC(const edata& E  /*Woodall: Monotonocity of single-seat preferential election rules, Discrete Applied Maths 77 (1997) 81-98.*/)
{
	uint WoodHashCount[3*MaxNumCands*MaxNumVoters]={0};
	uint WoodHashSet[3*MaxNumCands*MaxNumVoters]={0};
	uint WoodSetPerm[3*MaxNumCands*MaxNumVoters];
	bool leave;
	/* Hash Tab entries contain counter and set-code which is a single machine word. */
	int v,c,k,r;
	uint s, x, h, numsets;
	const uint64_t& numberOfCandidates = E.NumCands;
	const uint& numberOfVoters = E.NumVoters;
	const oneVoter (&allVoters)[MaxNumVoters] = E.Voters;
	if( (numberOfCandidates) > (4*sizeof(numberOfCandidates)) ) {
		printf("WoodallDAC: too many candidates %lld to use machine words(%d) to represent sets\n",
			numberOfCandidates,
			(int)(4*sizeof(numberOfCandidates)) );
		printf("You could rewrite the code to use uint128s to try allow up to 128 canddts\n");
		exit(EXIT_FAILURE);
	}
	for(v=ARTINPRIME-1; v>=0; v--) { WoodHashCount[v] = 0; WoodHashSet[v] = 0; }
	for(v=0; v<numberOfVoters; v++) {
		const oneVoter& theVoter = allVoters[v];
			s = 0;
			for(c=0; c < (int)numberOfCandidates; c++) {
				s |= (1U<<(theVoter.topDownPrefs[c]));
				h = s%ARTINPRIME;
				assert( !EmptySet(s) );
				leave = false;
				for(;;) {
					x = WoodHashSet[h];
					if( EmptySet(x) ) { /* insert new set s into hash table */
						WoodHashSet[h] = s;
						WoodHashCount[h] = 1;
						leave = true;
					}else if( x==s ) { /* already there so increment count */
						WoodHashCount[h]++;
						leave = true;
					}
					if (leave) {
						break;
					}
					h++; /* hash table collision so walk ("linear probing") */
				}
			}
	}
	numsets=0;
	for(v=ARTINPRIME-1; v>=0; v--) {
		if( !EmptySet(WoodHashSet[v]) ) {
			WoodSetPerm[numsets] = v;
			numsets++;
		}
	}
	assert(numsets <= (numberOfCandidates * numberOfVoters));
	PermShellSortDown<int>(numsets, (int*)WoodSetPerm, (int*)WoodHashCount);
	s = WoodHashSet[WoodSetPerm[0]];
	assert( !EmptySet(s) );
	if(SingletonSet(s) == false) {
		for(k=1; k<(int)numsets; k++) { /*decreasing order!*/
			h = WoodSetPerm[k];
			x = s & WoodHashSet[h];
			if(!EmptySet(x)) {
				s = x;
				if(SingletonSet(s)) {  break;  }
			}
		}
	}
	/*printf("C%d/%d\n", CardinalitySet(s), numberOfCandidates);*/
	/*It is extremely rare that s is not a singleton set.  In fact may never happen.*/
	RandomlyPermute( numberOfCandidates, RandCandPerm );
	for(k=(int)numberOfCandidates-1; k>=0; k--) {
		r = RandCandPerm[k];
		if( (s>>r)&1U ) {
			return r; /* return random set-element */
		}
	}
	return(-1); /*failure*/
}

/******** uber-routines which package many voting methods into one: **********/

/*	PrintMethName(WhichMeth, Padding):	prints an electoral method name and
 *						potentially some padding
 *	WhichMeth:	an integer indicating the electoral method
 *	Padding:	whether to output some padding
 */
void PrintMethName( int WhichMeth, bool Padding )
{
	const char *name;
	int spaces;
	switch(WhichMeth) {
	case(0) :
		name = "SociallyBest";
		spaces = 3;
		break;
	case(1) :
		name = "SociallyWorst";
		spaces = 2;
		break;
	case(2) :
		name = "RandomWinner";
		spaces = 3;
		break;
	case(3) :
		name = "Plurality";
		spaces = 5;
		break;
	case(4) :
		name = "Borda";
		spaces = 8;
		break;
	case(5) :
		name = "IRV";
		spaces = 10;
		break;
	case(6) :
		name = "Approval";
		spaces = 9;
		break;
	case(7) :
		name = "Range";
		spaces = 8;
		break;
	case(8) :
		name = "SmithSet";
		spaces = 8;
		break;
	case(9) :
		name = "SchwartzSet";
		spaces = 6;
		break;
	case(10) :
		name = "AntiPlurality";
		spaces = 3;
		break;
	case(11) :
		name = "Top2Runoff";
		spaces = 3;
		break;
		/****** above methods were "Core"; below are optional *****/
	case(12) :
		name = "CondorcetLR";
		spaces = 7;
		break;
	case(13) :
		name = "SimpsonKramer";
		spaces = 3;
		break;
	case(14) :
		name = "Bucklin";
		spaces = 7;
		break;
	case(15) :
		name = "Copeland";
		spaces = 6;
		break;
	case(16) :
		name = "SimmonsCond";
		spaces = 3;
		break;
	case(17) :
		name = "SmithIRV";
		spaces = 8;
		break;
	case(18) :
		name = "BTRIRV";
		spaces = 7;
		break;
	case(19) :
		name = "DMC";
		spaces = 10;
		break;
	case(20) :
		name = "Dabagh";
		spaces = 7;
		break;
	case(21) :
		name = "VtForAgainst";
		spaces = 3;
		break;
	case(22) :
		name = "SchulzeBeatpaths";
		spaces = 0;
		break;
	case(23) :
		name = "PlurIR";
		spaces = 9;
		break;
	case(24) :
		name = "Black";
		spaces = 8;
		break;
	case(25) :
		name = "RandomBallot";
		spaces = 1;
		break;
	case(26) :
		name = "RandomPair";
		spaces = 3;
		break;
	case(27) :
		name = "NansonBaldwin";
		spaces = 0;
		break;
	case(28) :
		name = "Nauru";
		spaces = 8;
		break;
	case(29) :
		name = "TopMedianRating";
		spaces = 0;
		break;
	case(30) :
		name = "LoMedianRank";
		spaces = 3;
		break;
	case(31) :
		name = "RaynaudElim";
		spaces = 4;
		break;
	case(32) :
		name = "ArrowRaynaud";
		spaces = 3;
		break;
	case(33) :
		name = "Sinkhorn";
		spaces = 7;
		break;
	case(34) :
		name = "KeenerEig";
		spaces = 7;
		break;
	case(35) :
		name = "MDDA";
		spaces = 12;
		break;
	case(36) :
		name = "VenzkeDisqPlur";
		spaces = 2;
		break;
	case(37) :
		name = "CondorcetApproval";
		spaces = 0;
		break;
	case(38) :
		name = "UncoveredSet";
		spaces = 4;
		break;
	case(39) :
		name = "BramsSanverPrAV";
		spaces = 3;
		break;
	case(40) :
		name = "Coombs";
		spaces = 12;
		break;
	case(41) :
		name = "Top3IRV";
		spaces = 11;
		break;
	case(42) :
		name = "ContinCumul";
		spaces = 7;
		break;
	case(43) :
		name = "IterCopeland";
		spaces = 6;
		break;
	case(44) :
		name = "HeitzigRiver";
		spaces = 6;
		break;
	case(45) :
		name = "MCA";
		spaces = 15;
		break;
	case(46) :
		name = "Range3";
		spaces = 12;
		break;
	case(47) :
		name = "Range10";
		spaces = 11;
		break;
	case(48) :
		name = "HeismanTrophy";
		spaces = 2;
		break;
	case(49) :
		name = "BaseballMVP";
		spaces = 4;
		break;
	case(50) :
		name = "App2Runoff";
		spaces = 5;
		break;
	case(51) :
		name = "Range2Runoff";
		spaces = 3;
		break;
	case(52) :
		name = "HeitzigDFC";
		spaces = 5;
		break;
	case(53) :
		name = "ArmytagePCSchulze";
		spaces = 0;
		break;
	case(54) :
		name = "Hay";
		spaces = 12;
		break;
	case(55) :
		name = "HeitzigLFC";
		spaces = 5;
		break;
	case(56) :
		name = "Benham2AppRunoff";
		spaces = 0;
		break;
	case(57) :
		name = "Benham2AppRunB";
		spaces = 0;
		break;
	case(58) :
		name = "WoodallDAC";
		spaces = 5;
		break;
	case(59) :
		name = "UncAAO";
		spaces = 9;
		break;
		/****** below methods are "Slow": *****/
	case(NumFastMethods+0) :
		name = "TidemanRankedPairs";
		spaces = 0;
		break;
	case(NumFastMethods+1) :
		name = "IRNR2";
		spaces = 10;
		break;
	case(NumFastMethods+2) :
		name = "IRNR1";
		spaces = 10;
		break;
	case(NumFastMethods+3) :
		name = "IRNR3";
		spaces = 10;
		break;
	case(NumFastMethods+4) :
		name = "IRNR9";
		spaces = 10;
		break;
	case(NumFastMethods+5) :
		name = "IRNRlinear";
		spaces = 10;
		break;
	case(NumFastMethods+6) :
		name = "IRNRmin-max";
		spaces = 10;
		break;
	case(NumFastMethods+7) :
		name = "Rouse";
		spaces = 13;
		break;
	default :
		printf("[Unsupported voting method %d]", WhichMeth);
		exit(EXIT_FAILURE);
	}
	printName(name, Padding, spaces);
}

void PrintAvailableVMethods(){
	int i;
	printf("\nAvailable Voting methods:\n");
	for(i=0; i<NumMethods; i++){
		printf("%d=",i); PrintMethName(i,false); printf("\n");
	}
}

voteVector traditionalVoteVectorNormalization(const oneVoter&, const CandidateSlate&, const uint64_t&);
voteVector linearVoteVectorNormalization(const oneVoter&, const CandidateSlate&, const uint64_t&);
voteVector minMaxStyleVoteVectorNormalization(const oneVoter&, const CandidateSlate&, const uint64_t&);

/*	GimmeWinner(E, WhichMeth):	returns the Winner of the election as determined
 *					by a specific method
 *	E:		the election data used to determine the Winner
 *	WhichMeth:	the method to use to determine the Winner
 */
int GimmeWinner( edata& E, int WhichMeth )
{
	int w;
	switch(WhichMeth) {
	case(0) : w=SociallyBest(E); break;
	case(1) : w=SociallyWorst(E); break;
	case(2) : w=RandomWinner(E); break;
	case(3) : w=Plurality(E); break;
	case(4) : w=Borda(E); break;
	case(5) : w=IRV(E); break;
	case(6) : w=Approval(E); break;
	case(7) : w=Range(E); break;
	case(8) : w=SmithSet(E); break;
	case(9) : w=SchwartzSet(E); break;
	case(10) : w=AntiPlurality(E); break;
	case(11) : w=Top2Runoff(E); break;
		/****** above methods were "Core"; below are optional *****/
	case(12) : w=CondorcetLR(E); break;
	case(13) : w=SimpsonKramer(E); break;
	case(14) : w=Bucklin(E); break;
	case(15) : w=Copeland(E); break;
	case(16) : w=SimmonsCond(E); break;
	case(17) : w=SmithIRV(E); break;
	case(18) : w=BTRIRV(E); break;
	case(19) : w=DMC(E); break;
	case(20) : w=Dabagh(E); break;
	case(21) : w=VtForAgainst(E); break;
	case(22) : w=SchulzeBeatpaths(E); break;
	case(23) : w=PlurIR(E); break;
	case(24) : w=Black(E); break;
	case(25) : w=RandomBallot(E); break;
	case(26) : w=RandomPair(E); break;
	case(27) : w=NansonBaldwin(E); break;
	case(28) : w=Nauru(E); break;
	case(29) : w=TopMedianRating(E); break;
	case(30) : w=LoMedianRank(E); break;
	case(31) : w=RaynaudElim(E); break;
	case(32) : w=ArrowRaynaud(E); break;
	case(33) : w=Sinkhorn(E); break;
	case(34) : w=KeenerEig(E); break;
	case(35) : w=MDDA(E); break;
	case(36) : w=VenzkeDisqPlur(E); break;
	case(37) : w=CondorcetApproval(E); break;
	case(38) : w=UncoveredSet(E); break;
	case(39) : w=BramsSanverPrAV(E); break;
	case(40) : w=Coombs(E); break;
	case(41) : w=Top3IRV(E); break;
	case(42) : w=ContinCumul(E); break;
	case(43) : w=IterCopeland(E); break;
	case(44) : w=HeitzigRiver(E); break;
	case(45) : w=MCA(E); break;
	case(46) : RangeGranul=3;  w=RangeN(E); break;
	case(47) : RangeGranul=10; w=RangeN(E); break;
	case(48) : w=HeismanTrophy(E); break;
	case(49) : w=BaseballMVP(E); break;
	case(50) : w=App2Runoff(E); break;
	case(51) : w=Range2Runoff(E); break;
	case(52) : w=HeitzigDFC(E); break;
	case(53) : w=ArmytagePCSchulze(E); break;
	case(54) : w=Hay(E); break;
	case(55) : w=HeitzigLFC(E); break;
	case(56) : w=Benham2AppRunoff(E, false); break;
	case(57) : w=Benham2AppRunoff(E, true); break;
	case(58) : w=WoodallDAC(E); break;
	case(59) : w=UncAAO(E); break;
		/****** below methods are "Slow": *****/
	case(NumFastMethods+0) : w=TidemanRankedPairs(E); break;
	case(NumFastMethods+1) :
			IRNRPOWER=2.0;
			w=IRNR(E, traditionalVoteVectorNormalization);
			break;
	case(NumFastMethods+2) :
			IRNRPOWER=1.0;
			w=IRNR(E, traditionalVoteVectorNormalization);
			break;
	case(NumFastMethods+3) :
			IRNRPOWER=3.0;
			w=IRNR(E, traditionalVoteVectorNormalization);
			break;
	case(NumFastMethods+4) :
			IRNRPOWER=9.0;
			w=IRNR(E, traditionalVoteVectorNormalization);
			break;
	case(NumFastMethods+5) :
			w=IRNR(E, linearVoteVectorNormalization);
			break;
	case(NumFastMethods+6) :
			w=IRNR(E, minMaxStyleVoteVectorNormalization);
			break;
	case(NumFastMethods+7) : w=Rouse(E); break;
	default :
		printf("Unsupported voting method %d\n", WhichMeth);
		exit(EXIT_FAILURE);
	} /*end switch*/
	if((w<0) || (w>=(int)E.NumCands)) {
		printf("Voting method %d=", WhichMeth); PrintMethName(WhichMeth,false);
		printf(" returned erroneous winner %d\n", w);
		exit(EXIT_FAILURE);
	}
	return(w);
}

typedef struct dum2 {
	uint NumVoters;
	uint64_t NumCands;
	uint NumElections;
	real IgnoranceAmplitude;
	real Honfrac;
	std::array<oneVotingMethod, NumMethods> votingMethods;
} brdata;

void EDataPrep(edata& E, const brdata& B);
void PrepareForBayesianRegretOutput(brdata& regretObject, const int &iglevel, bool (&VotMethods)[NumMethods]);
void PrintBROutput(const brdata& regretObject, uint &scenarios);

//	Function: updateMethodAgreements
//
//	updates the agreement counts of each voting method from
//	index 0 to 'm' which result in the same winning Candidate
//	as a given method
//
//	Parameters:
//		theGivenMethod   - the method against which all
//		                   other methods are compared
//		m                - the index of 'theGivenMethod'
//		                   into the array of all voting
//		                   methods
//		allVotingMethods - the array of all voting methods
//		methods          - an array indicating which voting
//		                   methods to compare
//		w                - a value indicating the winning
//		                   Candidate of 'theGivenMethod'
void updateMethodAgreements(oneVotingMethod& theGivenMethod,
                            const int& m,
                            std::array<oneVotingMethod, NumMethods>& allVotingMethods,
                            const bool methods[],
                            const int& w)
{
	for(int j=0; j<m; j++) {
		if(methods[j] || j<NumCoreMethods) {
			oneVotingMethod& methodJ = allVotingMethods[j];
			if( methodJ.Winner == w ) {
				uint& MAgreesWithJ = theGivenMethod.agreementCountWithMethod[j];
				uint& JAgreesWithM = methodJ.agreementCountWithMethod[m];
				MAgreesWithJ++;
				JAgreesWithM++;
				ensure( MAgreesWithJ == JAgreesWithM, 31 );
			}
		}
	}
}

//	Function: updateCondorcetAgreementOfOneMethod
//
//	updates the agreement count of a given voting method with
//	both the Condorcet Winner and the true Condorcet Winner
//	if the current method is actively tested for the current
//	election
//
//	Parameters:
//		theGivenMethod - the voting method to be analyzed
//		tested         - whether that method is actively
//		                 tested for the current election
//		core           - whether the given method is a "core
//		                 method"
void updateCondorcetAgreementOfOneMethod(oneVotingMethod& theGivenMethod,
                                         const bool& tested,
                                         const bool& core)
{
	if(tested || core) {
		const int& theWinner = theGivenMethod.Winner;
		if(theWinner==CondorcetWinner) {
			theGivenMethod.CondorcetAgreementCount++;
		}
		if(theWinner==TrueCW) {
			theGivenMethod.trueCondorcetAgreementCount++;
		}
	}
}

/* all arrays here have NumMethods entries */
//	Function: FindWinnersAndRegrets
//
//	Returns:
//		the Condorcet Winner of an election
//
//	Parameters:
//		E       - the election data to use for determining
//		          the Winner
//		B       - the Bayesian regret data for this election
//		Methods - an array indicating which voting methods
//		          to perform
int FindWinnersAndRegrets( edata& E,  brdata& B,  const bool Methods[] )
{
	std::array<oneVotingMethod, NumMethods>& allVotingMethods = B.votingMethods;
	const CandidateSlate& allCandidates = E.Candidates;
	const int& sociallyBestWinner = allVotingMethods[0].Winner;
	int m,w;
	real r;
	BuildDefeatsMatrix(E);
	InitCoreElState();
	for(m=0; m<NumMethods; m++) {
		if(Methods[m] || m<NumCoreMethods) { /* always run Core Methods */
			oneVotingMethod& methodM = allVotingMethods[m];
			w = GimmeWinner(E, m);
			methodM.Winner = w;
			r = allCandidates[BestWinner].utilitySum - allCandidates[w].utilitySum;
			ensure( BestWinner == sociallyBestWinner, 32 );
			ensure( r>=0.0, 52 );
			assert(BestWinner == sociallyBestWinner); /*can only fail if somebody overwrites array...*/
			assert(r>=0.0); /*can only fail if somebody overwrites array...*/
			WelfordUpdateMeanSD(r, methodM);
			updateMethodAgreements(methodM, m, allVotingMethods, Methods, w);
		}
	}
	if(CondorcetWinner >= 0) {
		for(m=0; m<NumMethods; m++) {
			oneVotingMethod& methodM = allVotingMethods[m];
			const bool& coreMethod = m<NumCoreMethods;
			updateCondorcetAgreementOfOneMethod(methodM, Methods[m], coreMethod);
		}
	}
	return(CondorcetWinner);
}


/*************************** VOTING STRATEGIES: ****
all are subroutines with a common format (embraced by the edata structure):

input:
uint NumVoters;
uint64_t NumCands;

real PerceivedUtility[NumVoters*NumCands];
Entry x*NumCands+y says the utility (a floating point real; greater means better candidate)
of the yth candidate (y=0..NumCands-1) according to voter x, x=0..E->NumVoters-1,
and note these are undistorted by strategies but are distorted by ignorance.

output:
same as the input for VOTING METHODS [other components of edata]
Note, all strategies assume (truthfully?) that the pre-election
polls are a statistical dead heat, i.e. all candidates equally likely to win.
WELL NO: BIASED 1,2,3...
That is done because pre-biased elections are exponentially well-predictable and result in
too little interesting data.
***************************/

//	Function: voteHonestly
//
//	determines the votes of an honest Voter; "voting honestly"
//	is defined as "approving all Candidates above mean utility,
//	approve2 Candidates have utility of at least the mean utility
//	of 'approved' Candidates, ranking all Candidates honestly,
//	and a linear transform of range scores so best == 1 and
//	worst == 0 and All Others are linearly interpolated"
//
//	Parameters:
//		theVoter           - the Voter voting honestly
//		numberOfCandidates - the number of Candidate in
//		                   - the current election
void voteHonestly(oneVoter& theVoter, const uint64_t& numberOfCandidates)
{
	oneCandidateToTheVoter (&allCandidates)[MaxNumCands] = theVoter.Candidates;
	uint (&preferences)[MaxNumCands] = theVoter.topDownPrefs;
	MakeIdentityPerm( numberOfCandidates, preferences );
	PermShellSortDown<real>( numberOfCandidates, (int*)preferences, allCandidates );
	ensure( IsPerm(numberOfCandidates, preferences), 33 );
	real MaxUtil = -HUGE;
	real MinUtil =  HUGE;
	real SumU = 0.0;
	real ThisU;
	for(int i=0; i<numberOfCandidates; i++) {
		allCandidates[preferences[i]].ranking = i;
		ThisU = allCandidates[i].perceivedUtility;
		if(MaxUtil < ThisU) {  MaxUtil = ThisU; }
		if(MinUtil > ThisU) {  MinUtil = ThisU; }
		SumU += ThisU;
	}
	assert(IsPerm(numberOfCandidates, theVoter.Candidates));
	real utilityRange = MaxUtil-MinUtil;
	real RecipDiffUtil;
	if(utilityRange != 0.0) {
		RecipDiffUtil = 1.0 / utilityRange;
	} else {
		RecipDiffUtil = 0.0;
	}
	real MeanU = SumU / numberOfCandidates;
	real Mean2U = 0.0;
	int ACT=0;
	for(int j=0; j<numberOfCandidates; j++) {
		oneCandidateToTheVoter &theCandidate = allCandidates[j];
		ThisU = theCandidate.perceivedUtility;
		theCandidate.score = ( ThisU-MinUtil ) * RecipDiffUtil;
		/* mean-based threshold (with coin toss if exactly at thresh) for approvals */
		if( ThisU > MeanU ) {
			theCandidate.approve = true;
			Mean2U += ThisU;
			ACT++;
		} else if( ThisU < MeanU ) {
			theCandidate.approve = false;
		} else {
			theCandidate.approve = RandBool();
		}
	}
	ensure((ACT!=0), 4);
	Mean2U /= ACT;
	for(int k=0; k<numberOfCandidates; k++) {
		oneCandidateToTheVoter &theCandidate = allCandidates[k];
		ThisU = theCandidate.perceivedUtility;
		if( ThisU >= Mean2U ) {
			theCandidate.approve2 = true;
		} else {
			theCandidate.approve2 = false;
		}
	}
}

//	Function: voteStrategically
//
//	determines the votes of a strategic Voter; "voting
//	strategically" is defined as "presuming the Candidates
//	are pre-ordered in order of decreasing likelihood of
//	winning and the chances decline very rapidly and
//	attempting to maximize One's vote's effect on
//	lower-numbered Candidtes"
//
//	Parameters:
//		theVoter           - the Voter voting honestly
//		numberOfCandidates - the number of Candidate in
//		                   - the current election
void voteStrategically(oneVoter& theVoter, const uint64_t& numberOfCandidates)
{
	real MovingAvg, tmp, Mean2U, ThisU;
	oneCandidateToTheVoter (&allCandidates)[MaxNumCands] = theVoter.Candidates;
	uint (&preferences)[MaxNumCands] = theVoter.topDownPrefs;
	int ACT = 0;
	Mean2U = 0.0;
	MovingAvg = 0.0;
	int nexti = -1;
	uint64_t hibd = numberOfCandidates-1;
	int lobd = 0;
	for(int i=0; i<(int)numberOfCandidates; i++) {
		oneCandidateToTheVoter &theCandidate = allCandidates[i];
		ThisU = theCandidate.perceivedUtility;
		if(i > nexti) {
			nexti++;
			assert(nexti >= 0);
			for( ; nexti<(int)numberOfCandidates; nexti++) {
				tmp = allCandidates[nexti].perceivedUtility;
				MovingAvg += (tmp-MovingAvg)/(nexti+1.0);
				if( fabs(tmp-MovingAvg) > 0.000000000001) {
					break;
				}
			}
		}
		assert(lobd >= 0);
		assert(hibd >= 0);
		assert(lobd < (int)numberOfCandidates);
		assert(hibd < (int)numberOfCandidates);
		if( ThisU > MovingAvg ) {
			theCandidate.approve = true;
			theCandidate.score = 1.0;
			Mean2U += ThisU; ACT++;
			theCandidate.ranking = lobd;
			lobd++;
		}
		else if( ThisU < MovingAvg ) {
			theCandidate.approve = false;
			theCandidate.score = 0.0;
			theCandidate.ranking = hibd;
			hibd--;
		}
		else{
			const bool& rb = RandBool();
			theCandidate.approve = rb;
			theCandidate.score = 0.5;
			if(rb) {
				theCandidate.ranking = lobd;
				lobd++;
			} else {
				theCandidate.ranking = hibd;
				hibd--;
			}
		}
	}
	ensure((ACT!=0), 6);
	Mean2U /= ACT;
	for(int j=0; j<numberOfCandidates; j++) {
		oneCandidateToTheVoter &theCandidate = allCandidates[j];
		const uint64_t &theRanking = theCandidate.ranking;
		ThisU = theCandidate.perceivedUtility;
		if( ThisU >= Mean2U ) {
			theCandidate.approve2 = true;
		} else {
			theCandidate.approve2 = false;
		}
		assert( theRanking < numberOfCandidates );
		preferences[theRanking] = j;
	}
}

//	Function: HonestyStrat
//
//	determines election results based on Voter honesty; the
//	probability of a Voter voting honestly is equal to the value
//	of 'honfrac'; the probability of a Voter voting strategically
//	is equal to the value of '1-honfrac'
//
//	Parameters:
//		E       - the election data to use for determining
//		          the Winner
//		honfrac - probability of a Voter voting honestly
void HonestyStrat( edata& E, real honfrac )
{
	int v;
	oneVoter (&allVoters)[MaxNumVoters] = E.Voters;
	const uint64_t& numberOfCandidates = E.NumCands;
	const uint& numberOfVoters = E.NumVoters;
	assert(numberOfVoters <= MaxNumVoters);
	assert(numberOfCandidates <= MaxNumCands);
	assert(honfrac >= 0.0);
	assert(honfrac <= 1.0);
	for(v=numberOfVoters -1; v>=0; v--) {
		oneVoter& theVoter = allVoters[v];
		if( Rand01() < honfrac ) { /*honest voter*/
			voteHonestly(theVoter, numberOfCandidates);
		} else { /*strategic voter*/
			voteStrategically(theVoter, numberOfCandidates);
		}/*end if(...honfrac) else clause*/
	}/*end for(v)*/
}

/*************************** VOTER IGNORANCE: ***********
input:
uint NumVoters;
uint64_t NumCands;
real IgnoranceAmplitude;

real Utility[NumVoters*NumCands];
Entry x*NumCands+y says the utility (a floating point real; greater means better candidate)
of the yth candidate (y=0..NumCands-1) according to voter x, x=0..NumVoters-1.

output:
real PerceivedUtility[NumVoters*NumCands];
Same thing but adjusted by adding voter ignorance

(everything is embraced inside the edata structure)
***************************/

/*	AddIgnorance(E, IgnoranceAmplitude):	adds a certain degree of ignorance to
 *						Voters
 *	E:			the election data to which ignorance is to be added
 *	IgnoranceAmplitude:	a scaling factor for the degree of ignorance added to
 *				the Voters; a negative value indicates a variable level
 *				of ignorance is to be added depending on the Voter
 *				(stratified); a positive value indicates a constant
 *				level of ignorance across all Voters
 */
void AddIgnorance( edata& E, real IgnoranceAmplitude )
{
	int i;
	real ignorance;
	const uint64_t& numberOfCandidates = E.NumCands;
	const uint& numberOfVoters = E.NumVoters;
	if(IgnoranceAmplitude < 0.0) {
		/* negative is flag to cause VARIABLE levels of ignorance depending on voter (stratified) */
		for(i=numberOfVoters*(int)numberOfCandidates-1; i>=0; i--) {
			uint64_t VoterIndex = i/numberOfCandidates;
			uint64_t CandidateIndex = i % numberOfCandidates;
			ignorance = ((IgnoranceAmplitude * ((2*i)+1)) / numberOfCandidates) * RandNormal();
			E.Voters[VoterIndex].Candidates[CandidateIndex].perceivedUtility = ignorance + E.Voters[VoterIndex].Candidates[CandidateIndex].actualUtility;
		}
	} else {
		/* positive is flag to cause CONSTANT level of ignorance across all voters */
		for(i=numberOfVoters*(int)numberOfCandidates-1; i>=0; i--) {
			uint64_t VoterIndex = i/numberOfCandidates;
			uint64_t CandidateIndex = i % numberOfCandidates;
			ignorance = IgnoranceAmplitude * RandNormal();
			E.Voters[VoterIndex].Candidates[CandidateIndex].perceivedUtility = ignorance + E.Voters[VoterIndex].Candidates[CandidateIndex].actualUtility;
		}
	}
	/* Both positive & negative modes have the same mean ignorance amplitude */
}

/*************************** UTILITY GENERATORS: ***********
input:
uint NumVoters;
uint64_t NumCands;

output:
real Utility[NumVoters*NumCands];
Entry x*NumCands+y says the utility (a floating point real; greater means better candidate)
of the yth candidate (y=0..NumCands-1) according to voter x, x=0..NumVoters-1.
***************************/

real VoterLocation[MaxNumVoters*MaxNumIssues];
real CandLocation[MaxNumCands*MaxNumIssues];

typedef void UTGEN;	/*allows fgrep UTGEN IEVS.c to find out what utility-generators now available*/

void GenNormalLocations( /*input:*/ uint NumVoters, uint64_t NumCands, uint Issues,
			 /*output:*/ real vLocation[], real cLocation[] )
{
	GenRandNormalArr(NumVoters*Issues, vLocation);
	GenRandNormalArr(NumCands*Issues, cLocation);
}

void GenWackyLocations( /*input:*/ uint NumVoters, uint64_t NumCands, uint Issues,
			 /*output:*/ real vLocation[], real cLocation[] )
{
	GenRandWackyArr(NumVoters*Issues, vLocation);
	GenRandNormalArr(NumCands*Issues, cLocation);
}

UTGEN GenNormalUtils( edata& E ){ /* simplest possible utility generator: random normal numbers: */
	const uint64_t& numberOfVoters = E.NumVoters;
	const uint64_t& numberOfCandidates = E.NumCands;
	for(int i=0; i<numberOfVoters; i++) {
		for(int j=0; j<numberOfCandidates; j++) {
			E.Voters[i].Candidates[j].actualUtility = RandNormal();
		}
	}
}

/* if Issues<0 then it uses wacky skew voter distribution instead of normal. Uses Lp distance: */
UTGEN GenIssueDistanceUtils( edata& E, int Issues, real Lp ){  /* utility = distance-based formula in Issue-space */
	uint off2, y, x;
	real KK;
	const uint64_t& numberOfCandidates = E.NumCands;
	const uint& numberOfVoters = E.NumVoters;
	oneVoter (&allVoters)[MaxNumVoters] = E.Voters;
	if(Issues<0){
		Issues = -Issues;
		GenWackyLocations( numberOfVoters, numberOfCandidates, Issues, VoterLocation, CandLocation );
	}else{
		GenNormalLocations( numberOfVoters, numberOfCandidates, Issues, VoterLocation, CandLocation );
	}
	KK = 0.6*Issues;
	for(x=0; x < numberOfVoters; x++){
		oneVoter& theVoter = allVoters[x];
		oneCandidateToTheVoter (&allCandidatesToTheVoter)[MaxNumCands] = theVoter.Candidates;
		off2   = x * Issues;
		for(y=0; y<numberOfCandidates; y++){
			allCandidatesToTheVoter[y].actualUtility = 1.0 / sqrt(KK +
			LpDistanceSquared(Issues, VoterLocation+off2, CandLocation+y*Issues, Lp));
		}
	}
}

UTGEN GenIssueDotprodUtils( edata& E, uint Issues ){  /* utility = canddt*voter vector dot-product in Issue-space */
	uint off2, y, x;
	real s;
	const uint64_t& numberOfCandidates = E.NumCands;
	const uint& numberOfVoters = E.NumVoters;
	oneVoter (&allVoters)[MaxNumVoters] = E.Voters;
	GenNormalLocations( numberOfVoters, numberOfCandidates, Issues, VoterLocation, CandLocation );
	assert(Issues>0);
	s = 1.0/sqrt((real)Issues);
	for(x=0; x < numberOfVoters; x++){
		oneVoter& theVoter = allVoters[x];
		oneCandidateToTheVoter (&allCandidatesToTheVoter)[MaxNumCands] = theVoter.Candidates;
		off2   = x * Issues;
		for(y=0; y < numberOfCandidates; y++){
			allCandidatesToTheVoter[y].actualUtility = s*DotProd(Issues, VoterLocation+off2, CandLocation+y*Issues);
		}
	}
}

const int NumHilFiles = 87;
const int NumDebFiles = 6;
const int MaxNumRanks = 339999;

uint NVotersData[NumHilFiles+NumDebFiles],  NCandsData[NumHilFiles+NumDebFiles];
uint8_t ElData[MaxNumRanks];
int NumElectionsLoaded = 0;

#define VERBOSELOAD 0

//	Function: LoadEldataFiles
//
//	loads real-world election data from a set files; the value
//	returned is the number of elections loaded and processed
int LoadEldataFiles()
{
	static const char*const electionDEBnames[NumDebFiles] = {
		"DB2001.DEB",   "DB2002.DEB",   "DB2003.DEB",   "DB2004.DEB",   "DB2005.DEB",   "DB2006.DEB" };

	static const char*const electionHILnames[NumHilFiles] = {
		"A1.HIL", "A2.HIL", "A3.HIL", "A4.HIL", "A5.HIL", "A6.HIL", "A7.HIL", "A8.HIL", "A9.HIL",
		"A10.HIL", "A11.HIL", "A12.HIL", "A13.HIL", "A14.HIL", "A15.HIL",
		"A16.HIL", "A17.HIL", "A18.HIL", "A19.HIL",
		"A20.HIL", "A21.HIL", "A22.HIL", "A23.HIL", "A24.HIL", "A25.HIL",
		"A26.HIL", "A27.HIL", "A28.HIL", "A29.HIL",
		"A30.HIL", "A31.HIL", "A32.HIL", "A33.HIL", "A34.HIL", "A35.HIL",
		"A48.HIL", "A49.HIL",
		"A50.HIL", "A51.HIL", "A52.HIL", "A53.HIL", "A54.HIL", "A55.HIL",
		"A56.HIL", "A57.HIL", "A58.HIL", "A59.HIL",
		"A60.HIL", "A61.HIL", "A62.HIL", "A63.HIL", "A64.HIL", "A65.HIL",
		"A66.HIL", "A67.HIL", "A68.HIL", "A69.HIL",
		"A70.HIL", "A71.HIL", "A72.HIL", "A73.HIL", "A74.HIL", "A75.HIL",
		"A76.HIL", "A77.HIL", "A78.HIL", "A79.HIL",
		"A80.HIL", "A81.HIL", "A82.HIL", "A83.HIL", "A84.HIL", "A85.HIL",
		"A86.HIL", "A87.HIL", "A88.HIL", "A89.HIL",
		"A90.HIL", "A91.HIL", "A92.HIL", "A93.HIL", "A94.HIL", "A95.HIL",
		"A96.HIL", "A97.HIL", "A98.HIL", "A99.HIL" };
	const int TOOMANYELVOTERS = 7000;
	char c;
	int i,j,v,x,y,elcount,votcount,prefcount,ncands,nvoters,nwinners,itcount;
	FILE *fp;
	elcount=0;  votcount=0;  prefcount=0;  itcount=0;
	printf("Loading %d HIL-format elections files...\n", NumHilFiles);
	for(i=0; i<NumHilFiles; i++) {
		printf("loading %s itcount=%d\n", electionHILnames[i], itcount);
		fp = fopen(electionHILnames[i], "r");
		if(fp==NULL) {
			printf("failure to open file %s for read - terminating\n", electionHILnames[i]);
			printf("Tideman election data files can be got from\n");
			printf("  http://rangevoting.org/TidemanData.html\n");
			exit(EXIT_FAILURE);
		}
		fscanf(fp, "%d %d", &ncands, &nwinners);
		if((ncands<3) || (ncands>MaxNumCands)) {
			printf("bad #candidates %d in %s - terminating\n", ncands, electionHILnames[i]);
			exit(EXIT_FAILURE);
		}
		if((nwinners<1) || (nwinners>=ncands)) {
			printf("bad #winners %d in %s - terminating\n", nwinners, electionHILnames[i]);
			exit(EXIT_FAILURE);
		}
		v = 0;
		do{
			for(j=0; j<ncands; j++) { ElData[itcount+j] = ncands-1; }
			for(y=0;;) {
				x = 0;
				do{ c = (char)getc(fp); }until(isdigit(c));
				do {
					x = (x*10) + (int)(c-'0');
					c = (char)getc(fp);
				} while( isdigit(c) );
#if VERBOSELOAD
				printf("%d ", x);
#endif
				if(x==0) {
					break;
				}
				/*Now do something with x>=1 which is preference y>=0 for voter v>=0*/
				if(x>ncands) {
					printf("bad vote %d in %s - terminating\n", x, electionHILnames[i]);
					exit(EXIT_FAILURE);
				}
				assert(x>0);
				ElData[(itcount+x)-1] = y;
				prefcount++;
				y++;
			}
			itcount += ncands;
			ensure((itcount+ncands) < MaxNumRanks, 54);
#if VERBOSELOAD
			printf("[%d]\n", v);
#endif
			votcount++;
			v++;
		}until(y==0);
		/*Now do something re the election that just ended with v votes*/
		if(v<3) {
			printf("bad #voters %d in %s - terminating\n", v, electionHILnames[i]);
			exit(EXIT_FAILURE);
		}
		NVotersData[elcount] = v;
		NCandsData[elcount] = ncands;
		elcount++;
		fclose(fp);
	}/*end for(i)*/
	printf("Loading %d DEB-format elections files...\n", NumDebFiles);
	for(i=0; i<NumDebFiles; i++) {
		printf("loading %s itcount=%d\n", electionDEBnames[i], itcount);
		fp = fopen(electionDEBnames[i], "r");
		if(fp==NULL) {
			printf("failure to open file %s for read - terminating\n", electionDEBnames[i]);
			printf("Tideman election data files can be got from\n");
			printf("  http://rangevoting.org/TidemanData.html\n");
			exit(EXIT_FAILURE);
		}
		fscanf(fp, "%d %d", &nvoters, &ncands);
		if((nvoters<4) || (nvoters>TOOMANYELVOTERS)) {
			printf("bad #voters %d in %s - terminating\n", nvoters, electionDEBnames[i]);
			exit(EXIT_FAILURE);
		}
		if((ncands<3) || (ncands>MaxNumCands) || (ncands>9)) {
			printf("bad #candidates %d in %s - terminating\n", ncands, electionDEBnames[i]);
			exit(EXIT_FAILURE);
		}
		for(v=0; v<nvoters; v++) {
			for(j=0; j<ncands; j++) { ElData[itcount+j] = ncands-1; }
			for(j=0; j<ncands; j++) {
				do{ c = (char)getc(fp); }while((c=='\n') || (c==' '));
				x = c-'0';  /* watch out for c=='-' */
				/*Now do something with j>=0 which is preference x>0 for voter v>=0*/
#if VERBOSELOAD
				putchar(c);
#endif
				if(c!='-') {
					ElData[itcount+j] = x-1;
				}
				prefcount++;
			}
			itcount += ncands;
			ensure((itcount+ncands) < MaxNumRanks, 53);
#if VERBOSELOAD
			printf("[%d]\n", v);
#endif
			votcount++;
		}
		/* Now do something re the election that just ended with v votes */
		NVotersData[elcount] = v;
		NCandsData[elcount] = ncands;
		elcount++;
		fclose(fp);
	}/*end for(i)*/
	printf("done loading files; loaded %d elections constituting %d prefs and %d votes and %d ranks in all\n",
		elcount, prefcount, votcount, itcount);
	NumElectionsLoaded = elcount;
	assert(NumElectionsLoaded>0);
	return(elcount);
}

//	Function: GenRealWorldUtils
//
//	loads real-world election data from a set files; the value
//	returned is the number of elections loaded and processed
UTGEN GenRealWorldUtils( edata& E ){  /** based on Tideman election dataset **/
	uint y, x,V,C;
	uint64_t VV;
	static int WhichElection=0, offset=0;
	real scalefac;
	uint64_t& numberOfCandidates = E.NumCands;
	uint& numberOfVoters = E.NumVoters;
	oneVoter (&allVoters)[MaxNumVoters] = E.Voters;
	if(WhichElection >= NumElectionsLoaded){
		WhichElection = 0; offset = 0;
	}
	V = NVotersData[WhichElection];
	C = NCandsData[WhichElection];
	assert(C>2);
	assert(V>2);
	numberOfCandidates = C;
	numberOfVoters = 53;  /* always will be 53 voters */
	scalefac = 1.0/sqrt((real)C);
	for(x=0; x < numberOfVoters; x++){
		oneVoter& theVoter = allVoters[x];
		oneCandidateToTheVoter (&allCandidatesToTheVoter)[MaxNumCands] = theVoter.Candidates;
		VV = RandInt(V);  /* choose random voter in the real world election */
		VV *= C;
		for(y=0; y < numberOfCandidates; y++){
			allCandidatesToTheVoter[y].actualUtility = ((numberOfCandidates - (real)ElData[offset+VV+y]) + RandNormal())*scalefac;
		}
	}
	offset += NVotersData[WhichElection]*NCandsData[WhichElection];
	WhichElection++;
}



void UtilDispatcher( edata& E, int WhichMeth ){   /*WhichMeth = -1 ==> real world utils*/
	switch(WhichMeth){
	case(-1) : GenRealWorldUtils(E);  break;
	case(0) : GenNormalUtils(E); break;
	case(1) : GenIssueDotprodUtils(E, 1); break;
	case(2) : GenIssueDotprodUtils(E, 2); break;
	case(3) : GenIssueDotprodUtils(E, 3); break;
	case(4) : GenIssueDotprodUtils(E, 4); break;
	case(5) : GenIssueDotprodUtils(E, 5); break;
	case(6) : GenIssueDistanceUtils(E, 1, 1.0);  break; /* Now using L1 distance */
	case(7) : GenIssueDistanceUtils(E, 2, 1.0);  break;
	case(8) : GenIssueDistanceUtils(E, 3, 1.0);  break;
	case(9) : GenIssueDistanceUtils(E, 4, 1.0);  break;
	case(10) : GenIssueDistanceUtils(E, 5, 1.0);  break;
	case(11)  : GenIssueDistanceUtils(E, -1, 2.0); break; /* These L2 distance */
	case(12)  : GenIssueDistanceUtils(E, -2, 2.0); break;
	case(13)  : GenIssueDistanceUtils(E, -3, 2.0); break;
	case(14)  : GenIssueDistanceUtils(E, -4, 2.0); break;
	case(15) : GenIssueDistanceUtils(E, -5, 2.0); break;

	default : printf("Unsupported util gen %d\n", WhichMeth); exit(EXIT_FAILURE);
	} /*end switch*/
}

/*	PrintUtilName(WhichMeth, Padding):	prints a utility method name and
 *						potentially some padding
 *	WhichMeth:	an integer indicating the electoral method
 *	Padding:	whether to output some padding
 */
void PrintUtilName( int WhichMeth, bool Padding )
{
	const char *name;
	int spaces;
	switch(WhichMeth) {
	case(-1) :
		name = "RealWorld";
		spaces = 9;
		break;
	case(0) :
		name = "RandomNormalUtils";
		spaces = 0;
		break;
	case(1) :
		name = "IssueDotProd[1]";
		spaces = 2;
		break;
	case(2) :
		name = "IssueDotProd[2]";
		spaces = 2;
		break;
	case(3) :
		name = "IssueDotProd[3]";
		spaces = 2;
		break;
	case(4) :
		name = "IssueDotProd[4]";
		spaces = 2;
		break;
	case(5) :
		name = "IssueDotProd[5]";
		spaces = 2;
		break;
	case(6) :
		name = "IssueDistance[1]";
		spaces = 0;
		break;
	case(7) :
		name = "IssueDistance[2]";
		spaces = 0;
		break;
	case(8) :
		name = "IssueDistance[3]";
		spaces = 0;
		break;
	case(9) :
		name = "IssueDistance[4]";
		spaces = 0;
		break;
	case(10) :
		name = "IssueDistance[5]";
		spaces = 0;
		break;
	case(11) :
		name = "IssueDistance[-1]";
		spaces = 1;
		break;
	case(12) :
		name = "IssueDistance[-2]";
		spaces = 1;
		break;
	case(13) :
		name = "IssueDistance[-3]";
		spaces = 1;
		break;
	case(14) :
		name = "IssueDistance[-4]";
		spaces = 1;
		break;
	case(15):
		name = "IssueDistance[-5]";
		spaces = 1;
		break;

	default : printf("UnsupportedUtilGen[%d]\n", WhichMeth); exit(EXIT_FAILURE);
	} /*end switch*/
	printName(name, Padding, spaces);
}

/************ useful IO stuff... ***********/
void PrintConsts()
{
	printf("\nConstants:\n");
	printf("sizeof(uint)=%d bytes\t", (int)sizeof(uint));
	printf("sizeof(uint32_t)=%d\t", (int)sizeof(uint32_t));
	printf("sizeof(uint64_t)=%d\t", (int)sizeof(uint64_t));
	printf("sizeof(real)=%d\n", (int)sizeof(real));
	printf("sizeof(edata)=%d\t", (int)sizeof(edata));
	printf("MaxNumCands=%d\t", MaxNumCands);
	printf("MaxNumVoters=%d\t", MaxNumVoters);
	printf("MaxNumIssues=%d\n", MaxNumIssues);
	printf("NumMethods=%d\t", NumMethods);
	printf("NumCoreMethods=%d\t", NumCoreMethods);
	printf("true=%d\t", (int)true);
	printf("false=%d\n", (int)false);
	ARTINPRIME = FindArtinPrime(MaxNumCands*3*MaxNumVoters);
	printf("ArtinPrime=%d\n", ARTINPRIME);

	printf("BROutputMode=%x\n", BROutputMode);

	ensure(sizeof(uint32_t)==4, 24);
	ensure(sizeof(uint64_t)==8, 25);
}

/************ Bayesian Regret ***********/
void ComputeBRs( brdata& B, const bool VotMethods[], int UtilMeth )
{
	const uint& numberOfElections = B.NumElections;
	uint elnum;
	edata E;

	reset(B.votingMethods);
	InitCoreElState();
	EDataPrep(E, B);
	for(elnum=0; elnum < numberOfElections; elnum++){
		UtilDispatcher(E, UtilMeth);
		AddIgnorance(E, B.IgnoranceAmplitude);
		HonestyStrat(E, B.Honfrac);
		FindWinnersAndRegrets(E, B, VotMethods);
	}
	B.NumVoters =  E.NumVoters;
	B.NumCands = E.NumCands;
	ScaleRegrets(B.votingMethods, 1.0/((B.NumElections - 1.0)*B.NumElections)); /*StdDev/sqrt(#) = StdErr.*/
}

void TestEDataStructs( const brdata& B )
{
	const uint& numberOfElections = B.NumElections;
	uint elnum;
	edata E;
	EDataPrep(E, B);
	for(elnum=0; elnum < numberOfElections; elnum++) {
		printf("GenNormalUtils:\n");
		GenNormalUtils(E);
		printf("AddIgnorance:\n");
		AddIgnorance(E, B.IgnoranceAmplitude);
		printf("HonestyStrat:\n");
		HonestyStrat(E, 1.0);
		printf("BuildDefetasMatrix:\n");
		BuildDefeatsMatrix(E);
		printf("PrintEdata:\n");
		PrintEdata(stdout, E);
	}
}

/*************************** BMP BITMAP GRAPHICS: ***************************/

/**************
Description of MS-Windows bmp format, especially 16 colors uncompressed mode.
All multibyte numbers are little-endian, meaning the
least significant byte (LSB) value is at the lowest address.

File begins with these bytes forming the BITMAPFILEHEADER:
#bytes
2        "BM" = 19778 decimal = 42 4D hex = 66 dec, 77 dec.  Note 19778 = 77*256 + 66.
4        size of file in bytes.  For 4 bits per pixel and 200x200 will be 20118 decimal.
2        zero
2        zero
4        offset from begin of file to bitmap data start in bytes.
               This is 118 decimal for 4 bits per pixel:
              2+4+2+2+4+4+4+4+2+2+4+4+4+4+4+4+64 = 118 bytes.

Then comes the BITMAPINFOHEADER:
4       40 decimal  (in windows) = 28 hex.
4       200 decimal  width of image in pixels (here for 200x200).  At most 32K.
4       200 decimal  height of image in pixels (here for 200x200)  At most 32K.
 The #BYTEs per image row must be a multiple of 4 or else effectively will be due to 0-padding.
 With 200 pixels and 4 bits per pixel that is 100 bytes which is indeed a multiple of 4.
2       1 decimal    #planes
2       4 decimal    bits per pixel, can be 1, 4, 8, or 24.  (maybe also 16 & 32.)
    With 4, the 16 colors are usually the "windows standard 16"
    (and 0 is black and 15 is white) but this is not necessary.
4       0     (but can be nonzero to indicate compression of various kinds; 2 is for
          run-length encoded 4-bit-per-pixel; compression not available in windows before windows 3.)
4       x     x=size of bitmap data in bytes, rounded to a multiple of 4 bytes.
              For 200x200, 16-color, will be 20000 decimal.
              But often left 0 which is legit for uncompressed images.
4       0     preferred horizontal resolutn in pixels/meter, but often left 0.
4       0     preferred vertical resolution in pixels/meter, but often left 0.
4       0     #colors actually used. For 4-bit per pixel will be <=16.  Generally 1<<bitsperpixel.
4       0     #important colors.  Make same as #colors used, or less.
The previous 6 fields often are left 0.  All 6 were added in windows 3 and used to be 0 earlier.

Then comes the RGBQUAD:
It for 16-color images will be a 16*4=64 bytes, each quadruple giving a color.
[For 2-color images it is 2 quads = 8 bytes; for 256-colors it is 256 quads = 1024 bytes.]
The quads are bl, gr, rd, 00 each.

Most commonly this table would be filled with the 16 windows standard colors in order:
color#  rrggbb
0 0x000000   black
1 0x800000   dark red
2 0x008000   dark green
3 0x808000   dark yellow
4 0x000080   dark blue
5 0x800080   dark purple
6 0x008080   dark aqua
7 0xC0C0C0     light grey
248 0x808080   darker grey
249 0xFF0000   red
250 0x00FF00   green
251 0xFFFF00   yellow
252 0x0000FF   blue
253 0xFF00FF   magenta
254 0x00FFFF   aqua
255 0xFFFFFF   white

It seems better, however, to use the best packing of 16 balls in a 3-cube
as the colorset.  The FCC packing of 14 balls would be eight (+-1, +-1, +-1)
plus six (0, 0, +-1).
Th. Gensane: Dense Packings of Equal Spheres in a Cube,
The electronic journal of combinatorics 11,1 (2004) #R33, gives the best known
packings for N balls, N<=32 and the coordinates are at
http://www.randomwalk.de/sequences/a084824.txt .

Then comes BYTE:
this is just the pixels, which for 4-bit-per-pixel will be 2 pixels per byte.
[Apparently the hi byte comes first in the row and rows move from left to right?]
Note each row is padded out if necessary to make it a multiple of 4 bytes long.

If we are using type-1or2 compression, then
each run of N pixels is specified by
1st byte: N
2nd byte: the color's index. In 4-bit-color mode, contains two indices in the two nybbles;
    the even#d pixels (starting with 0) are using the high-order nybble,
    the odd#d  pixels (starting with 1) are using the low-order nybble.
    That permits dithering.
If N=0 that means an escape whose meaning depends on the 2nd byte:
00=end of line,
01=end of file.
**************************************/

void OutputLittleEndianUint32( uint x, FILE *F ){
	putc( x%256, F);
	x /= 256;
	putc( x%256, F);
	x /= 256;
	putc( x%256, F);
	x /= 256;
	putc( x%256, F);
}

void OutputLittleEndianUint16( uint x, FILE *F ){
	putc( x%256, F);
	x /= 256;
	putc( x%256, F);
}

/*	OutputBarray(imgsize, Barray, F):	outputs the byte array, 'Barray[]', to
 *						the file, 'F'
 *	imgsize:	the image size in bytes
 *	Barray:		the byte array
 *	F:		the output file
 */
void OutputBarray( uint imgsize, const uint8_t Barray[], FILE *F )
{
	int i;
	assert(imgsize == 20000);
	for(i=0; i<20000; i++) {
		putc( Barray[i], F );
	}
}

uint ReadPixel( uint x, uint y, const uint8_t Barray[] ){ /*assumes 200x200, 4bits per pixel*/
	int addr;
	uint q;
	addr = 100*y + x/2;
	q = Barray[addr];
	if(!(x%2)){ q >>= 4; }
	return(q&15U);
}

uint OutputCompressedBarray( uint imgsize, const uint8_t Barray[], bool really, FILE *F )
{
	int j,N;
	uint bc, pix, oldpix;
	extern uint shiftForCompressedBMP(const uint8_t [], const int &);
	assert(imgsize <= 20000);
	j=0; N=1; bc=0;
	oldpix = shiftForCompressedBMP(Barray, j);
	for(j=1; j<40000; j++) {
		pix = shiftForCompressedBMP(Barray, j);
		assert(pix<16);
		if( pix == oldpix) {
			N++;
		}else{
			assert(N<256);
			if(really) {
				putc(N, F); /*run length*/
			}
			assert(oldpix<16);
			if(really) {
				putc(oldpix | (oldpix<<4), F); /*color in both nybbles*/
			}
			N=1; bc+=2;
		}
		if(j%200==199) {
			assert(N<256);
			if(really) {
				putc(N, F); /*run length*/
			}
			if(really) {
				putc(pix | (pix<<4), F); /*color in both nybbles*/
			}
			N=1; bc+=4;
			if(really) {
				putc(0, F);
			}
			if(j>=39999) {
				if(really) {
					putc(1, F); /*end of file*/
				}
			} else {
				if(really) {
					putc(0, F); /*end of line*/
				}
				j++;
				pix = shiftForCompressedBMP(Barray, j);
			}
		}
		oldpix=pix;
	} /*end for(j)*/
	return bc;
}

uint OutputBMPHead( uint width, uint height, uint bitsperpixel, uint compression, const uint8_t Barray[], FILE *F ){
	uint sf, roundedwidth, colors, colortabsize, offset, imgsize, compressedsize;
	assert(bitsperpixel==1 || bitsperpixel==4 || bitsperpixel==8 || bitsperpixel==24);
	putc(66, F); /*B*/
	putc(77, F); /*M*/
	roundedwidth = (width/2) + 2*(width%2);
	if(roundedwidth%4 != 0){ roundedwidth = (roundedwidth/4 + 1)*4; }
	colors = 1U<<bitsperpixel;
	colortabsize = 4*colors;
	if(colortabsize > 1024) {
		colortabsize=0;
	}
	offset =54 + colortabsize;
	imgsize = roundedwidth*height;
	if(imgsize==20000 && compression){
		compressedsize = OutputCompressedBarray( imgsize, Barray, false, F );
		if(compressedsize < imgsize){ imgsize=compressedsize; compression=2; }
		else{ compression=0; }
	}
	sf = imgsize + offset;
	OutputLittleEndianUint32(sf, F);
	putc(0,F);  putc(0,F);
	putc(0,F);  putc(0,F);
	OutputLittleEndianUint32(offset, F);
	OutputLittleEndianUint32(40, F);
	OutputLittleEndianUint32(width, F);
	OutputLittleEndianUint32(height, F);
	OutputLittleEndianUint16(1, F); /*1 plane*/
	OutputLittleEndianUint16(bitsperpixel, F);
	OutputLittleEndianUint32(compression, F); /*0=uncompressed*/
	OutputLittleEndianUint32(imgsize, F);
	OutputLittleEndianUint32(0, F);
	OutputLittleEndianUint32(0, F);
	OutputLittleEndianUint32(colors, F);
	OutputLittleEndianUint32(colors, F);
	/*The user will now want to
	 *1. output colortabsize worth of palette data.
	 *2. output imgsize bytes worth of bitmap data.
	 *3. close the file F.
	 *But neither is not done by this routine.*/
	return(imgsize);
}

uint8_t PaletteColorArray[64];

void BogoPutc(uint8_t x, FILE *F){
	static uint i=0;
	if(i>=64) {
		i=0;
	}
	PaletteColorArray[i] = x; i++;
	if(F!=NULL) {
		putc(x, F);
	}
}

void OutputFCC16ColorPalette( FILE *F ){
	BogoPutc(255, F);   BogoPutc(  0, F);   BogoPutc(  0, F);   BogoPutc(0, F); /*red*/
	BogoPutc(  0, F);   BogoPutc(  0, F);   BogoPutc(255, F);   BogoPutc(0, F); /*blue*/
	BogoPutc(  0, F);   BogoPutc(255, F);   BogoPutc(  0, F);   BogoPutc(0, F); /*green*/

	BogoPutc(  0, F);   BogoPutc(255, F);   BogoPutc(255, F);   BogoPutc(0, F); /*bluegreen*/
	BogoPutc(255, F);   BogoPutc(  0, F);   BogoPutc(255, F);   BogoPutc(0, F); /*bluered*/
	BogoPutc(255, F);   BogoPutc(255, F);   BogoPutc(  0, F);   BogoPutc(0, F); /*redgreen*/

	BogoPutc(255, F);   BogoPutc(127, F);   BogoPutc(127, F);   BogoPutc(0, F);
	BogoPutc(127, F);   BogoPutc(255, F);   BogoPutc(177, F);   BogoPutc(0, F); /*added red to try to make yellower*/
	BogoPutc(127, F);   BogoPutc(127, F);   BogoPutc(255, F);   BogoPutc(0, F);

	BogoPutc(  0, F);   BogoPutc(127, F);   BogoPutc(127, F);   BogoPutc(0, F);
	BogoPutc(127, F);   BogoPutc(  0, F);   BogoPutc(127, F);   BogoPutc(0, F);
	BogoPutc(127, F);   BogoPutc(127, F);   BogoPutc(  0, F);   BogoPutc(0, F);

	BogoPutc(127, F);   BogoPutc(127, F);   BogoPutc(127, F);   BogoPutc(0, F); /*add on dark grey*/
	BogoPutc(191, F);   BogoPutc(191, F);   BogoPutc(191, F);   BogoPutc(0, F); /*add on light grey*/

	BogoPutc(255, F);   BogoPutc(255, F);   BogoPutc(255, F);   BogoPutc(0, F); /*white*/
	BogoPutc(  0, F);   BogoPutc(  0, F);   BogoPutc(  0, F);   BogoPutc(0, F); /*black*/
}

void CreatePixel(int x, int y, uint color, uint8_t Barray[])
{ /*Create 16 color pixel at x,y*/
	uint q, addr;
	if(x>=200 || y>=200 || x<0 || y<0) {
		return;  /*auto-clipping to 200x200 region*/
	}
	color &= 15U; /*at most 16 colors, like it or not*/
	addr = 100*y + x/2;
	if(x%2) {  /* use of !(x%2) would be wrong, experimentally verified. */
		q = Barray[addr];
		q &= 0xF0U;
		q |= color;
		Barray[addr] = q;
	}else{
		q = Barray[addr];
		q &= 0x0FU;
		q |= color<<4;
		Barray[addr] = q;
	}
}

void DrawCircle(int x, int y, uint radius, uint BorderColor, uint FillColor, uint8_t Barray[])
{ /*Bresenham alg.*/
	int xd, yd, i, d, mode;
	for(mode=0; mode<=1; mode++){
		xd = 0;
		yd = radius;
		d = 3 - (2 * yd);
		while(yd >= xd){
			if(mode==0){ /*fill:*/
				for(i=(x-xd)+1; i<x+xd; i++){ CreatePixel(i, y + yd, FillColor, Barray); }
				for(i=(x-xd)+1; i<x+xd; i++){ CreatePixel(i, y - yd, FillColor, Barray); }
				for(i=(x-yd)+1; i<x+yd; i++){ CreatePixel(i, y + xd, FillColor, Barray); }
				for(i=(x-yd)+1; i<x+yd; i++){ CreatePixel(i, y - xd, FillColor, Barray); }
			}else{	/*border points in 8 octants:*/
				CreatePixel(x - xd, y + yd, BorderColor, Barray);
				CreatePixel(x + xd, y + yd, BorderColor, Barray);
				CreatePixel(x - xd, y - yd, BorderColor, Barray);
				CreatePixel(x + xd, y - yd, BorderColor, Barray);
				CreatePixel(x - yd, y + xd, BorderColor, Barray);
				CreatePixel(x + yd, y + xd, BorderColor, Barray);
				CreatePixel(x - yd, y - xd, BorderColor, Barray);
				CreatePixel(x + yd, y - xd, BorderColor, Barray);
			}
			/*update coords:*/
			if (d < 0){
				d += 4*xd + 6;
			}else{
				d += 10 + 4*(xd-yd);
				yd--;
			}
			xd++;
		} /*end while*/
	}
}

void DrawVoronoi(uint NumSites, const int xx[], const int yy[], uint8_t Barray[20000], int LpPow)
{
	int x,y,ds,min;
	uint col,i;
	if(NumSites>16) {
		NumSites=16;
	}
	for(x=0; x<200; x++) {
		for(y=0; y<200; y++) {
			min = 9999999; col=BIGINT;
			for(i=0; i<NumSites; i++) {
				if(LpPow==2) {
					ds = (int)(pow((x-xx[i]) + Square(y-yy[i]), 2));
				} else {/*1*/
					ds = (int)(fabs(x-xx[i]) + fabs(y-yy[i]));
				}
				if(ds < min) { min=ds; col=i; }
			}
			assert(col<NumSites);
			CreatePixel(x,y,col,Barray);
		}
	}
	for(i=0; i<NumSites; i++) {
		DrawCircle(xx[i], yy[i], 2, ((i!=15)?15:14), i, Barray);
	}
}

void DrawFPvor(uint NumSites, const int xx[], const int yy[], uint8_t Barray[20000], int LpPow)
{
	uint x,y,ds,mx,col,i;
	if(NumSites>16) {
		NumSites=16;
	}
	for(x=0; x<200; x++) {
		for(y=0; y<200; y++) {
			mx = 0; col=BIGINT;
			for(i=0; i<NumSites; i++) {
				if(LpPow==2) {
					ds = (int)(pow((x-xx[i]) + Square(y-yy[i]), 2));
				} else {/*1*/
					ds = (int)(fabs(x-xx[i]) + fabs(y-yy[i]));
				}
				if(ds >= mx) { mx=ds; col=i; }
			}
			assert(col<NumSites);
			CreatePixel(x,y,col,Barray);
		}
	}
	for(i=0; i<NumSites; i++) {
		DrawCircle(xx[i], yy[i], 2, ((i!=15)?15:14), i, Barray);
	}
}

/***
Do equiangular spacing in [0, 180)  and
on each line place voters exactly centrosymmetically.
Use reals.  RandNormalRadius.
Then add 1 voter at the center.

Autopilot:
increase voters 10% each time until too many voters, or until
some weight exceeds runnerup by factor 5 or more.
CreatePixel for winner with maxweight.
***/
void YeePicture( uint NumSites, int MaxK, const int xx[], const int yy[], int WhichMeth,
		 uint8_t Barray[], uint GaussStdDev, real honfrac, int LpPow )
{
	edata E;
	uint x;
	uint y;
	int k,v,j,ja,i,s,maxw,col,w;
	uint pass;
	int x0,x1,y_0,y_1,PreColor;
	uint p0,p1,p2,p3;
	int weight[16];
	uint RandPerm[16];
	real xt, yt, xto, yto, th, ds, ut, rr;
	bool leave;
	oneVoter (&allVoters)[MaxNumVoters] = E.Voters;
	assert(honfrac >= 0.0);
	assert(honfrac <= 1.0);
	if(NumSites>16) {
		NumSites=16;
	}
	E.NumCands = NumSites;
	MakeIdentityPerm(NumSites, RandPerm);
	for(pass=8; pass>=1; pass /= 2) {
		for(y=0; y<200; y+=pass) {
			printf("[%d]", y);
			if((pass==1 && y%10==9) ||
			   (pass==2 && y%20==18) ||
			   (pass==4 && y%40==36) ||
			   (pass==8 && y%80==72) ) {
				printf("\n");
			}
			for(x=0; x<200; x+=pass) {
				PreColor = -1;
				if( pass>=8 || ((x&pass) || (y&pass)) ) {
					if(x>0 && y>0 && x<200-pass && y<200-pass) {
						if(pass<=4) {
							/*speedup hack: examine previously-computed 4 neighbors*/
							if(x&pass) { x0 = x-pass; x1= x+pass; } else { x0=x-pass; x1=x; }
							if(y&pass) { y_0 = y-pass; y_1= y+pass; } else { y_0=y-pass; y_1=y; }
							p0 = ReadPixel(x0,y_0,Barray);
							p1 = ReadPixel(x0,y_1,Barray);
							p2 = ReadPixel(x1,y_0,Barray);
							p3 = ReadPixel(x1,y_1,Barray);
							if(p0==p1 && p0==p2 && p0==p3) { /*if all agree then use common color*/
								PreColor = p0;
							}
						}/*end if(pass)*/
					}/*end if(x>0 &&..)*/
					ZeroArray( NumSites, weight );
					for(k=10; k<MaxK; k = (k*11)/10) { /*try k-voter election*/
						v = k+k+1;
						E.NumVoters = v;
						j=0; ja=0;
						while(ja<k) { /* generate voter locations and their candidate-utilities */
							th = (ja*PI)/k;
							xt = cos(th);
							yt = sin(th);
							rr = RandRadialNormal()*GaussStdDev;
							xt *= rr; yt *= rr;  /*xt, yt are offsets of voters from central pixel*/
							for(s = -1; s<=1; s++) {
								if(s!=0 || ja==0) { /*centro-symmetric: s=sign*/
									/*printf("k=%d v=%d j=%d ja=%d s=%d x=%f y=%f\n",k,v,j,ja,s,xt*s,yt*s);*/
									oneVoter& theVoter = allVoters[j];
									oneCandidateToTheVoter (&allCandidatesToTheVoter)[MaxNumCands] = theVoter.Candidates;
									xto = xt*s + x;
									yto = yt*s + y;
									for(i=0; i<(int)NumSites; i++) { /*go thru canddts generating utils for voter j*/
										oneCandidateToTheVoter& theCandidateToTheVoter = allCandidatesToTheVoter[i];
										if(LpPow==2) {
											ds = pow((xto-xx[i]) + Square(yto-yy[i]), 2);
										} else {/*1*/
											ds = 0.7*pow(( fabs(xto-xx[i]) + fabs(yto-yy[i]) ), 2);
										}
										ut = 1.0 / sqrt(12000.0 + ds);
										theCandidateToTheVoter.actualUtility = ut;
										theCandidateToTheVoter.perceivedUtility = ut;
									}/*end for(i)*/
									j++;
								}
							}/*end for(s)*/
							ja++;
						}/*end while(ja)*/
						InitCoreElState();
						HonestyStrat(E, honfrac);
						BuildDefeatsMatrix(E);
						SmithIRVwinner = 0;
						w = GimmeWinner( E,  WhichMeth );
						assert(w>=0);
						weight[w] += v;
						maxw = -BIGINT;
						for(i=0; i<(int)NumSites; i++) {
							if(i!=w) {
								if( maxw<weight[i] ) { maxw=weight[i]; }
							}
						}
						leave = false;
						if (PreColor==w && k==10) { /*early break; precolor agree with first election pass*/
							leave = true;
						} else if (weight[w] >= 5*maxw && k>16) { /*early break - good confidence*/
							leave = true;
						} else if (weight[w] - maxw > MaxK * log(((real)MaxK)/k)*10.1) { /*early break - futility*/
							leave = true;
						} else {
							/* do nothing */
						}
						if (leave) {
							break;
						}
					} /*end for(k)*/
					col = ArgMaxUIntArr( NumSites, (uint*)weight, (int*)RandPerm );
					CreatePixel(x,y,col,Barray);
				} /*end if(pass)*/
			} /*end for(x)*/
		} /*end for(y)*/
	} /*end for(pass)*/
	for(i=0; i<(int)NumSites; i++) {
		DrawCircle(xx[i], yy[i], 2, ((i!=15)?15:14), i, Barray);
	}
}

void MakeYeePict( char filename[], const int xx[], const int yy[], int NumSites, int WhichMeth,
		  uint TopYeeVoters, uint GaussStdDev, real honfrac, int LpPow ){
	FILE *F;
	uint8_t Barray[20000];
	int i;
	uint imgsize;
#if defined(MSWINDOWS) && MSWINDOWS
	F = fopen(filename, "wb");
#else
	F = fopen(filename, "w");
#endif
	if(F==NULL){
		printf("failed to open %s\n", filename);
		exit(EXIT_FAILURE);
	}
	if(WhichMeth==0 && LpPow==2) {
		DrawVoronoi( NumSites, (const int*)xx, (const int*)yy, Barray, LpPow );
	} else if(WhichMeth==1 && LpPow==2) {
		DrawFPvor( NumSites, (const int*)xx, (const int*)yy, Barray, LpPow );
	} else {
		YeePicture( NumSites, (TopYeeVoters-1)/2, xx, yy, WhichMeth, Barray, GaussStdDev, honfrac, LpPow );
	}
	imgsize = OutputBMPHead(200, 200, 4, true, Barray, F);
	OutputFCC16ColorPalette(F);
	if(imgsize>=20000) {
		OutputBarray(imgsize, Barray, F);
	} else {
		imgsize=OutputCompressedBarray(imgsize, Barray, true, F);
	}
	printf("%s: %d bytes\n", filename, imgsize);
	fclose(F);
	printf("coordinates:\n");
	for(i=0; i<NumSites; i++){
		printf("(%d,%d)", xx[i], yy[i]);
		if(i<NumSites-1) {
	printf(", ");
		}
		if(i%8==7) {
	printf("\n");
		}
	}
	printf("\n");
}

real ColorContrastScore( uint NumSites, const int xx[], const int yy[] )
{
	int i,j;
	real s;
	assert(NumSites <= 16);
	s = 0;
	for(i=0; i<(int)NumSites; i++) {
		for(j=NumSites-1; j>i; j--) {
			s += (
				abs( PaletteColorArray[i*4]-PaletteColorArray[j*4] ) +
				abs( PaletteColorArray[i*4+1]-PaletteColorArray[j*4+1] ) +
				abs( PaletteColorArray[i*4+2]-PaletteColorArray[j*4+2] )
				) / (abs(xx[i]-xx[j])+abs(yy[j]-yy[j])+3.0);
		}
	}
	assert(s>0.0);
	return s;
}

real ReorderForColorContrast( uint NumSites, int xx[], int yy[] ){
	int i;
	int temp;
	real cs2,cscore;
	uint64_t s1;
	uint64_t s2;
	OutputFCC16ColorPalette( NULL );
	cscore = ColorContrastScore( NumSites, xx, yy );
	for(i=0; i<9000; i++){
		s1 = RandInt(NumSites);
		s2 = RandInt(NumSites);
		if(s1 != s2){
			temp=xx[s1]; xx[s1]=xx[s2]; xx[s2]=temp;
			temp=yy[s1]; yy[s1]=yy[s2]; yy[s2]=temp;
			cs2 = ColorContrastScore( NumSites, xx, yy );
			if(cs2<cscore){
				temp=xx[s1]; xx[s1]=xx[s2]; xx[s2]=temp;
				temp=yy[s1]; yy[s1]=yy[s2]; yy[s2]=temp;

			}else{ cscore = cs2; }
		}
	}
	return cscore;
}

/******************** BR driver and output: *************/

int      honfraclower=0, honfracupper=100;
uint     candnumlower=2, candnumupper=7;
int      votnumlower=2, votnumupper=MaxNumVoters;
int      numelections2try = 59;
int      utilnumlower=0,  utilnumupper = NumUtilGens;
const real HonLevels[] = {1.0, 0.5, 0.0, 0.75, 0.25};
int MethPerm[NumMethods];
real RegretData[MaxScenarios*NumMethods];

struct PopulaceState_t
{
	bool realWorld;
	int numberOfVoters;
	int ignoranceLevel;
	int utilityGeneratorMethod;
};

void PrintTheVotersBayesianRegret(brdata& regretObject, const PopulaceState_t&populaceState, uint &ScenarioCount);

/*In IEVS 2.59 with NumElections=2999 and MaxNumVoters=3000,
 *this driver runs for 80-200 hours
 *on a 2003-era computer, producing several 100 Mbytes output.
 *In IEVS 2.59 and above we include a summarizer so that if
 *you ignore the voluminous output, you still get a nice summary of it at
 *the end:*/
/*	BRDriver():	produces and outputs Bayesian Regret data
 */
void BRDriver()
{
	static const int Pow2Primes[] = {2, 3, 7, 13, 31, 61, 127, 251, 509, 1021, 2039, 4093, 8191, 16381};
	/** Greatest prime <=2^n. **/

	int prind;
	int whichhonlevel;
	int UtilMeth;
	int iglevel;
	uint ScenarioCount=0;
	brdata B;
	PopulaceState_t P;

	P.realWorld = false;
	for(iglevel=0; iglevel<5; iglevel++) {
		P.ignoranceLevel = iglevel;
		for(UtilMeth=0; UtilMeth<NumUtilGens; UtilMeth++) {
			if(UtilMeth>=utilnumlower && UtilMeth<=utilnumupper) {
				P.utilityGeneratorMethod = UtilMeth;
				for(whichhonlevel=0; whichhonlevel<5; whichhonlevel++) {
					B.Honfrac = HonLevels[whichhonlevel];
					if(B.Honfrac*100 < honfracupper + 0.0001 &&
						B.Honfrac*100 > honfraclower - 0.0001 ) {
							for(prind=0; Pow2Primes[prind]<MaxNumVoters; prind++) {
								P.numberOfVoters = Pow2Primes[prind];
								PrintTheVotersBayesianRegret(B, P, ScenarioCount);
							}
					} /*end for(whichhonlevel)*/
				}
			}
		}
	}
	PrintSummaryOfNormalizedRegretData(ScenarioCount);
}


/* Like BRDriver only based on the real world election dataset: */
/*	RWBRDriver():	produces and outputs Bayesian Regret data based upon the real
 *			world election dataset
 */
void RWBRDriver()
{
	int whichhonlevel;
	int iglevel;
	uint ScenarioCount=0;
	brdata B;
	PopulaceState_t P;

	P.realWorld = true;
	for(iglevel=0; iglevel<4; iglevel++) {
		P.ignoranceLevel = iglevel;
		for(whichhonlevel=0; whichhonlevel<5; whichhonlevel++) {
			B.Honfrac = HonLevels[whichhonlevel];
			if(B.Honfrac*100 < honfracupper + 0.0001 && B.Honfrac*100 > honfraclower - 0.0001 ) {
				PrintTheVotersBayesianRegret(B, P, ScenarioCount);
			}
		} /*end for(whichhonlevel)*/
	}/*end for(ignlevel)*/
	PrintSummaryOfNormalizedRegretData(ScenarioCount);
}

/*************************** MAIN CODE: ***************************/

int main(int argc, const char *const *argv)
{
	extern void ErectMainMenu(uint seed);
	uint InitializeRandomNumberGenerator(void);
	void PrintOpeningCredits(void);
	void RunBasicAssertionTests(void);
	uint seed;
	if (argc > 1) {
		if (!strcmp(argv[1], "--test")) {
			extern void runTests(void);

			runTests();
			exit(EXIT_SUCCESS);
		}
	}
	PrintOpeningCredits();
	PrintConsts();
	seed = InitializeRandomNumberGenerator();

	BuildLCMfact();
	RunBasicAssertionTests();

	ErectMainMenu(seed);
	exit(EXIT_SUCCESS);
}



/****Chris Benham's FBC-satisfying 3-slot method:
* http://lists.electorama.com/pipermail/election-methods-electorama.com/2007-January/019160.html
* 1. Voters give each candidate a top rating , a middle rating or no
* rating.
*
* 2. Fix the winning threshold T at 50% of the total valid ballots. Give
* each candidate a score equal to
* the number of ballots on which it is top-rated. If the candidate X
* with the highest score has a score
* equal or greater than  T, elect  X.
*
* 3. If not, eliminate the (remaining) candidate which is given a top or
* middle rating on the fewest ballots, and
* on ballots that now top-rate none of the remaining candidates promote
* all the middle-rated candidates to "top-rated"
* and accordingly amend the scores.
*
* 4. Again, if the now highest scoring candidate X has a score of at
* least T then elect X. (T does not shrink as ballots 'exhaust').
*
* 5. Repeat steps 3 and 4 until there is a winner. If  no candidate ever
* reaches a score of T, elect the candidate
* that is top or middle rated on the most ballots (i.e. the Approval winner).
* NOTES:
* we use top-rank as"top", "approval" as middle, and disapproval as unrated.
*
* Note: later Benham abandoned this method and no longer recommends it; it does not
* actually obey FBC.
*****/

/*	Additions by Me		*/

//	Function: adjustYeeCoordinates
//
//	adjusts (I think) the 'x' and 'y' values created for a Yee
//	picture
//
//	Parameters:
//		numSites   - the number of point sites for the Yee
//		             picture
//		xx         - the set of x-coordinates to adjust
//		yy         - the set of y-coordinates to adjust
//		subsqsideX - the length of one side of the Yee picture
//		subsqsideY - the length of the other side of the
//		             Yee picture
void adjustYeeCoordinates(const int &numSites, int (&xx)[16], int (&yy)[16], const int &subsqsideX, const int &subsqsideY)
{
	int a;
	int b;
	bool advance;
	for(a=0; a<numSites; a += (advance ? 1 : 0)) {
		const auto& xxa = (int)((100 - (subsqsideX*0.5)) + (Rand01()*(subsqsideX-0.001)));
		const auto& yya = (int)((100 - (subsqsideY*0.5)) + (Rand01()*(subsqsideY-0.001)));
		xx[a] = xxa;
		yy[a] = yya;
		advance = true;
		for(b=0; b<a && advance; b++) {
			const auto& xxb = xx[b];
			const auto& yyb = yy[b];
			const auto& xxDifference = xxa-xxb;
			const auto& yyDifference = yya-yyb;
			const auto& absXXDifference = abs(xxDifference);
			const auto& absYYDifference = abs(yyDifference);
			const auto& absDifferenceSum = absXXDifference + absYYDifference;
			const auto& threshold = (((7*subsqsideX)+(7*subsqsideY))/400);
			advance = (absDifferenceSum > threshold);
		}
	}
}

//	Function: ArgMaxArr
//
//	Returns:
//		index of a randomly selected entry from all entries
//		in 'Arr[0..N-1]' with the maximum value
//
//	Parameters:
//		N        - the expected number of elements in 'Arr'
//		           and 'RandPerm'
//		Arr      - array of values to examine
//		RandPerm - array to aid in the random selection
template< class T >
int ArgMaxArr(uint64_t N, const T Arr[], int RandPerm[])
{
	T maxc;
	int a;
	int r;
	int winner;
	winner = -1;
	maxc = (typeid(T)==typeid(int)) ? (T)(-BIGINT) : (T)(-HUGE);
	RandomlyPermute( N, (uint*)RandPerm );
	for(a=0; a<(int)N; a++) {
		r = RandPerm[a];
		if(Arr[r] > maxc) {
			maxc=Arr[r];
			winner=r;
		}
	}
	assert(winner>=0);
	return(winner);
}

/*	ArgMaxArr(N, Candidates[], RandPerm[]):	returns index of random
 *						max perceive utility
 *						entry of Arr[0..N-1]
 *	N:		the expected number of elements in 'RandPerm'
 *	Candidates:	reference to an array of Candidates to provide
 *			information to examine
 *	RandPerm:	array of perm.
 */
template< class T >
int ArgMaxArr(uint64_t N, const oneCandidateToTheVoter (&Candidates)[MaxNumCands], int RandPerm[])
{
	T maxc;
	int a;
	int r;
	int winner;
	winner = -1;
	if(typeid(T)==typeid(int)) {
		maxc = (T)(-BIGINT);
	} else if(typeid(T)==typeid(uint64_t)) {
		maxc = (T) 0;
	} else {
		maxc = (T)(-HUGE);
	}
	RandomlyPermute( N, (uint*)RandPerm );
	for(a=0; a<(int)N; a++) {
		r = RandPerm[a];
		if(Candidates[r].perceivedUtility > maxc) {
			maxc=Candidates[r].perceivedUtility;
			winner=r;
		}
	}
	assert(winner>=0);
	return(winner);
}

//	Function: ArgMinArr
//
//	Returns:
//		the index of a Candidate randomly selected from
//		all Candidates with the minimum value of 'Arr[0..N-1]'
//
//	Parameters:
//		N                - the expected number of elements in
//		                   'Arr' and 'RandPerm'
//		Arr              - array of values to examine
//		RandPerm         - array of randomly ordered indices
//		                   into 'Arr'
template< class T >
int ArgMinArr(uint64_t N, const T Arr[], int RandPerm[])
{
	T minc;
	int a;
	int r;
	int winner;
	winner = -1;
	minc = (typeid(T)==typeid(int)) ? (T)BIGINT : (T)HUGE;
	RandomlyPermute( N, (uint*)RandPerm );
	for(a=0; a<(int)N; a++) {
		r = RandPerm[a];
		if(Arr[r]<minc) {
			minc=Arr[r];
			winner=r;
		}
	}
	assert(winner>=0);
	assert( Arr[winner] <= Arr[0] );
	assert( Arr[winner] <= Arr[N-1] );
	return(winner);
}

//	Function: Maximum
//
//	Returns:
//		the index of a Candidate randomly selected from
//		all Candidates with the maximum value of 'allCandidates[0..N-1].member'
//
//	Parameters:
//		N                - the expected number of elements
//		                   in 'allCandidates'
//		allCandidates    - a slate of Candidates to examine
//		member           - the member of Each Candidate
//		                   to use for comparison
//		permute          - whether to randomly permute 'RandCandPerm'
//		                   before analyzing; default is 'true'
//		checkElimination - check if a Candidate is eliminated
//		                   before testing its 'member';
//		                   default is 'false'
template< class T >
int Maximum(uint64_t N,
	    const CandidateSlate& allCandidates,
	    T oneCandidate::*member,
            const bool& permute,
            const bool& checkElimination)
{
	T maxc;
	int a;
	bool test;
	int r;
	bool typeIsUnsigned = typeid(T) == typeid(uint);
	int winner;
	winner = -1;
	if(typeid(T)==typeid(uint64_t)) {
		maxc = (T) 0;
	} else if(typeid(T)==typeid(int64_t)) {
		maxc = (T)(LLONG_MIN);
	} else {
		maxc = (T)(-HUGE);
	}
	if(permute) {
		RandomlyPermute( N, (uint*)RandCandPerm );
	}
	for(a=0; a<(int)N; a++) {
		r = RandCandPerm[a];
		const oneCandidate& theCandidate = allCandidates[r];
		if(checkElimination) {
			test = !theCandidate.eliminated;
		} else {
			test = true;
		}
		if(test &&
		   ((typeIsUnsigned && allCandidates[r].*member>=maxc) || allCandidates[r].*member>maxc)) {
			maxc=allCandidates[r].*member;
			winner=r;
		}
	}
	assert(winner>=0);
	return(winner);
}

//	Function: Minimum
//
//	Returns:
//		the index of a Candidate randomly selected from
//		all Candidates with the minimum value of 'allCandidates[0..N-1].member'
//
//	Parameters:
//		N                - the expected number of elements
//		                   in 'allCandidates'
//		allCandidates    - a slate of Candidates to examine
//		member           - the member of Each Candidate
//		                   to use for comparison
//		permute          - whether to randomly permute 'RandCandPerm'
//		                   before analyzing; default is 'true'
//		checkElimination - check if a Candidate is eliminated
//		                   before testing its 'member';
//		                   default is 'false'
template< class T >
int Minimum(uint64_t N,
            const CandidateSlate& allCandidates,
            T oneCandidate::*member,
            const bool& permute,
            const bool& checkElimination)
{
	T minc;
	int a;
	int r;
	int winner;
	bool test;
	winner = -1;
	if(typeid(T)==typeid(uint64_t)) {
		minc = (T)MAXUINT64;
	} else if(typeid(T)==typeid(int64_t)) {
		minc = (T)BIGINT;
	} else {
		minc = (T)HUGE;
	}
	if(permute) {
		RandomlyPermute( N, (uint*)RandCandPerm );
	}
	for(a=0; a<(int)N; a++) {
		r = RandCandPerm[a];
		const oneCandidate& theCandidate = allCandidates[r];
		if(checkElimination) {
			test = !theCandidate.eliminated;
		} else {
			test = true;
		}
		if(test && theCandidate.*member<minc) {
			minc=theCandidate.*member;
			winner=r;
		}
	}
	assert(winner>=0);
	if(false==checkElimination) {
		assert( allCandidates[winner].*member <= allCandidates[0].*member );
		assert( allCandidates[winner].*member <= allCandidates[N-1].*member );
	}
	return(winner);
}

/*	calculateForRunoff(E, first, second):	uses 'E' to perform a 'runoff' between 'first' and
 *						'second'; the Winner is returned;
 *						'first' and 'second' are expected to be
 *						different
 *	E:	the election data used to determine the Winner
 *	first:	an index representing one particular Candidate, such as the plurality
 *		winner or the approval Winner for example
 *	second:	an index representing an alternate Candidate with respect to 'first'
 */
int calculateForRunoff(const edata& E, int first, int second)
{
	uint a;
	uint pwct=0;
	uint wct=0;
	const uint numberOfVoters = E.NumVoters;
	real perceivedUtilityOfFirst;
	real perceivedUtilityOfSecond;
	for(a=0; a<numberOfVoters; a++) {
		const oneVoter& theVoter = E.Voters[a];
		const oneCandidateToTheVoter (&allCandidatesToTheVoter)[MaxNumCands] = theVoter.Candidates;
		perceivedUtilityOfFirst = allCandidatesToTheVoter[first].perceivedUtility;
		perceivedUtilityOfSecond = allCandidatesToTheVoter[second].perceivedUtility;
		if( perceivedUtilityOfFirst > perceivedUtilityOfSecond ) {
			pwct++;
		}else if( perceivedUtilityOfFirst < perceivedUtilityOfSecond ) {
			wct++;
		} else {
			/* do nothing */
		}
	}
	if(pwct > wct) {
		return first;
	} else if(pwct < wct) {
		return second;
	} else {
		return flipACoin(first, second);
	}
}

/*	FloydRivestSelect(L, R, K, A[]):	rearranges A[L..R] so A[i]<=A[K]<= A[j]
 *						if L<=i<K<j<=R. Then returns A[K].
 *
 *						If R-L+1=N, then expected runtime for
 *						randomly ordered A[] is 1.5*N compares
 *						asymptotically to find median (and less
 *						for other K, finds max in only 1.0*N).
 *
 *						In particular, if L=0, R=N-1, and
 *						K=(N-1)/2 with N odd, then the function
 *						returns the median. If N is even, then
 *						the function uses K=N/2 to find the
 *						upper bimedian, and then the max of
 *						A[0..K-1] is the lower bimedian.
 *
 *						Note: the magic constants 601, 0.5, 0.5,
 *						20, 5, and 5 in the below should be
 *						tuned to optimize speed. They are
 *						probably improveable.
 *	L:	the 'leftmost' index, assuming indices increase from left to right
 *	R:	the 'rightmost' index, assuming indices increase from left to right
 *	K:	the median index of 'A'
 *	A:	the array of values to rearrange
 *
 *	Note:	the descriptions of 'L', 'R', 'K', and/or 'A' may not be accurate; the
 *		description of what the function does, however, is.
 */
template< class T1 > T1 FloydRivestSelect(uint L, uint R, uint K, T1 A[])
{
	int I, J, N, S, SD, LL, RR;
	real Z;
	T1 T,W;
	while( R > L ) {
		N = (R - L) + 1;
		if( N > 601 ) { /* big enough so worth doing FR-type subsampling to find strong splitter */
			I = (K - L) + 1;
			Z = log((real)(N));
			S = (int)(0.5 * exp((Z * 2.0)/3));
			SD = (int)(0.5 * sqrt((Z * S * (N - S))/N) * Sign<int>((2*I) - N));
			LL = MaxInt(L, (K - ((I * S) / N)) + SD);
			RR = std::min(R, K + (((N - I) * S) / N) + SD);
			/* Recursively select inside small subsample to find an element A[K] usually very
			 * near the desired one: */
			FloydRivestSelect<T1>( LL, RR, K, A );
		}else if( (N > 20) && ((int)(5*(K-L)) > N) && ((int)(5*(R-K)) > N)) {
				/* not big enough to be worth evaluating those expensive logs etc;
				 * but big enough so random splitter is poor; so use median-of-3 */
			I = K-1; if(A[K] < A[I]) { W = A[I]; A[I] = A[K]; A[K] = W; }
			I = K+1; if(A[I] < A[K]) { W = A[I]; A[I] = A[K]; A[K] = W; }
			I = K-1; if(A[K] < A[I]) { W = A[I]; A[I] = A[K]; A[K] = W; }
		} /* otherwise using random splitter (i.e. current value of A[K]) */
		/* now use A[K] to split A[L..R] into two parts... */
		T = A[K]; I = L; J = R;
		W = A[L]; A[L] = A[K]; A[K] = W;
		if( A[R] > T ) { W = A[R]; A[R] = A[L]; A[L] = W; }
		while( I < J ) {
			W = A[I]; A[I] = A[J]; A[J] = W; I++; J--;
			while(A[I] < T) { I++; }
			while(A[J] > T) { J--; }
		}
		if( A[L] == T ) { W = A[L]; A[L] = A[J]; A[J] = W; }
		else{ J++; W = A[J]; A[J] = A[R]; A[R] = W; }
		if( J <= (int)K ) { L = J + 1; }
		if( (int)K <= J ) { R = J - 1; }
		/* Now continue on using contracted [L..R] interval... */
	}
	return( A[K] );
}

/*	runSelfTests():		perfomrs a test of the various
 *				PRNGs and edata structures
 */
void runSelfTests()
{
	brdata B;

	printf("Test of randgen & other tests\n");
	TestsOfRand();
	printf("\nTest edata structure:\n");
	B.NumVoters=6;
	B.NumCands=5;
	B.NumElections=1;
	B.IgnoranceAmplitude=0.001;
	TestEDataStructs(B);
}

/*	runSingleTest(aSeed):	simulates User interation of a single
 *				scenario to help ensure output remains
 *				the same
 *	aSeed:			the seed value for randomization
 */
void runSingleTest(uint aSeed)
{
	extern void runSingleYeeTest(uint aSeed);

	InitRand(aSeed);
	BROutputMode = SORTMODE|ALLMETHS;
	candnumupper=16;//MaxNumCands-1;
	votnumupper=50;//MaxNumVoters;
	numelections2try=1;//99999999;
	utilnumupper=15;
	BRDriver();
	runSingleYeeTest(aSeed);
	runSelfTests();
	LoadEldataFiles();
	RWBRDriver();
}

/*	runTests():	simulates User interaction of multiple
 *			scenarios to help ensure output remains the
 *			same
 */
void runTests()
{
	PrintConsts();
	runSingleTest(1);
	runSingleTest(2);
}

/*	ensure(good, number):	assertion function making sure 'good' is
 *				true before continuing and issuing a
 *				diagnostic if not
 *	good:			the condition to test
 *	number:			the error number to report
 */
void ensure(bool good, int number)
{
	if(good) { /* do nothing */ }
	else {
		printf("TOLCINSE #%d\n", number);
		exit(EXIT_FAILURE);
	}
}

/*	EDataPrep(E, B):	prepares 'E' for testing with data
 *				from 'B'
 *	E:			the edata object which will be
 *				tested
 *	B:			the information for initializing
 *				'E'
 */
void EDataPrep(edata& E, const brdata& B)
{
	const uint& numberOfElections = B.NumElections;
	E.NumVoters = B.NumVoters;
	E.NumCands = B.NumCands;
	if(numberOfElections < 1){
		printf("NumElections=%d<1, error\n", numberOfElections);
		exit(EXIT_FAILURE);
	}
}

//	Function: flipACoin
//
//	Returns:
//		a psuedo-random choice between the two arguments given;
//		the probability of selecting either choice is the same
//		as selecting the other choice
//
//	Parameters:
//		choice1 - one possible choice
//		choice2 - the other possible choice
int flipACoin(int choice1, int choice2)
{
	int selected;
	if(RandBool()) {
		selected = choice1;
	} else {
		selected = choice2;
	}
	return selected;
}

//	Function: swapRandomCandidatesAtIntervals
//
//	swaps all entries in 'RandCandPerm[]' from one index down to
//	another at a certain interval as long as a particular member
//	of the Candidates corresponding to said entires remains below
//	the value of said member of a particular "threshold Candidate"
//
//	Parameters:
//		b          - the index of the "threshold Candidate"
//		a          - the interval at which to swap
//		Candidates - a collection of All Candidates in the current
//		             election
//		member     - the data member of Each Candidate to examine
template<class T>
void swapRandomCandidatesAtIntervals(const int& b,
				     const int& a,
				     const CandidateSlate& Candidates,
				     const T oneCandidate::*member)
{
	const int& x=RandCandPerm[b];
	const T& threshold = Candidates[x].*member;
	int c = b-a;
	for(; (c>=0) && (Candidates[RandCandPerm[c]].*member<threshold); c -= a) {
		RandCandPerm[c+a]=RandCandPerm[c];
	}
	RandCandPerm[c+a]=x;
}

//	Function: swapRandomCandidatesAtIntervals
//
//	swaps all entries in 'Perm[]' from one index down to another
//	at a certain interval as long as the corresponding value in 'Key[]'
//	remains below the corresponding value of a particular "threshold
//	entry"
//
//	Parameters:
//		b    - the index of the "threshold entry"
//		a    - the interval at which to swap
//		Perm - a collection of All Candidates in the current
//		       election
//		Key  - the data member of Each Candidate to examine
template<class T>
void swapRandomCandidatesAtIntervals(const int& b,
				     const int& a,
				     int Perm[],
				     const T Key[])
{
	const int x=Perm[b];
	const T& threshold = Key[x];
	int c = b-a;
	for(; (c>=0) && (Key[Perm[c]]<threshold); c -= a) {
		Perm[c+a]=Perm[c];
	}
	Perm[c+a]=x;
}

//	Function: swapRandomCandidatesAtIntervals
//
//	swaps all entries in 'Perm[]' from one index down to another
//	at a certain interval as long as the corresponding value in 'Candidates[Perm[0..N-1]].perceivedUtility'
//	remains below the corresponding value of a particular "threshold
//	entry"
//
//	Parameters:
//		b          - the index of the "threshold entry"
//		a          - the interval at which to swap
//		Perm       - a collection of All Candidates in the current
//		             election
//		Candidates - a reference to a set of constant Candidates
//		             to guide sorting
void swapRandomCandidatesAtIntervals(const int& b,
				     const int& a,
				     int Perm[],
				     const oneCandidateToTheVoter (&Candidates)[MaxNumCands])
{
	const int x=Perm[b];
	const auto& threshold = Candidates[x].perceivedUtility;
	int c = b-a;
	for(; (c>=0) && (Candidates[Perm[c]].perceivedUtility<threshold); c -= a) {
		Perm[c+a]=Perm[c];
	}
	Perm[c+a]=x;
}

//	Function: PermShellSortDown
//
//	rearranges 'RandCandPerm[0..N-1]' so 'Candidates[RandCandPerm[0..N-1]].member'
//	is in decreasing order
//
//	Parameters:
//		N          - the number of expected entries in 'Candidates'
//		Candidates - a slate of Candidate upon which to base
//		             the sorting
//		member     - the member of Each Candidate to guide sorting
template<class T>
void PermShellSortDown(uint64_t N, const CandidateSlate& Candidates, const T oneCandidate::*member)
{
	int a;
	int b;
	int d;
	for(d=0; (a=ShellIncs[d])>0; d++) {
		for(b=a; b<(int)N; b++) {
			swapRandomCandidatesAtIntervals(b, a, Candidates, member);
		}
	}
	assert(SortedKey<T>(N,Candidates,member)<=0);
}

//	Function: PermShellSortDown
//
//	rearranges 'Perm[0..N-1]' so 'Key[Perm[0..N-1]]' is in decreasing
//	order
//
//	Parameters:
//		N    - the number of expected entries in 'Perm' and 'Key'
//		Perm - a set of values to rearrange
//		Key  - a set of values to guide sorting
template<class T>
void PermShellSortDown(uint64_t N, int Perm[], const T Key[])
{
	int a;
	int b;
	int d;
	for(d=0; (a=ShellIncs[d])>0; d++) {
		for(b=a; b<(int)N; b++) {
			swapRandomCandidatesAtIntervals(b, a, Perm, Key);
		}
	}
	assert(SortedKey<T>(N,Perm,Key)<=0);
}

//	Function: PermShellSortDown
//
//	rearranges 'Perm[0..N-1]' so 'Candidates[Perm[0..N-1]].perceivedUtility'
//	is in decreasing order
//
//	Parameters:
//		N          - the number of expected entries in 'Perm'
//		             and 'Candidates'
//		Perm       - a set of values to rearrange
//		Candidates - a reference to a set of constant Candidates
//		             to guide sorting
template<class T>
void PermShellSortDown(uint64_t N, int Perm[], const oneCandidateToTheVoter (&Candidates)[MaxNumCands])
{
	int a;
	int b;
	int d;
	for(d=0; (a=ShellIncs[d])>0; d++) {
		for(b=a; b<(int)N; b++) {
			swapRandomCandidatesAtIntervals(b, a, Perm, Candidates);
		}
	}
	assert(SortedKey<T>(N,Perm,Candidates)<=0);
}

/*	printName(name, padding, spaces):	prints the given name of a utility
 *						generator or electoral method and
 *						potentially some padding blanks
 *	name:		the name to print
 *	padding:	whether to print padding blanks
 *	spaces:		number of padding blanks to print if 'padding' is 'true'
 */
void printName(const char *name, bool padding, int spaces)
{
	printf("%s", name);
	if(padding) {
		PrintNSpaces(spaces);
	}
}

/*	RandomTest(s, mn, mx, v, ct, func1, func2):	performs a
 *			loop of 100000 times to test if the randgen
 *			calls represented by 'func1' and 'func2'
 *			perform as expected.
 *	s:		the sum of all randomly created values during
 *			the loop
 *	mn:		the minimum value created during the loop
 *	mx:		the maximum value created during the loop
 *	v:		the sum of variances created during the loop
 *	ct:		the counts of value created in blocks of 0.1
 *	func1:		the first function used to create a random value
 *	func2:		the second function used to create a random value
 */
void RandomTest(real &s, real &mn, real &mx, real &v, int (&ct) [10], real (*func1)(void), real (*func2)(void))
{
	int a;
	int y;
	real w;
	real x;

	for(a=0; a<100000; a++){
		x = func1();
		if( func2 != NullFunction ) {
			w = func2();
			x = sqrt((x*x)+(w*w));
		}
		s += x;
		if(mx<x) {
			mx=x;
		}
		if(x<mn) {
			mn=x;
		}
		v += x*x;
		y = (int)(x*10);
		if((x>=0) && (y<10)) {
			ct[y]++;
		}
	}
}

/*	runSingleYeeTest(aSeed):	simulates User interaction to
 *					help ensure Yee output remains
 *					the same
 *	aSeed:				the seed value for
 *					randomization
 */
void runSingleYeeTest(uint aSeed)
{
	const int LpPow=1;
	const int ihonfrac=50;
	const int WhichMeth=2;
	char fname[100];
	const int NumSites=16;
	const int subsqsideX=200;
	const int subsqsideY=200;
	int i;
	int xx[16], yy[16];
	real cscore;
	const int TopYeeVoters=26;
	const int GaussStdDev=200;

	adjustYeeCoordinates(NumSites, xx, yy, subsqsideX, subsqsideY);
	cscore = ReorderForColorContrast(  NumSites, xx, yy );
	printf("Color score %f (big=more constrast); Your coords are:\n", cscore);
	for(i=0; i<NumSites; i++) {
		printf("(%d, %d)\n", xx[i], yy[i]);
	}
	printf("Reordering...\n");
	cscore = ReorderForColorContrast(  NumSites, xx, yy );
	printf("Color score %f (big=more constrast); Your (reordered) coords are:\n", cscore);
	for(i=0; i<NumSites; i++) {
		printf("(%d, %d)\n", xx[i], yy[i]);
	}
	sprintf(fname,"IEVSbmp%u",aSeed);
	MakeYeePict( fname, xx, yy, NumSites, WhichMeth, TopYeeVoters, GaussStdDev, ihonfrac*0.01, LpPow );
	printf("seed=%d\n", aSeed);
}

/*	RandomTestReport(mean_str, meansq_str, s, mn, mx, v, class T):
 *			outputs the results of a call to 'RandomTest()'
 *	mean_str:	a string showing the expected arithmetic mean
 *	meansq_str:     a string showing the expected mean of the squares
 *	s:		the sum of all randomly created values during
 *			the call to 'RandomTest()'
 *	mn:		the minimum value created during the call to
 *			    'RandomTest()'
 *	mx:		the maximum value created during the call to
 *			    'RandomTest()'
 *	v:		the sum of variances created during the call
 *			    to 'RandomTest()'
 *	ct:		the counts of value created in blocks of 0.1
 *			    during the call to 'RandomTest()'
 */
void RandomTestReport(const char *mean_str, const char *meansq_str, real s, real mn, real mx, real v, int (&ct)[10])
{
	int i;
	printf("mean=%g(should be %s)  min=%f  max=%g   meansq=%g(should be %s)\n",
	       s/100000.0,
	       mean_str,
	       mn,
	       mx,
	       v/100000.0,
	       meansq_str);
	printf("counts in 10 bins 0-0.1, 0.1-0.2, etc: ");
	for(i=0; i<10; i++) {
		printf(" %d", ct[i]);
	}
	printf("\n");
}

/*	runoffForApprovalVoting(E):	returns the approval voting runoff Winner
 *
 *	E:				the election data
 */
EMETH runoffForApprovalVoting(const edata& E)
{
	int ASecond;
	ASecond = SecondMaximum(E.NumCands, E.Candidates, &oneCandidate::approvals, ApprovalWinner);
	assert(ASecond>=0);
	return calculateForRunoff(E, ApprovalWinner, ASecond);
}

/*	SelectedRight(L, R, K, A[]):	verifies all entries in 'A[]' from index 'L'
 *					thru 'R' are sorted as expected
 *	L:	the lowest index of entries to examine
 *	R:	the highest index of entries to examine
 *	K:	an index of entries around which the others are to be sorted; a "pivot
 *			index", so to speak
 *	A:	an array of values to examine
 */
template< class T > bool SelectedRight( uint L, uint R, uint K, const T A[] )
{
	uint a;
	for(a=L; a<K; a++) {
		if( A[a]>A[K] ) {
			return false;
		}
	}
	for(a=R; a>K; a--) {
		if( A[K]>A[a] ) {
			return false;
		}
	}
	return true;
}

/*	shiftForCompressedBMP(array, j):	performs an adjustment of an entry in
 *						'array' based upon the value of 'j' and
 *						returns the adjusted value; the value
 *						stored in the array is not modified
 *	array:	an array of values to be adjusted
 *	j:	a 'pseudo-index' into 'array'; 'pseudo' in the sense the function
 *		accesses the 'j/2'-th element of 'array', as oppose to the 'j'-th
 *		element
 */
uint shiftForCompressedBMP(const uint8_t array[], const int &j)
{
	const int shift = (4*(1-(j%2)));
	const int psuedoIndex = j/2;
	const uint8_t rawValue = array[psuedoIndex];
	const uint8_t shiftedValue = rawValue >> shift;
	const uint returnValue = shiftedValue & 15U;
	return returnValue;
}

/*	Sign(x):	returns an indication value to signify whether 'x' is positive,
 *			negative, or 0
 *	x:	the value to compare
 */
template<class T> int Sign(T x)
{
	int rv;
	const T zero = T(0);
	if(x>zero) {
		rv = 1;
	}
	else if(x<zero) {
		rv = -1;
	}
	else {
		rv = 0;
	}
	return rv;
}

/*	Sign(start, finish):	returns an indication value to signify
 *				whether the trend from 'start' to 'finish'
 *				is positive, negative, or 0
 *	start:	the value at which the trend begins
 *	finish:	the value at which the trend ends
 */
template<class T> int Sign(const T& start, const T& finish)
{
	int rv;
	if(start<finish) {
		rv = 1;
	}
	else if(start>finish) {
		rv = -1;
	}
	else {
		rv = 0;
	}
	return rv;
}

/*	SortedKey(N, Candidates, member):	helps to verify 'RandCandPerm'
 *						is sorted in a manner
 *						which is expected; returns
 *						1 if 'Candidates[RandCandPerm[i]].member'
 *						is in increasing order
 *						for all 'i' from 0 to
 *						'N-1', -1 if 'Candidates[RandCandPerm[i]].member'
 *						is in decreasing order,
 *						0 if 'Candidates[RandCandPerm[i]].member'
 *						remains the same, or
 *						2 for any other circumstance
 *	N:		number of elements expected in 'Candidates'
 *	Candidates:	a slate of Candidates used for sorting 'RandCandPerm[]'
 *	member:		a member of Each Candidate used for guiding sorting
 */
template<class T> int SortedKey(uint64_t N, const CandidateSlate& Candidates, const T oneCandidate::*member)
{
	int a;
	int overallTrend;
	T currentValue;
	T nextValue;
	overallTrend = Sign<T>(Candidates[RandCandPerm[0]].*member, Candidates[RandCandPerm[N-1]].*member);
	for(a=1; a<(int)N; a++) {
		currentValue = Candidates[RandCandPerm[a-1]].*member;
		nextValue = Candidates[RandCandPerm[a]].*member;
		switch( overallTrend ) {
			case 1:
				if( nextValue < currentValue ) {
					return 2;
				}
				break;
			case 0:
				if( nextValue != currentValue ) {
					return 2;
				}
				break;
			case -1:
				if( nextValue > currentValue ) {
					return 2;
				}
				break;
			default:
				ensure(false, 28);
		}
	}
	return overallTrend;
}

/*	SortedKey(N, Arr, Key):	helps to verify 'Arr' is sorted in
 *				a manner which is expected;
 *				returns 1 if 'Key[Arr[i]]' is in
 *				increasing order for all 'i' from
 *				0 to 'N-1', -1 if 'Key[Arr[i]]' is
 *				in decreasing order, 0 if
 *				'Key[Arr[i]]' remains the same, or
 *				2 for any other circumstance
 *	N:	number of elements expected in Arr
 *	Arr:	a set of values to have been sorted
 *	Key:	a set of values which are expected to have guided the sorting
 */
template<class T> int SortedKey(uint64_t N, const int Arr[], const T Key[])
{
	int a;
	int overallTrend;
	T currentValue;
	T nextValue;
	overallTrend = Sign<T>( Key[Arr[N-1]]-Key[Arr[0]] );
	for(a=1; a<(int)N; a++) {
		currentValue = Key[Arr[a-1]];
		nextValue = Key[Arr[a]];
		switch( overallTrend ) {
		case 1:
			if( nextValue < currentValue ) {
				return 2;
			}
			break;
		case 0:
			if( nextValue != currentValue ) {
				return 2;
			}
			break;
		case -1:
			if( nextValue > currentValue ) {
				return 2;
			}
			break;
		default:
			ensure(false, 28);
		}
	}
	return overallTrend;
}

/*	SortedKey(N, Arr, Key):	helps to verify 'Arr' is sorted in
 *				a manner which is expected;
 *				returns 1 if 'Key[Arr[i]]' is in
 *				increasing order for all 'i' from
 *				0 to 'N-1', -1 if 'Key[Arr[i]]' is
 *				in decreasing order, 0 if
 *				'Key[Arr[i]]' remains the same, or
 *				2 for any other circumstance
 *	N:		number of elements expected in Arr
 *	Arr:		a set of values to have been sorted
 *	Candidates:	a reference to a set of Candidates with
 *			perceived utility values which are expected
 *			to have guided the sorting
 */
template<class T> int SortedKey(uint64_t N, const int Arr[], const oneCandidateToTheVoter (&Candidates)[MaxNumCands])
{
	int a;
	int overallTrend;
	T currentValue;
	T nextValue;
	overallTrend = Sign<T>( Candidates[Arr[N-1]].perceivedUtility-Candidates[Arr[0]].perceivedUtility );
	for(a=1; a<(int)N; a++) {
		currentValue = Candidates[Arr[a-1]].perceivedUtility;
		nextValue = Candidates[Arr[a]].perceivedUtility;
		switch( overallTrend ) {
			case 1:
				if( nextValue < currentValue ) {
					return 2;
				}
				break;
			case 0:
				if( nextValue != currentValue ) {
					return 2;
				}
				break;
			case -1:
				if( nextValue > currentValue ) {
					return 2;
				}
				break;
			default:
				ensure(false, 29);
		}
	}
	return overallTrend;
}

/*	SortedKey(Arr, Methods):	helps to verify 'Arr' is sorted
 *					in a manner which is expected;
 *					returns 1 if 'Key[Arr[i]]' is
 *					in increasing order for all
 *					'i' from 0 to 'N-1', -1 if
 *					'Key[Arr[i]]' is in decreasing
 *					order, 0 if 'Key[Arr[i]]' remains
 *					the same, or 2 for any other
 *					circumstance
 *	Arr:		a set of values to have been sorted
 *	methods:	a reference to a set of voting methods with
 *			perceived mean regret values which are expected
 *			to have guided the sorting
 */
int SortedKey(const int Arr[], const std::array<oneVotingMethod, NumMethods>& methods)
{
	size_t size = methods.size();
	const range& theRange = range(size, 1);
	int overallTrend = Sign<real>( methods[Arr[size-1]].meanRegret-methods[Arr[0]].meanRegret );
	for(int64_t a : theRange ) {
		const real& currentValue = methods[Arr[a-1]].meanRegret;
		const real& nextValue = methods[Arr[a]].meanRegret;
		switch( overallTrend ) {
			case 1:
				if( nextValue < currentValue ) {
					return 2;
				}
				break;
			case 0:
				if( nextValue != currentValue ) {
					return 2;
				}
				break;
			case -1:
				if( nextValue > currentValue ) {
					return 2;
				}
				break;
			default:
				ensure(false, 29);
		}
	}
	return overallTrend;
}

/*	Test(name, direction, func1, func2, mean_str, meansq_str):
 *				runs a 'random' number generation
 *				test involving 'func1' and 'func2'
 *				to ensure the functions are
 *				performing as expected
 *	name:		the name of the test
 *	direction:	a 'direction' or behavior
 *	func1:		the first function used to create a random value
 *	func2:		the second function used to create a random value
 *	mean_str:	a string showing the expected arithmetic mean
 *	meansq_str:     a string showing the expected mean of the squares
 */
void Test(const char *name, const char *direction, real (*func1)(void), real (*func2)(void), const char *mean_str, const char *meansq_str)
{
	int a;
	int ct[10];
	real s,mx,mn,v;
	s=0.0; v=0.0;
	mn = HUGE;
	mx = -HUGE;
	for(a=0; a<10; a++) {
		ct[a]=0;
	}
	printf("Performing 100000 randgen calls to test that %s behaving ok%s:\n", name, direction);
	RandomTest(s, mn, mx, v, ct, func1, func2);
	RandomTestReport(mean_str, meansq_str, s, mn, mx, v, ct);
}

/*	TwiceMedian(N, A[] ):	returns twice the median of 'A[0..N-1]' or sum of the
 *				bimedians if 'N' is even
 *	N:	the expected number of elements in 'A'
 *	A:	the set of values to examine
 */
template< class T1 >
		T1 TwiceMedian(uint N, T1 A[] )
{
	T1 M,T;
	int b;
	assert(N>0);
	M = FloydRivestSelect<T1>( 0, N-1, N/2, A );
	assert( SelectedRight<T1>( 0, N-1, N/2, A ) );
	if((N%2)==0) { /*N is even*/
		T = A[(N/2) - 1];
		for(b=(N/2) - 2; b>=0; b--) {
			if(A[b] > T) {
				T=A[b];
			}
		}
	}else{ T=M; }
	return T+M;
}

/*	PrintSummaryOfNormalizedRegretData(scenarios):	does what it 'says on the tin'
 *
 *	scenarios:		the number of scenarios
 */
void PrintSummaryOfNormalizedRegretData(uint scenarios)
{
	real BPStrength[NumMethods*NumMethods];
	int CoombCt[NumMethods];
	bool coombElimination[NumMethods];
	int coombrd;
	int i;
	int j;
	int k;
	int r;
	int minc;
	real maxc;
	real RegretSum[NumMethods];
	int VMPerm[NumMethods];
	printf("==================SUMMARY OF NORMALIZED REGRET DATA FROM %d SCENARIOS=============\n",
		scenarios);
	/* regret-sum, Coombs, and Schulze beatpaths used as summarizers
	 * since are good for honest voters and cloneproof. */
	printf("1. voting methods sorted by sum of (normalized so RandWinner=1) regrets (best first):\n");
	for(i=0; i<NumMethods; i++) { RegretSum[i] = 0.0; }
	for(j=0; j<(int)scenarios; j++) {
		r = j*NumMethods;
		for(i=0; i<NumMethods; i++) {
			RegretSum[i] += RegretData[r+i];
		}
	}
	for(i=0; i<NumMethods; i++) { VMPerm[i] = i; MethPerm[i] = i; }
	RealPermShellSortUp( NumMethods, VMPerm, RegretSum );
	for(i=0; i<NumMethods; i++) {
		r = VMPerm[i];
		printf("%d=", r); PrintMethName(r,true);
		printf("\t %g\n", RegretSum[r]);
	}

	printf("\n2. in order of elimination by the Coombs method (worst first):\n");
	for(i=NumMethods -1; i>=0; i--) { coombElimination[i] = false; }
	for(coombrd=NumMethods-2; coombrd>=0; coombrd--) {
		for(i=NumMethods -1; i>=0; i--) { CoombCt[i] = 0; }
		for(j=0; j<(int)scenarios; j++) {
			k = -1;
			r = j*NumMethods;
			maxc = -HUGE;
			for(i=0; i<NumMethods; i++) {
				if(!coombElimination[i]) {
					if(RegretData[r+i] >= maxc) { maxc=RegretData[r+i]; k=i; }
				}
			}
			assert(k>=0);
			ensure(k>=0, 20);
			CoombCt[k]++;
		}
		k = -1; j = -1;
		for(i=0; i<NumMethods; i++) {
			if(CoombCt[i] > k) { k=CoombCt[i]; j=i; }
		}
		assert(j>=0);
		ensure(j>=0, 21);
		assert(!coombElimination[j]);
		coombElimination[j] = true;
		printf("%d=",j); PrintMethName(j,true);
		printf("\n");
	}
	for(i=0; i<NumMethods; i++) {
		if(!coombElimination[i]) {
			printf("Coombs Winner: %d=",i);
			PrintMethName(i,true);
			printf("\n");
			break;
		}
	}

	printf("\n3. voting methods sorted via Schulze beatpaths ordering (best first):\n");
	for(i=NumMethods -1; i>=0; i--) {
		for(j=NumMethods -1; j>=0; j--) {
			BPStrength[i*NumMethods +j]=0;
		}
	}
	for(r=0; r<(int)scenarios; r++) {
		k = r*NumMethods;
		for(i=NumMethods -1; i>=0; i--) {
			for(j=NumMethods -1; j>=0; j--) {
				if(i != j) {
					BPStrength[i*NumMethods +j] += (RegretData[k+i] < RegretData[k+j]) ? 1 : -1;
				}
			}
		}
	}

	for(i=NumMethods -1; i>=0; i--) {
		for(j=NumMethods -1; j>=0; j--) {
			if(i != j) {
				for(k=0; k<NumMethods; k++) {
					if(k != j && k != i) {
						minc = (int)(BPStrength[j*NumMethods+i]);
						if( BPStrength[i*NumMethods +k] < minc ) minc = (int)BPStrength[i*NumMethods +k];
						if( BPStrength[j*NumMethods +k] < minc ) BPStrength[j*NumMethods +k] = minc;
					}
				}
			}
		}
	}

	for(i=0; i<NumMethods; i++) {
		for(j=i+1; j<NumMethods; j++) {
			if( BPStrength[MethPerm[j]*NumMethods +MethPerm[i]] >
					BPStrength[MethPerm[i]*NumMethods +MethPerm[j]] ) {
				/*i is not as good as j, so swap:*/
				r = MethPerm[i]; MethPerm[i] = MethPerm[j]; MethPerm[j] = r;
			}
		}
		printf("%d=",MethPerm[i]); PrintMethName(MethPerm[i], true);
		printf("\n");
	}
	printf("==========end of summary============\n");
}

/*	PrintBROutput(regretObject, scenarios):	produces prints Bayesian regret
 *						information stored in 'regretObject',
 *						incrementing 'scenarios' each time the
 *						function is called
 *	regretObject:	the Bayesian regret data to output
 *	scenarios:	the count of scenarios for which Bayesian regret output has been
 *			created as of the time this function has been called
 */
void PrintBROutput(const brdata& regretObject, uint &scenarios)
{
	int i;
	int j;
	int r;
	real reb;
	real scalefac;
	int TopMeth = 0;
	const bool allMethods = !!(BROutputMode&ALLMETHS);
	const bool top10Methods = !!(BROutputMode&TOP10METHS);
	const bool htmlMode = !!(BROutputMode&HTMLMODE);
	const bool texMode = !!(BROutputMode&TEXMODE);
	const bool normalizeRegrets = !!(BROutputMode&NORMALIZEREGRETS);
	const bool normalizeRegretsByShentrup = !!(BROutputMode&SHENTRUPVSR);
	const bool omitErrorBars = !!(BROutputMode&OMITERRORBARS);
	const real meanRegretBase = regretObject.votingMethods[2].meanRegret;
	real currentMethodsMeanRegret;
	if(allMethods) TopMeth = NumMethods;
	if(top10Methods) TopMeth = 10;
	for(i=0; i<TopMeth; i++) {
		r=i;
		if(BROutputMode&SORTMODE) r=MethPerm[i];
		currentMethodsMeanRegret = regretObject.votingMethods[r].meanRegret;
		if(htmlMode) printf("<tr><td>");
		printf("%d=",r); PrintMethName(r,true);
		if(htmlMode) printf("</td><td>");
		else if(texMode) printf(" & ");
		if(normalizeRegrets) printf(" \t %8.5g", currentMethodsMeanRegret/meanRegretBase);
		else if(normalizeRegretsByShentrup) printf(" \t %8.5g", 100.0*(1.0 - currentMethodsMeanRegret/meanRegretBase));
		else printf(" \t %8.5g", currentMethodsMeanRegret);
		if(!omitErrorBars) {
			const real currentMethodsSRegret = regretObject.votingMethods[r].sRegret;
			if(htmlMode) printf("&plusmn;");
			else if(texMode) printf("\\pm");
			else printf("+-");
			if(normalizeRegrets) reb = sqrt(currentMethodsSRegret)/meanRegretBase;
			else if(normalizeRegretsByShentrup) reb = 100.0*sqrt(currentMethodsSRegret)/meanRegretBase;
			else reb = sqrt(currentMethodsSRegret);
			printf("%5.2g", reb);
		}
		if(htmlMode) printf("</td><td>");
		else if(texMode) printf(" & ");
		if(BROutputMode&VBCONDMODE) printf(" \t  %7d", regretObject.votingMethods[r].CondorcetAgreementCount);
		else printf(" \t  %7d", regretObject.votingMethods[r].trueCondorcetAgreementCount);
		if(htmlMode) printf("</td></tr>\n");
		else if(texMode) printf(" \\\\ \n");
		else printf("\n");
	}/*end for(i)*/
	for(i=0; i<NumMethods; i++) {  /*accumulate regret data for later summary*/
		RegretData[scenarios*NumMethods + i] =
			(regretObject.votingMethods[i].meanRegret + 0.00000000001*Rand01())/( meanRegretBase+0.00000000001 );
	}
	scenarios++;
	if(scenarios > MaxScenarios) {
		printf("ScenarioCount=%d exceeded upper limit; terminating\n", scenarios);
		exit(EXIT_FAILURE);
	}
	if(htmlMode) printf("</td></tr>");
	else if(texMode) printf(" \\\\ ");
	printf("\n");
	if( omitErrorBars && (allMethods || top10Methods) ) {
		const real SRegretBase = regretObject.votingMethods[2].sRegret;
		if(normalizeRegrets) reb = sqrt(SRegretBase)/meanRegretBase;
		else if(normalizeRegretsByShentrup) reb = 100.0*sqrt(SRegretBase)/meanRegretBase;
		else reb = sqrt(SRegretBase);
		printf("ErrorBar for RandomWinner's regret=");
		if(htmlMode) printf("&plusmn;");
		else if(texMode) printf("\\pm");
		else printf("+-");
		printf("%5.2g;\n", reb);
		printf("This (experimentally always?) upperbounds the error bar for every other regret.\n");
	}

	if(BROutputMode&DOAGREETABLES) {
		scalefac = 1.0;
		if(regretObject.NumElections > 999) {
			scalefac = 999.5/regretObject.NumElections;
			printf("\nScaling agreementCountWithMethod values into 0-999.");
		}
		printf("\nAGREE 0 ");
		for(i=1; i<NumMethods; i++) { printf(" %3d ", i); }
		for(i=0; i<NumMethods; i++) {
			printf("\n%2d ", i);
			for(j=0; j<NumMethods; j++) {
				if(j==i) printf("  *  ");
				else printf(" %4d", (int)(0.4999 + regretObject.votingMethods[i].agreementCountWithMethod[j] * scalefac));
			}
		}
		printf("\n");
	}
}

/*	PrepareForBayesianRegretOutput(regretObject, iglevel, VotMethods):	prepares 'regretObject'
 *										and 'VotMethods' for the
 *										outputting of Bayesian regret
 *										information, based in
 *										part on 'iglevel'
 *	regretObject:	the Bayesian regret object to be prepared
 *	iglevel:	the ignorance level
 *	VotMethods:	an array to prepare which will show which voting methods to
 *			perform
 */
void PrepareForBayesianRegretOutput(brdata& regretObject, const int &iglevel, bool (&VotMethods)[NumMethods])
{
	static const real IgnLevels[] = {0.001, 0.01, 0.1, 1.0, -1.0};
	regretObject.NumElections=numelections2try;
	/*1299999=good enough to get all BRs accurate to at least 3 significant digits*/
	/*2999=good enough for usually 2 sig figs, and is 400X faster*/
	regretObject.IgnoranceAmplitude = IgnLevels[iglevel];
	FillArray(NumMethods, VotMethods, true); /*might want to only do a subset... ??*/
	printf("\n");
}

/*	PrintBRPreamble():	prints some 'preambular' text for Bayesian regret output
 */
void PrintBRPreamble()
{
	if(BROutputMode&(ALLMETHS|TOP10METHS)) {
		const bool htmlMode = !!(BROutputMode&HTMLMODE);
		if(htmlMode) {
			printf("<tr><th>Voting Method</th><th>Regrets</th><th>#Agreements with ");
		} else {
			printf("Voting Method & Regrets & #Agreements with ");
		}
		if(BROutputMode&VBCONDMODE) {
			printf("(vote-based) ");
		} else {
			printf("(true-utility-based) ");
		}
		printf("Condorcet Winner (when CW exists)");
		if(htmlMode) {
			printf("</th></tr>");
		}
		printf("\n");
	}
}

/*	PrintTheVotersBayesianRegret(regretObject, populaceState, ScenarioCount):	prints the
 *											Bayesian regret
 *											information for
 *											Voters in either
 *											a 'realWorld'
 *											election of for
 *											a given number
 *											of hypothetical
 *											Voters
 *	regretObject:	the Bayesian regret data to analyze and print
 *	populaceState:	the state of the polucate as a whole, ignorance, number of
 *			Voters, etc.
 *	ScenarioCount:	the number of election scenarios processed
 */
void PrintTheVotersBayesianRegret(brdata& regretObject, const PopulaceState_t &populaceState, uint &ScenarioCount)
{
	bool VotMethods[NumMethods];
	const bool &realWorld = populaceState.realWorld;
	const int &numberOfVoters = populaceState.numberOfVoters;
	const int &iglevel = populaceState.ignoranceLevel;
	const int &UtilMeth = populaceState.utilityGeneratorMethod;
	if(realWorld || (numberOfVoters<=votnumupper && numberOfVoters>=votnumlower)) {
		uint64_t completedLoops;
		uint64_t maximumLoopCount;
		if(not realWorld) {
			regretObject.NumVoters=numberOfVoters;
			maximumLoopCount = (candnumupper - candnumlower) + 1;
		} else {
			maximumLoopCount = 1;
		}
		for(completedLoops = 0; completedLoops < maximumLoopCount; completedLoops++) {
			if(not realWorld) {
				regretObject.NumCands = candnumlower + completedLoops;
			}
			PrepareForBayesianRegretOutput(regretObject, iglevel, VotMethods);
			if(realWorld) {
				MakeIdentityPerm(NumMethods, (uint*)MethPerm);
				ComputeBRs(regretObject, VotMethods, -1);
				RealPermShellSortUp(MethPerm, regretObject.votingMethods);
			}
			printf("(Scenario#%d:", ScenarioCount);
			if(not realWorld) {
				printf(" UtilMeth=");
				PrintUtilName(UtilMeth, false);
			}
			printf(" Honfrac=%.2f, NumVoters=%d, NumCands=%lld, NumElections=%d, IgnoranceAmplitude=%f)\n",
				regretObject.Honfrac, regretObject.NumVoters, regretObject.NumCands,
				regretObject.NumElections, regretObject.IgnoranceAmplitude);
			PrintBRPreamble();
			if(not realWorld) {
				ComputeBRs(regretObject, VotMethods, UtilMeth);
				MakeIdentityPerm(NumMethods, (uint*)MethPerm);
				RealPermShellSortUp(MethPerm, regretObject.votingMethods);
			}
			PrintBROutput(regretObject, ScenarioCount);
		} /*end for(prind)*/
	}
}

/*	InitializeRandomNumberGenerator():	initializes the random number generator
 *						with User input and returns the value
 *						used to seed said generator
 */
uint InitializeRandomNumberGenerator()
{
	uint seed;
	printf("\nPlease enter random seed (0 causes machine to auto-generate from TimeOfDay)\n");
	scanf("%u", &seed);
	InitRand(seed);
	return seed;
}

/*	RunBasicAssertionTests():	does what it says on the tin
 */
void RunBasicAssertionTests()
{
	assert(SingletonSet(8));
	assert(SingletonSet(256));
	assert(!SingletonSet(256+8));
	assert(!SingletonSet(256+512));
	assert(!SingletonSet(3));
	assert(!SingletonSet(7));
	assert(!SingletonSet(5));
	assert(!SingletonSet(10));
	assert(!EmptySet(5));
	assert(EmptySet(0));
}

/*	PrintOpeningCredits():	does what it says on the tin
 */
void PrintOpeningCredits()
{
	const real VERSION = 3.24;
	const int VERSIONYEAR = 2007;
	const int VERSIONMONTH = 2;
	printf("IEVS (Warren D. Smith's infinitely extendible voting system comparator) at your service!\n");
	printf("Version=%f  Year=%d  Month=%d\n", VERSION, VERSIONYEAR, VERSIONMONTH);
}

enum parameters { unsetParameters, defaultParameters, customParameters };

/*	ErectMainMenu(seed):	presents the main menu for the User
 *
 *	seed:	the seed value used for the random number generator
 */
void ErectMainMenu(uint seed)
{
	uint choice;
	extern void PerformBayesianRegretAnalysis(void);
	extern void DisplayMainQuestion(void);
	extern void MakeYeePicture(uint);
	extern void runSelfTests(void);
	DisplayMainQuestion();
	do{
		scanf("%u", &choice);
		switch(choice) {
		case(1) :
			PerformBayesianRegretAnalysis();
			break;
		case(2) :
			MakeYeePicture(seed);
			break;
		case(3):
			runSelfTests();
			break;
		case(4):
			printf("Tally an election with votes you enter (NOT IMPLEMENTED HERE - try\n");
			printf("http://RangeVoting.org/VoteCalc.html)\n");
			break;
		default : printf("Wrong choice %d, moron - try again\n", choice);
			continue;
		}
	} while(false); /* end switch */
}

/*	DisplayMainQuestion(void):	asks the User what They want to do
 */
void DisplayMainQuestion(void)
{
	printf("What do you want to do?\n1=BayesianRegrets\n2=YeePicture\n");
	printf("3=Test RandGen (and other self-tests)\n");
	printf("4=Tally an election with votes you enter\n");
}

/*	DisplayBayesianRegretMenuIntroduction(void):	displays the opening text of the
 *							Bayesian regret analysis menu
 */
void DisplayBayesianRegretMenuIntroduction(void)
{
	printf("Answer a sequence of questions indicating what output format you want for\n");
	printf("the regret tables:\n");
}

/*	RequestBayesianRegretSortingMethod(void):	asks the User to choose a
 *							sorting method for the outputted
 *							results
 */
void RequestBayesianRegretSortingMethod(void)
{
	bool finished;
	uint u;
	printf("I. voting methods (1) sorted-by-regret or (2) by voting-method-number?\n");
	printf("[The latter, while arbitrary, has the advantage of invariance throughout the run.]\n");
	for( finished = false; not finished; ) {
		scanf("%u", &u);
		switch(u) {
		case(1) :
			printf("sorted by regrets.\n");
			BROutputMode |= SORTMODE;
			finished = true;
			break;
		case(2) :
			printf("sorting by voting method number.\n");
			finished = true;
			break;
		default :
			printf("Wrong choice %u, moron - try again\n", u);
			break;
		}
	}
}

/*	RequestOutputFormat():	asks the User to choose an output format
 */
void RequestOutputFormat()
{
	bool finished;
	uint u;
	printf("II. output (1) plain ASCII (2) TeX table formatting (3) HTML table formatting?\n");
	for( finished = false; not finished; ) {
		scanf("%u", &u);
		switch(u) {
		case(1) :
			printf("plain ASCII.\n");
			finished = true;
			break;
		case(2) :
			printf("TeX.\n");
			BROutputMode |= TEXMODE;
			finished = true;
			break;
		case(3) :
			printf("HTML.\n");
			BROutputMode |= HTMLMODE;
			finished = true;
			break;
		default :
			printf("Wrong choice %u, moron - try again\n", u);
			break;
		}
	}
}

/*	RequestBayesianRegretNormalizationMethod():	asks the User to select a method
 *							of 'normalizing' Bayesian regret
 *							values
 */
void RequestBayesianRegretNormalizationMethod(void)
{
	bool finished;
	uint u;
	printf("III. BRs (1) plain (2) normalized so SociallyBest=0, RandomWinner=1\n");
	printf("     (3) normalized so SociallyBest=100, RandomWinner=0, WorseThanRandom<0?\n");
	for( finished = false; not finished; ) {
		scanf("%u", &u);
		switch(u) {
		case(1) :
			printf("plain.\n");
			finished = true;
			break;
		case(2) :
			printf("Best=0, Random=1.\n");
			BROutputMode |= NORMALIZEREGRETS;
			finished = true;
			break;
		case(3) :
			printf("Best=100, Random=0.\n");
			BROutputMode |= SHENTRUPVSR;
			finished = true;
			break;
		default :
			printf("Wrong choice %u, moron - try again\n", u);
			break;
		}
	}
}

/*	RequestErrorBarStatus():	asks the User to specify when to produce error
 *					bars
 */
void RequestErrorBarStatus()
{
	bool finished;
	uint u;
	printf("IV. Error bars (1) on every BR value (2) omit & only compute for RandomWinner\n");
	for( finished = false; not finished; ) {
		scanf("%u", &u);
		switch(u) {
		case(1) :
			printf("all error bars.\n");
			finished = true;
			break;
		case(2) :
			printf("omit error bars.\n");
			BROutputMode |= OMITERRORBARS;
			finished = true;
			break;
		default :
			printf("Wrong choice %u, moron - try again\n", u);
			break;
		}
	}
}

/*	RequestAgreementCountStatus():	asks the User to specify how to print Condorcet
 *					Winner agreement counts
 */
void RequestAgreementCountStatus()
{
	bool finished;
	uint u;
	printf("V. Print Agreement counts with (1) true-utility(undistorted) Condorcet Winners, (2) vote-based CWs\n");
	for( finished = false; not finished; ) {
		scanf("%u", &u);
		switch(u) {
		case(1) :
			printf("true-utility CWs.\n");
			finished = true;
			break;
		case(2) :
			printf("vote-based CWs.\n");
			BROutputMode |= VBCONDMODE;
			finished = true;
			break;
		default :
			printf("Wrong choice %u, moron - try again\n", u);
			break;
		}
	}
}

/*	RequestIntermethodWinnerAgreementCountDisplayStatus():	asks the User if
 *								intermethod
 *								Winner-agreement-count
 *								tables should be
 *								produced
 */
void RequestIntermethodWinnerAgreementCountDisplayStatus(void)
{
	bool finished;
	uint u;
	printf("VI. Print out intermethod winner-agreement-count tables (1) no, (2) yes\n");
	for( finished = false; not finished; ) {
		scanf("%u", &u);
		switch(u) {
		case(1) :
			printf("NO agree-count tables.\n");
			finished = true;
			break;
		case(2) :
			printf("Yes agree-count tables.\n");
			BROutputMode |= DOAGREETABLES;
			finished = true;
			break;
		default :
			printf("Wrong choice %u, moron - try again\n", u);
			break;
		}
	}
}

//	Function: RequestRegretOutput
//
//	asks the User which Bayesian regrets should be outputted
void RequestRegretOutput(void)
{
	bool finished;
	uint u;
	printf("VII. Print out regrets for (1) no, (2) only best 10, (3) all methods\n");
	for( finished = false; not finished; ) {
		scanf("%u", &u);
		switch(u) {
		case(1) :
			printf("No regrets printed (minimum verbosity).\n");
			finished = true;
			break;
		case(2) :
			printf("Top10 methods regrets only printed.\n");
			BROutputMode |= TOP10METHS;
			finished = true;
			break;
		case(3) :
			printf("All regrets printed (maximum verbosity).\n");
			BROutputMode |= ALLMETHS;
			finished = true;
			break;
		default :
			printf("Wrong choice %u, moron - try again\n", u);
			break;
		}
	}
}

//	Function: RequestParameters
//
//	Returns:
//		a value indicating whether the default parameters
//		are to be used in the analysis or a set of custom
//		parameters
parameters RequestParameters(void)
{
	uint u;
	parameters rv;
	printf("VIII. (1) All parameter knob-settings, or (2) restricted ranges?\n");
	honfraclower=0;
	honfracupper=100;
	candnumlower=2;
	candnumupper=7;
	votnumlower=2;
	votnumupper=MaxNumVoters;
	numelections2try = 59;
	for( rv = unsetParameters; unsetParameters == rv; ) {
		scanf("%u", &u);
		switch(u) {
		case(1) :
			printf("All settings.\n");
			rv = defaultParameters;
			break;
		case(2) :
			printf("Restricted Ranges...\n");
			printf("Honesty fraction range - default is 0 100:\n");
			scanf("%d %d", &honfraclower, &honfracupper);
			honfraclower = std::max(honfraclower,0);
			honfracupper = std::min(honfracupper,100);
			printf("Honesty fraction range [%d, %d] chosen.\n", honfraclower, honfracupper);
			printf("Candidate Number range - default is 2 7 [but this range ignored if real-world dataset]:\n");
			scanf("%d %d", &candnumlower, &candnumupper);
			candnumlower = std::max(candnumlower,2U);
			candnumupper = std::max(candnumupper,(unsigned)(MaxNumCands-1));
			printf("Candidate number range [%d, %d] chosen.\n", candnumlower, candnumupper);
			printf("Voter Number range - default is 2 %d [but this range ignored if real-world dataset:\n",
				votnumupper);
			scanf("%d %d", &votnumlower, &votnumupper);
			votnumlower = std::max(votnumlower,0);
			votnumupper = std::min(votnumupper,MaxNumVoters);
			printf("Voter number range [%d, %d] chosen.\n", votnumlower, votnumupper);
			printf("Number of elections to try per scenario - default is %d\n", numelections2try);
			scanf("%d", &numelections2try);
			numelections2try = std::max(numelections2try,29);
			numelections2try = std::min(numelections2try,99999999);
			printf("Trying %d elections per scenario.\n", numelections2try);
			rv = customParameters;
			break;
		default :
			printf("Wrong choice %u, moron - try again\n", u);
			break;
		}
	}
	ensure(rv not_eq unsetParameters, 22);
	return rv;
}

//	Function: requestCustomUtilityGenerators
//
//	requests the specific utility generators to use if
//	custom parameters have been requested
//
//	Parameters:
//		parameterType - the type of the various parameters
//		                requested by the User
void requestCustomUtilityGenerators(parameters parameterType)
{
	if(parameterType==customParameters) {
		printf("Select which utility-generators you want (default 0 thru 15):\n");
		for(int i=0; i<16; i++) {
			printf("%2d: ", i);
			PrintUtilName(i,true);
			printf("\n");
		}
		scanf("%d %d", &utilnumlower, &utilnumupper);
		utilnumlower = std::max(utilnumlower, 0);
		utilnumupper = std::min(utilnumupper, 15);
		printf("Utility gens t com [%d, %d] chosen.\n", utilnumlower, utilnumupper);
		/**** if ???
		printf("Select LPpow???d):\n");
		scanf("%d", &LPpow);
		if(LPpow<1) LPpow=1;
		if(LPpow>5) LPpow=5;
		printf("Using L%d distances.\n", LPpow);
		*****/
	}
}

//	Function: RequestUtilityGenerators
//
//	Returns:
//		a a pointer to the relevant utility driver
//
//	Parameters:
//		parameterType - the type of the various parameters
//		                requested by the User
driver_t RequestUtilityGenerators(parameters parameterType)
{
	bool finished;
	driver_t driver = NULL;
	uint u;
	printf("IX. (1) Machine or (2) Real-world-based utilities?\n");
	for(finished = false; not finished; ) {
		scanf("%u", &u);
		switch(u) {
		case(1) :
			printf("Machine.\n");
			requestCustomUtilityGenerators(parameterType);
			driver = BRDriver;
			finished = true;
			break;
		case(2) :
			printf("Real-world-based.\n");
			LoadEldataFiles();
			driver = RWBRDriver;
			finished = true;
			break;
		default :
			printf("Wrong choice %u, moron - try again\n", u);
			break;
		}
	}
	ensure( driver not_eq NULL, 23);
	return driver;
}

/*	RequestVotingMethod(void):	asks the User to select a voting method for
 *					creating a Yee Picture; the procedure returns a
 *					value representing the chosen voting method
 */
int RequestVotingMethod(void)
{
	int chosenMethod;
	printf("Which voting method? Your choices:");
	PrintAvailableVMethods();
	scanf("%d", &chosenMethod);
	printf("using %d=", chosenMethod);
	PrintMethName(chosenMethod, false);
	printf(".\n");
	return chosenMethod;
}

/*	RequestNameForBMP(name):	asks the User for a name of a file into which to
 *					output the Yee picture; the procedure
 *					automatically adds ".bmp" to the end of the file
 *					name
 *
 *	name:	a reference to a character array which receives the file name
 */
void RequestNameForBMP(char (&name)[100])
{
	size_t i;
	printf("What filename [.bmp suffix will be auto-added for you]?\n");
	scanf("%s", name);
	i = strlen(name);
	if(i>30) {
		printf("filename too long, moron\n");
		exit(EXIT_FAILURE);
	}
	strcat(name, ".bmp");
}

/*	RequestPointSiteCount(void):	asks the User for a number of 'point-sites' for
 *					the Yee picture; the procedure returns the
 *					number entered by the User as long as it is at
 *					least 1 and at most 16; otherwise, the procedure
 *					returns 16.
 */
int RequestPointSiteCount(void)
{
	int numberOfSites;
	printf("how many point-sites do you want [1 to 16]?\n");
	scanf("%d", &numberOfSites);
	if((numberOfSites<1) || (numberOfSites>16)) {
		printf("out of bounds value %d moron, using 16 instead\n",numberOfSites);
		numberOfSites=16;
	}
	return numberOfSites;
}

/*	RequestCoordPairs(numberOfSites, xarray, yarray):	asks the User if the
 *								coord-pairs should be
 *								entered manually or
 *								created for the User
 *
 *	numberOfSites:	the number of coord-pairs to enter/create
 *	xarray:		a reference to an array to receive x coordinates
 *	yarray:		a reference to an array to receive y coordinates
 */
void RequestCoordPairs(int numberOfSites, int (&xarray)[16], int (&yarray)[16])
{
	int i;
	uint u;
	bool finished;
	printf("Do you want to:\n1. enter the %d coord-pairs yourself;\n",numberOfSites);
	printf("2. random coordinate auto-generation?\n");
	for(finished = false; not finished; ) {
		scanf("%u", &u);
		switch(u) {
		case(1) :
			printf("Enter coord pairs X Y with space (not comma) between X & Y, newline between pairs\n");
			printf("(0,0) is in the lower left.  For example an equilateral triangle would be\n");
			printf("99 197\n186 47\n12 47\nCoords outside of the [[0,199] range are permitted.\nYour coords:\n");
			for(i=0; i<numberOfSites; i++) {
				scanf("%d %d", &(xarray[i]), &(yarray[i]));
			}
			printf("Your coords are:\n");
			for(i=0; i<numberOfSites; i++) {
				printf("(%d, %d)\n", xarray[i], yarray[i]);
			}
			finished = true;
			break;
		case(2) :
			extern void adjustYeeCoordinates(const int &numSites, int (&xx)[16], int (&yy)[16], const int &subsqsideX, const int &subsqsideY);
			real cscore;
			int subsqsideX;
			int subsqsideY;
			printf("X Y sidelengths of subsquare in which you want the random points (200 for full square):\n");
			scanf("%d %d", &subsqsideX,  &subsqsideY);
			if((subsqsideX <=0) || (subsqsideX>=200)) {
				subsqsideX = 200;
			}
			if((subsqsideY <=0) || (subsqsideY>=200)) {
				subsqsideY = 200;
			}
			printf("using %dx%d centered subrectangle\n", subsqsideX, subsqsideY);
			adjustYeeCoordinates(numberOfSites, xarray, yarray, subsqsideX, subsqsideY);
			cscore = ReorderForColorContrast(  numberOfSites, xarray, yarray );
			printf("Color score %f (big=more constrast); Your coords are:\n", cscore);
			for(i=0; i<numberOfSites; i++) {
				printf("(%d, %d)\n", xarray[i], yarray[i]);
			}
			finished = true;
			break;
		default :
			printf("Wrong choice %u, moron - try again\n", u);
			break;
		}
	}
}

/*	RequestCoordPairReordering(numberOfSites, xarray, yarray):	asks the User if
 *									the coord-pairs
 *									should be
 *									reordered to try
 *									for maximum
 *									color-contrast
 *
 *	numberOfSites:	the number of coord-pairs to possibly reorder
 *	xarray:		a reference to an array of x coordinates
 *	yarray:		a reference to an array of y coordinates
 */
void RequestCoordPairReordering(int numberOfSites, int (&xarray)[16], int (&yarray)[16])
{
	uint u;
	bool finished;
	printf("Do you want IEVS to re-order the points to try for maximum color-contrast? (1) yes (2) no\n");
	for(finished = false; not finished; ) {
		scanf("%u", &u);
		switch(u) {
		case(1) :
			real cscore;
			int i;
			printf("Reordering...\n");
			cscore = ReorderForColorContrast(  numberOfSites, xarray, yarray );
			printf("Color score %f (big=more constrast); Your (reordered) coords are:\n", cscore);
			for(i=0; i<numberOfSites; i++) {
				printf("(%d, %d)\n", xarray[i], yarray[i]);
			}
			finished = true;
			break;
		case(2) :
			printf("OK, leaving points ordered as is.\n");
			finished = true;
			break;
		default :
			printf("Wrong choice %u, moron - try again\n", u);
			break;
		}
	}
}

/*	RequestElectorateSize(void):	asks the User for the number of Voters to use
 *					when making the Yee picture; the procedure
 *					returns the number inputted, subject to
 *					limitation
 */
int RequestElectorateSize(void)
{
	int Voters;
	printf("What max election size (#voters) would you like?\n");
	printf("256 recommended as good compromise between speed and randomness.\n");
	printf("You're allowed to go as high as %d (for slowest speed).\n", MaxNumVoters);
	printf("Algorithm keeps redoing elections with 10%% more voters each time until\n");
	printf("either confident know the winner, or reach this voter# bound.\n");
	for(;;) {
		scanf("%d", &Voters);
		if((Voters<=0) || (Voters>MaxNumVoters)) {
			printf("%d out of range, moron - try again\n", Voters);
		} else {
			printf("Using TopYeeVoters=%d.\n", Voters);
			break;
		}
	}
	return Voters;
}

/*	RequestStandardDeviation(void):	asks the User for a size for the standard
 *					deviation in the 1D gaussian
 */
int RequestStandardDeviation(void)
{
	int size;
	printf("What standard deviation on the 1D gaussian do you want? (Whole picture width is 200.)\n");
	for(;;) {
		scanf("%d", &size);
		if((size<=0) || (size>999)) {
			printf("%d out of range, moron - try again\n", size);
		} else {
			printf("Using GaussStdDevX=%d.\n", size);
			break;
		}
	}
	return size;
}

/*	RequestHonestyPercentage(void):	asks the User for a Voter honesty percentage;
 *					the procedure returns that percentage, subject
 *					to limitations
 */
int RequestHonestyPercentage(void)
{
	int percentage;
	printf("What honesty-percentage do you want? (0 to 100.)\n");
	for(;;) {
		scanf("%d", &percentage);
		if((percentage<=0) || (percentage>100)) {
			printf("%d out of range, moron - try again\n", percentage);
		} else {
			printf("Using honfrac=%d%%.\n", percentage);
			break;
		}
	}
	return percentage;
}

/*	RequestUtilityDistancePow(void):	asks the User to choose utilities to be
 *						based on either the L1 or L2 distance;
 *						the procedure returns the chosen
 *						distance type
 */
int RequestUtilityDistancePow(void)
{
	int LpPow;
	printf("Utilities based on (1) L1 or (2) L2 distance?\n");
	for(;;) {
		scanf("%d", &LpPow);
		if((LpPow<=0) || (LpPow>2)) {
			printf("%d out of range, moron - try again\n", LpPow);
		} else {
			printf("Using LpPow=%d.\n", LpPow);
			break;
		}
	}
	return LpPow;
}

/*	PerformBayesianRegretAnalysis(void):	asks the User a series of questions and
 *						uses those answers to perform
 *						Bayesian regret analysis and to output
 *						the results
 */
void PerformBayesianRegretAnalysis(void)
{
	extern void DisplayBayesianRegretMenuIntroduction(void);
	driver_t driver;
	parameters parameterType;
	extern void RequestAgreementCountStatus(void);
	extern void RequestBayesianRegretNormalizationMethod(void);
	extern void RequestBayesianRegretSortingMethod(void);
	extern void RequestErrorBarStatus(void);
	extern void RequestIntermethodWinnerAgreementCountDisplayStatus(void);
	extern void RequestOutputFormat(void);
	extern parameters RequestParameters(void);
	extern void RequestRegretOutput(void);
	extern driver_t RequestUtilityGenerators(parameters parameterType);
	DisplayBayesianRegretMenuIntroduction();
	RequestBayesianRegretSortingMethod();
	RequestOutputFormat();
	RequestBayesianRegretNormalizationMethod();
	RequestErrorBarStatus();
	RequestAgreementCountStatus();
	RequestIntermethodWinnerAgreementCountDisplayStatus();
	RequestRegretOutput();
	parameterType = RequestParameters();
	driver = RequestUtilityGenerators(parameterType);
	driver();
}

/*	MakeYeePicture(seed):	asks the USer a series of questions and uses those
 *				answers to create a Yee picture and to output the result
 *
 *	seed:	the seed value used for the random number generator
 */
void MakeYeePicture(uint seed)
{
	char fname[100];
	int ihonfrac;
	int GaussStdDev;
	int LpPow;
	int NumSites;
	int TopYeeVoters;
	int WhichMeth;
	int xx[16];
	int yy[16];
	extern void RequestCoordPairs(int numberOfSites, int (&xarray)[16], int (&yarray)[16]);
	extern void RequestCoordPairReordering(int numberOfSites, int (&xarray)[16], int (&yarray)[16]);
	extern int RequestElectorateSize(void);
	extern int RequestHonestyPercentage(void);
	extern void RequestNameForBMP(char (&name)[100]);
	extern 	int RequestPointSiteCount(void);
	extern int RequestStandardDeviation(void);
	extern int RequestUtilityDistancePow(void);
	extern int RequestVotingMethod(void);
	WhichMeth = RequestVotingMethod();
	RequestNameForBMP(fname);
	NumSites = RequestPointSiteCount();
	RequestCoordPairs(NumSites, xx, yy);
	RequestCoordPairReordering(NumSites, xx, yy);
	TopYeeVoters = RequestElectorateSize();
	GaussStdDev = RequestStandardDeviation();
	ihonfrac = RequestHonestyPercentage();
	LpPow = RequestUtilityDistancePow();
	printf("grinding...\n");
	MakeYeePict( fname, xx, yy, NumSites, WhichMeth, TopYeeVoters, GaussStdDev, ihonfrac*0.01, LpPow );
	printf("seed=%d\n", seed);
}

/*	resetFavorites(Voters, Candidate):	resets the favorite
 *						Candidate of each Voter
 *						to the given Candidate
 *	Voters:		a set of Voters to have favorites reset
 *	Candidate:	the Candidate to which the favorite is to be
 *			reset
 */
void resetFavorites(oneVoter (&Voters)[MaxNumVoters], const int64_t& Candidate) {
	for( oneVoter& eachVoter : Voters ) {
		eachVoter.favoriteCandidate = Candidate;
	}
}

/*	Zero(number, allCandidates, member):	sets the 'member' of
 *						the first 'number' of
 *						Candidates in 'allCandidates'
 *						to 0
 *	number:		the number of Candidates to set a given member
 *			to 0
 *	allCandidates:	the slate of Candidates with a member to set
 *	member:		the member of Each Candidate to set to 0
 */
template <class T>
	void Zero(const uint64_t& number, CandidateSlate& allCandidates, T oneCandidate::*member) {
		for(uint64_t count = 0; count < number; count++) {
			allCandidates[count].*member = 0;
		}
}

/*	traditionalVoteVectorNormalization(theVoter, Candidates, count):	normalizes
 *				the vote vector of non-eliminated Candidates provided by
 *				'theVoter' and returns it; the normalization process is as
 *				follows: (1) each vote value of non-eleiminated Candidates
 *				is raised to the IRNRPOWER (2) the resulting values are
 *				summed (3) the IRNRPOWERth root of the resulting sum is
 *				calcualted (4) each vote value is divided by the resulting
 *				root (5) the collection of quotients constitutes the
 *				normalized vote vector
 *	theVoter:	the Voter with the vote vector to
 *			normalize
 *	Candidates:	the slate of Candidates to consider
 *	count:		the number of Candidates
 *
 *	NOTE #1:	As of this point in time, I am not sure
 *			why each vote is first reduced by 0.5. If
 *			I learn the reason, I will update this
 *			comment.
 *	Note #2:	Is there a good reason to require the sum
 *			exceed 0? If IRNRPOWER is even, the sum is
 *			guaranteed to be non-negative. If
 *			IRNRPOWER is negative, a negative sum is
 *			acceptable for taking the root.
 */
voteVector traditionalVoteVectorNormalization(const oneVoter& theVoter, const CandidateSlate& Candidates, const uint64_t& count)
{
	real s = 0.0;
	real t;
	voteVector normalizedVoteVector;
	const oneCandidateToTheVoter (&allCandidatesToTheVoter)[MaxNumCands] = theVoter.Candidates;
	for(int j=0; j<count; j++) {
		if(not Candidates[j].eliminated) {
			t = allCandidatesToTheVoter[j].score - 0.5;
			t = std::abs(t);
			s += pow(t, IRNRPOWER);
		}
	}
	if(s>0.0) {
		s = pow(s, -1.0/IRNRPOWER);
		for(int j=0; j<count; j++) {
			if(not Candidates[j].eliminated) {
				normalizedVoteVector[j] = allCandidatesToTheVoter[j].score * s;
			}
		}
	}
	return normalizedVoteVector;
}

/*	addToNormalizedRatingSum(normalizedVotes, Candidates, count):	adds normalized vote values
 *															to the normalized rating sums
 *															of Each non-eliminated Candidate
 *	normalizedVotes:	the normalized vote vector
 *	Candidates:		the set of Candidates to add the normalized votes
 *	count:			the number of Candidates to consider
 */
void addToNormalizedRatingSum(const voteVector& normalizedVotes, CandidateSlate& allCandidates, const uint64_t& count)
{
	for(uint64_t j=0; j<count; j++) {
		oneCandidate& EachCandidate = allCandidates[j];
		if(not EachCandidate.eliminated) {
			EachCandidate.normalizedRatingSum += normalizedVotes[j];
		}
	}
}

/*	linearVoteVectorNormalization(theVoter, Candidates, count):	normalizes
 *				the vote vector of non-eliminated Candidates provided
 *				by 'theVoter' and returns it; the normalization
 *				process is just like that described for
 *				'traditionalVoteVectorNormalization()' but it
 *				performs a TWO-parameter linear transformation
 *				so mean=0 and variance=1.
 *	theVoter:	the Voter with the vote vector to
 *			normalize
 *	Candidates:	the slate of Candidates to consider
 *	count:		the number of Candidates
 *
 *	NOTE #1:	As of this point in time, I am not sure
 *			why each vote is reduced by the average. If
 *			I learn the reason, I will update this
 *			comment.
 *	Note #2:	There seems to be conditions in the test output
 *			when All Candidates have a score of 1. How should
 *			these cases be handled?
 */
voteVector linearVoteVectorNormalization(const oneVoter& theVoter, const CandidateSlate& Candidates, const uint64_t& count)
{
	uint64_t countedThisRound = 0;
	real mean;
	real s = 0.0;
	real t;
	voteVector normalizedVoteVector;
	const oneCandidateToTheVoter (&allCandidatesToTheVoter)[MaxNumCands] = theVoter.Candidates;
	for(int j=0; j<count; j++) {
		if(not Candidates[j].eliminated) {
			s += allCandidatesToTheVoter[j].score;
			countedThisRound++;
		}
	}
	assert(countedThisRound>0);
	ensure(countedThisRound>0, 5);
	mean = s/countedThisRound;
	s = 0.0;
	normalizedVoteVector.fill(0);
	for(int64_t j=(count-1); j>=0; j--) {
		if(not Candidates[j].eliminated) {
			t = allCandidatesToTheVoter[j].score - mean;
			s += t*t;
		}
	}
	if(s>0.0) {
		s = 1.0/sqrt(s);
		for(int j=0; j<count; j++) {
			if(not Candidates[j].eliminated) {
				const real& adjustedScore = (allCandidatesToTheVoter[j].score - mean);
				normalizedVoteVector[j] = adjustedScore * s;
			}
		}
	}
	return normalizedVoteVector;
}

/*	minMaxStyleVoteVectorNormalization(theVoter, Candidates, count):	normalizes
 *				the vote vector of non-eliminated Candidates
 *				provided by 'theVoter' and returns it;
 *				the normalization process is a 2-parameter
 *				linear transformation renormalization
 *				so the maximum equals 1 and the minimum
 *				equals 0; the procedure returns the
 *				index of the Winner or -1 if an error
 *				occurs
 *	theVoter:	the Voter with the vote vector to
 *			normalize
 *	Candidates:	the slate of Candidates to consider
 *	count:		the number of Candidates
 */
voteVector minMaxStyleVoteVectorNormalization(const oneVoter& theVoter, const CandidateSlate& Candidates, const uint64_t& count)
{
	real s = 0.0;
	voteVector normalizedVoteVector;
	const oneCandidateToTheVoter (&allCandidatesToTheVoter)[MaxNumCands] = theVoter.Candidates;
	real theMaximum = -HUGE;
	real theMinimum = HUGE;
	for(int j=0; j<count; j++) {
		if(not Candidates[j].eliminated) {
			const real& t = allCandidatesToTheVoter[j].score;
			theMaximum = std::max(t, theMaximum);
			theMinimum = std::min(t, theMinimum);
		}
	}
	normalizedVoteVector.fill(0);
	s = theMaximum - theMinimum;
	if(s>0.0) {
		s = 1.0/s;
		for(int j=0; j<count; j++) {
			if(not Candidates[j].eliminated) {
				const real& adjustedScore = allCandidatesToTheVoter[j].score - theMinimum;
				normalizedVoteVector[j] = adjustedScore * s;
			}
		}
	}
	return normalizedVoteVector;
}
