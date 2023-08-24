#ifndef _applic_h_
#define _applic_h_


#include <stdlib.h>
#include "parallel.h"
#include "lodepng.h"
#define DEBUG

 
#define SMALLEST_TIME_LOGARITMIZE	(1.0e-9)
#define UCHAR_RANGE_MIN				0
#define UCHAR_RANGE_MAX				255
#define INTEGER_RANGE_MAX			2147483647
#define SHORT_MAX					32767
#define RED							0
#define GREEN						1
#define BLUE						2
#define ALPHA						3
#define RGB							3
#define WHITE_RGB					(255)

//lodepng
#define REDCHAN						0
#define GREENCHAN					1
#define BLUECHAN					2
#define ALPHACHAN					3
#define	DEFAULT_ZSECTIONING_ZPOS	(0.5)
#define BLACK						0
#define WHITE						1
#define DEFMODUS					2


#define _STR(x) _VAL(x)
#define _VAL(x) #x
#define ASSERT(cond) if(!(cond)) exitus("failed assertion:"__FILE__"line"_STR(__LINE__)":"#cond)
#define TOLFP(x) ((1.0)-(0.99999/((double) x))) //tolerance
#define ERRTXT(text) (text" : file: "__FILE__" line:"_STR(__LINE__)"\n")
#define FACC8(x) (float)((int)(100000000*x))/100000000.

#define NOW true
#define NOTYET false

#define NEVER true


#ifndef DEBUG
        #define QUICKASSERT(cond)       ((void)0)
#else
#define QUICKASSERT(cond)       ASSERT(cond)
#endif

void exitus (const char *s);
typedef double Real;

//class trigger for time dependent events

class trigger
{
public:
	trigger( long freq = 0 ) : frequency(freq) { never = false; whenOut = 10000000; } 		//constructor of trigger class
	~trigger( void ); //destructor of trigger class
	bool now(void); //Boolean method now for time dependent events

public:
	long frequency;
	bool never;
	long whenOut;
};

#endif