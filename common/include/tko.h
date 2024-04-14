#ifndef _TKO_H
#define _TKO_H

// TKO Header Length
#define TKO_HEADER_LEN 64

#define STAR 0
#define CMR  1
#define SNCR 2
#define DSR  3
#define ENR  4
#define BCR  5
#define ISIZE 6
#define CID  8
#define SEC  9
#define USEC 10

#define NORMAL 1
#define DRT 0

//****Please modify these params to suit the data structure ***//

#define MAXSMP 11
#define NSCALER 40

//crate type now defined in the counter map man

//static int crate_type[MAXSMP]={NORMAL, DRT, NORMAL, NORMAL, DRT, DRT};
//static int crate_type[MAXSMP]={NORMAL, DRT, NORMAL, NORMAL, DRT, DRT, DRT};
//static int crate_type[MAXSMP]={NORMAL, DRT, NORMAL, NORMAL, DRT, DRT, DRT, NORMAL};
//static int crate_type[MAXSMP]={DRT, DRT, NORMAL, NORMAL, DRT, DRT, DRT, NORMAL};
//static int crate_type[MAXSMP]={DRT, DRT, NORMAL, DRT, DRT, DRT, NORMAL}; // modified by k.tsukada 2010/9/27
//static int crate_type[MAXSMP]={DRT, DRT, NORMAL, DRT, DRT, DRT, NORMAL,NORMAL,NORMAL}; // modified by sada 2012/4/2
//static int crate_type[MAXSMP]={DRT, DRT, NORMAL, DRT, DRT, DRT, NORMAL,NORMAL,NORMAL,NORMAL,DRT}; // FDC1 added 2012/12/16

#define NUMBER_OF_RUNTYPE 3
//static char file_prefixes[NUMBER_OF_RUNTYPE][6]={"e471-","e5490","e5700"};
//static char *file_prefix;
//static char file_prefix[]=
//static char file_prefix[]=

//*************************************************************//

#define OVFL 0x0000FFFF // Overflow = 0xFFFF

#define DRT_OVFL 0x7FF // Overflow of DR.T (11bit)

#define NPAWC  8100000 // <- kinko ok in hlimap
//#define NPAWC   50370000  <- itapc
//#define NPAWC   300000000 <- kinko x
//#define NPAWC  200000000 // <- kinko ok in hlimit

#define BSIZE 0x3FFFFF // Buffer size for serial buffer
#define MAXCH 0xFFFF   // Maximum number of channel in a single event

#define NA -1

#define UNIDAQ_HEADER_SIZE 6
#define HEADER_LENGTH 64
#define SCALER_LENGTH 16

#define TDC 0
#define ADC 1

typedef struct {
  // A struct for raw data
  int c; // crate number
  int m; // module address=module slot number
  int s; // sub address=channel number
  int h; // multiple hit
  int v; // value of the channel
} s_raw;

typedef struct {
  // A struct for smp information
  int star;
  int cmr;
  int sncr;
  int dsr;
  int enr;
  int bcr;
  int isize;
  int cid;
  int sec;
  int usec;
} s_header;

typedef struct {
  int isize;
  int cid;
  int *beginning;
} s_buffer;

typedef struct {
  // A struct for UNIDAQ information
  int b_event_length;
  int b_event_ID;
  int run_no;
  int evt_no;
  int hdr_mode;
  int reserved;
} s_unidaq_header;

typedef struct {
  s_header header;
  int length;
  int serial_number;
  int *data;
} s_sds;

typedef struct{
  int number;
  char fname[80];
  int type;
  char title[64];
  char comment[256];
} s_run;

#else

#endif













