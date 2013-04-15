/* Written by Shamsundar as a poor man's replacement for
   the vararglist function error(status,errnum,fmt,...) */
int error_message_count=0;
#define ERR0(status,errnum,fmt) \
   {fprintf(stderr,"Program %s, File %s, Line %d: ", \
    argv[0],__FILE__,__LINE__); \
   fprintf(stderr,fmt); \
   if(errnum)fprintf(stderr,": Err.Num. %d\n",errnum); \
   if(status)exit(status); \
   error_message_count++; \
   }
#define ERR1(status,errnum,fmt,p1) \
   {fprintf(stderr,"Program %s, File %s, Line %d: ", \
    argv[0],__FILE__,__LINE__); \
   fprintf(stderr,fmt,p1); \
   if(errnum)fprintf(stderr,": Err.Num. %d\n",errnum); \
   if(status)exit(status); \
   error_message_count++; \
   }

