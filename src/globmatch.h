
#ifndef GLOBMATCH_H
#define GLOBMATCH_H 1

#ifdef  __cplusplus
extern "C" {
#endif

  
#ifndef GLOBMATCH_NEGATE
#define GLOBMATCH_NEGATE '^'       /* std char set negation char */
#endif


int globmatch(char *string, char *pattern);


#ifdef  __cplusplus
}
#endif

#endif /* globmatch.h  */
