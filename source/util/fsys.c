/* $Id: fsys.c,v 1.1.1.1 2003/06/17 10:15:23 hal Exp $ */
/* 

   fsys_: emulate system(3) call for Fortran

*/


#ifdef DEBUG
#include <stdio.h>
#endif /* DEBUG */

#include <stdlib.h>
#include <errno.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/times.h>

#define SHELL "sh"
#define SHELLPATH "/bin/sh"

typedef int f77_int;		/* enhance portability? */

extern char **environ;



void fsys_ (const char *command, f77_int *status, f77_int *error, f77_int len)
{
  pid_t pid;
  char *buff, *p;



  *status = 0;
  *error = errno = 0;

  /* allocate enough memory for command + \0 */
  if ( (p = (char *) malloc (len+1)) == NULL) {
    *status = -1;
    *error = errno;

    return;
  }

  buff = p;

  /* copy command into temporary buffer and terminate with \0 */
  while (len > 0) {
    *p++ = *command++;
    len--;
  }

  *p = 0;

  /* if no command was given return immediately with an error */
  if (buff == 0) {
    free (buff);

    *status = 1;

    return;
  }


  /* fork and exec child */
  switch (pid = fork ()) {

  case -1:			/* an error occured while forking */
    free (buff);

    *status = -1;
    *error = errno;

    return;

  case 0:			/* child process */
    {
      char *argv[4];

      /* create command line */
      argv[0] = SHELL;
      argv[1] = "-c";
      argv[2] = buff;
      argv[3] = 0;

      *status = execve (SHELLPATH, argv, environ);

      exit (127);		/* if anything goes wrong */
    }
  default:			/* parent process */
    waitpid (pid, (int *) status, 0);
  }

  free (buff);

  *error = errno;


  return;
}
