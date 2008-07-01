#include <stdio.h>

void flush_hack_(void)
{
  fflush(stdout);
  fflush(stderr);
}

