#include <stdio.h>
// #include <sys/types.h> // for ftruncate...

int
main ()
{
  fputs("output1\n",stdout);
  fputs("output2\n",stdout);
  fputs("\033[A\033[2K\033[A\033[2K",stdout);
  rewind(stdout);
  // ftruncate(1,0);
  fputs("output3\n",stdout);
  fputs("output4\n",stdout);
  return 0;
}
