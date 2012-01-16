#include <stdio.h>
#include <stdlib.h>

#include "include.h"

int main()
{
   bool x = true;

   if (X::intVar != MAGIC)
   {
      puts ("\"static int\" in library is broken!\n");
      x = false;
   }
   if (X::constIntVar != MAGIC)
   {
      puts ("\"static const int\" in library is broken!\n");
      x = false;
   }
   if (X::a.get() != MAGIC)
   {
      puts ("static object in library is broken!\n");
      x = false;
   }

   // 110   ok
   // 111   failed

   exit (x ? 110 : 111);
}

