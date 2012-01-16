#!/bin/sh

./main
x=$?

if test $x -eq 111; then

cat>result.h << EOF

/**
 *  Class-static variables are NOT initialized properly in libraries.
 */

#define HINTLIB_STATIC_WORKS 0
EOF

else
if test $x -eq 110; then

cat>result.h <<EOF

// created by   libtest/create

/**
 *  Class-static variables are initialized properly in libraries.
 */

#define HINTLIB_STATIC_WORKS 1
EOF

  else
cat>result.h <<EOF

// created by   libtest/create

/**
 *  Unable to determin if class-static variables wor in libraries.
 *
 *  Assuming fail-save defaults and hoping for the best.
 */

#define HINTLIB_STATIC_WORKS 0
EOF

cat <<EOF

// created by   libtest/create

***************************************************

Unable to create testlib!

Assuming failsave defaults and hoping for the best.

***************************************************

EOF

  fi
fi

