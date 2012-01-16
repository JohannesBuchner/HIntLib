
#ifndef INCLUDE_H
#define INCLUDE_H

#define MAGIC 0xaffe

class A
{
public:
   A();
   int get() const  { return a; }

private:

   int a;
};

int f();

class X
{
public:

   static const int    constIntVar;
   static const double constDoubleVar;
   static int    intVar;
   static double doubleVar;

   static A a;
};

#endif

