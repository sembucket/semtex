#ifndef STACK_H
#define STACK_H
///////////////////////////////////////////////////////////////////////////////
// Stack.h: templated operations for LIFO stack.
//
// Summary:
//   creator, destructor, push, pop, isEmpty.
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include <iostream.h>

template<class T> class StackNode {
public:
  StackNode (T x) : next (0), datum (x) {}
  StackNode<T>* next;
  T     datum;
};

template<class T> class Stack {
private:
  StackNode<T>* top;

  Stack(const Stack<T>&);                // Prohibit, since not implemented.
  Stack<T>& operator=(const Stack<T>&);  // Prohibit, since not implemented.

public:
  int   isEmpty () { return top == 0; }

  Stack  () : top (0) { }
  ~Stack (){ while (!isEmpty()) {
    StackNode<T>* p = top->next;
    delete top;
    top = p;} 
  }

  void  push    (T x) {
    if (isEmpty())
      top = new StackNode<T> (x);
    else {
      StackNode<T>* p = new StackNode<T> (x);
      p -> next = top;
      top = p;
    }
  }

  T     pop     () {
    if (isEmpty()) {
      return 0;
    } else {
      StackNode<T>* p     = top;
      T             value = top -> datum;
      top = top -> next;
      delete p;
      return value;
    }
  }

};

#endif
