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


template<class T>
class Stack {
public:
  int   isEmpty () { return top == 0; }

  Stack  () : top (0) { }
  ~Stack (){ while (!isEmpty()) {Node* p = top -> next; delete top; top = p;} }

  void  push    (T x) {
    if (isEmpty())
      top = new Node (x);
    else {
      Node *p = new Node (x);
      p -> next = top;
      top = p;
    }
  }

  T     pop     () {
    if (isEmpty()) {
      return 0;
    } else {
      Node *p     = top;
      T     value = top -> datum;
      top = top -> next;
      delete p;
      return value;
    }
  }

private:
  class Node {
  public:
    Node (T x) : next (0), datum (x) {}
    Node* next;
    T     datum;
  };

  Node* top;

  Stack(const Stack<T>&);                // Prohibit, since not implemented.
  Stack<T>& operator=(const Stack<T>&);  // Prohibit, since not implemented.
};

#endif
