#ifndef STACK_H
#define STACK_H
///////////////////////////////////////////////////////////////////////////////
// Stack.h: templated operations for LIFO stack.
//
// Summary: creator, destructor, push, pop, depth.
//
// $Id$
///////////////////////////////////////////////////////////////////////////////


template<class T>
class Stack {
public:
  Stack () : top (0), stack_depth (0) {}
 ~Stack ();

  void  push    (T x);
  T     pop     ();
  int   depth   () const { return stack_depth; }

private:
  class Node {
  public:
    Node (T x) : next (0), datum (x) {}
    Node* next;
    T     datum;
  };

  Node* top;
  int   stack_depth;

  Stack(const Stack<T>&);                // Prohibit, since not implemented.
  Stack<T>& operator=(const Stack<T>&);  // Prohibit, since not implemented.
};


template<class T>
inline void Stack<T>::push(T x) {
  if (stack_depth) {
    Node* p   = new Node (x);
    p -> next = top;
    top       = p;
  } else {
    top = new Node (x);
  }
  stack_depth++;
}


template<class T>
inline T Stack<T>::pop() {
  if (stack_depth) {
    Node* p     = top;
    T     value = top -> datum;
    top         = top -> next;
    delete p;
    stack_depth--;
    return value;
  } else {
    return 0;
  }

}


template<class T>
inline Stack<T>::~Stack() {
  while (stack_depth--) {
    Node* p = top -> next;
    delete top;
    top = p;
  }
}

#endif
