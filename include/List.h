///////////////////////////////////////////////////////////////////////////////
// List.h:  template list operations.
//
// Reference: Barton & Nackman, "Scientific & Engineering C++".
///////////////////////////////////////////////////////////////////////////////

// $Id$

#ifndef ListH
#define ListH

template<class T> class ListIterator;


template<class T> class List {
friend class ListIterator<T>;
private:
  class Node {
  public:
    Node(T x) : link(0), datum(x) {}
    Node* link;
    T     datum;
  };

  Node* head;
  Node* tail;
  int   nel;

  List(const List<T>&);                // -- Prohibit, since not implemented. 
  List<T>& operator=(const List<T>&);  // -- Prohibit, since not implemented. 
public:
  List() : head(0), tail(0), nel(0) {}

  ~List() {
    while (head != 0) {
      Node* p = head -> link;
      delete head;
      head    = p;
    }
    nel = 0;
  }

  void add (T x) {		// -- Unconditional insertion.
    if (head == 0) {
      head = new Node (x);
      tail = head;
    } else
      tail = tail -> link = new Node (x);
    nel++;
  }

  int xadd (T x) {		// -- Conditional insertion.
    register int   found = 0;
    register Node* ptr;

    for (ptr = head; !found && ptr; ptr = ptr -> link) found = x == ptr->datum;
    if   (found)  {          return 0; }
    else          { add (x); return 1; }
  }


  void remove (T x) {
    Node* prev = 0;
    Node* curr = head;
    while (curr != 0) {
      if (curr -> datum == x) {
	if (prev == 0) {
	  head = curr -> link;
	  delete curr;
	  curr = head;
	} else {
	  prev -> link = curr -> link;
	  delete curr;
	  curr = prev -> link;
	}
	nel--;
      } else {
	prev = curr;
	curr = curr -> link;
      }
    }
  }

  void clear  () { head = tail = 0; nel = 0; }

  int  length () const { return nel;                       }
  T    first  () const { return (nel) ? head -> datum : 0; }

};





template<class T>
class ListIterator {
public: 
  ListIterator (const List<T>& list) : cur(list.head), top(list.head) { }
  
  int   more    () const  { return cur != 0;           }
  T     current () const  { return cur -> datum;       }
  void  next    ()        {        cur =  cur -> link; }
  void  reset   ()        {        cur =  top;         }

private:
  List<T>::Node* cur;
  List<T>::Node* top;
};


#endif
