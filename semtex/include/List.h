/*****************************************************************************
 * List.h:  template list operations.
 *
 * Reference: Barton & Nackman, "Scientific & Engineering C++".
 *****************************************************************************/

// $Id$

#ifndef ListH
#define ListH

template<class T>
class ListIterator;





template<class T>
class List {
friend class ListIterator<T>;

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

  void add (T x) {
    if (head == 0)
      head = tail = new Node(x);
    else
      tail = tail -> link = new Node(x);
    nel++;
  }

  void remove (T x) {
    Node* prev = 0;
    Node* cur  = head;
    while (cur != 0) {
      if (cur -> datum == x) {
	if (prev == 0) {
	  head = cur -> next;
	  delete cur;
	  cur = head;
	} else {
	  prev -> next = cur -> link;
	  delete cur;
	  cur = prev -> link;
	}
	nel--;
      } else {
	prev = cur;
	cur  = cur -> link;
      }
    }
  }

  void clear  () { head = tail = 0; nel = 0; }

  int  length () const { return nel;                       }
  T    first  () const { return (nel) ? head -> datum : 0; }

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

  List(const List<T>&);                // Prohibit, since not implemented 
  List<T>& operator=(const List<T>&);  // Prohibit, since not implemented 

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
