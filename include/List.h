/*****************************************************************************
 * List.h:  template list operations.                                        *
 *****************************************************************************/

/* $Id$ */

#ifndef ListH
#define ListH

template<class T>
class ListIterator;


template<class T>
class List {
public:
  List() : head(0), tail(0) {}
  ~List();


  void add    (T x);
  void remove (T x);

  friend class ListIterator<T>;

private:
  class Node {
  public:
    Node(T x) : next(0), datum(x) {}
    Node* next;
    T     datum;
  };

  Node* head;
  Node* tail;
  List(const List<T>&);                     // Prohibit, since not implemented.
  List<T>& operator=(const List<T>&);       // Prohibit, since not implemented.
};



template<class T>
class ListIterator {
public:
  ListIterator(const List<T>& list) : cur(list.head) {}
  
  int     more    () const { return cur != 0;   }
  T       current () const { return cur->datum; }
  void    advance ()       { cur = cur->next;   }

private:
  List<T>::Node* cur;
};




template<class T>
inline List<T>::~List() {
  while (head != 0) {
    Node* p = head->next;
    delete head;
    head = p;
  }
};


template<class T>
inline void List<T>::add(T x) {
  if (head == 0)
    head = tail = new Node(x);
  else
    tail = tail->next = new Node(x);
};


template<class T>
inline void List<T>::remove(T x) {
  Node* prev = 0;
  Node* cur = head;
  while (cur != 0) {
    if (cur->datum == x) {
      if (prev == 0) {
        head = cur->next;
        delete cur;
        cur = head;
      } else {
        prev->next = cur->next;
        delete cur;
        cur = prev->next;
      }
    } else {
      prev = cur;
      cur = cur->next;
    }
  }
};


#endif
