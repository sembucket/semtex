/*****************************************************************************
 * Domain.C:  implement Domain class.                                        *
 *****************************************************************************/

// $Id$


#include "Fem.h"


Domain::Domain() {
  u  = 0;
  uB = 0;

  domain_name      = new char[STR_MAX];
  domain_field_tag = new char[STR_MAX];
  domain_edge_tag  = new char[STR_MAX];
  domain_name[0]   = domain_field_tag[0] = domain_edge_tag[0] = '\0';
  domain_step      = 0;
  domain_time      = 0.0;
  domain_nfield    = 0;
  domain_nedge     = 0;
  domain_nel       = 0;
}

Domain::~Domain() {
  delete [] domain_name;
  delete [] domain_field_tag;
  delete [] domain_edge_tag;
  domain_state_file.  close();
  domain_history_file.close();
  domain_field_file.  close();
}

const char *Domain::name() const          { return domain_name;      }
void        Domain::name(const char* s)   { strcpy (domain_name, s); }
const char* Domain::fTag() const          { return domain_field_tag; }
const char* Domain::eTag() const          { return domain_edge_tag;  }

int&        Domain::step()                { return domain_step;   }
double&     Domain::time()                { return domain_time;   }
int         Domain::nField() const        { return domain_nfield; }
int         Domain::nBedge() const        { return domain_nedge;  }
int&        Domain::nEl()                 { return domain_nel;    }

ofstream& Domain::stateout()              {
  return domain_state_file; }
void      Domain::stateout(const char* s) {
  domain_state_file.close(); domain_state_file.open(s); }

ofstream& Domain::histout()               {
  return domain_history_file; }
void      Domain::histout(const char* s)  {
  domain_history_file.close(); domain_history_file.open(s); }

ofstream& Domain::fieldout()              {
  return domain_field_file; }
void      Domain::fieldout(const char* s) {
  domain_field_file.close(); domain_field_file.open(s); }


void Domain::addField (Element *E) {
  Element **X = new Element*[domain_nfield+1];

  for (int i = 0; i < domain_nfield; i++) X[i] = u[i];
  X[domain_nfield] = E;

  delete [] u;
  u = X;

  strncat (domain_field_tag, &(E -> name), 1);
  domain_nfield++;
}


void Domain::addBedge (Bedge *B) {
  Bedge **X = new Bedge*[domain_nedge+1];
  
  for (int i = 0; i < domain_nedge; i++) X[i] = uB[i];
  X[domain_nedge] = B;

  delete [] uB;
  uB = X;

  strncat (domain_edge_tag, &(B -> name), 1);
  domain_nedge++;
}

#include <time.h>

ostream& operator << (ostream& os, Domain& D)
/* ========================================================================= *
 * Output all Domain field variables on ostream in prism-compatible form.    *
 * ========================================================================= */
{
  char *hdr_fmt[] = { 
    "%-25s "    "Session\n",
    "%-25s "    "Created\n",
    "%-25s "    "Nr, Ns, Nz, Elements\n",
    "%-25d "    "Step\n",
    "%-25.6g "  "Time\n",
    "%-25.6g "  "Time step\n",
    "%-25.6g "  "Kinvis\n",
    "%-25.6g "  "Beta\n",
    "%-25s "    "Fields Written\n",
    "%-25s "    "Format\n"
  };

  char      routine[] = "Domain << operator";
  char      buf1[STR_MAX], buf2[STR_MAX];
  int       np   = D.u[0] -> np;
  int       ntot = sqr(np) * D.domain_nel;
  time_t    tp   = ::time(0);

  sprintf (buf1, hdr_fmt[0], D.domain_name);
  os << buf1;

  strftime (buf2, 25, "%a %b %d %H:%M:%S %Y", localtime(&tp));
  sprintf  (buf1, hdr_fmt[1], buf2);
  os << buf1;

  sprintf (buf2, "%1d %1d %1d %1d", np, np, 1, D.domain_nel);
  sprintf (buf1, hdr_fmt[2], buf2);
  os << buf1;

  sprintf (buf1, hdr_fmt[3], D.domain_step);
  os << buf1;

  sprintf (buf1, hdr_fmt[4], D.domain_time);
  os << buf1;

  sprintf (buf1, hdr_fmt[5], dparam("DELTAT"));
  os << buf1;

  sprintf (buf1, hdr_fmt[6], dparam("KINVIS"));
  os << buf1;

  sprintf (buf1, hdr_fmt[7], dparam("BETA"));
  os << buf1;

  sprintf (buf1, hdr_fmt[8], D.domain_field_tag);
  os << buf1;

  sprintf (buf1, hdr_fmt[9], (option("BINARY")) ? "binary" : "ASCII");
  os << buf1;

  if (option ("BINARY")) {
    register int n;
    for (n = 0; n < D.domain_nfield; n++)
      os.write((char*) *D.u[n] -> value, ntot*sizeof(double));
  } else {
    register int  i, n;
    os.setf (ios::scientific, ios::floatfield);
    os.setf (ios::uppercase);
    for (i = 0; i < ntot; i++) {
      for (n = 0; n < D.domain_nfield; n++)
        os << setw(14) << D.u[n] -> value[0][i];
      os << endl;
    }
  }

  os << flush;
  if (!os) message (routine, "failed writing field file", ERROR);

  return os;
}
