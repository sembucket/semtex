///////////////////////////////////////////////////////////////////////////////
// PSplot package introduced in Numerical Recipes 3e Ch. 22, described in
// webnote 26:
//   http://www.nr.com/webnotes?26
//
// We have replaced Doub by double, Int with int, VecDoub with vector<double>
// (assuming templated vector library).
///////////////////////////////////////////////////////////////////////////////

#ifndef PSPLOT_H
#define PSPLOT_H

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

///////////////////////////////////////////////////////////////////////////////

struct PSpage {
  // Object for generating a PostScript page.
  
  static FILE* PLT;
  static char* file;
  char         fontname[128];
  double       fontsize;

  PSpage(char* filnam) {
    // Constructor. Argument is name of PostScript file to be created.
    file = new char[128];
    strcpy(file,filnam);
    PLT = fopen(file,"wb");
    if (!PLT) throw("failure opening output file for plot");
    fprintf(PLT,"%%!\n/mt{moveto}def /lt{lineto}def /np{newpath}def\n");
    fprintf(PLT,"/st{stroke}def /cp{closepath}def /fi{fill}def\n");
    fprintf(PLT,"/zp {gsave /ZapfDingbats findfont exch ");
    fprintf(PLT,"scalefont setfont moveto show grestore} def\n");
    setfont("Times-Roman",12.);
    setlinewidth(0.5);
  }

  PSpage(FILE* fp)
  // Constructor that takes an (already-initialised/open) FILE*.
  {
    PLT = fp;
    fprintf(PLT,"%%!\n/mt{moveto}def /lt{lineto}def /np{newpath}def\n");
    fprintf(PLT,"/st{stroke}def /cp{closepath}def /fi{fill}def\n");
    fprintf(PLT,"/zp {gsave /ZapfDingbats findfont exch ");
    fprintf(PLT,"scalefont setfont moveto show grestore} def\n");
    setfont("Times-Roman",12.);
    setlinewidth(0.5);
  }  
  PSpage() {}
  // Alternative contructor used internally, binds a PSplot to a PSpage.
  ~PSpage() {if (PLT) close();}

  void setfont(const char *fontnam, double size) {
    // Change font from default (12pt Times-Roman)
    strcpy(fontname,fontnam);
    fontsize = size;
    fprintf(PLT,"/%s findfont %g scalefont setfont\n",fontnam,size);
  }

  void setcolor(int r, int g, int b) {
    // Change colour from default (0,0,0=black).  Range is 0 to 255.
    fprintf(PLT,"%g %g %g setrgbcolor\n",r/255.,g/255.,b/255.);}

  void setdash(char *patt, int phase=0) {
    // Change line pattern from default (solid). Call with "" to reset.
    fprintf(PLT,"[%s] %d setdash\n",patt,phase);}

  void setlinewidth(double w) {fprintf(PLT,"%g setlinewidth\n",w);}
  // Change linewidth from default (0.5pt).

  void setgray(double w) {fprintf(PLT,"%g setgray\n",w);}
  // Change gray level from default (0.0=black). Range is 0 to 1.

  void gsave() {fprintf(PLT,"gsave\n");}
  // Do a PostScript gsave.

  void grestore() {fprintf(PLT,"grestore\n");}
  // Do a PostScript grestore.  Restores settings prior to last gsave.

  void rawps(char *text) {fprintf(PLT,"%s\n",text);}
  // Emit a string to the plotfile as raw PostScript.

  void addtext(char *text) { fprintf(PLT,"(%s) show ",text); }
  // Plot text at current location.

  void puttext(char *text, double x, double y, double rot=0.0) {
    // Plot text at page location x, y (in pts) at angle rot.
    fprintf(PLT,"gsave %g %g translate %g rotate 0 0 mt ",x,y,rot);
    addtext(text);
    fprintf(PLT,"grestore \n");
  }

  void putctext(char *text, double x, double y, double rot=0.0) {
    // Plot horizontally-centred text at page location
    // x, y (in pts) at angle rot.
    fprintf(PLT,"gsave %g %g translate %g rotate 0 0 mt (%s) ",x,y,rot,text);
    fprintf(PLT,"dup stringwidth pop 2 div neg 0 rmoveto show grestore\n");
  }

  void putCtext(char *text, double x, double y, double rot=0.0) {
    // Plot horizontally and vertically centred text at page location
    // x, y (in pts) at angle rot.
    fprintf(PLT,"gsave %g %g translate %g rotate 0 0 mt (%-s) ",x,y,rot,text);
#if 0    // -- Don't know why this won't work well:
    fprintf(PLT,"dup true charpath pathbbox "
	        "3 -1 roll sub 2 div neg "
	        "3  1 roll sub 2 div exch "
	        "rmoveto show grestore\n");
#else    // -- This is also imperfect, but better:
    fprintf(PLT,"dup stringwidth 2 div neg exch 2 div neg exch "
	        "rmoveto show grestore\n");
#endif    
  }

  void putrtext(char *text, double x, double y, double rot=0.0) {
    // Plot right-justified  text at page location x, y (in pts) at angle rot.
    fprintf(PLT,"gsave %g %g translate %g rotate 0 0 mt (%s) ",x,y,rot,text);
    fprintf(PLT,"dup stringwidth pop neg 0 rmoveto show grestore\n");
  }

  void close() {fprintf(PLT,"showpage\n"); fclose(PLT); PLT = NULL;}
  // Close the plot file.  Called automatically by destructor.
  
  void display() {
    // Start external process to display the plot file.
    // Here we assume gv (ghostview) is in your path.
    char cmd[128];
    if (PLT) close();
    strcpy(cmd,"gv ");
    strcat(cmd,file);
    system(cmd);
  }
  
  void display(char* prog) {
    // Start external process to display the plot file,
    // assuming prog (e.g. gv) is in your path.
    char cmd[128];
    if (PLT) close();
    strcpy(cmd, prog);
    strcat(cmd," ");
    strcat(cmd,file);
    system(cmd);
  }
  
  void display(char* prog, char* fname) {
    // Start external process to display the plot file,
    // assuming prog (e.g. gv) is in your path.
    // It is also assumed you already closed file fname.
    char cmd[128];
    strcpy(cmd, prog);
    strcat(cmd," ");
    strcat(cmd,fname);
    system(cmd);
  }

  void pointsymbol(double x, double y, int num, double size) {
    // Plot Zapf Dingbat symbol num in page coordinates with specified size.
    fprintf(PLT,"(\\%03o) %g %g %g zp\n",num,x-0.394*size,y-0.343*size,size);
  }

  void lineseg(double xs, double ys,
	       double xf, double yf) {
    // Draw a line segment in page coordinates (pts).
    fprintf(PLT,"np %g %g mt %g %g lt st\n",xs,ys,xf,yf);
  }

  void polyline(vector<double>& x            ,
		vector<double>& y            ,
		bool            close = false,
		bool            fill  = false,
		bool            clip  = false)
  // Draw connected line segments in page coordinates (pts), with options
  // to close and/or fill the curve, or set the curve as a clip area.
  {

    int       i;
    const int n = MIN (x.size(), y.size());

    fprintf (PLT,"np %g %g mt\n", x[0], y[0]);
    for (i = 1; i < n; i++)
      fprintf (PLT,"%g %g lt\n",  x[i], y[i]);

    if (close || fill || clip) fprintf(PLT,"cp ");
    
    if      (fill) fprintf (PLT,"fi\n");
    else if (clip) fprintf (PLT,"clip\n");
    else           fprintf (PLT,"st\n");
  }

};

///////////////////////////////////////////////////////////////////////////////

struct PSplot : PSpage {
  // Object that represents an x, y plot box on the page.  Note that a
  // PSplot object can call all the methods of its PSpage.  It overloads
  // many of these methods with versions taking x, y, user coordinates
  // instead of p, q page coordinates.

  double         pll,qll,pur,qur;
  double         xll,yll,xur,yur;
  vector<double> xbox,ybox;
  double         majticsz,minticsz;

  PSplot(PSpage &page, double ppll, double ppur, double qqll, double qqur)
    // Constructor.  Bind to page, with the specified p, q, page coordinates
    // (measured in pts) for the lower left and upper-right corners.
    : pll(ppll), qll(qqll), pur(ppur), qur(qqur),
    xll(ppll), yll(qqll), xur(ppur), yur(qqur), xbox(4), ybox(4),
    majticsz(8.), minticsz(4.)
    {
      strcpy(fontname,page.fontname);
      fontsize = page.fontsize;
      setlimits(xll,xur,yll,yur);
    }

  double p (double x) { return pll + (pur-pll)*(x-xll)/(xur-xll); }
  double q (double y) { return qll + (qur-qll)*(y-yll)/(yur-yll); }
  // Functions returning page coordinates p, q (in pts) from
  // user plot coordinates x, y.
	
  void setlimits(double xxll, double xxur, double yyll, double yyur)
  // Set x and y user values for the lower-left and upper-right
  // corners of the plot object.  Always required.
  {
    xbox[0] = xbox[3] = xll = xxll; ybox[0] = ybox[1] = yll = yyll;
    xbox[1] = xbox[2] = xur = xxur; ybox[2] = ybox[3] = yur = yyur;
  }
	
  void lineseg(double xs, double ys, double xf, double yf)
  // Draw line segment using user coordinates.
  { PSpage::lineseg(p(xs),q(ys),p(xf),q(yf)); }

  void polyline(vector<double>& x            ,
		vector<double>& y            ,
		bool            close = false,
		bool            fill  = false,
		bool            clip  = false)
  // Draw connected line segments using user coordinates.
  // See PSpage::polyline for meaning of options.
  {
    int             i;
    vector<double> xx(x), yy(y);
    for (i=0;i<x.size();i++) xx[i] = p(x[i]);
    for (i=0;i<y.size();i++) yy[i] = q(y[i]);
    PSpage::polyline(xx,yy,close,fill,clip);
  }

  void polyline(int     n            ,
		double* x            ,
		int     incx         ,
		double* y            ,
		int     incy         ,
		bool    close = false,
		bool    fill  = false,
		bool    clip  = false)
  // Draw connected line segments using user coordinates.
  // Incx/y work the same as in BLAS (if negative we start at high end).    
  // See PSpage::polyline for meaning of trailing options.
  {
    int            i;
    vector<double> xx(n), yy(n);
    x += (incx<0) ? (-n+1)*incx : 0;
    y += (incy<0) ? (-n+1)*incy : 0;
    
    for (i = 0; i < n; i++) {
      xx[i] = p(x[i*incx]);
      yy[i] = q(y[i*incy]);
    }
    PSpage::polyline(xx,yy,close,fill,clip);
  }
	
  void lineplot(vector<double>& x,
		vector<double>& y)
  // Plot a curve from x and y vectors of points.
  { polyline(x,y); }
	
  void dot(double x, double y, double size=2.)
  // Plot a filled circle at the specified location, by default small.
  { PSpage::pointsymbol(p(x),q(y),108,size); }
	
  void pointsymbol(double x, double y, int num, double size)
  // Plot a Zapf Dingbat symbol at the specified location and size (in pts).
  { PSpage::pointsymbol(p(x),q(y),num,size); }
	
  void frame()
  // Draw a frame around this plot.
  { polyline(xbox,ybox,true); }
	
  void clear()
  // Erase the interior of this plot.
  { gsave(); setgray(1.); polyline(xbox,ybox,true,true); grestore(); }
	
  void clip()
  // Set interior of this plot as a clip area.
  { gsave(); polyline(xbox,ybox,true,false,true); }
	
  void clip(vector<double>& x,
	    vector<double>& y)
  // Set a clip area from x and y vectors of points.
  { gsave(); polyline(x,y,true,false,true); }
	
  void unclip()
  // Undo previous clip area (or anything else set subsequently.
  { grestore(); }
	
  void xlabel(char* text)
  // Put text label on the x axis.
  { putctext(text,0.5*(pll+pur),qll-2.*fontsize-8.); }
	
  void ylabel(char* text)
  // Put text label on the y axis.
  { putctext(text,pll-3.*fontsize-8.,0.5*(qll+qur),90.); }
	
  void label(char* text, double x, double y, double rot=0.)
  // Put text label at arbitrary location and rotation.
  { puttext(text,p(x),q(y),rot); }
	
  void clabel(char* text, double x, double y, double rot=0.)
  // Put text label horizontally centred
  // at arbitrary location and rotation.
  { putctext(text,p(x),q(y),rot); }
	
  void Clabel(char* text, double x, double y, double rot=0.)
  // Put text label horizontally and vertically centred
  // at arbitrary location and rotation.
  { putCtext(text,p(x),q(y),rot); }
	
  void scalestr(char *str, double x)
  // Format a string for the axis labels.  Used internally.
  { if (abs(x) < 1.e-15) x = 0.;  sprintf(str,"%g",x);
  }
	
  void scales(double xmajd, double xmind,
	      double ymajd, double ymind,
	      int dox=2, int doy=2, int doxx=1, int doyy=1)
  // Draw scales (tick marks) on the plot.  The x and y major and minor
  // division intervals are specified by the first four arguments.
  // The "do" arguments have values 0 (no ticks), 1 (ticks only),
  // 2 (ticks and numbers).  x and xx are the bottom and top sides,
  // y and yy are the left and right sides.
  {
    char   str[128];
    double x,y,xlo,ylo;
     if (dox || doxx) {
      xlo = ceil(MIN(xll,xur)/xmajd)*xmajd;
      for (x=xlo;x<=MAX(xll,xur);x+=xmajd) {
	scalestr(str,x);
	if (dox>1) putctext(str,p(x),qll-fontsize-2.);
	if (dox)  PSpage::lineseg(p(x),qll,p(x),qll+majticsz);
	if (doxx) PSpage::lineseg(p(x),qur,p(x),qur-majticsz);
      }
      xlo = ceil(MIN(xll,xur)/xmind)*xmind;
      for (x=xlo;x<=MAX(xll,xur);x+=xmind) {
	if (dox)  PSpage::lineseg(p(x),qll,p(x),qll+minticsz);
	if (doxx) PSpage::lineseg(p(x),qur,p(x),qur-minticsz);
      }
    }
    if (doy || doyy) {
      ylo = ceil(MIN(yll,yur)/ymajd)*ymajd;
      for (y=ylo;y<=MAX(yll,yur);y+=ymajd) {
	scalestr(str,y);
	if (doy>1) putrtext(str,pll-4.,q(y)-0.3*fontsize);
	if (doy)  PSpage::lineseg(pll,q(y),pll+majticsz,q(y));
	if (doyy) PSpage::lineseg(pur,q(y),pur-majticsz,q(y));
      }
      ylo = ceil(MIN(yll,yur)/ymind)*ymind;
      for (y=ylo;y<=MAX(yll,yur);y+=ymind) {
	if (doy)  PSpage::lineseg(pll,q(y),pll+minticsz,q(y));
	if (doyy) PSpage::lineseg(pur,q(y),pur-minticsz,q(y));
      }
    }
  }
	
  void autoscales() {
    // Draw scales making reasonable default choices.
    double xmajd, xmind, ymajd, ymind;
    xmajd = pow(10.,((int)(log10(abs(xur-xll))-1.1)));
    xmind = xmajd/5.;
    ymajd = pow(10.,((int)(log10(abs(yur-yll))-1.1)));
    ymind = ymajd/5.;
    scales(xmajd,xmind,ymajd,ymind);
  }
};
FILE *PSpage::PLT;
char *PSpage::file;

#endif
