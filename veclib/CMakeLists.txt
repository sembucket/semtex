set (vec_lib_src
  util
  memory
  xcopy   xfill   xneg    xvneg   xsadd   xvadd   xssub
  xvsub   xsmul   xvmul   xsdiv   xvrecp  xvdiv   xzero
  xspow   xvsgn   xznan
  xvabs   xvamax  xvpow   xvexp   xvlg10  xvlog   xvatan
  xvatn2  xvcos   xvsin   xvsqrt  xrand   xvhypot xvtanh
  xsvmvt  xsvpvt  xsvtsp  xsvtvm  xsvtvp  xsvvmt  xsvvpt
  xsvvtm  xsvvtp  xvvmvt  xvvpvt  xvvtvm  xvvtvp  xvvvtm
  xsvvtt  xvvvtt
  xseq    xsge    xsle    xslt    xsne
  icount  ifirst  ixmax   ixmin   xsum    lany    lxsame
  vdble       xvfloa         vsngl            xbrev
  xscatr      xgathr         xgathr_scatr
  xscatr_sum  xgathr_sum     xgathr_scatr_sum
  xramp       xcndst         xmask
  xvpoly      xpolint        xclip
  xspline     xsplquad       xiclip
  xmxv        xmxva
  xvvtvvtp    xsvvttvp       xvvtvvtm
  xmxm.f
)

add_library (vec STATIC ${vec_lib_src})
