#!/usr/bin/env gnuplot
f1(x)=a1*exp(-0.5*x**2*b1**2);
f2(x)=a1*exp(-0.5*x**2*b1**2)+\
      a2*exp(-0.5*x**2*b2**2);
f3(x)=a1*exp(-0.5*x**2*b1**2)+\
      a2*exp(-0.5*x**2*b2**2)+\
      a3*exp(-0.5*x**2*b3**2);
f4(x)=a1*exp(-0.5*x**2*b1**2)+\
      a2*exp(-0.5*x**2*b2**2)+\
      a3*exp(-0.5*x**2*b3**2)+\
      a4*exp(-0.5*x**2*b4**2);
set fit quiet;
set fit logfile "/dev/null";
set print "-";
set fit errorvariables;
b1=10.5**0.5*2;
a1=0.5*sqrt(2*pi/b1**2);
b2=26**0.5*2;
a2=0.3*sqrt(2*pi/b2**2);
b3=3.1**0.5*2;
a3=0.14*sqrt(2*pi/b3**2);
b4=58**0.5*2;
a4=0.04*sqrt(2*pi/b4**2);
FIT_LIMIT=1e-8;
do for [i=2:212]{
#a1=i;
#a2=i/2.0;
#a3=i/4.0;
#a4=i/8.0;
#b1=1;
#b2=1;
#b3=1;
#b4=1;
  fit f1(x) "<tail -n+2 diffraction_table2" u 1:i via a1,b1; 
  print a1,abs(b1),FIT_STDFIT;
#fit f2(x) "<tail -n+2 diffraction_table2" u 1:i via a1,b1,a2,b2; 
#print a1,abs(b1),a2,abs(b2),FIT_STDFIT;
#fit f3(x) "<tail -n+2 diffraction_table2" u 1:i via a1,b1,a2,b2,a3,b3; 
#print a1,abs(b1),a2,abs(b2),a3,abs(b3),FIT_STDFIT;
#fit f4(x) "<tail -n+2 diffraction_table2" u 1:i via a1,b1,a2,b2,a3,b3,a4,b4;
#print a1,abs(b1),a2,abs(b2),a3,abs(b3),a4,abs(b4),FIT_STDFIT;
}
