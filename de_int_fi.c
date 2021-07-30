#include "de_int_fi.h"

void deft_ctf(double *phi,double *dphi,double *psi,double t,double a,double b);


double deintd_fts(double (*f)(double,void *),void *pa,double a,double omega,double veps,int *erc)
{
  double M,h,alpha,beta,ofa,oa,t,cx,cw,xl,wl,phi,dphi,psi,sig,il,di;
  int l,lb;

  M=-M_PI/eerc*log(veps);
  h=M_PI/M;
  oa=omega*a;
  ofa=oa/M;
  cx=M/omega;
  cw=M_PI/omega;
  beta=1.0/4.0;
  alpha=beta/sqrt(1.0+M*log(1.0+M)/(4.0*M_PI)); 
  
  lb=(int)ceil(ofa/h);
  l=lb;
  sig=pow(-1.0,lb);
  t=(double)l*h-ofa;
  deft_ctf(&phi,&dphi,&psi,t,alpha,beta);
  xl=cx*phi+a;
  wl=(psi>1.0)? cw*sin(M*phi+oa)*dphi : cw*sig*sin(M*psi*phi)*dphi;
  il=f(xl,pa)*wl;

  for(l=lb+1;l<LMAX+lb;l+=1){ // t>0
    sig*=-1.0;
    t=(double)l*h-ofa;
    deft_ctf(&phi,&dphi,&psi,t,alpha,beta);
    xl=cx*phi+a;
    wl=(psi>1.0)? cw*sin(M*phi+oa)*dphi : cw*sig*sin(M*psi*phi)*dphi;
    di=f(xl,pa)*wl;
    il+=di;
    if(fabs(di)<fabs(il)*veps*omega) break;
  }
  if(l==LMAX+lb){
    *erc=-1;
    return 0.0;
  }

  for(l=lb-1;l>-LMAX-lb;l-=1){ // t<0
    t=(double)l*h-ofa;
    deft_ctf(&phi,&dphi,&psi,t,alpha,beta);
    xl=cx*phi+a;
    wl=cw*sin(M*phi+oa)*dphi;
    di=f(xl,pa)*wl;
    il+=di;
    if(fabs(di)<fabs(il)*veps*omega) break;
  }
  if(l==-LMAX-lb){
    *erc=-2;
    return 0.0;
  }
  *erc=0;
  return il;
}
double deintd_ffs(double (*f)(double,void *),void *pa,double a,double omega,double veps,int *erc)
{
  double M,h,alpha,beta,ofa,oa,t,cx,cw,xl,wl,phi,dphi,psi,il,di;
  int l,lb;

  M=-M_PI/eerc*log(veps);
  h=M_PI/M;
  oa=omega*a;
  ofa=oa/M;
  cx=M/omega;
  cw=M_PI/omega;
  beta=1.0/4.0;
  alpha=beta/sqrt(1.0+M*log(1.0+M)/(4.0*M_PI)); 
  
  lb=(int)ceil(ofa/h);
  l=lb;
  t=(double)l*h-ofa;
  deft_ctf(&phi,&dphi,&psi,t,alpha,beta);
  xl=cx*phi+a;
  wl=cw*dphi;
  il=f(xl,pa)*wl;

  for(l=lb+1;l<LMAX+lb;l+=1){ // t>0
    t=(double)l*h-ofa;
    deft_ctf(&phi,&dphi,&psi,t,alpha,beta);
    xl=cx*phi+a;
    wl=cw*dphi;
    di=f(xl,pa)*wl;
    il+=di;
    if(fabs(di)<fabs(il)*veps*omega) break;
  }
  if(l==LMAX+lb){
    *erc=-1;
    return 0.0;
  }

  for(l=lb-1;l>-LMAX-lb;l-=1){ // t<0
    t=(double)l*h-ofa;
    deft_ctf(&phi,&dphi,&psi,t,alpha,beta);
    xl=cx*phi+a;
    wl=cw*dphi;
    di=f(xl,pa)*wl;
    il+=di;
    if(fabs(di)<fabs(il)*veps*omega) break;
  }
  if(l==-LMAX-lb){
    *erc=-2;
    return 0.0;
  }
  *erc=0;
  return il;
}

double deintd_ftc(double (*f)(double,void *),void *pa,double a,double omega,double veps,int *erc)
{
  double M,h,alpha,beta,ofa,oa,t,cx,cw,xl,wl,phi,dphi,psi,sig,il,di;
  int l,lb;

  M=-M_PI/eerc*log(veps);
  h=M_PI/M;
  oa=omega*a;
  ofa=oa/M+M_PI/(2.0*M);
  cx=M/omega;
  cw=M_PI/omega;
  beta=1.0/4.0;
  alpha=beta/sqrt(1.0+M*log(1.0+M)/(4.0*M_PI));

  lb=(int)ceil(ofa/h);
  l=lb;
  sig=pow(-1.0,lb);
  t=(double)l*h-ofa;
  deft_ctf(&phi,&dphi,&psi,t,alpha,beta);
  xl=cx*phi+a;
  wl=(psi>1.0)? cw*cos(M*phi+oa)*dphi : cw*sig*sin(M*psi*phi)*dphi;
  il=f(xl,pa)*wl;

  for(l=lb+1;l<LMAX+lb;l+=1){ // t>0
    sig*=-1.0;
    t=(double)l*h-ofa;
    deft_ctf(&phi,&dphi,&psi,t,alpha,beta);
    xl=cx*phi+a;
    wl=(psi>1.0)? cw*cos(M*phi+oa)*dphi : cw*sig*sin(M*psi*phi)*dphi;
    di=f(xl,pa)*wl;
    il+=di;
    if(fabs(di)<fabs(il)*veps*omega) break;
  }
  if(l==LMAX+lb){
    *erc=-1;
    return 0.0;
  }

  for(l=lb-1;l>-LMAX-lb;l-=1){ // t<0
    t=(double)l*h-ofa;
    deft_ctf(&phi,&dphi,&psi,t,alpha,beta);
    xl=cx*phi+a;
    wl=cw*cos(M*phi+oa)*dphi;
    di=f(xl,pa)*wl;
    il+=di;
    if(fabs(di)<fabs(il)*veps*omega) break;
  }
  if(l==-LMAX-lb){
    *erc=-2;
    return 0.0;
  }
  *erc=0;
  return il;
}

double deintd_ffc(double (*f)(double,void *),void *pa,double a,double omega,double veps,int *erc)
{
  double M,h,alpha,beta,ofa,oa,t,cx,cw,xl,wl,phi,dphi,psi,il,di;
  int l,lb;

  M=-M_PI/eerc*log(veps);
  h=M_PI/M;
  oa=omega*a;
  ofa=oa/M+M_PI/(2.0*M);
  cx=M/omega;
  cw=M_PI/omega;
  beta=1.0/4.0;
  alpha=beta/sqrt(1.0+M*log(1.0+M)/(4.0*M_PI));

  lb=(int)ceil(ofa/h);
  l=lb;
  t=(double)l*h-ofa;
  deft_ctf(&phi,&dphi,&psi,t,alpha,beta);
  xl=cx*phi+a;
  wl=cw*dphi;
  il=f(xl,pa)*wl;

  for(l=lb+1;l<LMAX+lb;l+=1){ // t>0
    t=(double)l*h-ofa;
    deft_ctf(&phi,&dphi,&psi,t,alpha,beta);
    xl=cx*phi+a;
    wl=cw*dphi;
    di=f(xl,pa)*wl;
    il+=di;
    if(fabs(di)<fabs(il)*veps*omega) break;
  }
  if(l==LMAX+lb){
    *erc=-1;
    return 0.0;
  }
  
  for(l=lb-1;l>-LMAX-lb;l-=1){ // t<0
    t=(double)l*h-ofa;
    deft_ctf(&phi,&dphi,&psi,t,alpha,beta);
    xl=cx*phi+a;
    wl=cw*dphi;
    di=f(xl,pa)*wl;
    il+=di;
    if(fabs(di)<fabs(il)*veps*omega) break;
  }
  if(l==-LMAX-lb){
    *erc=-2;
    return 0.0;
  }
  *erc=0;
  return il;
}

double complex deintz_ffs(double complex (*f)(double,void *),void *pa,double a,double omega,double veps,int *erc)
{
  double complex il,di;
  double M,h,alpha,beta,ofa,oa,t,cx,cw,xl,wl,phi,dphi,psi;
  int l,lb;

  M=-M_PI/eerc*log(veps);
  h=M_PI/M;
  oa=omega*a;
  ofa=oa/M;
  cx=M/omega;
  cw=M_PI/omega;
  beta=1.0/4.0;
  alpha=beta/sqrt(1.0+M*log(1.0+M)/(4.0*M_PI)); 
  
  lb=(int)ceil(ofa/h);
  l=lb;
  t=(double)l*h-ofa;
  deft_ctf(&phi,&dphi,&psi,t,alpha,beta);
  xl=cx*phi+a;
  wl=cw*dphi;
  il=f(xl,pa)*wl;

  for(l=lb+1;l<LMAX+lb;l+=1){ // t>0
    t=(double)l*h-ofa;
    deft_ctf(&phi,&dphi,&psi,t,alpha,beta);
    xl=cx*phi+a;
    wl=cw*dphi;
    di=f(xl,pa)*wl;
    il+=di;
    if(cabs(di)<cabs(il)*veps*omega) break;
  }
  if(l==LMAX+lb){
    *erc=-1;
    return 0.0;
  }

  for(l=lb-1;l>-LMAX-lb;l-=1){ // t<0
    t=(double)l*h-ofa;
    deft_ctf(&phi,&dphi,&psi,t,alpha,beta);
    xl=cx*phi+a;
    wl=cw*dphi;
    di=f(xl,pa)*wl;
    il+=di;
    if(cabs(di)<cabs(il)*veps*omega) break;
  }
  if(l==-LMAX-lb){
    *erc=-2;
    return 0.0;
  }
  *erc=0;
  return il;
}

double complex deintz_ffc(double complex (*f)(double,void *),void *pa,double a,double omega,double veps,int *erc)
{
  double complex il,di;
  double M,h,alpha,beta,ofa,oa,t,cx,cw,xl,wl,phi,dphi,psi;
  int l,lb;

  M=-M_PI/eerc*log(veps);
  h=M_PI/M;
  oa=omega*a;
  ofa=oa/M+M_PI/(2.0*M);
  cx=M/omega;
  cw=M_PI/omega;
  beta=1.0/4.0;
  alpha=beta/sqrt(1.0+M*log(1.0+M)/(4.0*M_PI));

  lb=(int)ceil(ofa/h);
  l=lb;
  t=(double)l*h-ofa;
  deft_ctf(&phi,&dphi,&psi,t,alpha,beta);
  xl=cx*phi+a;
  wl=cw*dphi;
  il=f(xl,pa)*wl;

  for(l=lb+1;l<LMAX+lb;l+=1){ // t>0
    t=(double)l*h-ofa;
    deft_ctf(&phi,&dphi,&psi,t,alpha,beta);
    xl=cx*phi+a;
    wl=cw*dphi;
    di=f(xl,pa)*wl;
    il+=di;
    if(cabs(di)<cabs(il)*veps*omega) break;
  }
  if(l==LMAX+lb){
    *erc=-1;
    return 0.0;
  }
  
  for(l=lb-1;l>-LMAX-lb;l-=1){ // t<0
    t=(double)l*h-ofa;
    deft_ctf(&phi,&dphi,&psi,t,alpha,beta);
    xl=cx*phi+a;
    wl=cw*dphi;
    di=f(xl,pa)*wl;
    il+=di;
    if(cabs(di)<cabs(il)*veps*omega) break;
  }
  if(l==-LMAX-lb){
    *erc=-2;
    return 0.0;
  }
  *erc=0;
  return il;
}

//////////////////////////////////////////////////////////////////////
void deft_ctf(double *phi,double *dphi,double *psi,double t,double a,double b)
{
  double et,iet;
  double t1,t2,t3,t4;

  if(fabs(t)>1.0e-9){
    et=exp(t);
    iet=1.0/et;
    
    t1=1.0+t*(2.0+a*iet+b*et);
    t2=exp(-2.0*t-a*(1.0-iet)-b*(et-1.0));
    t3=1.0-t2;

    *psi=t2;
    *phi=t/t3;
    if(t2<1.0e308) *dphi=(1.0-t1*t2)/(t3*t3);
    else *dphi=0.0;
  }
  else {
    et=exp(t);
    iet=1.0/et;
    *psi=exp(-2.0*t-a*(1.0-iet)-b*(et-1.0));
    
    t1=2.0+a+b;
    t2=1.0/t1;
    t3=(a+1.0)*t2*t2-0.5*t2+0.5;
    t4=2.0*(a+1.0)*(a+1.0)*t2*t2*t2-2.0/3.0*(3.0*a+2.0)*t2*t2+1.0/6.0*t2+t1/6.0;
    *phi=t2+t3*t;
    *dphi=t3+t4*t;
  }
}
