'open AGCM-Template.ctl'
'reset'

tstep=1200
long=262.5
lati=36.5
kmax=28

'set lon 'long
'set lat 'lati
'set z 1 'kmax
'set time 0Z18JUN1997 0Z17JUL1997'

'q dims'
lmaxlin=sublin(result,5)
lmax=subwrd(lmaxlin,13)
l=subwrd(lmaxlin,11)
*kmaxlin=sublin(result,4)
*kmax=subwrd(kmaxlin,13)
*kmax=28

'set t 1'
'set z 1'
'd TOPO'
topo=subwrd(result,4)
'd LSMK'
lsmk=subwrd(result,4)

* model parameters and constants
write('SEMIPROG_IN',kmax'	'tstep'	'long'	'lati'	'topo'	'lsmk)

* model time varying fields
while(l<=lmax)
   'set lon 'long
   'set lat 'lati
   'set z 1'
   'set t 'l
   'q dims'
   tlevlin=sublin(result,5)
   tlev=subwrd(tlevlin,6)
   'd PSLC'
   pslc=subwrd(result,4)
   'd USST'
   usst=subwrd(result,4)
   'd VSST'
   vsst=subwrd(result,4)
   'd CSSF'
   cssf=subwrd(result,4)
   'd CLSF'
   clsf=subwrd(result,4)
   'd OCIS'
   ocis=subwrd(result,4)
   'd OCES'
   oces=subwrd(result,4)
   'd ISWF'
   iswf=subwrd(result,4)
   'd ROCE'
   roce=subwrd(result,4)
   'd OLIS'
   olis=subwrd(result,4)
   'd OLES'
   oles=subwrd(result,4)
   'd ROLE'
   role=subwrd(result,4)
   'd SWTC'
   swtc=subwrd(result,4)
   'd OCIC'
   ocic=subwrd(result,4)
   'd LWTC'
   lwtc=subwrd(result,4)
   'd LWBC'
   lwbc=subwrd(result,4)
   write('SEMIPROG_IN',tlev'	'pslc'	'usst'	'vsst'	'cssf'	'clsf'	'ocis'	'oces'	'iswf,append)
   write('SEMIPROG_IN',roce'	'olis'	'oles'	'role'	'swtc'	'ocic'	'lwtc'	'lwbc,append)
* model time-vertical varying fields
   k=1
   while(k<=kmax)
      'set t 'l
      'set lon 'long
      'set lat 'lati
      'set z 'k
      'q dims'
      zlevlin=sublin(result,4)
      zlev=subwrd(zlevlin,6)
      'd TEMP'
      temp=subwrd(result,4)
      'd UMES'
      umes=subwrd(result,4)
*     'd LIQM'
*     liqm=subwrd(result,4)
      liqm=0
*     'd ICEM'
*     icem=subwrd(result,4)
      icem=0
      'd UVEL'
      uvel=subwrd(result,4)
      'd VVEL'
      vvel=subwrd(result,4)
      'd SWRH'
      swrh=subwrd(result,4)
      'd LWRH'
      lwrh=subwrd(result,4)
      write('SEMIPROG_IN',zlev'	'temp'	'umes'	'liqm'	'icem'	'uvel'	'vvel'	'swrh'	'lwrh,append)
      k=k+1
   endwhile
   l=l+1
endwhile

'quit'

***
