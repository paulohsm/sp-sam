'open AGCM-Template.ctl'
'reset'

tstep=1200
long=262.5
lati=36.5

'set lon 'long
'set lat 'lati
'set time 0Z18JUN1997 0Z17JUL1997'

'q dims'
lmaxlin=sublin(result,5)
lmax=subwrd(lmaxlin,13)
l=subwrd(lmaxlin,11)
kmaxlin=sublin(result,4)
kmax=subwrd(kmaxlin,13)


'set t 'lmax
'd TOPO'

while(l<=lmax)
   k=1
   while(k<=kmax)
      'set t 'l
      'set z 'k
      'd TEMP'
      'q dims'
      k=k+1
   endwhile
   l=l+1
endwhile

'quit'

***
