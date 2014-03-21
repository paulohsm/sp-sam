'open Exparm_20/GFCTNMC1997061823F.fct.T062L28.ctl'
'reset'

tstep=900

'q dims'
longlin=sublin(result,2)
latilin=sublin(result,3)
kmaxlin=sublin(result,4)
long=subwrd(longlin,6)
lati=subwrd(latilin,6)
kmax=subwrd(kmaxlin,13)

'set t 2'
'set x 1'
'set y 1'
'set z 1'
'd topo'
topo=subwrd(result,4)
'd pslc'
pslc=subwrd(result,4)
'd ocis'
ocis=subwrd(result,4)
'd oces'
oces=subwrd(result,4)
'd iswf'
iswf=subwrd(result,4)
'd roce'
roce=subwrd(result,4)
'd role'
role=subwrd(result,4)
'd swtc'
swtc=subwrd(result,4)
'd ocic'
ocic=subwrd(result,4)
'd lwbc'
lwbc=subwrd(result,4)
'd lsmk'
lsmk=subwrd(result,4)
'd cssf'
cssf=subwrd(result,4)
'd clsf'
clsf=subwrd(result,4)
'd usst'
usst=subwrd(result,4)
'd vsst'
vsst=subwrd(result,4)
'd tsfc'
tsfc=subwrd(result,4)

write('SEMIPROG_IN',kmax'	'tstep'	'long'	'lati'	'topo'	'pslc'	'ocis'	'oces'	'iswf'	'roce)
write('SEMIPROG_IN',role'	'swtc'	'ocic'	'lwbc'	'lsmk'	'cssf'	'clsf'	'usst'	'vsst'	'tsfc,append)

*write('var_prof.txt',zmax)

k=1
while(k<=kmax)
   'set z 'k
   'set t 2'
   'set x 1'
   'set y 1'
   'q dims'
   zlevlin=sublin(result,4)
   zlev=subwrd(zlevlin,6)
   'd TEMP'
   temp=subwrd(result,4)
   'd UMES'
   umes=subwrd(result,4)
   'd LIQM'
   liqm=subwrd(result,4)
   'd ICEM'
   icem=subwrd(result,4)
   'd UVEL'
   uvel=subwrd(result,4)
   'd VVEL'
   vvel=subwrd(result,4)
   'd SWRH'
   swrh=subwrd(result,4)
   'd LWRH'
   lwrh=subwrd(result,4)
   vars=write('SEMIPROG_IN',zlev'	'temp'	'umes'	'liqm'	'icem'	'uvel'	'vvel'	'swrh'	'lwrh,append)
   k=k+1
endwhile

***
