PROGRAM transf_variables
IMPLICIT NONE

INTEGER :: k, m, plev, pplev

pplev = 28
plev = pplev - 1

WRITE(*,*) "k, m, m+1"
DO k=1, plev
   m = plev - k + 1
   WRITE(*,*) k, m, m+1
END DO

END PROGRAM transf_variables
