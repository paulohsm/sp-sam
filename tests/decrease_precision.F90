 PROGRAM decrease_precision
 USE shr_kind_mod, ONLY: r8 => shr_kind_r8, r4 => shr_kind_r4
 IMPLICIT NONE

!integer,parameter :: r8 = shr_kind_r8 !selected_real_kind(12) ! 8 byte real
!integer,parameter :: r4 = shr_kind_r4 !selected_real_kind( 6) ! 4 byte
 integer,parameter :: nx = 32, ny = 1, nz = 27
 integer,parameter :: rl = nx * ny * nz ! record length
 real(kind=r8) :: a(nx,ny,nz)
 real(kind=r4) :: b(nx,ny,nz)
 real(kind=r8),dimension(nx,ny,nz) :: u, v, w

!a = 2.0_r8
!b = a

!print*, a, b

!open(11,file='SEMIPROG_OUT',access='direct',form='unformatted',recl=32*27*4)
!read(11,rec=3) b

!write(*,*) b(1:4,1,1:4)

 a(:,:,:) = sqrt(2.0_r8)
 b(:,:,:) = sqrt(2.0_r4)

 open(21,file='test',access='direct',form='unformatted',recl=rl*8)
 write(21,rec=1)  a
 close(21)

 open(13,file='SEMIPROG_OUT',access='direct',form='unformatted',recl=rl*8)
 read(13,rec=1) u
 read(13,rec=2) v
 read(13,rec=3) w
 write(*,*) w(1:4,1,1:4)

!open(17,file='teste.bin',access='direct',form='unformatted',
!a(:,:,:) = 36.
!write(*,*) a(1:4,1,1:4)


 END PROGRAM
