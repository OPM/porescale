!  Copyright 2015 Pierre Adler
!
!  This file is part of the Open Porous Media project (OPM).
!
!  OPM is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  OPM is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Publi! License
!  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

      program latticespringmodel
c___________________________________________________
c    LSM - Lattice spring model
c___________________________________________________
c    include '/usr/include/fexcp.h'      ! for debug
!___________________________________________________
! Son of lsm_elastic_element_shift_10.05.10.f
!    and lsm_elastic_element_shift_03.05.10.f
! just the complete version of both  

c Input files:
! INPUT_PARAMETERS - program parameters 
! ind - geometry of the medium (1 - pore, 0 - solid)
! velocities19 - discrete velocities
c Output files:
! lsm_settings - input program settings and parameters
! coord_field - coordinates of solid points
! displ_field - displacements of solid points
! velocity_field - velocities of solid points
! force_field - forces acting on solid points

      implicit  real*8 (a-h,o-z)
      implicit  integer*8 (i-n)

c      parameter ( nthread = 1)      ! NUMBER OF THREADS
      integer   nthread
      integer   nthread2            ! = nthread*2

      integer*4 ncx ! x unit cell size
      integer*4 ncy ! y unit cell size
      integer*4 ncz ! z unit cell size

      integer*4 ns  ! number of solid points
c      integer*4 nnc ! number of connected solid points
c      integer*4 nnb ! number of bounday solid points
      integer*4 nsc ! dimension of isc
      integer*4 nsb ! dimension of isb

      real*8    vv(3,19) ! discrete velocities (real numbers)
      real*8    vn(3,19) ! discrete velocities normalized
      real*8    vx(3,19) ! xn/|rn|
      real*8    vnr(19) ! discrete velocities norms
      real*8    vni(19) ! discrete velocities norms inverse
      integer   iv(3,19) ! discrete velocities (integer numbers)
      integer   niv(19) ! indices of the opposite velocities

      integer   ib(3,9) ! bond vectors (9 instead of 18 for uniqueness) 
      real*8    vb(3,9) ! bond vectors (real numbers)
      real*8    vbn(9) ! norms
      real*8    vbi(9) ! norms inverse
      real*8    bn(3,9) ! normalized vectors
      real*8    bx(3,9) ! xn/|rn|

      real*8 rs ! solid body density
      real*8 rsi ! inverse of solid body density
      real*8 alf1 ! alpha linear spring constant of 1st solid component 
      real*8 alf2 ! alpha linear spring constant of 1st solid component           
      real*8 bet1 ! beta angular spring constant
      real*8 bet2 ! beta angular spring constant
      real*8 Em1 ! Young modulus
      real*8 nu1 ! Poisson's ratio
      real*8 Em2 ! Young modulus
      real*8 nu2 ! Poisson's ratio
      integer*4 itm ! simulation duration
      integer*4 its ! save parameter

      integer*4, dimension(:,:), allocatable :: bon ! bonds
      integer*4, dimension(:,:), allocatable :: ang ! angles

      integer*4, dimension(:,:), allocatable :: icors ! coordinates of the solid points
      integer*4, dimension(:,:,:), allocatable :: mijks ! point numbers
      real*8,    dimension(:,:,:), allocatable :: sigma   
      integer*4, dimension(:,:), allocatable :: isc ! solid medium connectivity table
      integer*4, dimension(:,:), allocatable :: isb ! solid medium boundary table
      integer*4, dimension(:), allocatable :: nsn ! intervals for number of solid neighbors
      integer*4, dimension(:), allocatable :: nbn ! intervals for number of boundary neighbors
  
      real*8, dimension(:,:), allocatable :: adat ! elastic solid points acceleration
      real*8, dimension(:,:), allocatable :: vdat ! elastic solid points velocity
      real*8, dimension(:,:), allocatable :: rdat ! elastic solid points coordinates

      integer, dimension(:,:,:), allocatable :: ind  ! porous medium structure
      integer, dimension(:,:,:), allocatable :: inda ! all solid points indicator matrix
      integer, dimension(:,:,:), allocatable :: inds ! solid point indicator matrix
      integer, dimension(:), allocatable :: ipcx ! x-axis periodic boundary conditions
      integer, dimension(:), allocatable :: ipcy ! y-axis periodic boundary conditions
      integer, dimension(:), allocatable :: ipcz ! z-axis periodic boundary conditions
      integer, dimension(:), allocatable :: iicx ! x-axis periodic boundary conditions
      integer, dimension(:), allocatable :: iicy ! y-axis periodic boundary conditions
      integer, dimension(:), allocatable :: iicz ! z-axis periodic boundary conditions
      real*8, dimension(:), allocatable :: sbond1 ! linear bond constants
      real*8, dimension(:), allocatable :: sbond2 ! linear bond constants
      real*8, dimension(:), allocatable :: sbond  ! linear bond constants
      real*8, dimension(:), allocatable :: sbondin ! linear bond constants in cube
      real*8, dimension(:), allocatable :: abond1 ! angular bond constants
      real*8, dimension(:), allocatable :: abond2 ! angular bond constants
      real*8, dimension(:), allocatable :: abond  ! angular bond constants
      real*8, dimension(:), allocatable :: abondin ! angular bond constants in cube
     
      real*8 gstrain
c      real*8, dimension(:,:), allocatable :: forces ! global boundary conditions

c      real*8  frc(3)       ! local force
c      real*8, dimension(:,:), allocatable :: fmacro ! global forces
c      integer*4, dimension(:), allocatable :: fmind ! global forces index
c      integer nfpoint ! number of force points
c      real*8 vlayer ! volume of elemmentary cubes in the layer for global conditions
c___________________________________________________
c     technical variables
c___________________________________________________      
      character*8  date ! current date
      character*10 time ! current time
      character*30 format211 ! format for pgf95

      integer checkstat ! check the success of allocate command
      integer*4 nx,ny,nz ! variables for point coordinates
      integer*4 iz(1:3) ! technical variable for comparation 
      integer*4 npc(1:3) ! coordinates of adjacent point 
      integer*4 npc2(1:3) ! coordinates of adjacent point
      integer*4 npc3(1:3) ! coordinates of adjacent point
      real*8 xn(3),xr(3),rn(3),rt(3) ! vectors in linear spring force calculation
      real*8 anf(9,48) ! angle plane normal, angular force vectors j k; for angular springs
      real*8 c1ds2,c1ds3,c1ds6,cs2d3 ! numerical constants
      real*8 c1d2,c1d2s3
      real*8 RM(3,3) ! rotation matrix
      real*8 w ! frequency
      integer*8 nbo,nna ! bonds counter
      integer nta(2,48) ! angles table, angles numeration is the same as in anf
      integer*4 npj(3),npk(3),npl(3) ! angles counter, angleside endpoints coordinates
      real*8 flin(3) ! linear spring force
      real*8 dfi ! delta phi, angle change for angular force
      real*8 xnj(3),xrj(3),rnj(3),rtj(3),xnk(3),xrk(3),rnk(3),rtk(3) ! force calculation vectors
      real*8 RTO(3),VTO(3),FTO(3)
      real*8 sigm(3,3),sigmt(3,3) 
      real*8 FX(3),FY(3),FZ(3),FXJ(3),FXK(3)
      real*8 F20(1:3,20)
      real*8 F20Y(1:3,20)
      real*8 F20Z(1:3,20)
      real*8 scr(3),scrj(3),scrk(3)
      integer ich100(3,9) ! check 100 type elastic elements to calculate spring constant
      integer ich110(6,12) ! check 110 type elastic elements to calculate spring constant
      integer*4 jv(3),kv(3),ip0(3),mp(3)
      integer itab199(19) ! correspondence table 9-19

c      integer layers(2*nthread)  ! layers width in ncxXncyXncz domain 
c      integer*4 layb(2*nthread+1)  ! layers boundaries in 1,ns range
c      integer*4 lib(2*nthread+1),lia(2*nthread+1) ! layers boundaries in bon and ang massives
      integer, dimension(:), allocatable :: layers
      integer*8, dimension(:), allocatable :: layb
      integer*8, dimension(:), allocatable :: lib
      integer*8, dimension(:), allocatable :: lia

      real*8 sumloc(3)
      integer icon
c      include 'omp_lib.h'
c      integer(OMP_LOCK_KIND) LCK
c      integer(OMP_LOCK_KIND) MCK
c      real*8 gsigma(3,3) ! global sigma tensor
c___________________________________________________
c read geometrical structure, initialize all the data
c___________________________________________________
c      gsigma(1,1:3) = (/0.01D0/27.D0,0.D0         ,0.D0/) ! global stress tensor
c      gsigma(2,1:3) = (/   0.D0     ,0.005D0/54.D0,0.D0/)
c      gsigma(3,1:3) = (/   0.D0     ,0.D0         ,0.005D0/54.D0/)
c      nthread2=2*nthread
      itab199=(/0,1,1,2,2,3,3,4,5,5,4,6,7,7,6,8,9,9,8/)
      ich100(1,1:9)=(/5,19,7,18,4,16,6,17,5/)
      ich100(2,1:9)=(/7,14,2,12,6,13,3,15,7/)
      ich100(3,1:9)=(/2,8,4,9,3,11,5,10,2/)

      ich110(1,1:12)=(/6,12,2,3,13,6,7,14,2,3,15,7/)
      ich110(2,1:12)=(/7,14,2,7,15,3,2,12,6,6,13,3/)
      ich110(3,1:12)=(/4,16,6,4,18,7,6,17,5,5,19,7/)
      ich110(4,1:12)=(/4,18,7,4,16,6,5,19,7,5,17,6/)
      ich110(5,1:12)=(/6,12,2,2,14,7,6,13,3,3,15,7/)
      ich110(6,1:12)=(/7,14,2,2,12,6,3,15,7,3,13,6/)

      c1ds2=1.d0/sqrt(2.d0) ! some numerical constants
      c1ds3=1.d0/sqrt(3.d0)
      c1ds6=1.d0/sqrt(6.d0)
      cs2d3=sqrt(2.d0/3.d0)
      c1d2=1.d0/2.d0
      c1d2s3=1.d0/(2.d0*sqrt(3.d0))
! 1:3 plane normal n, 4:6 fj, 7:9 fk
! xy plane
      anf(1:9,1)=(/0.d0,0.d0,1.d0,0.d0,1.d0,0.d0
     &                           ,c1d2,-c1d2,0.d0/) ! 2 8
      anf(1:9,2)=(/0.d0,0.d0,1.d0,-c1d2,c1d2,0.d0
     &                           ,1.d0,0.d0,0.d0/) ! 8 4
      anf(1:9,3)=(/0.d0,0.d0,1.d0,-1.d0,0.d0,0.d0
     &                           ,c1d2,c1d2,0.d0/) ! 4 9
      anf(1:9,4)=(/0.d0,0.d0,1.d0,-c1d2,-c1d2,0.d0
     &                           ,0.d0,1.d0,0.d0/) ! 9 3
      anf(1:9,5)=(/0.d0,0.d0,1.d0,0.d0,-1.d0,0.d0
     &                           ,-c1d2,c1d2,0.d0/) ! 3 11
      anf(1:9,6)=(/0.d0,0.d0,1.d0,c1d2,-c1d2,0.d0
     &                           ,-1.d0,0.d0,0.d0/) ! 11 5
      anf(1:9,7)=(/0.d0,0.d0,1.d0,1.d0,0.d0,0.d0
     &                           ,-c1d2,-c1d2,0.d0/) ! 5 10
      anf(1:9,8)=(/0.d0,0.d0,1.d0,c1d2,c1d2,0.d0
     &                           ,0.d0,-1.d0,0.d0/) ! 10 2
! xz plane
      anf(1:9,9)=(/0.d0,1.d0,0.d0,0.d0,0.d0,-1.d0
     &                           ,-c1d2,0.d0,c1d2/) ! 2 12
      anf(1:9,10)=(/0.d0,1.d0,0.d0,c1d2,0.d0,-c1d2
     &                           ,-1.d0,0.d0,0.d0/) ! 12 6
      anf(1:9,11)=(/0.d0,1.d0,0.d0,1.d0,0.d0,0.d0
     &                           ,-c1d2,0.d0,-c1d2/) ! 6 13
      anf(1:9,12)=(/0.d0,1.d0,0.d0,c1d2,0.d0,c1d2
     &                           ,0.d0,0.d0,-1.d0/) ! 13 3
      anf(1:9,13)=(/0.d0,1.d0,0.d0,0.d0,0.d0,1.d0
     &                           ,c1d2,0.d0,-c1d2/) ! 3 15
      anf(1:9,14)=(/0.d0,1.d0,0.d0,-c1d2,0.d0,c1d2
     &                           ,1.d0,0.d0,0.d0/) ! 15 7
      anf(1:9,15)=(/0.d0,1.d0,0.d0,-1.d0,0.d0,0.d0
     &                           ,c1d2,0.d0,c1d2/) ! 7 14
      anf(1:9,16)=(/0.d0,1.d0,0.d0,-c1d2,0.d0,-c1d2
     &                           ,0.d0,0.d0,1.d0/) ! 14 2
! yz plane
      anf(1:9,17)=(/1.d0,0.d0,0.d0,0.d0,0.d0,1.d0
     &                           ,0.d0,c1d2,-c1d2/) ! 4 16
      anf(1:9,18)=(/1.d0,0.d0,0.d0,0.d0,-c1d2,c1d2
     &                           ,0.d0,1.d0,0.d0/) ! 16 6
      anf(1:9,19)=(/1.d0,0.d0,0.d0,0.d0,-1.d0,0.d0
     &                           ,0.d0,c1d2,c1d2/) ! 6 17
      anf(1:9,20)=(/1.d0,0.d0,0.d0,0.d0,-c1d2,-c1d2
     &                           ,0.d0,0.d0,1.d0/) ! 17 5
      anf(1:9,21)=(/1.d0,0.d0,0.d0,0.d0,0.d0,-1.d0
     &                           ,0.d0,-c1d2,c1d2/) ! 5 19
      anf(1:9,22)=(/1.d0,0.d0,0.d0,0.d0,c1d2,-c1d2
     &                           ,0.d0,-1.d0,0.d0/) ! 19 7
      anf(1:9,23)=(/1.d0,0.d0,0.d0,0.d0,1.d0,0.d0
     &                           ,0.d0,-c1d2,-c1d2/) ! 7 18
      anf(1:9,24)=(/1.d0,0.d0,0.d0,0.d0,c1d2,c1d2
     &                           ,0.d0,0.d0,-1.d0/) ! 18 4
! xy1
      anf(1:9,25)=(/-c1ds3,-c1ds3,c1ds3,-c1d2s3,c1ds3,c1d2s3
     &                           ,c1ds3,-c1d2s3,c1d2s3/) ! 12 16
      anf(1:9,26)=(/c1ds3,-c1ds3,c1ds3,-c1ds3,-c1d2s3,c1d2s3
     &                           ,c1d2s3,c1ds3,c1d2s3/) ! 16 13
      anf(1:9,27)=(/c1ds3,c1ds3,c1ds3,c1d2s3,-c1ds3,c1d2s3
     &                           ,-c1ds3,c1d2s3,c1d2s3/) ! 13 17
      anf(1:9,28)=(/-c1ds3,c1ds3,c1ds3,c1ds3,c1d2s3,c1d2s3
     &                           ,-c1d2s3,-c1ds3,c1d2s3/) ! 17 12
! xy-1
      anf(1:9,29)=(/c1ds3,c1ds3,c1ds3,-c1d2s3,c1ds3,-c1d2s3
     &                           ,c1ds3,-c1d2s3,-c1d2s3/) ! 14 18
      anf(1:9,30)=(/-c1ds3,c1ds3,c1ds3,-c1ds3,-c1d2s3,-c1d2s3
     &                           ,c1d2s3,c1ds3,-c1d2s3/) ! 18 15
      anf(1:9,31)=(/-c1ds3,-c1ds3,c1ds3,c1d2s3,-c1ds3,-c1d2s3
     &                           ,-c1ds3,c1d2s3,-c1d2s3/) ! 15 19
      anf(1:9,32)=(/c1ds3,-c1ds3,c1ds3,c1ds3,c1d2s3,-c1d2s3
     &                           ,-c1d2s3,-c1ds3,-c1d2s3/) ! 19 14
! yz1
      anf(1:9,33)=(/c1ds3,-c1ds3,-c1ds3,c1d2s3,-c1d2s3,c1ds3
     &                           ,c1d2s3,c1ds3,-c1d2s3/) ! 8 12
      anf(1:9,34)=(/c1ds3,c1ds3,-c1ds3,c1d2s3,-c1ds3,-c1d2s3
     &                           ,c1d2s3,c1d2s3,c1ds3/) ! 12 10
      anf(1:9,35)=(/c1ds3,c1ds3,c1ds3,c1d2s3,c1d2s3,-c1ds3
     &                           ,c1d2s3,-c1ds3,c1d2s3/) ! 10 14
      anf(1:9,36)=(/c1ds3,-c1ds3,c1ds3,c1d2s3,c1ds3,c1d2s3
     &                           ,c1d2s3,-c1d2s3,-c1ds3/) ! 14 8
! yz-1
      anf(1:9,37)=(/c1ds3,c1ds3,c1ds3,-c1d2s3,-c1d2s3,c1ds3
     &                           ,-c1d2s3,c1ds3,-c1d2s3/) ! 9 13
      anf(1:9,38)=(/c1ds3,-c1ds3,c1ds3,-c1d2s3,-c1ds3,-c1d2s3
     &                           ,-c1d2s3,c1d2s3,c1ds3/) ! 13 11
      anf(1:9,39)=(/c1ds3,-c1ds3,-c1ds3,-c1d2s3,c1d2s3,-c1ds3
     &                           ,-c1d2s3,-c1ds3,c1d2s3/) ! 11 15
      anf(1:9,40)=(/c1ds3,c1ds3,-c1ds3,-c1d2s3,c1ds3,c1d2s3
     &                           ,-c1d2s3,-c1d2s3,-c1ds3/) ! 15 9
! xz1
      anf(1:9,41)=(/c1ds3,-c1ds3,c1ds3,-c1d2s3,c1d2s3,c1ds3
     &                           ,c1ds3,c1d2s3,-c1d2s3/) ! 8 16
      anf(1:9,42)=(/-c1ds3,-c1ds3,c1ds3,-c1ds3,c1d2s3,-c1d2s3
     &                           ,c1d2s3,c1d2s3,c1ds3/) ! 16 9
      anf(1:9,43)=(/-c1ds3,-c1ds3,-c1ds3,c1d2s3,c1d2s3,-c1ds3
     &                           ,-c1ds3,c1d2s3,c1d2s3/) ! 9 18
      anf(1:9,44)=(/c1ds3,-c1ds3,-c1ds3,c1ds3,c1d2s3,c1d2s3
     &                           ,-c1d2s3,c1d2s3,-c1ds3/) ! 18 8
! xz-1
      anf(1:9,45)=(/-c1ds3,-c1ds3,-c1ds3,-c1d2s3,-c1d2s3,c1ds3
     &                           ,c1ds3,-c1d2s3,-c1d2s3/) ! 10 17
      anf(1:9,46)=(/c1ds3,-c1ds3,-c1ds3,-c1ds3,-c1d2s3,-c1d2s3
     &                           ,c1d2s3,-c1d2s3,c1ds3/) ! 17 11
      anf(1:9,47)=(/c1ds3,-c1ds3,c1ds3,c1d2s3,-c1d2s3,-c1ds3
     &                           ,-c1ds3,-c1d2s3,c1d2s3/) ! 11 19
      anf(1:9,48)=(/-c1ds3,-c1ds3,c1ds3,c1ds3,-c1d2s3,c1d2s3
     &                           ,-c1d2s3,-c1d2s3,-c1ds3/) ! 19 10

! pi/4 angles clock counterwise
! xy plane
      nta(1:2,1)=(/2,8/) ! xy
      nta(1:2,2)=(/8,4/) ! xy
      nta(1:2,3)=(/4,9/) ! xy
      nta(1:2,4)=(/9,3/) ! xy
      nta(1:2,5)=(/3,11/) ! xy
      nta(1:2,6)=(/11,5/) ! xy
      nta(1:2,7)=(/5,10/) ! xy
      nta(1:2,8)=(/10,2/) ! xy
! xz plane
      nta(1:2,9)=(/2,12/) ! xz
      nta(1:2,10)=(/12,6/) ! xz
      nta(1:2,11)=(/6,13/) ! xz
      nta(1:2,12)=(/13,3/) ! xz
      nta(1:2,13)=(/3,15/) ! xz
      nta(1:2,14)=(/15,7/) ! xz
      nta(1:2,15)=(/7,14/) ! xz
      nta(1:2,16)=(/14,2/) ! xz
! yz plane
      nta(1:2,17)=(/4,16/) ! yz
      nta(1:2,18)=(/16,6/) ! yz
      nta(1:2,19)=(/6,17/) ! yz
      nta(1:2,20)=(/17,5/) ! yz
      nta(1:2,21)=(/5,19/) ! yz
      nta(1:2,22)=(/19,7/) ! yz
      nta(1:2,23)=(/7,18/) ! yz
      nta(1:2,24)=(/18,4/) ! yz
! pi/3 angles
! comments - sides
      nta(1:2,25)=(/12,16/) ! xy1
      nta(1:2,26)=(/16,13/) ! xy1
      nta(1:2,27)=(/13,17/) ! xy1
      nta(1:2,28)=(/17,12/) ! xy1
      nta(1:2,29)=(/14,18/) ! xy-1
      nta(1:2,30)=(/18,15/) ! xy-1
      nta(1:2,31)=(/15,19/) ! xy-1
      nta(1:2,32)=(/19,14/) ! xy-1

      nta(1:2,33)=(/8,12/) ! yz1
      nta(1:2,34)=(/12,10/) ! yz1
      nta(1:2,35)=(/10,14/) ! yz1
      nta(1:2,36)=(/14,8/) ! yz1
      nta(1:2,37)=(/9,13/) ! yz-1
      nta(1:2,38)=(/13,11/) ! yz-1
      nta(1:2,39)=(/11,15/) ! yz-1
      nta(1:2,40)=(/15,9/) ! yz-1

      nta(1:2,41)=(/8,16/) ! xz1
      nta(1:2,42)=(/16,9/) ! xz1
      nta(1:2,43)=(/9,18/) ! xz1
      nta(1:2,44)=(/18,8/) ! xz1
      nta(1:2,45)=(/10,17/) ! xz-1
      nta(1:2,46)=(/17,11/) ! xz-1
      nta(1:2,47)=(/11,19/) ! xz-1
      nta(1:2,48)=(/19,10/) ! xz-1

      open(1,file='INPUT_PARAMETERS') ! read input parameters
      read(1,*)ncx ! x unit cell size
      read(1,*)ncy ! y unit cell size
      read(1,*)ncz ! z unit cell size
      read(1,*)Em1 ! Young modulus of 1st solid component
      read(1,*)nu1 ! Poisson's ratio of 1st solid component 
      read(1,*)Em2 ! Young modulus of 2nd solid conpomnent
      read(1,*)nu2 ! Poisson's ratio of 2nd solid component
      read(1,*)rs ! solid medium density
      read(1,*)dt ! integration time step
      read(1,*)dl ! lattice unit step (distance between two lattice nodes)
      read(1,*)itm ! simulation duration
      read(1,*)its ! simulation save parameter
      read(1,*)gstrain !read strain
      read(1,*)nthread !read nthread
      read(1,*)icon
!      read(1,*)(bfr(i),i=1,3) ! red fluid body force
!      read(1,*)(bfb(i),i=1,3) ! blue fluid body force
      close(1)

      nthread2=2*nthread

      allocate(layers(2*nthread),stat=checkstat) ! create medium structure array
      if(checkstat>0) then
       write(*,*)'ERROR! impossible to allocate memory for ind'
       stop
      endif

      allocate(layb(2*nthread+1),stat=checkstat) ! create medium structure array
      if(checkstat>0) then
       write(*,*)'ERROR! impossible to allocate memory for ind'
       stop
      endif

      allocate(lib(2*nthread+1),stat=checkstat) ! create medium structure array
      if(checkstat>0) then
       write(*,*)'ERROR! impossible to allocate memory for ind'
       stop
      endif

      allocate(lia(2*nthread+1),stat=checkstat) ! create medium structure array
      if(checkstat>0) then
       write(*,*)'ERROR! impossible to allocate memory for ind'
       stop
      endif

c      ncx=ncx+2
      alf1=dl*Em1/(5.d0-10.d0*nu1) ! linear spring constant
      bet1=dl*dl*dl*Em1*(4*nu1-1)/(20.d0*(2*nu1*nu1+nu1-1)) ! angular spring constant
      alf2=dl*Em2/(5.d0-10.d0*nu2) ! linear spring constant
      bet2=dl*dl*dl*Em2*(4*nu2-1)/(20.d0*(2*nu2*nu2+nu2-1)) ! angular spring constant

      rsi=1.d0/rs

c      alf=0.d0 !!!!!!!!!!!!!!! just angular springs are considered

      open(11,file='velocities9') ! read bond vectors
      do i=1,3
       read(11,129) (ib(i,j),j=1,9) ! read bond vectors from the file
      enddo
      vb = ib ! real numbers version of bond vectors 
129   format(9i3)

      open(11,file='velocities19') ! read discrete velocities
      do i=1,3
       read(11,119) (iv(i,j),j=1,19) ! read discrete velocities from the file 
      enddo
      vv = iv ! real numbers version of discrete velocities
119   format(19i3)

      niv=(/(0,i=1,19)/) ! initialize niv
      do i=1,19 ! define opposite velocities numbers
        do j=1,19
         iz(1:3)=iv(1:3,i)+iv(1:3,j)
         if(iz(1).eq.0.and.iz(2).eq.0.and.iz(3).eq.0) then
          niv(i)=j
         endif 
        enddo
      enddo
      close(11)

      allocate(ind(ncx,ncy+1,ncz),stat=checkstat) ! create medium structure array
      if(checkstat>0) then
       write(*,*)'ERROR! impossible to allocate memory for ind'
       stop
      endif
      ind=1

      open(12,file='ind') ! open medium structure index file (1-pore, 0-solid)  
      write(format211,*)'(" ",',ncx,'i1)' ! create format to read from 'ind'
      do k=1,ncz
       do j=1,ncy
c        read(12,211)(ind(i,j,k),i=1,ncx) ! read structure data from the file
        read(12,fmt=format211)(ind(i,j,k),i=1,ncx) ! read structure data from the file
       enddo
      enddo
      close(12)
c211   format(1i2,2i1)

c      do k=1,ncz
c       do j=1,ncy+1
c        do i=1,ncx
c         if(ind(i,j,k).gt.1)then ! make all pores ind = 1
c          ind(i,j,k)=1
c         endif
c        enddo
c       enddo
c      enddo

      allocate(inda(ncx+1,ncy+2,ncz+1),stat=checkstat) ! create solid point indicator matrix
      if(checkstat>0) then
       write(*,*)'ERROR! impossible to allocate memory for inda'
       stop
      endif

      inda=1 ! initialize inds with pores
      do k=1,ncz
       do j=1,ncy+1
        do i=1,ncx
         if(ind(i,j,k).eq.0.or.ind(i,j,k).eq.3)then ! solid element is detected
          inda(i,j,k)=0
          inda(i+1,j,k)=0
          inda(i,j+1,k)=0
          inda(i,j,k+1)=0
          inda(i+1,j+1,k)=0
          inda(i+1,j,k+1)=0
          inda(i,j+1,k+1)=0
          inda(i+1,j+1,k+1)=0
         endif
        enddo
       enddo
      enddo

      allocate(inds(ncx,ncy+2,ncz),stat=checkstat) ! create solid point indicator matrix
      if(checkstat>0) then
       write(*,*)'ERROR! impossible to allocate memory for inds'
       stop
      endif

      do j=1,ncy+2
       do k=1,ncz+1
        inda(1,j,k)=min(inda(1,j,k),inda(ncx+1,j,k)) ! stick together x sides
       enddo
      enddo

      do i=1,ncx+1
       do j=1,ncy+2
        inda(i,j,1)=min(inda(i,j,1),inda(i,j,ncz+1)) ! stick together z sides
       enddo
      enddo

      inds=1 ! initialize with pores
      inds(1:ncx,1:ncy+2,1:ncz)=inda(1:ncx,1:ncy+2,1:ncz) ! initialize from inda

      deallocate(inda) ! dont need it
      ncx=ncx   ! new bounds of the domain
      ncy=ncy+2
      ncz=ncz

!!! Each layer must contain solid points 
      layers = INT(ncz/(nthread2)) ! initial layers width (2*number of threads)
      do i=1,mod(ncz,nthread2)     ! adjust width by reminder redistribution
       layers(i)=layers(i)+1
      enddo

      ns      = 0 ! initialize number of solid points
      kend    = 0
      layb(1) = 0 
      do l=1,nthread2 ! loop over the layers
       kstart = kend+1
       kend   = kstart+layers(l)-1
! ATTENTION ! here must be the same order for k-j-i as in the loop where ICORS defined
       do k=kstart,kend ! count the number of solid points 
        do j=1,ncy
         do i=1,ncx
          if(inds(i,j,k).eq.0) then ! 1 - liquid or pore point, 0 - solid point
           ns = ns + 1
          endif
         enddo 
        enddo
       enddo
       layb(l+1)=ns
      enddo
      rns=ns

      allocate(icors(3,ns),stat=checkstat) ! create icors (solid points coordinates)
      if(checkstat>0) then
       write(*,*)'ERROR! impossible to allocate memory for icors'
       stop
      endif

      allocate(mijks(ncx,ncy,ncz),stat=checkstat) ! create mijks (solid points numbers)
      if(checkstat>0) then
       write(*,*)'ERROR! impossible to allocate memory for mijk'
       stop
      endif

      mijks = 0 ! initialize mijks (0 - pore, not zero - solid)
      ns = 0 ! initialize number of solid points
      do k=1,ncz ! define coordinates and numbers of solid points 
       do j=1,ncy
        do i=1,ncx
         if(inds(i,j,k).eq.0) then ! 0 - solid point
          ns = ns + 1
! to take into account dl icors should be changed
          icors(1,ns) = i ! define coordinate of the solid point
          icors(2,ns) = j
          icors(3,ns) = k
          mijks(i,j,k) = ns ! assign the number to the solid point 
         endif
        enddo 
       enddo
      enddo

      allocate(ipcx(ncx+2),stat=checkstat) ! create x-axis periodic boundary conditions
      if(checkstat>0) then
       write(*,*)'ERROR! impossible to allocate memory for ipcx'
       stop
      endif

      allocate(ipcy(ncy+2),stat=checkstat) ! create y-axis periodic boundary conditions
      if(checkstat>0) then
       write(*,*)'ERROR! impossible to allocate memory for ipcy'
       stop
      endif

      allocate(ipcz(ncz+2),stat=checkstat) ! create z-axis periodic boundary conditions
      if(checkstat>0) then
       write(*,*)'ERROR! impossible to allocate memory for ipcz'
       stop
      endif

      allocate(iicx(ncx+2),stat=checkstat) ! create x-axis periodic boundary conditions
      if(checkstat>0) then
       write(*,*)'ERROR! impossible to allocate memory for iicx'
       stop
      endif
      allocate(iicy(ncy+1),stat=checkstat) ! create y-axis periodic boundary conditions
      if(checkstat>0) then
       write(*,*)'ERROR! impossible to allocate memory for iicy'
       stop
      endif
      allocate(iicz(ncz+2),stat=checkstat) ! create z-axis periodic boundary conditions
      if(checkstat>0) then
       write(*,*)'ERROR! impossible to allocate memory for iicz'
       stop
      endif

      iicx=(/ncx,(i,i=1,ncx),1/) ! s.p. bound.cond. ncx+1 - max, ncx+2 - imposible !
      iicy=(/ncy-1,(i,i=1,ncy-1),1/)
      iicz=(/ncz,(i,i=1,ncz),1/)

      ipcx=(/ncx,(i,i=1,ncx),1/) ! periodic boundary conditions vectors
      ipcy=(/ncy,(i,i=1,ncy),1/)
      ipcz=(/ncz,(i,i=1,ncz),1/)

      nbo = 0 ! initialize the number of bonds counter
      do i=1,ns ! solid points loop
       do j=1,9 ! bond vectors loop (velocities9 file)
        npc(1:3)=icors(1:3,i)+ib(1:3,j) ! coordinates of the point adjacent to i by bond j
        nx1=icors(1,i)
        ny1=icors(2,i)
        nz1=icors(3,i)
        nx=ipcx(npc(1)+1) ! apply periodic conditions
        ny=ipcy(npc(2)+1) ! apply periodic conditions
        nz=ipcz(npc(3)+1) ! apply periodic conditions
        if(inds(nx,ny,nz).eq.0) then ! the next point is solid (may be there is a bond)
         ny1=iicy(ny1+1) ! s.p. for F1-Fn=E
         ncube1=0
         ncube2=0
         if(j.le.3)then ! first three cases
          if(j.eq.1)then ! 1 0 0 bond
c           ncube=ind(nx1,ny1,nz1)+ind(nx1,ny1,iicz(nz1))
c     &     +ind(nx1,iicy(ny1),nz1)+ind(nx1,iicy(ny1),iicz(nz1))
              if (ind(nx1,ny1,nz1).eq.0) then 
                                              ncube1=ncube1+1
              elseif (ind(nx1,ny1,nz1).eq.3) then 
                                              ncube2=ncube2+1
              endif
              if (ind(nx1,ny1,iicz(nz1)).eq.0) then
                                              ncube1=ncube1+1
              elseif (ind(nx1,ny1,iicz(nz1)).eq.3) then 
                                              ncube2=ncube2+1
              endif    
              if (ind(nx1,iicy(ny1),nz1).eq.0) then 
                                              ncube1=ncube1+1
              elseif (ind(nx1,iicy(ny1),nz1).eq.3) then
                                              ncube2=ncube2+1
              endif
              if (ind(nx1,iicy(ny1),iicz(nz1)).eq.0) then 
                                                     ncube1=ncube1+1
              elseif (ind(nx1,iicy(ny1),iicz(nz1)).eq.3) then
                                                      ncube2=ncube2+1
              endif
          elseif(j.eq.2)then ! 0 1 0 bond
c           ncube=ind(nx1,ny1,nz1)+ind(iicx(nx1),ny1,nz1)
c     &     +ind(nx1,ny1,iicz(nz1))+ind(iicx(nx1),ny1,iicz(nz1))
              if (ind(nx1,ny1,nz1).eq.0) then
                                            ncube1=ncube1+1
              elseif (ind(nx1,ny1,nz1).eq.3) then 
                                            ncube2=ncube2+1
              endif
              if (ind(iicx(nx1),ny1,nz1).eq.0) then
                                             ncube1=ncube1+1
              elseif (ind(iicx(nx1),ny1,nz1).eq.3) then
                                             ncube2=ncube2+1
              endif
              if (ind(nx1,ny1,iicz(nz1)).eq.0) then
                                             ncube1=ncube1+1
              elseif (ind(nx1,ny1,iicz(nz1)).eq.3) then 
                                             ncube2=ncube2+1
              endif
              if (ind(iicx(nx1),ny1,iicz(nz1)).eq.0) then 
                                                     ncube1=ncube1+1
              elseif (ind(iicx(nx1),ny1,iicz(nz1)).eq.3) then 
                                                      ncube2=ncube2+1
              endif

               
          elseif(j.eq.3)then ! 0 0 1 bond
c           ncube=ind(nx1,ny1,nz1)+ind(iicx(nx1),ny1,nz1)
c     &     +ind(nx1,iicy(ny1),nz1)+ind(iicx(nx1),iicy(ny1),nz1)
              if (ind(nx1,ny1,nz1).eq.0) then
                                             ncube1=ncube1+1
              elseif (ind(nx1,ny1,nz1).eq.3) then 
                                             ncube2=ncube2+1
              endif
              if (ind(iicx(nx1),ny1,nz1).eq.0) then
                                             ncube1=ncube1+1
              elseif (ind(iicx(nx1),ny1,nz1).eq.3) then
                                             ncube2=ncube2+1
              endif
              if (ind(nx1,iicy(ny1),nz1).eq.0) then
                                             ncube1=ncube1+1
              elseif (ind(nx1,iicy(ny1),nz1).eq.3) then
                                             ncube2=ncube2+1
              endif
              if (ind(iicx(nx1),iicy(ny1),nz1).eq.0) then 
                                                     ncube1=ncube1+1
              elseif (ind(iicx(nx1),iicy(ny1),nz1).eq.3) then 
                                                      ncube2=ncube2+1
              endif
           endif
          nucube=0
          ncube=ncube1+ncube2
          if(ncube.gt.0)then
           nbo=nbo+1 ! there is a bond
          endif
         else ! last six cases
          if(j.eq.4)then ! 1 1 0
c           ncube=ind(nx1,ny1,nz1)+ind(nx1,ny1,iicz(nz1))
              if (ind(nx1,ny1,nz1).eq.0) then
                                               ncube1=ncube1+1
              elseif (ind(nx1,ny1,nz1).eq.3) then 
                                               ncube2=ncube2+1
              endif
              if (ind(nx1,ny1,iicz(nz1)).eq.0) then 
                                               ncube1=ncube1+1
              elseif (ind(nx1,ny1,iicz(nz1)).eq.3) then
                                               ncube2=ncube2+1
              endif
          elseif(j.eq.5)then ! 1 -1 0
c           ncube=ind(nx1,iicy(ny1),nz1)+ind(nx1,iicy(ny1),iicz(nz1))
         
              if (ind(nx1,iicy(ny1),nz1).eq.0) then 
                                               ncube1=ncube1+1
              elseif (ind(nx1,iicy(ny1),nz1).eq.3) then
                                               ncube2=ncube2+1
              endif
              if (ind(nx1,iicy(ny1),iicz(nz1)).eq.0) then 
                                                         ncube1=ncube1+1
              elseif (ind(nx1,iicy(ny1),iicz(nz1)).eq.3) then 
                                                         ncube2=ncube2+1
              endif

          elseif(j.eq.6)then ! 1 0 1
c           ncube=ind(nx1,ny1,nz1)+ind(nx1,iicy(ny1),nz1)
              if (ind(nx1,ny1,nz1).eq.0) then 
                                              ncube1=ncube1+1
              elseif (ind(nx1,ny1,nz1).eq.3) then
                                              ncube2=ncube2+1
              endif
              if (ind(nx1,iicy(ny1),nz1).eq.0) then
                                              ncube1=ncube1+1
              elseif (ind(nx1,iicy(ny1),nz1).eq.3) then 
                                              ncube2=ncube2+1
              endif
 
          elseif(j.eq.7)then ! 1 0 -1
c           ncube=ind(nx1,ny1,iicz(nz1))+ind(nx1,iicy(ny1),iicz(nz1))
              if (ind(nx1,ny1,iicz(nz1)).eq.0) then 
                                              ncube1=ncube1+1
              elseif (ind(nx1,ny1,iicz(nz1)).eq.3) then
                                              ncube2=ncube2+1
              endif
              if (ind(nx1,iicy(ny1),iicz(nz1)).eq.0) then
                                                         ncube1=ncube1+1
              elseif (ind(nx1,iicy(ny1),iicz(nz1)).eq.3) then 
                                                         ncube2=ncube2+1
              endif
     
          elseif(j.eq.8)then ! 0 1 1
c           ncube=ind(nx1,ny1,nz1)+ind(iicx(nx1),ny1,nz1)
              if (ind(nx1,ny1,nz1).eq.0) then 
                                             ncube1=ncube1+1
              elseif (ind(nx1,ny1,nz1).eq.3) then
                                             ncube2=ncube2+1
              endif
              if (ind(iicx(nx1),ny1,nz1).eq.0) then 
                                                         ncube1=ncube1+1
              elseif (ind(iicx(nx1),ny1,nz1).eq.3) then 
                                                         ncube2=ncube2+1
              endif
          
          elseif(j.eq.9)then ! 0 1 -1
c           ncube=ind(nx1,ny1,iicz(nz1))+ind(iicx(nx1),ny1,iicz(nz1))
              if (ind(nx1,ny1,iicz(nz1)).eq.0) then
                                              ncube1=ncube1+1
              elseif (ind(nx1,ny1,iicz(nz1)).eq.3) then 
                                              ncube2=ncube2+1
              endif
              if (ind(iicx(nx1),ny1,iicz(nz1)).eq.0) then 
                                                         ncube1=ncube1+1
              elseif (ind(iicx(nx1),ny1,iicz(nz1)).eq.3) then 
                                                         ncube2=ncube2+1
              endif
          endif
          ncube=ncube1+ncube2
          if(ncube.gt.0)then
           nbo=nbo+1 ! there is a bond
          endif
         endif
        endif
       enddo
      enddo

      allocate(bon(3,nbo),stat=checkstat) ! create bon - bonds array
      if(checkstat>0) then
       write(*,*)'ERROR! impossible to allocate memory for bon'
       stop
      endif

      allocate(sbond(nbo),stat=checkstat) ! create sbond - surface bonds indicator array
      if(checkstat>0) then
       write(*,*)'ERROR! impossible to allocate memory for sbond'
       stop
      endif
      allocate(sbond1(nbo),stat=checkstat) ! create sbond - surface bonds indicator array
      if(checkstat>0) then
       write(*,*)'ERROR! impossible to allocate memory for sbond'
       stop
      endif
      allocate(sbond2(nbo),stat=checkstat) ! create sbond - surface bonds indicator array
      if(checkstat>0) then
       write(*,*)'ERROR! impossible to allocate memory for sbond'
       stop
      endif
      allocate(sbondin(nbo),stat=checkstat) ! create sbond - surface bonds indicator array
      if(checkstat>0) then
       write(*,*)'ERROR! impossible to allocate memory for sbond'
       stop
      endif
    
c THE start
      nbo    = 0 ! initialize the number of bonds counter
      lib(1) = 0
      do ii=1,nthread2
      istart = layb(ii)+1 ! layer boundaries
      iend   = layb(ii+1) ! layer boundaries
      do i=istart,iend ! solid points loop
       do j=1,9 ! bond vectors loop (velocities9 file)
        npc(1:3)=icors(1:3,i)+ib(1:3,j) ! coordinates of the point adjacent to i by bond j
        nx1=icors(1,i)
        ny1=icors(2,i)
        nz1=icors(3,i)
        nx=ipcx(npc(1)+1) ! apply periodic conditions
        ny=ipcy(npc(2)+1) ! apply periodic conditions
        nz=ipcz(npc(3)+1) ! apply periodic conditions
         if(inds(nx,ny,nz).eq.0) then ! the next point is solid (may be there is a bond)
         ny1=iicy(ny1+1) ! s.p. for F1-Fn=E
         ncube1=0
         ncube2=0
          if(j.le.3)then ! first three cases
           if(j.eq.1)then ! 1 0 0 bond
              if (ind(nx1,ny1,nz1).eq.0) then 
                                             ncube1=ncube1+1
              elseif (ind(nx1,ny1,nz1).eq.3) then 
                                             ncube2=ncube2+1
              endif
              if (ind(nx1,ny1,iicz(nz1)).eq.0) then
                                             ncube1=ncube1+1
              elseif (ind(nx1,ny1,iicz(nz1)).eq.3) then
                                             ncube2=ncube2+1
              endif
              if (ind(nx1,iicy(ny1),nz1).eq.0) then 
                                             ncube1=ncube1+1
              elseif (ind(nx1,iicy(ny1),nz1).eq.3) then
                                             ncube2=ncube2+1
              endif
              if (ind(nx1,iicy(ny1),iicz(nz1)).eq.0) then
                                             ncube1=ncube1+1
              elseif (ind(nx1,iicy(ny1),iicz(nz1)).eq.3) then
                                             ncube2=ncube2+1
              endif
          elseif(j.eq.2)then ! 0 1 0 bond
              if (ind(nx1,ny1,nz1).eq.0) then
                                             ncube1=ncube1+1
              elseif (ind(nx1,ny1,nz1).eq.3) then
                                             ncube2=ncube2+1
              endif
              if (ind(iicx(nx1),ny1,nz1).eq.0) then 
                                             ncube1=ncube1+1
              elseif (ind(iicx(nx1),ny1,nz1).eq.3) then 
                                             ncube2=ncube2+1
              endif
              if (ind(nx1,ny1,iicz(nz1)).eq.0) then 
                                             ncube1=ncube1+1
              elseif (ind(nx1,ny1,iicz(nz1)).eq.3) then 
                                             ncube2=ncube2+1
              endif
              if (ind(iicx(nx1),ny1,iicz(nz1)).eq.0) then
                                                     ncube1=ncube1+1
              elseif (ind(iicx(nx1),ny1,iicz(nz1)).eq.3) then
                                                      ncube2=ncube2+1
              endif
           elseif(j.eq.3)then ! 0 0 1 bond
              if (ind(nx1,ny1,nz1).eq.0) then
                                              ncube1=ncube1+1
              elseif (ind(nx1,ny1,nz1).eq.3) then
                                              ncube2=ncube2+1
              endif
              if (ind(iicx(nx1),ny1,nz1).eq.0) then
                                              ncube1=ncube1+1
              elseif (ind(iicx(nx1),ny1,nz1).eq.3) then
                                              ncube2=ncube2+1
              endif
              if (ind(nx1,iicy(ny1),nz1).eq.0) then
                                              ncube1=ncube1+1
              elseif (ind(nx1,iicy(ny1),nz1).eq.3) then
                                              ncube2=ncube2+1
              endif
              if (ind(iicx(nx1),iicy(ny1),nz1).eq.0) then
                                                     ncube1=ncube1+1
              elseif (ind(iicx(nx1),iicy(ny1),nz1).eq.3) then
                                                      ncube2=ncube2+1
              endif
           endif
           ncube=ncube1+ncube2
           if(ncube.gt.0)then
           bconst1 = ncube1*0.25D0 ! bond constant component 1
           bconst2 = ncube2*0.25D0 ! bon constant compnent 2
           nbo=nbo+1 ! there is a bond
           bon(1,nbo)=i ! point 'from'
           bon(2,nbo)=mijks(nx,ny,nz) ! point 'to'
           bon(3,nbo)=j ! bond vector index 'from'->'to'
           sbond1(nbo)=bconst1 ! bond elastic constant
           sbond2(nbo)=bconst2
          endif
         else ! last six cases
           if(j.eq.4)then ! 1 1 0
c           ncube=ind(nx1,ny1,nz1)+ind(nx1,ny1,iicz(nz1))
              if (ind(nx1,ny1,nz1).eq.0) then
                                              ncube1=ncube1+1
              elseif (ind(nx1,ny1,nz1).eq.3) then 
                                              ncube2=ncube2+1
              endif
              if (ind(nx1,ny1,iicz(nz1)).eq.0) then 
                                              ncube1=ncube1+1
              elseif (ind(nx1,ny1,iicz(nz1)).eq.3) then
                                              ncube2=ncube2+1
              endif
          elseif(j.eq.5)then ! 1 -1 0
c           ncube=ind(nx1,iicy(ny1),nz1)+ind(nx1,iicy(ny1),iicz(nz1))

              if (ind(nx1,iicy(ny1),nz1).eq.0) then 
                                              ncube1=ncube1+1
              elseif (ind(nx1,iicy(ny1),nz1).eq.3) then 
                                              ncube2=ncube2+1
              endif
              if (ind(nx1,iicy(ny1),iicz(nz1)).eq.0) then
                                                         ncube1=ncube1+1
              elseif (ind(nx1,iicy(ny1),iicz(nz1)).eq.3) then
                                                         ncube2=ncube2+1
              endif
           elseif(j.eq.6)then ! 1 0 1
c           ncube=ind(nx1,ny1,nz1)+ind(nx1,iicy(ny1),nz1)
              if (ind(nx1,ny1,nz1).eq.0) then 
                                              ncube1=ncube1+1
              elseif (ind(nx1,ny1,nz1).eq.3) then 
                                              ncube2=ncube2+1
              endif
              if (ind(nx1,iicy(ny1),nz1).eq.0) then
                                              ncube1=ncube1+1
              elseif (ind(nx1,iicy(ny1),nz1).eq.3) then
                                              ncube2=ncube2+1
              endif

          elseif(j.eq.7)then ! 1 0 -1
c           ncube=ind(nx1,ny1,iicz(nz1))+ind(nx1,iicy(ny1),iicz(nz1))
              if (ind(nx1,ny1,iicz(nz1)).eq.0) then 
                                              ncube1=ncube1+1
              elseif (ind(nx1,ny1,iicz(nz1)).eq.3) then
                                              ncube2=ncube2+1
              endif
              if (ind(nx1,iicy(ny1),iicz(nz1)).eq.0) then
                                                         ncube1=ncube1+1
              elseif (ind(nx1,iicy(ny1),iicz(nz1)).eq.3) then
                                                         ncube2=ncube2+1
              endif

          elseif(j.eq.8)then ! 0 1 1
c           ncube=ind(nx1,ny1,nz1)+ind(iicx(nx1),ny1,nz1)
              if (ind(nx1,ny1,nz1).eq.0) then
                                               ncube1=ncube1+1
              elseif (ind(nx1,ny1,nz1).eq.3) then
                                               ncube2=ncube2+1
              endif
              if (ind(iicx(nx1),ny1,nz1).eq.0) then
                                                         ncube1=ncube1+1
              elseif (ind(iicx(nx1),ny1,nz1).eq.3) then
                                                         ncube2=ncube2+1
              endif
          elseif(j.eq.9)then ! 0 1 -1
c           ncube=ind(nx1,ny1,iicz(nz1))+ind(iicx(nx1),ny1,iicz(nz1))
              if (ind(nx1,ny1,iicz(nz1)).eq.0) then
                                                 ncube1=ncube1+1
              elseif (ind(nx1,ny1,iicz(nz1)).eq.3) then
                                                 ncube2=ncube2+1
              endif
              if (ind(iicx(nx1),ny1,iicz(nz1)).eq.0) then
                                                         ncube1=ncube1+1
              elseif (ind(iicx(nx1),ny1,iicz(nz1)).eq.3) then
                                                         ncube2=ncube2+1
              endif
          endif
          ncube=ncube1+ncube2
          if(ncube.gt.0)then
           bconst1 = ncube1*0.5D0 ! bond constant
           bconst2 = ncube2*0.5D0 ! bond constant component 2
           nbo=nbo+1 ! there is a bond
           bon(1,nbo)=i ! point 'from'
           bon(2,nbo)=mijks(nx,ny,nz) ! point 'to'
           bon(3,nbo)=j ! bond vector index 'from'->'to'
           sbond1(nbo)=bconst1 ! bond elastic constant
           sbond2(nbo)=bconst2 ! compo 2
          endif
         endif
        endif
       enddo
      enddo
      lib(ii+1)=nbo
      enddo

      nna = 0 ! initialize the number of angles counter
      nn1=0
      do i=1,ns ! solid points loop
       do j=1,48 ! angles loop (angle i-j-k)
        npj(1:3)=icors(1:3,i)+iv(1:3,nta(1,j)) ! coordinates of the point j (i-j)
        nxj=ipcx(npj(1)+1) ! apply periodic conditions
        nyj=ipcy(npj(2)+1) ! apply periodic conditions
        nzj=ipcz(npj(3)+1) ! apply periodic conditions

        npk(1:3)=icors(1:3,i)+iv(1:3,nta(2,j)) ! coordinates of the point k (i-k)
        nxk=ipcx(npk(1)+1) ! apply periodic conditions
        nyk=ipcy(npk(2)+1) ! apply periodic conditions
        nzk=ipcz(npk(3)+1) ! apply periodic conditions

        nx1=icors(1,i) ! top point coordinates
        ny1=icors(2,i)
        nz1=icors(3,i)

        if(inds(nxj,nyj,nzj).eq.0.and.inds(nxk,nyk,nzk).eq.0) then ! may be there is an angle
         if(nta(1,j).gt.7.and.nta(2,j).gt.7)then ! only one cube possible
          inx1=iicx(min(nx1,npj(1),npk(1))+1)
          iny1=iicy(min(ny1,npj(2),npk(2))+1)
          inz1=iicz(min(nz1,npj(3),npk(3))+1)
          iscube = ind(inx1,iny1,inz1)
          if(iscube.eq.0.or.iscube.eq.3)then ! there is a cube and so there is an angle
           nna=nna+1 ! there is an angle i-j-k
           nn1=nn1+1
          endif
         else ! two cubes possible
          jv=iv(1:3,nta(1,j))
          kv=iv(1:3,nta(2,j))
          iscube1=0
          iscube2=0
          if(jv(1).eq.0.and.kv(1).eq.0)then ! yz plane
           inx1=iicx(min(nx1,npj(1),npk(1))+1)
           iny1=iicy(min(ny1,npj(2),npk(2))+1)
           inz1=iicz(min(nz1,npj(3),npk(3))+1)
c           iscube = ind(inx1,iny1,inz1)+ind(iicx(inx1),iny1,inz1)
              if (ind(inx1,iny1,inz1).eq.0) then
                                             iscube1=iscube1+1
              elseif (ind(inx1,iny1,inz1).eq.3) then 
                                             iscube2=iscube2+1
              endif
              if (ind(iicx(inx1),iny1,inz1).eq.0) then
                                                     iscube1=iscube1+1
              elseif (ind(iicx(inx1),iny1,inz1).eq.3) then
                                                     iscube2=iscube2+1
              endif
          elseif(jv(2).eq.0.and.kv(2).eq.0)then ! xz plane
           inx1=iicx(min(nx1,npj(1),npk(1))+1)
           iny1=iicy(min(ny1,npj(2),npk(2))+1)
           inz1=iicz(min(nz1,npj(3),npk(3))+1)
c           iscube = ind(inx1,iny1,inz1)+ind(inx1,iicy(iny1),inz1)
              if (ind(inx1,iny1,inz1).eq.0) then
                                             iscube1=iscube1+1
              elseif (ind(inx1,iny1,inz1).eq.3) then 
                                             iscube2=iscube2+1
              endif
              if (ind(inx1,iicy(iny1),inz1).eq.0) then
                                                     iscube1=iscube1+1
              elseif (ind(inx1,iicy(iny1),inz1).eq.3) then
                                                     iscube2=iscube2+1
              endif
          elseif(jv(3).eq.0.and.kv(3).eq.0)then ! xy plane
           inx1=iicx(min(nx1,npj(1),npk(1))+1)
           iny1=iicy(min(ny1,npj(2),npk(2))+1)
           inz1=iicz(min(nz1,npj(3),npk(3))+1)
c           iscube = ind(inx1,iny1,inz1)+ind(inx1,iny1,iicz(inz1))
              if (ind(inx1,iny1,inz1).eq.0) then
                                             iscube1=iscube1+1
              elseif (ind(inx1,iny1,inz1).eq.3) then 
                                             iscube2=iscube2+1
              endif
              if (ind(inx1,iny1,iicz(inz1)).eq.0) then
                                                     iscube1=iscube1+1
              elseif (ind(inx1,iny1,iicz(inz1)).eq.3) then
                                                     iscube2=iscube2+1
              endif
          endif
          iscube=0
          iscube=iscube1+iscube2
          if(iscube.gt.0)then ! there is an angle
           nna=nna+1 ! there is an angle i-j-k
          endif
         endif
        endif
       enddo
      enddo
c      write(*,*)nna,nn1
c THE end
 
      allocate(ang(6,nna),stat=checkstat) ! create ang - angles array
      if(checkstat>0) then
       write(*,*)'ERROR! impossible to allocate memory for ang'
       stop
      endif

      allocate(abond(nna),stat=checkstat) ! create abond - surface bonds indicator array
      if(checkstat>0) then
       write(*,*)'ERROR! impossible to allocate memory for abond'
       stop
      endif
      allocate(abond1(nna),stat=checkstat) ! create abond - surface bonds indicator array
      if(checkstat>0) then
       write(*,*)'ERROR! impossible to allocate memory for abond'
       stop
      endif
      allocate(abond2(nna),stat=checkstat) ! create abond - surface bonds indicator array
      if(checkstat>0) then
       write(*,*)'ERROR! impossible to allocate memory for abond'
       stop
      endif
      allocate(abondin(nna),stat=checkstat) ! create sbond - surface bonds indicator array
      if(checkstat>0) then
       write(*,*)'ERROR! impossible to allocate memory for sbond'
       stop
      endif
c THE start
      nna = 0 ! initialize the number of angles counter
      nn1=0
      lia(1) = 0
      do ii=1,nthread2
      istart = layb(ii)+1 ! layer boundaries
      iend   = layb(ii+1) ! layer boundaries
      do i=istart,iend ! solid points loop
       do j=1,48 ! angles loop (angle i-j-k)
        npj(1:3)=icors(1:3,i)+iv(1:3,nta(1,j)) ! coordinates of the point j (i-j)
        nxj=ipcx(npj(1)+1) ! apply periodic conditions
        nyj=ipcy(npj(2)+1) ! apply periodic conditions
        nzj=ipcz(npj(3)+1) ! apply periodic conditions

        npk(1:3)=icors(1:3,i)+iv(1:3,nta(2,j)) ! coordinates of the point k (i-k)
        nxk=ipcx(npk(1)+1) ! apply periodic conditions
        nyk=ipcy(npk(2)+1) ! apply periodic conditions
        nzk=ipcz(npk(3)+1) ! apply periodic conditions

        nx1=icors(1,i) ! top point coordinates
        ny1=icors(2,i)
        nz1=icors(3,i)

        if(inds(nxj,nyj,nzj).eq.0.and.inds(nxk,nyk,nzk).eq.0) then ! may be there is an angle
         if(nta(1,j).gt.7.and.nta(2,j).gt.7)then ! only one cube possible
          inx1=iicx(min(nx1,npj(1),npk(1))+1)
          iny1=iicy(min(ny1,npj(2),npk(2))+1)
          inz1=iicz(min(nz1,npj(3),npk(3))+1)
          iscube = ind(inx1,iny1,inz1)
          if(iscube.eq.0)then ! there is a cube and so there is an angle
           nna=nna+1 ! there is an angle i-j-k
           nn1=nn1+1
           ang(1,nna)=i ! angle top point i
           ang(2,nna)=mijks(nxj,nyj,nzj) ! angleside endpoint j
           ang(3,nna)=mijks(nxk,nyk,nzk) ! angleside endpoint k
           ang(4,nna)=nta(1,j) ! i-j vector index
           ang(5,nna)=nta(2,j) ! i-k vector index
           ang(6,nna)=j ! angle indentificator
           abond1(nna)=1.D0 ! elastic constant
           abond2(nna)=0.D0 ! elastic constan of angular spring compo 2
          endif
           if(iscube.eq.3)then ! there is a cube and so there is an angle
           nna=nna+1 ! there is an angle i-j-k
           nn1=nn1+1
           ang(1,nna)=i ! angle top point i
           ang(2,nna)=mijks(nxj,nyj,nzj) ! angleside endpoint j
           ang(3,nna)=mijks(nxk,nyk,nzk) ! angleside endpoint k
           ang(4,nna)=nta(1,j) ! i-j vector index
           ang(5,nna)=nta(2,j) ! i-k vector index
           ang(6,nna)=j ! angle indentificator
           abond1(nna)=0.D0 ! elastic constant
           abond2(nna)=1.D0 ! elastic constan of angular spring compo 2
          endif
         else ! two cubes possible
          jv=iv(1:3,nta(1,j))
          kv=iv(1:3,nta(2,j))
          iscube1=0
          iscube2=0
          if(jv(1).eq.0.and.kv(1).eq.0)then ! yz plane
           inx1=iicx(min(nx1,npj(1),npk(1))+1)
           iny1=iicy(min(ny1,npj(2),npk(2))+1)
           inz1=iicz(min(nz1,npj(3),npk(3))+1)
c           iscube = ind(inx1,iny1,inz1)+ind(iicx(inx1),iny1,inz1)
              if (ind(inx1,iny1,inz1).eq.0) then
                                             iscube1=iscube1+1
              elseif (ind(inx1,iny1,inz1).eq.3) then
                                             iscube2=iscube2+1
              endif
              if (ind(iicx(inx1),iny1,inz1).eq.0) then
                                                     iscube1=iscube1+1
              elseif (ind(iicx(inx1),iny1,inz1).eq.3) then
                                                     iscube2=iscube2+1
              endif
           elseif(jv(2).eq.0.and.kv(2).eq.0)then ! xz plane
           inx1=iicx(min(nx1,npj(1),npk(1))+1)
           iny1=iicy(min(ny1,npj(2),npk(2))+1)
           inz1=iicz(min(nz1,npj(3),npk(3))+1)
c           iscube = ind(inx1,iny1,inz1)+ind(inx1,iicy(iny1),inz1)
              if (ind(inx1,iny1,inz1).eq.0) then
                                             iscube1=iscube1+1
              elseif (ind(inx1,iny1,inz1).eq.3) then
                                             iscube2=iscube2+1
              endif
              if (ind(inx1,iicy(iny1),inz1).eq.0) then
                                                     iscube1=iscube1+1
              elseif (ind(inx1,iicy(iny1),inz1).eq.3) then
                                                     iscube2=iscube2+1
              endif
          elseif(jv(3).eq.0.and.kv(3).eq.0)then ! xy plane
           inx1=iicx(min(nx1,npj(1),npk(1))+1)
           iny1=iicy(min(ny1,npj(2),npk(2))+1)
           inz1=iicz(min(nz1,npj(3),npk(3))+1)
c           iscube = ind(inx1,iny1,inz1)+ind(inx1,iny1,iicz(inz1))
              if (ind(inx1,iny1,inz1).eq.0) then
                                             iscube1=iscube1+1
              elseif (ind(inx1,iny1,inz1).eq.3) then
                                             iscube2=iscube2+1
              endif
              if (ind(inx1,iny1,iicz(inz1)).eq.0) then
                                                     iscube1=iscube1+1
              elseif (ind(inx1,iny1,iicz(inz1)).eq.3) then
                                                     iscube2=iscube2+1
              endif
          endif
          iscube=iscube1+iscube2
          if(iscube.gt.0)then ! there is an angle
           nna=nna+1 ! there is an angle i-j-k
           nn1=nn1+1
           ang(1,nna)=i ! angle top point i
           ang(2,nna)=mijks(nxj,nyj,nzj) ! angleside endpoint j
           ang(3,nna)=mijks(nxk,nyk,nzk) ! angleside endpoint k
           ang(4,nna)=nta(1,j) ! i-j vector index
           ang(5,nna)=nta(2,j) ! i-k vector index
           ang(6,nna)=j ! angle indentificator
           abond1(nna)=iscube1*0.5D0 ! elastic constant
           abond2(nna)=iscube2*0.5D0 ! elastic constan of angular spring compo 2
          endif
         endif
        endif
       enddo
      enddo
      lia(ii+1)=nna
      enddo

      allocate(adat(3,ns),stat=checkstat) ! create adat
      if(checkstat>0) then
       write(*,*)'ERROR! impossible to allocate memory for adat'
       stop
      endif

      allocate(vdat(3,ns),stat=checkstat) ! create vdat
      if(checkstat>0) then
       write(*,*)'ERROR! impossible to allocate memory for vdat'
       stop
      endif

      allocate(rdat(3,ns),stat=checkstat) ! create rdat
      if(checkstat>0) then
       write(*,*)'ERROR! impossible to allocate memory for rdat'
       stop
      endif

! Initialize the adat,vdat,rdat arrays
      adat = 0.d0 ! initialize acceleration, velocity with zeroes
      vdat = 0.d0
      rdat = 0.D0
      if(icon.eq.1)then
       open(1,file='iteration2')
       read(1,*)((rdat(j,i),j=1,3),i=1,ns)
       read(1,*)((vdat(j,i),j=1,3),i=1,ns)
       read(1,*)((adat(j,i),j=1,3),i=1,ns)
       close(1)
      endif

! Calculate unit cell directions vectors
      vn(1:3,1) = vv(1:3,1) ! first is zero vector
      vx(1:3,1) = vv(1:3,1) ! first is zero vector
      vnr(1) = 0.d0 ! norm is zero
      vni(1) = 0.d0 ! zero (this one is not used anyway)
      do i=2,19 ! normalization of direction vectors and norms calculation
       vnr(i) = sqrt(dot_product(vv(1:3,i),vv(1:3,i))) ! norms
       vni(i) = 1.d0/vnr(i) ! norms inverse
       vn(1:3,i) = vv(1:3,i)/vnr(i) ! normalized vectors
       vx(1:3,i) = vv(1:3,i)/(dot_product(vv(1:3,i),vv(1:3,i))) ! xn/|rn|
      enddo

      do i=1,9 ! normalization of bond vectors and norms calculation
       vbn(i) = sqrt(dot_product(vb(1:3,i),vb(1:3,i))) ! norms
       vbi(i) = 1.d0/vbn(i) ! norms inverse
       bn(1:3,i) = vb(1:3,i)/vbn(i) ! normalized vectors
       bx(1:3,i) = vb(1:3,i)/(dot_product(vb(1:3,i),vb(1:3,i))) ! xn/|rn|
      enddo
! Save program parameters in lsm_settings

      open(19,file='lsm_settings') ! save program parameters

      write(19,*)'Initial program settings'
      write(19,*)'Physical parameters'
      write(19,*)'Young modulus Em1 = ',Em1
      write(19,*)'Poissons ratio nu1 = ',nu1
      write(19,*)'linear spring parameter alf1 = ',alf1
      write(19,*)'angular spring parameter bet1 = ',bet1
      write(19,*)'Young modulus Em2 = ',Em2
      write(19,*)'Poissons ratio nu2 = ',nu2
      write(19,*)'linear spring parameter alf2 = ',alf2
      write(19,*)'angular spring parameter bet2 = ',bet2
      write(19,*)'solid medium density rs = ',rs
      write(19,*)'integration time step dt = ',dt
      write(19,*)'lattice unit step dl = ',dl

      write(19,*)'Computational domain properties'
      write(19,*)'x unit cell size ncx = ',ncx
      write(19,*)'y unit cell size ncy = ',ncy
      write(19,*)'z unit cell size ncz = ',ncz
      write(19,*)'number of solid points ns = ',ns

      write(19,*)'LSM algorithm parameters'
      write(19,*)'simulation duration itm = ',itm
      write(19,*)'save parameter its = ',its
      write(19,*)'nbo = ',nbo
      write(19,*)'nna = ',nna
      write(19,*)'layers = ',(layers(i),i=1,2*nthread)
      write(19,*)'layb = ',(layb(i),i=1,2*nthread+1)
      write(19,*)'lib = ',(lib(i),i=1,2*nthread+1)
      write(19,*)'lia = ',(lia(i),i=1,2*nthread+1)
  
      close(19)

c      do i=1,nbo
c       if(icors(3,bon(2,i)).eq.1.and.icors(3,bon(1,i)).eq.1)then
c        write(111,*)bon(1:3,i),sbond(i)
c       endif
c      enddo
c      do i=1,nna
c       write(222,*)ang(1:6,i),abond(i)
c      enddo

c      sbond = sbond*alf
c      abond = abond*bet
       sbond = sbond1*alf1+sbond2*alf2
       abond = abond1*bet1+abond2*bet2

c      gstrain=0.0004D0*0.5D0
      do i=1,ncx
       do k=1,ncz
        l1=mijks(i,1,k)
        l2=mijks(i,ncy-1,k)
        if(l1.gt.0.and.l2.gt.0)then
         rdat(1:3,l1)=(/0.D0,0.D0,-gstrain/)
         rdat(1:3,l2)=(/0.D0,0.D0,gstrain/)
        endif
       enddo
      enddo

      do iti=1,itm ! main time simulations circle
      do itj=1,its ! time simulation subcircle (after each subcircle program check and saves the data)

! Integrate LSM - velocity Verlet algorithm: Start
! Note: now in LSM dt = 1, as well as in LBM dt = 1, so the relations can be simplified

! Second step: calculate new velocities v_i(t+dt/2)=v_i(t)+a_i(t)dt/2
!$OMP  PARALLEL DO DEFAULT(NONE)
!$OMP& PRIVATE(i)
!$OMP& SHARED(ns,vdat,adat,dt)
      do i=1,ns
       vdat(1:3,i)=vdat(1:3,i)+adat(1:3,i)*dt*0.5d0
      enddo
!$OMP END PARALLEL DO

! First step: calculate new nodes positions r_i(t+dt)=r_i(t)+v_i(t)dt+(a_i(t)dt^2)/2
! First step: calculate new nodes positions r_i(t+dt)=r_i(t)+v_i(t+dt/2)dt

!$OMP  PARALLEL DO DEFAULT(NONE)
!$OMP& PRIVATE(i)
!$OMP& SHARED(ns,rdat,vdat,dt)
      do i=1,ns
       rdat(1:3,i)=rdat(1:3,i)+vdat(1:3,i)*dt
      enddo
!$OMP END PARALLEL DO

! Third step: calculate new force F_i(t+dt)  or , a_i(t+dt)=F_i(t+dt) /M_i
      adat = 0.d0 ! initialize all the forces

!$OMP  PARALLEL DEFAULT(NONE)
!$OMP& PRIVATE(i,l,xn,rn,rt,flin,lstart,lend)
!$OMP& SHARED(nbo,bn,bon,vb,rdat,icors,sbond,adat,lib,nthread2)
!$OMP DO
      do l=1,nthread2,2
       lstart = lib(l)+1
       lend   = lib(l+1)
       do i=lstart,lend ! bonds loop to calculate linear springs force
        xn=bn(1:3,bon(3,i)) ! normalized equilibrium i-j vector
        rn=vb(1:3,bon(3,i)) ! not normalized equilibrium i-j vector
        rt=rdat(1:3,bon(2,i))-rdat(1:3,bon(1,i))+rn ! deformed i-j vector

        flin=sbond(i)*((rt(1)*xn(1)+rt(2)*xn(2)+rt(3)*xn(3))*xn-rn) ! linear spring force
        adat(1:3,bon(1,i))=adat(1:3,bon(1,i))+flin ! force acting on the first end of the bond
        adat(1:3,bon(2,i))=adat(1:3,bon(2,i))-flin ! force acting on the second end of the bond
       enddo
      enddo
!$OMP END DO
!$OMP BARRIER
!$OMP DO
      do l=2,nthread2,2
       lstart = lib(l)+1
       lend   = lib(l+1)
       do i=lstart,lend ! bonds loop to calculate linear springs force
        xn=bn(1:3,bon(3,i)) ! normalized equilibrium i-j vector
        rn=vb(1:3,bon(3,i)) ! not normalized equilibrium i-j vector
        rt=rdat(1:3,bon(2,i))-rdat(1:3,bon(1,i))+rn ! deformed i-j vector

        flin=sbond(i)*((rt(1)*xn(1)+rt(2)*xn(2)+rt(3)*xn(3))*xn-rn) ! linear spring force
        adat(1:3,bon(1,i))=adat(1:3,bon(1,i))+flin ! force acting on the first end of the bond
        adat(1:3,bon(2,i))=adat(1:3,bon(2,i))-flin ! force acting on the second end of the bond
       enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL 

!$OMP  PARALLEL DEFAULT(NONE)
!$OMP& PRIVATE(i,xrj,rnj,rtj,xrk,rnk,rtk,dfi,l,lstart,lend)
!$OMP& SHARED(nna,vn,vx,vv,ang,rdat,icors,anf,abond,adat,lia,nthread2)
!$OMP DO
      do l=1,nthread2,2
       lstart = lia(l)+1
       lend   = lia(l+1)
       do i=lstart,lend ! angles loop
        xrj=vx(1:3,ang(4,i)) ! xn/|rn|
        rnj=vv(1:3,ang(4,i)) ! not normalized equilibrium i-j vector
        rtj=rdat(1:3,ang(2,i))-rdat(1:3,ang(1,i))+rnj ! deformed i-j vector
 
        xrk=vx(1:3,ang(5,i)) ! xn/|rn|
        rnk=vv(1:3,ang(5,i)) ! not normalized equilibrium i-k vector
        rtk=rdat(1:3,ang(3,i))-rdat(1:3,ang(1,i))+rnk ! deformed i-k vector
        
        dfi=((xrk(2)*rtk(3)-xrk(3)*rtk(2)-xrj(2)*rtj(3)+xrj(3)*rtj(2)) ! angle change delta phi
     & *anf(1,ang(6,i))
     & -(xrk(1)*rtk(3)-xrk(3)*rtk(1)-xrj(1)*rtj(3)+xrj(3)*rtj(1))
     & *anf(2,ang(6,i))
     & +(xrk(1)*rtk(2)-xrk(2)*rtk(1)-xrj(1)*rtj(2)+xrj(2)*rtj(1))
     & *anf(3,ang(6,i)))*abond(i)

        adat(1:3,ang(2,i))=adat(1:3,ang(2,i))+dfi*anf(4:6,ang(6,i)) ! angular force for point j
        adat(1:3,ang(3,i))=adat(1:3,ang(3,i))+dfi*anf(7:9,ang(6,i)) ! angular force for point k
        adat(1:3,ang(1,i))=adat(1:3,ang(1,i))
     & -(dfi*anf(4:6,ang(6,i))+dfi*anf(7:9,ang(6,i)))
       enddo
      enddo
!$OMP END DO
!$OMP BARRIER
!$OMP DO
      do l=2,nthread2,2
       lstart = lia(l)+1
       lend   = lia(l+1)
       do i=lstart,lend ! angles loop
        xrj=vx(1:3,ang(4,i)) ! xn/|rn|
        rnj=vv(1:3,ang(4,i)) ! not normalized equilibrium i-j vector
        rtj=rdat(1:3,ang(2,i))-rdat(1:3,ang(1,i))+rnj ! deformed i-j vector

        xrk=vx(1:3,ang(5,i)) ! xn/|rn|
        rnk=vv(1:3,ang(5,i)) ! not normalized equilibrium i-k vector
        rtk=rdat(1:3,ang(3,i))-rdat(1:3,ang(1,i))+rnk ! deformed i-k vector

        dfi=((xrk(2)*rtk(3)-xrk(3)*rtk(2)-xrj(2)*rtj(3)+xrj(3)*rtj(2)) ! angle change delta phi
     & *anf(1,ang(6,i))
     & -(xrk(1)*rtk(3)-xrk(3)*rtk(1)-xrj(1)*rtj(3)+xrj(3)*rtj(1))
     & *anf(2,ang(6,i))
     & +(xrk(1)*rtk(2)-xrk(2)*rtk(1)-xrj(1)*rtj(2)+xrj(2)*rtj(1))
     & *anf(3,ang(6,i)))*abond(i)

        adat(1:3,ang(2,i))=adat(1:3,ang(2,i))+dfi*anf(4:6,ang(6,i)) ! angular force for point j
        adat(1:3,ang(3,i))=adat(1:3,ang(3,i))+dfi*anf(7:9,ang(6,i)) ! angular force for point k
        adat(1:3,ang(1,i))=adat(1:3,ang(1,i))
     & -(dfi*anf(4:6,ang(6,i))+dfi*anf(7:9,ang(6,i)))
       enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL 

c!$OMP  PARALLEL DO DEFAULT(NONE)
c!$OMP& PRIVATE(i)
c!$OMP& SHARED(nfpoint,adat,fmind,fmacro)
c      do i=1,nfpoint ! impose global force conditions
c        adat(1:3,fmind(i))=adat(1:3,fmind(i))+fmacro(1:3,i)
c      enddo
c!$OMP END PARALLEL DO


      do i=1,ncx
       do k=1,ncz
        l1=mijks(i,1,k)
        l2=mijks(i,ncy-1,k)
        if(l1.gt.0.and.l2.gt.0)then
         adat(1:3,l1)=adat(1:3,l1)+adat(1:3,l2)
         adat(1:3,l2)=adat(1:3,l1)
        endif
       enddo
      enddo

      adat = adat*rsi ! calculate acceleration

! Fourth step: calculate new velocities v_i(t+dt)=v_i(t+dt/2)+a_i(t+dt)dt/2

c      if(iti.le.5)then
!$OMP  PARALLEL DO DEFAULT(NONE)
!$OMP& PRIVATE(i)
!$OMP& SHARED(ns,adat,vdat)
       do i=1,ns
        adat(1:3,i)=adat(1:3,i)-0.5d0*vdat(1:3,i) ! viscous effects
       enddo
!$OMP END PARALLEL DO
c      endif

!$OMP  PARALLEL DO DEFAULT(NONE)
!$OMP& PRIVATE(i)
!$OMP& SHARED(ns,vdat,adat,dt)
      do i=1,ns
       vdat(1:3,i)=vdat(1:3,i)+adat(1:3,i)*dt*0.5d0
      enddo
!$OMP END PARALLEL DO

! Integrate LSM - velocity Verlet algorithm: End
      enddo ! end of time subcircle
! Save output data: coordinates, displacements, velocities, forces, strain and stress tensors

      open(1,file='LSM_coord_field')
      open(2,file='LSM_displ_field')
      open(3,file='LSM_velocity_field')
      open(4,file='LSM_force_field')

      do k=1,ncz
       do j=1,ncy
        do i=1,ncx
         l=mijks(i,j,k) ! solid point number
         if(l.ne.0.d0)then
          write(1,*)rdat(1:3,l) ! save coordinates
          write(2,*)rdat(1:3,l)-icors(1:3,l) ! save displacements
          write(3,*)vdat(1:3,l) ! save velocities
          write(4,*)adat(1:3,l)*rs ! save forces
         else
          write(3,*)(/0.D0,0.D0,0.D0/) ! save velocities
          write(4,*)(/0.D0,0.D0,0.D0/) ! save forces
          write(1,*)(/i,j,k/)          ! save coordinates
          write(2,*)(/0.D0,0.D0,0.D0/) ! save displacements
         endif
        enddo
       enddo
      enddo

      close(1)
      close(2)
      close(3)
      close(4)
! End saving data


      open(11,file='write')
      FTO=(/0.D0,0.D0,0.D0/)
      VTO=(/0.D0,0.D0,0.D0/)
      RTO=(/0.D0,0.D0,0.D0/)
      write(11,*)'iti = ',iti,'itj = ',itj
      do l1 = 1,ns
       FTO=FTO+adat(1:3,l1)*rs
       VTO=VTO+vdat(1:3,l1)
       RTO=RTO+rdat(1:3,l1)
      enddo
      write(11,*)'FTO = ',FTO
      write(11,*)'VTO = ',VTO
      write(11,*)'RTO = ',RTO
      close(11)

      call calc_stress(ns,rdat,icors,bn,vb,vn,vx,vv,anf,
     &                 ncx,ncy,ncz,mijks,ipcx,ipcy,ipcz,
     &                 ind,inds,alf1,bet1,alf2,bet2,Em,nu,gstrain*2.D0)

c save intermediate data
       open(1,file='iteration')
       write(1,*)((rdat(j,i),j=1,3),i=1,ns)
       write(1,*)((vdat(j,i),j=1,3),i=1,ns)
       write(1,*)((adat(j,i),j=1,3),i=1,ns)
       close(1)
       open(1,file='iteration2')
       write(1,*)((rdat(j,i),j=1,3),i=1,ns)
       write(1,*)((vdat(j,i),j=1,3),i=1,ns)
       write(1,*)((adat(j,i),j=1,3),i=1,ns)
       close(1)
       call system('rm iteration')
c end save

      enddo ! end of main time circle

      open(4,file='icors')
      do i=1,ns
       write(4,*)icors(1:3,i)
      enddo
      close(4)

      deallocate(bon,ang,icors,mijks,
     &           adat,vdat,rdat,ind,ipcx,ipcy,ipcz)
 
      stop
      end
c----------------------------------------------------------
c---------------------- subroutine ------------------------
c----------------------------------------------------------
      subroutine calc_stress(ns,rdat,icors,bn,vb,vn,vx,vv,anf,
     &                    ncx,ncy,ncz,mijks,ipcx,ipcy,ipcz,
     &                    ind,inds,alf1,bet1,alf2,bet2,Em,nu,eyz0)

      implicit  real*8 (a-h,o-z)
      implicit  integer*4 (i-n)

      integer checkstat ! check the success of allocate command
      integer*4 ns ! number of solid points
      real*8    rdat(3,ns) ! elastic solid points coordinates
c      real*8    rvirt(3,ns) ! virtual displacement
      integer*4 icors(3,ns) ! coordinates of the solid points
      real*8    bn(3,9) ! normalized vectors
      real*8    vb(3,9) ! bond vectors (real numbers)
      real*8    vn(3,19) ! discrete velocities normalized
      real*8    vx(3,19) ! xn/|rn|
      real*8    vv(3,19) ! discrete velocities (real numbers)
      real*8    anf(9,48) ! angle plane normal, angular force vectors j k; for angular springs

      real*8 flin(3) ! linear spring force
      real*8 xnj(3),xrj(3),rnj(3),rtj(3),xnk(3),xrk(3),rnk(3),rtk(3) ! force calculation vectors
      real*8 xn(3),rn(3),rt(3) ! vectors in linear spring force calculation
      real*8 dfi ! delta phi, angle change for angular force

      integer*4 ncx ! x unit cell size
      integer*4 ncy ! y unit cell size
      integer*4 ncz ! z unit cell size
      integer*4 mijks(ncx,ncy,ncz)

      integer icv(3,8) ! cubic vectors
      integer icp(3,8) ! cube points 

      integer ipcx(ncx+2),ipcy(ncy+2),ipcz(ncz+2) ! periodic boundary conditions
      integer isc ! check solid cube indicator
      integer ind(ncx,ncy-1,ncz) ! index file
      integer inds(ncx,ncy,ncz) ! index file

      integer*4 bon(3,24), ang(6,72) ! bonds and angles for a cube
      real*8    sbond(24), abond(72) ! elastic constants
      real*8    sbondin(24), abondin(72) ! elastic constants
      real*8    alf1,bet1,alf2,bet2              ! elastic constants

      real*8,    dimension(:,:,:), allocatable :: sigmap ! (+->-) plus stress tensors
      real*8,    dimension(:,:,:), allocatable :: sigmam ! (+->-) minus stress tensors
      real*8,    dimension(:,:,:), allocatable :: strain ! local strain tensor 
      real*8 asigmap(3,3),asigmam(3,3),astrain(3,3) ! average tensors
      real*8 sv(3,8) ! cube strain vectors
      real*8 up(3) ! cube point displacement

      real*8 vlayer ! layer volume
      real*8 Em ! Young modulus
      real*8 nu ! Poisson's ratio
      real*8 Cxxxx,Cyyxx,Cxyxy ! stiffness matrix coefficients
      real*8 Exxxx,Eyyxx,Exyxy,Eyzyz ! effective stiffness coefficients
      real*8 frc(3)            ! global force condition
      real*8 sxx0,syy0,eyz0    ! global stress and strains in the layer
      real*8 vcell             ! cell volume    

      Cxxxx = Em*(nu-1.D0)/(2.D0*nu*nu+nu-1.D0) ! stiffness coefficients
      Cyyxx = Em*nu/((1.D0+nu)*(1.D0-2.D0*nu))
      Cxyxy = Em/(2.D0*(1.D0+nu))
c      exx0  = 0.0004D0 
c      vcell = (ncx-2)*ncy*ncz
      vcell = (ncx)*(ncy-2)*(ncz)

      icv(1,1:8)=(/0,1,0,1,0,1,0,1/) ! cubic vectors
      icv(2,1:8)=(/0,0,1,1,0,0,1,1/)
      icv(3,1:8)=(/0,0,0,0,1,1,1,1/)

      sv(1,1:8)=(/-0.5D0,0.5D0,-0.5D0,0.5D0,-0.5D0,0.5D0,-0.5D0,0.5D0/) ! cube strain vectors
      sv(2,1:8)=(/-0.5D0,-0.5D0,0.5D0,0.5D0,-0.5D0,-0.5D0,0.5D0,0.5D0/)
      sv(3,1:8)=(/-0.5D0,-0.5D0,-0.5D0,-0.5D0,0.5D0,0.5D0,0.5D0,0.5D0/)

      open(1,file='cube_bon')        ! read elementary cube data
      read(1,*)((bon(i,j),i=1,3),j=1,24)
      close(1)
      open(1,file='cube_ang')
      read(1,*)((ang(i,j),i=1,6),j=1,72)
      close(1)
      open(1,file='cube_sbond')
      read(1,*)(sbond(i),i=1,24)
      close(1)
      open(1,file='cube_abond')
      read(1,*)(abond(i),i=1,72)
      close(1)
      sbondin = sbond
      abondin = abond

      ncub = 0
      do i=1,ncx
       do j=1,ncy-2
        do k=1,ncz

         isc=0
         do l=1,8  ! extract cube points
          icp(1,l)=ipcx(i+icv(1,l)+1)
          icp(2,l)=ipcy(j+icv(2,l)+1)
          icp(3,l)=ipcz(k+icv(3,l)+1)
         enddo
         isc=ind(i,j,k)

         if(isc.eq.0.or.isc.eq.3)then ! there is a cube, start calculate stress tensor
          ncub=ncub+1
         endif
        enddo
       enddo
      enddo

      allocate(sigmap(3,3,ncub),stat=checkstat) ! create sigmap
      if(checkstat>0) then
       write(*,*)'ERROR! impossible to allocate memory for sigmap'
       stop
      endif

      allocate(sigmam(3,3,ncub),stat=checkstat) ! create sigmam
      if(checkstat>0) then
       write(*,*)'ERROR! impossible to allocate memory for sigmam'
       stop
      endif

      allocate(strain(3,3,ncub),stat=checkstat) ! create strain
      if(checkstat>0) then
       write(*,*)'ERROR! impossible to allocate memory for strain'
       stop
      endif

      open(113,file='stress_tensor',position='append')
      open(333,file='strain')
      open(666,file='stress')
      ncub = 0
      sigmap = 0.D0
      sigmam = 0.D0
      strain = 0.D0
      do i=1,ncx   
       do j=1,ncy-2
        do k=1,ncz

         isc=0
         do l=1,8  ! extract cube points
          icp(1,l)=ipcx(i+icv(1,l)+1) 
          icp(2,l)=ipcy(j+icv(2,l)+1)
          icp(3,l)=ipcz(k+icv(3,l)+1)
         enddo
         isc=ind(i,j,k)
         if (isc.eq.0) then
                           sbond=sbondin*alf1
                           abond=abondin*bet1
         elseif (isc.eq.3) then
                           sbond=sbondin*alf2
                           abond=abondin*bet2
         endif
         if(isc.eq.0.or.isc.eq.3)then ! there is a cube, start calculate strain and stress tensors
          ncub=ncub+1

          do l=1,8 ! calculate strain tensor
           do m=1,3
            do n=1,3
             np=mijks(icp(1,l),icp(2,l),icp(3,l)) ! extract point number
             up=rdat(1:3,np)        ! displacement in the point
             strain(m,n,ncub)=strain(m,n,ncub)+sv(m,l)*up(n) ! calculate strain tensor
     &                                        +sv(n,l)*up(m)
            enddo
           enddo
          enddo 

          do l=1,24 ! bonds loop to calculate linear springs force
           n1=mijks(icp(1,bon(1,l)),icp(2,bon(1,l)),icp(3,bon(1,l)))
           n2=mijks(icp(1,bon(2,l)),icp(2,bon(2,l)),icp(3,bon(2,l)))
           xn=bn(1:3,bon(3,l)) ! normalized equilibrium i-j vector
           rn=vb(1:3,bon(3,l)) ! not normalized equilibrium i-j vector
           rt=rdat(1:3,n2)-rdat(1:3,n1)+rn ! deformed i-j vector
           flin=sbond(l)*((rt(1)*xn(1)+rt(2)*xn(2)+rt(3)*xn(3))*xn-rn) ! linear spring force

           if(xn(1).gt.0)then
            sigmap(1:3,1,ncub)=sigmap(1:3,1,ncub)+flin
            sigmam(1:3,1,ncub)=sigmam(1:3,1,ncub)-flin
           endif
           if(xn(1).lt.0)then
            sigmap(1:3,1,ncub)=sigmap(1:3,1,ncub)-flin
            sigmam(1:3,1,ncub)=sigmam(1:3,1,ncub)+flin
           endif
           if(xn(2).gt.0)then
            sigmap(1:3,2,ncub)=sigmap(1:3,2,ncub)+flin
            sigmam(1:3,2,ncub)=sigmam(1:3,2,ncub)-flin
           endif
           if(xn(2).lt.0)then
            sigmap(1:3,2,ncub)=sigmap(1:3,2,ncub)-flin
            sigmam(1:3,2,ncub)=sigmam(1:3,2,ncub)+flin
           endif
           if(xn(3).gt.0)then
            sigmap(1:3,3,ncub)=sigmap(1:3,3,ncub)+flin
            sigmam(1:3,3,ncub)=sigmam(1:3,3,ncub)-flin
           endif
           if(xn(3).lt.0)then
            sigmap(1:3,3,ncub)=sigmap(1:3,3,ncub)-flin
            sigmam(1:3,3,ncub)=sigmam(1:3,3,ncub)+flin
           endif
          enddo

          do l=1,72 ! angles loop
           n1=mijks(icp(1,ang(1,l)),icp(2,ang(1,l)),icp(3,ang(1,l)))
           n2=mijks(icp(1,ang(2,l)),icp(2,ang(2,l)),icp(3,ang(2,l)))
           n3=mijks(icp(1,ang(3,l)),icp(2,ang(3,l)),icp(3,ang(3,l)))
           xnj=vn(1:3,ang(4,l)) ! normalized equilibrium i-j vector
           xrj=vx(1:3,ang(4,l)) ! xn/|rn|
           rnj=vv(1:3,ang(4,l)) ! not normalized equilibrium i-j vector
           rtj=rdat(1:3,n2)-rdat(1:3,n1)+rnj ! deformed i-j vector

           xnk=vn(1:3,ang(5,l)) ! normalized equilibrium i-k vector
           xrk=vx(1:3,ang(5,l)) ! xn/|rn|
           rnk=vv(1:3,ang(5,l)) ! not normalized equilibrium i-k vector
           rtk=rdat(1:3,n3)-rdat(1:3,n1)+rnk ! deformed i-k vector

           dfi=((xrk(2)*rtk(3)-xrk(3)*rtk(2)
     &          -xrj(2)*rtj(3)+xrj(3)*rtj(2)) ! angle change delta phi
     &     *anf(1,ang(6,l))
     &     -(xrk(1)*rtk(3)-xrk(3)*rtk(1)-xrj(1)*rtj(3)+xrj(3)*rtj(1))
     &     *anf(2,ang(6,l))
     &     +(xrk(1)*rtk(2)-xrk(2)*rtk(1)-xrj(1)*rtj(2)+xrj(2)*rtj(1))
     &     *anf(3,ang(6,l)))*abond(l)


           if(xnj(1).gt.0)then
            sigmap(1:3,1,ncub)=sigmap(1:3,1,ncub)-dfi*anf(4:6,ang(6,l))
            sigmam(1:3,1,ncub)=sigmam(1:3,1,ncub)+dfi*anf(4:6,ang(6,l))
           endif
           if(xnj(1).lt.0)then
            sigmap(1:3,1,ncub)=sigmap(1:3,1,ncub)+dfi*anf(4:6,ang(6,l))
            sigmam(1:3,1,ncub)=sigmam(1:3,1,ncub)-dfi*anf(4:6,ang(6,l))
           endif
           if(xnj(2).gt.0)then
            sigmap(1:3,2,ncub)=sigmap(1:3,2,ncub)-dfi*anf(4:6,ang(6,l))
            sigmam(1:3,2,ncub)=sigmam(1:3,2,ncub)+dfi*anf(4:6,ang(6,l))
           endif
           if(xnj(2).lt.0)then
            sigmap(1:3,2,ncub)=sigmap(1:3,2,ncub)+dfi*anf(4:6,ang(6,l))
            sigmam(1:3,2,ncub)=sigmam(1:3,2,ncub)-dfi*anf(4:6,ang(6,l))
           endif
           if(xnj(3).gt.0)then
            sigmap(1:3,3,ncub)=sigmap(1:3,3,ncub)-dfi*anf(4:6,ang(6,l))
            sigmam(1:3,3,ncub)=sigmam(1:3,3,ncub)+dfi*anf(4:6,ang(6,l))
           endif
           if(xnj(3).lt.0)then
            sigmap(1:3,3,ncub)=sigmap(1:3,3,ncub)+dfi*anf(4:6,ang(6,l))
            sigmam(1:3,3,ncub)=sigmam(1:3,3,ncub)-dfi*anf(4:6,ang(6,l))
           endif

           if(xnk(1).gt.0)then
            sigmap(1:3,1,ncub)=sigmap(1:3,1,ncub)-dfi*anf(7:9,ang(6,l))
            sigmam(1:3,1,ncub)=sigmam(1:3,1,ncub)+dfi*anf(7:9,ang(6,l))
           endif
           if(xnk(1).lt.0)then
            sigmap(1:3,1,ncub)=sigmap(1:3,1,ncub)+dfi*anf(7:9,ang(6,l))
            sigmam(1:3,1,ncub)=sigmam(1:3,1,ncub)-dfi*anf(7:9,ang(6,l))
           endif
           if(xnk(2).gt.0)then
            sigmap(1:3,2,ncub)=sigmap(1:3,2,ncub)-dfi*anf(7:9,ang(6,l))
            sigmam(1:3,2,ncub)=sigmam(1:3,2,ncub)+dfi*anf(7:9,ang(6,l))
           endif
           if(xnk(2).lt.0)then
            sigmap(1:3,2,ncub)=sigmap(1:3,2,ncub)+dfi*anf(7:9,ang(6,l))
            sigmam(1:3,2,ncub)=sigmam(1:3,2,ncub)-dfi*anf(7:9,ang(6,l))
           endif
           if(xnk(3).gt.0)then
            sigmap(1:3,3,ncub)=sigmap(1:3,3,ncub)-dfi*anf(7:9,ang(6,l))
            sigmam(1:3,3,ncub)=sigmam(1:3,3,ncub)+dfi*anf(7:9,ang(6,l))
           endif
           if(xnk(3).lt.0)then
            sigmap(1:3,3,ncub)=sigmap(1:3,3,ncub)+dfi*anf(7:9,ang(6,l))
            sigmam(1:3,3,ncub)=sigmam(1:3,3,ncub)-dfi*anf(7:9,ang(6,l))
           endif
          enddo

c          write(113,*)'ncub   = ',ncub
c          write(113,*)'ijk = ',i,j,k
c          write(113,*)'sigmap = '
c          write(113,*)sigmap(1:3,1:3,ncub)
c          write(113,*)'sigmam = '
c          write(113,*)sigmam(1:3,1:3,ncub)
c          write(113,*)'strain = '
c          write(113,*)strain(1:3,1:3,ncub)*0.25D0

c          write(333,1984)strain(1:3,1,ncub)*0.25D0,
c     &                strain(1:3,2,ncub)*0.25D0,
c     &                strain(1:3,3,ncub)*0.25D0

c          write(666,1984)sigmap(1:3,1,ncub),
c     &                sigmap(1:3,2,ncub),
c     &                sigmap(1:3,3,ncub)


         else
c          write(333,1984)0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0
c          write(666,1984)0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0
         endif

        enddo
       enddo
      enddo
      strain=strain*0.25D0

      close(333)
      close(666)
1984  format(9(1X,E15.8))

      asigmap = 0.D0
      asigmam = 0.D0
      astrain = 0.D0
      do i=1,ncub
       asigmap = asigmap + sigmap(1:3,1:3,i)
       asigmam = asigmam + sigmam(1:3,1:3,i)
       astrain = astrain + strain(1:3,1:3,i)
      enddo
      rncub = ncub

      rncx=ncy-2
      Eyzyz=asigmap(2,3)/(vcell*eyz0/rncx)

      write(113,*)'Elastic medium properties'
      write(113,*)'Ypung modulus Em = ',Em
      write(113,*)'Poisson ratio nu = ',nu
      write(113,*)'Stiffness coefficients'
      write(113,*)'C_xxxx = ',Cxxxx
      write(113,*)'C_yyxx = ',Cyyxx
      write(113,*)'C_xyxy = ',Cxyxy
      write(113,*)'Simulation parameters'
      write(113,*)'Cell volume vcell = ',vcell
      write(113,*)'Solid part volume(elementary cubes) rncub = ',rncub
      write(113,*)'Layer volume (elementary cubes) vlayer = ',vlayer
      write(113,*)'Force per elementary cube vertex frc = ',frc
      write(113,*)'Global strain exx0 = ',eyz0
      write(113,*)'Stress imposed on the layer sxx0 = ',sxx0
      write(113,*)'Stress imposed on the layer syy0 = ',syy0
      write(113,*)'Simulation results'
      write(113,*)'averaged stress tensor over ALL the cell'
      write(113,*)'asigmap = '
      write(113,*)asigmap/rncub
      write(113,*)'asigmam = '
      write(113,*)asigmam/rncub
      write(113,*)'averaged strain'
      write(113,*)'astrain = '
      write(113,*)astrain/rncub
      write(113,*)'Effective stiffness coefficients'
      write(113,*)'Eyzyz = ',Eyzyz
      close(113)

      deallocate(sigmap,sigmam,strain)

      return
      end
c----------------------------------------------------------------------
c apply global stress/strain(by stress) conditions at the domain boundaries
c----------------------------------------------------------------------
      subroutine global_stress_cond(ns,ncx,ncy,ncz,mijks,ipcx,ipcy,ipcz,
     &                    ind,forces,gsigma)

      implicit  real*8 (a-h,o-z)
      implicit  integer*4 (i-n)

      integer*4 ns ! number of solid points

      integer*4 ncx ! x unit cell size
      integer*4 ncy ! y unit cell size
      integer*4 ncz ! z unit cell size
      integer*4 mijks(ncx,ncy,ncz)

      integer icv(3,8) ! cubic vectors
      integer icp(3,8) ! cube points

      integer ipcx(ncx+2),ipcy(ncy+2),ipcz(ncz+2) ! periodic boundary conditions
      integer isc ! check solid cube indicator
      integer ind(ncx,ncy,ncz) ! index file

      integer nrmls(3,6)   ! normals
      real*8  rmls(3,6)    ! normals
      real*8  forces(3,ns) ! surface forces
      integer nsds(6,4)    ! cube sides
      real*8  frc(3)       ! local force
      integer nspc(3,4)    ! side points coordinates
      real*8  gsigma(3,3)  ! global sigma tensor
      real*8  stot(3,3)    ! total stress

      icv(1,1:8)=(/0,1,0,1,0,1,0,1/) ! cubic vectors
      icv(2,1:8)=(/0,0,1,1,0,0,1,1/)
      icv(3,1:8)=(/0,0,0,0,1,1,1,1/)

      nrmls(1,1:6)=(/1,-1,0, 0,0, 0/) ! cube sides normals
      nrmls(2,1:6)=(/0, 0,1,-1,0, 0/)
      nrmls(3,1:6)=(/0, 0,0, 0,1,-1/)

      rmls=nrmls ! real*8

      nsds(1,1:4)=(/2,4,6,8/) ! cube sides
      nsds(2,1:4)=(/1,3,5,7/)
      nsds(3,1:4)=(/3,4,7,8/)
      nsds(4,1:4)=(/1,2,5,6/)
      nsds(5,1:4)=(/5,6,7,8/)
      nsds(6,1:4)=(/1,2,3,4/)

      forces=0.D0
      nsides=0
      stot=0.D0
      do i=1,ncx
       do j=1,ncy
        do k=1,ncz

         isc=0
         do l=1,8  ! extract cube points
          lifg=0
          icp(1,l)=ipcx(i+icv(1,l)+1)
          icp(2,l)=ipcy(j+icv(2,l)+1)
          icp(3,l)=ipcz(k+icv(3,l)+1)
          if (ind(icp(1,l),icp(2,l),icp(3,l)).eq.1) then
             lifg=1
          endif
          isc=isc+lifg  
        enddo

         if(isc.eq.0)then ! there is a cube
          do l=1,6  !(6 possible sides) check if cube side is a surfase
           isg=0
           do m=1,4 !(4 points to check)
            nspc(1,m)=ipcx(icp(1,nsds(l,m))+nrmls(1,l)+1)
            nspc(2,m)=ipcy(icp(2,nsds(l,m))+nrmls(2,l)+1)
            nspc(3,m)=ipcz(icp(3,nsds(l,m))+nrmls(3,l)+1)
c            isg=isg+ind(nspc(1,m),nspc(2,m),nspc(3,m))
           lisg=0
            if(ind(nspc(1,m),nspc(2,m),nspc(3,m)).eq.1) then
              lisg=1
            endif
            isg=isg+lisg
           enddo
           if(isg.ne.0)then ! there is a side
            nsides=nsides+1
            frc=MATMUL(gsigma,nrmls(1:3,l))*0.25D0

            do m=1,4
             nfsp =
     &       mijks(icp(1,nsds(l,m)),icp(2,nsds(l,m)),icp(3,nsds(l,m))) ! number of solid point
             forces(1:3,nfsp)=forces(1:3,nfsp)+frc       ! assign surface forces
              do ii=1,3
               do jj=1,3
                stot(ii,jj)=stot(ii,jj)+frc(ii)*icp(jj,nsds(l,m))
               enddo
              enddo
            enddo

           endif
          enddo
         endif

        enddo
       enddo
      enddo
      write(*,*)'nsides = ',nsides
      write(*,*)'stot   = '
      write(*,*)stot

      return
      end
c----------------------------------------------------------------------
c get nfpoint
c----------------------------------------------------------------------
      subroutine x_get_nfpoint(ns,ncx,ncy,ncz,mijks,ipcx,ipcy,ipcz,
     &                    ind,frc,nfpoint)

      implicit  real*8 (a-h,o-z)
      implicit  integer*4 (i-n)

      integer*4 ns ! number of solid points

      integer*4 ncx ! x unit cell size
      integer*4 ncy ! y unit cell size
      integer*4 ncz ! z unit cell size
      integer*4 mijks(ncx,ncy,ncz)

      integer icv(3,8) ! cubic vectors
      integer icp(3,8) ! cube points

      integer ipcx(ncx+2),ipcy(ncy+2),ipcz(ncz+2) ! periodic boundary conditions
      integer isc ! check solid cube indicator
      integer ind(ncx,ncy,ncz) ! index file

      real*8  forces(3,ns) ! surface forces
      real*8  frc(3)       ! local force
      integer nfpoint ! number of force points

      icv(1,1:8)=(/0,1,0,1,0,1,0,1/) ! cubic vectors
      icv(2,1:8)=(/0,0,1,1,0,0,1,1/)
      icv(3,1:8)=(/0,0,0,0,1,1,1,1/)

      forces=0.D0

      do i=1,1     ! apply stretching to the first layer
       do j=1,ncy
        do k=1,ncz

         isc=0
         do l=1,8  ! extract cube points
          icp(1,l)=ipcx(i+icv(1,l)+1)
          icp(2,l)=ipcy(j+icv(2,l)+1)
          icp(3,l)=ipcz(k+icv(3,l)+1)
          lifg=0
          if (ind(icp(1,l),icp(2,l),icp(3,l)).eq.1) then
             lifg=1
          endif
          isc=isc+lifg
c          isc=isc+ind(icp(1,l),icp(2,l),icp(3,l))
         enddo

         if(isc.eq.0)then ! there is a cube
          nfsp = mijks(icp(1,1),icp(2,1),icp(3,1)) ! stretching along x axis
          forces(1:3,nfsp)=forces(1:3,nfsp)-frc
          nfsp = mijks(icp(1,2),icp(2,2),icp(3,2))
          forces(1:3,nfsp)=forces(1:3,nfsp)+frc
          nfsp = mijks(icp(1,3),icp(2,3),icp(3,3))
          forces(1:3,nfsp)=forces(1:3,nfsp)-frc
          nfsp = mijks(icp(1,4),icp(2,4),icp(3,4))
          forces(1:3,nfsp)=forces(1:3,nfsp)+frc
          nfsp = mijks(icp(1,5),icp(2,5),icp(3,5))
          forces(1:3,nfsp)=forces(1:3,nfsp)-frc
          nfsp = mijks(icp(1,6),icp(2,6),icp(3,6))
          forces(1:3,nfsp)=forces(1:3,nfsp)+frc
          nfsp = mijks(icp(1,7),icp(2,7),icp(3,7))
          forces(1:3,nfsp)=forces(1:3,nfsp)-frc
          nfsp = mijks(icp(1,8),icp(2,8),icp(3,8))
          forces(1:3,nfsp)=forces(1:3,nfsp)+frc
         endif

        enddo
       enddo
      enddo

      nfpoint=0
      do i=1,ns ! calculate the number of non-zero force points
       if(forces(1,i).ne.0.D0)then   ! here check is only for the x-component !!!
        nfpoint = nfpoint+1
       endif
      enddo

      return
      end
c----------------------------------------------------------------------
c calculate global boundary conditions (forces)
c----------------------------------------------------------------------
      subroutine x_stretching(ns,ncx,ncy,ncz,mijks,ipcx,ipcy,ipcz,
     &                    ind,fmacro,fmind,frc,nfpoint,vlayer)

      implicit  real*8 (a-h,o-z)
      implicit  integer*4 (i-n)

      integer*4 ns ! number of solid points

      integer*4 ncx ! x unit cell size
      integer*4 ncy ! y unit cell size
      integer*4 ncz ! z unit cell size
      integer*4 mijks(ncx,ncy,ncz)

      integer icv(3,8) ! cubic vectors
      integer icp(3,8) ! cube points

      integer ipcx(ncx+2),ipcy(ncy+2),ipcz(ncz+2) ! periodic boundary conditions
      integer isc ! check solid cube indicator
      integer ind(ncx,ncy,ncz) ! index file

      real*8  forces(3,ns) ! surface forces
      real*8  frc(3)       ! local force
      integer nfpoint ! number of force points
      real*8  fmacro(3,nfpoint) ! global forces
      integer*4 fmind(nfpoint)  ! global forces index 
      real*8 vlayer             ! layer volume

      icv(1,1:8)=(/0,1,0,1,0,1,0,1/) ! cubic vectors
      icv(2,1:8)=(/0,0,1,1,0,0,1,1/)
      icv(3,1:8)=(/0,0,0,0,1,1,1,1/)

      forces=0.D0
      vlayer=0.D0

      do i=1,1     ! apply stretching to the first layer
       do j=1,ncy
        do k=1,ncz

         isc=0
         do l=1,8  ! extract cube points
          icp(1,l)=ipcx(i+icv(1,l)+1)
          icp(2,l)=ipcy(j+icv(2,l)+1)
          icp(3,l)=ipcz(k+icv(3,l)+1)
          lifg=0
          if (ind(icp(1,l),icp(2,l),icp(3,l)).eq.1) then
             lifg=1
          endif
          ics=ics+lifg
c          isc=isc+ind(icp(1,l),icp(2,l),icp(3,l))
         enddo

         if(isc.eq.0)then ! there is a cube
          vlayer=vlayer+1.D0
          nfsp = mijks(icp(1,1),icp(2,1),icp(3,1)) ! stretching along x axis
          forces(1:3,nfsp)=forces(1:3,nfsp)-frc
          nfsp = mijks(icp(1,2),icp(2,2),icp(3,2))
          forces(1:3,nfsp)=forces(1:3,nfsp)+frc
          nfsp = mijks(icp(1,3),icp(2,3),icp(3,3))
          forces(1:3,nfsp)=forces(1:3,nfsp)-frc
          nfsp = mijks(icp(1,4),icp(2,4),icp(3,4))
          forces(1:3,nfsp)=forces(1:3,nfsp)+frc
          nfsp = mijks(icp(1,5),icp(2,5),icp(3,5))
          forces(1:3,nfsp)=forces(1:3,nfsp)-frc
          nfsp = mijks(icp(1,6),icp(2,6),icp(3,6))
          forces(1:3,nfsp)=forces(1:3,nfsp)+frc
          nfsp = mijks(icp(1,7),icp(2,7),icp(3,7))
          forces(1:3,nfsp)=forces(1:3,nfsp)-frc
          nfsp = mijks(icp(1,8),icp(2,8),icp(3,8))
          forces(1:3,nfsp)=forces(1:3,nfsp)+frc
         endif

        enddo
       enddo
      enddo

      nfind=0
      do i=1,ns ! calculate the non-zero force points
       if(forces(1,i).ne.0.D0)then   ! here check is only for the x-component !!!
        nfind            = nfind + 1
        fmacro(1:3,nfind)= forces(1:3,i)
        fmind(nfind)     = i
       endif
      enddo

      return
      end
