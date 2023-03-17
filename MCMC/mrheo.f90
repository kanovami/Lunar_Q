MODULE mtype
!==================================================
! Data type module
!==================================================
IMPLICIT NONE
INTEGER, PARAMETER  :: DP  = KIND(1.0D0)
INTEGER, PARAMETER  :: DPC = KIND((1.0D0,1.0D0))
INTEGER, PARAMETER  :: I4B = SELECTED_INT_KIND(9)
END MODULE mtype
!==================================================

MODULE mconst
!==================================================
! Internal constants module
!==================================================
USE mtype
IMPLICIT NONE

! Physical constants - overwritten by Scipy
REAL(DP) :: kappa=6.674d-11,pi=4*atan(1.d0)

! Internal structure of planet 1 - can be changed by config file
INTEGER  :: iellid=0  ! 0-VE lithosphere, 1-elastic lihosphere
INTEGER  :: nlit=1    ! # of layers the litosphere consists of

! Common variables
!==================
! Numerical Recipes
INTEGER                  :: neq,np
REAL(DP)                 :: d
INTEGER,ALLOCATABLE      :: indx(:)
COMPLEX(DPC),ALLOCATABLE :: a(:,:),b(:)

! From rheo.in
INTEGER                  :: iost,nl
INTEGER,ALLOCATABLE      :: ire(:)
REAL(DP),ALLOCATABLE     :: rho(:),arad(:),eta(:),mu(:),alpha(:),zeta(:)
REAL(DP),ALLOCATABLE     :: eta_kv_rel(:),cdef(:)

! Complex rigidity
COMPLEX(DPC),ALLOCATABLE :: muC(:)
REAL(DP)                 :: ksi,grav

END MODULE mconst
!==================================================

MODULE mrheo
!===================================================
! Rheology module. Normal mode theory calculation
! of the deformation of a layered viscoelastic sphere
! + Love numbers for both differentiated and homo-
! geneous case.
!----------------------------
! Subroutines and functions:
!----------------------------
! peltier, compliance, fImk2, fRek2,
! gammln, profgen
!===================================================
USE mtype
USE mconst
IMPLICIT NONE

CONTAINS

!===============================================
! Love number k_j of layered viscoelastic sphere
! + complex rigidity muC(il)
!===============================================
SUBROUTINE peltier(j,freq,kj,hj)
USE mtype
USE mconst
IMPLICIT NONE
!-----------------------------------------------
! Variables
!-----------------------------------------------
INTEGER, INTENT(IN)      :: j    ! SH degree
REAL(DP),INTENT(IN)      :: freq ! loading frequency
COMPLEX(DP),INTENT(OUT)  :: kj ,hj  ! tidal Love number

INTEGER               :: ivr,ieq,i
COMPLEX(DPC)          :: JC
!-----------------------------------------------

neq = 6*nl-3
np  = neq
ALLOCATE(a(np,np),b(neq),indx(neq))
ALLOCATE(muC(size(mu)))

kj = dcmplx(0)
hj = dcmplx(0)

! Fill the array of complex rigidities
do i=1,nl
  if (abs(eta(i)).le.epsilon(eta(i))) then
    muC(i)=dcmplx(mu(i))
  else
    call compliance(JC,abs(freq),i)
    if (freq.gt.0._DP) then
      muC(i)  = 1._DP/JC
    else
      muC(i)  = dconjg(1._DP/JC)
    end if
  end if
end do

! Elastic lithosphere
if (iellid==1) then
  do i=nl,nl-nlit+1,-1
    muC(i) = dcmplx(mu(i))
  end do
end if

a=0._DP
b=0._DP

!--------------------------------------------------------------------------
! FILLING THE MATRIX
!--------------------------------------------------------------------------
IF (nl==1) THEN ! homogeneous case
  call fill_matrix_hom(j)
ELSE
  ! i=1
  !-----
  call fill_matrix1(j)

  ! other boundaries
  !------------------
  do i=2,nl-1
    grav = (arad(i-1)/arad(i))**2*grav &
    &    + kappa*4._DP*pi*rho(i)*(arad(i)**3-arad(i-1)**3)/3./arad(i)**2
    ksi  = grav/arad(i)
  
    call fill_matrix(j,i,0) ! lower block
    call fill_matrix(j,i,1) ! upper block
  end do

  ! i=nl
  !------
  call fill_matrixN(j)
END IF
!--------------------------------------------------------------------------

! Fill RHS
!----------
b = 0._DP
b(neq) = -(2*j+1)/arad(nl)

! LU decomposition
!------------------
call ludcmp(a,neq,np,indx,d)
call lubksb(a,neq,np,indx,b)

! Love numbers k_j and h_j
!--------------------------
if (nl==1) then ! homogeneous case
  kj = -b(3)*arad(nl)**j-dcmplx(1._DP)
  
  hj = b(1)*j*arad(nl)**(j+1)/2./dfloat(2*j+3) + &
  &    b(2)*arad(nl)**(j-1)
  hj = hj*grav
else ! multilayered case
  kj = -(b(neq-3)*arad(nl)**j+b(neq)*arad(nl)**(-j-1))-dcmplx(1._DP)
  
  hj = b(neq-5)*j*arad(nl)**(j+1)/2./dfloat(2*j+3) + &
  &    b(neq-4)*arad(nl)**(j-1) + &
  &    b(neq-2)*(j+1)*arad(nl)**(-j)/2./dfloat(2*j-1) + &
  &    b(neq-1)*arad(nl)**(-j-2)
  hj = hj*grav
end if

DEALLOCATE(a,b,indx,muC)

END SUBROUTINE peltier
!===============================================

!======================================
! Compute compliance
!======================================
SUBROUTINE compliance(Jc,freq,ilayer)
USE mtype
USE mconst
IMPLICIT NONE

INTEGER   :: ilayer
REAL(DP)  :: freq
REAL(DP)  :: sina,cosa,taum,denom
REAL(DP)  :: rej,imj
COMPLEX(DPC) :: Jc

if (abs(mu(ilayer)).le.epsilon(mu(ilayer))) then ! viscous limit
  Jc = dcmplx(0._DP,-1._DP/freq/eta(ilayer))
else if (abs(freq).le.epsilon(freq)) then ! fluid limit
  !Jc = dcmplx(1.d10,0._DP)
  Jc = dcmplx(1._DP/mu(ilayer),-1.d30)
else
  taum=eta(ilayer)/mu(ilayer)
  IF (ire(ilayer).eq.0) THEN      ! Maxwell
    rej=1._DP
    imj=-1._DP/(freq*taum)
    Jc = dcmplx(rej,imj)/mu(ilayer)
  ELSE IF (ire(ilayer).eq.1) THEN ! Andrade
    sina=sin(0.5_DP*alpha(ilayer)*pi)
    cosa=cos(0.5_DP*alpha(ilayer)*pi)
    imj=-1._DP/(freq*taum)-(freq*zeta(ilayer)*taum)**(-alpha(ilayer))*sina*exp(gammln(alpha(ilayer)+1._DP))
    rej=1._DP+(freq*zeta(ilayer)*taum)**(-alpha(ilayer))*cosa*exp(gammln(alpha(ilayer)+1._DP))
    Jc = dcmplx(rej,imj)/mu(ilayer)
  ELSE IF (ire(ilayer).eq.2) THEN ! Sundberg-Cooper
    sina=sin(0.5_DP*alpha(ilayer)*pi)
    cosa=cos(0.5_DP*alpha(ilayer)*pi)
    denom = 1._DP+(cdef(ilayer)/mu(ilayer)*eta_kv_rel(ilayer)*eta(ilayer)*freq)**2
    imj=-1._DP/(freq*taum)-(freq*zeta(ilayer)*taum)**(-alpha(ilayer))*sina*exp(gammln(alpha(ilayer)+1._DP)) &
      &      - cdef(ilayer)**2/mu(ilayer)*eta_kv_rel(ilayer)*eta(ilayer)*freq/denom
    rej=1._DP+(freq*zeta(ilayer)*taum)**(-alpha(ilayer))*cosa*exp(gammln(alpha(ilayer)+1._DP)) &
      &      + cdef(ilayer)/denom
    Jc = dcmplx(rej,imj)/mu(ilayer)
  ELSE
    write(*,*) "Unknown rheology in subroutine 'compliance'."
  END IF
end if

END SUBROUTINE compliance
!======================================

!======================================
! NR: logarithm of gamma function
!======================================
FUNCTION gammln(xx)
USE mtype
USE mconst
IMPLICIT NONE

REAL(DP)  :: gammln,xx
INTEGER   :: j
REAL(DP)  :: ser,stp,tmp,x,y,cof(6)
SAVE cof,stp
DATA cof,stp/76.18009172947146_DP,-86.50532032941677_DP, &
  24.01409824083091_DP,-1.231739572450155_DP,.1208650973866179d-2, &
  -.5395239384953d-5,2.5066282746310005_DP/
x=xx
y=x
tmp=x+5.5_DP
tmp=(x+0.5_DP)*log(tmp)-tmp
ser=1.000000000190015_DP
do j=1,6
  y=y+1.d0
  ser=ser+cof(j)/y
end do
gammln=tmp+log(stp*ser/x)
END FUNCTION gammln
!======================================

!========================================
! Fill matrix A:
!========================================
! ilr   - layer identifier
! icomp - 0=incompressible model
!       - 1=compressible model
! iupp  - 0=lower layer (ilr)
!       - 1=upper layer (ilr+1)
!========================================
SUBROUTINE fill_matrix(j,ilr,iupp)
USE mtype
USE mconst
IMPLICIT NONE

INTEGER :: j,ilr,iupp
INTEGER :: sgn,ivr,i,ieq

if (iupp==1) then
  sgn=1
else if (iupp==0) then
  sgn=-1
else
  stop "incorrect 'iupp' in 'fill_matrix'"
end if

ivr  = 6*(ilr-1)-3 + 6*iupp
i    = ilr

! Continuity of radial displacement
ieq  = 6*(i-1)+1
a(ieq,ivr+1)  = sgn*j*arad(i)**(j+1)/(2._DP*(2*j+3))
a(ieq,ivr+2)  = sgn*arad(i)**(j-1)
a(ieq,ivr+4)  = sgn*(j+1)/arad(i)**j/(2._DP*(2*j-1))
a(ieq,ivr+5)  = sgn*1._DP/arad(i)**(j+2)

! Continuity of tangential displacement
ieq  = 6*(i-1)+2
a(ieq,ivr+1)  = sgn*(j+3)*arad(i)**(j+1)/(2._DP*(2*j+3)*(j+1))
a(ieq,ivr+2)  = sgn*arad(i)**(j-1)/dfloat(j)
a(ieq,ivr+4)  = sgn*dfloat(2-j)/arad(i)**j/(2._DP*j*(2*j-1))
a(ieq,ivr+5)  = -sgn*1._DP/arad(i)**(j+2)/dfloat(j+1)
  
! Continuity of radial traction
ieq  = 6*(i-1)+3
a(ieq,ivr+1)  = rho(i+iupp)*ksi*j*arad(i)**(j+2) + muC(i+iupp)*2._DP*(j**2-j-3)*arad(i)**j
a(ieq,ivr+1)  = sgn*a(ieq,ivr+1)/(2._DP*(2*j+3))
a(ieq,ivr+2)  = rho(i+iupp)*ksi*arad(i)**j + 2._DP*muC(i+iupp)*(j-1)*arad(i)**(j-2)
a(ieq,ivr+2)  = sgn*a(ieq,ivr+2)
a(ieq,ivr+3)  = sgn*rho(i+iupp)*arad(i)**j
a(ieq,ivr+4)  = rho(i+iupp)*ksi*(j+1)/arad(i)**(j-1) - muC(i+iupp)*2._DP*(j**2+3*j-1)/arad(i)**(j+1)
a(ieq,ivr+4)  = sgn*a(ieq,ivr+4)/(2._DP*(2*j-1))
a(ieq,ivr+5)  = rho(i+iupp)*ksi/arad(i)**(j+1) - 2._DP*muC(i+iupp)*(j+2)/arad(i)**(j+3)
a(ieq,ivr+5)  = sgn*a(ieq,ivr+5)
a(ieq,ivr+6)  = sgn*rho(i+iupp)/arad(i)**(j+1)

! Continuity of tangential traction
ieq  = 6*(i-1)+4
a(ieq,ivr+1)  = sgn*muC(i+iupp)*j*(j+2)*arad(i)**j/dfloat((2*j+3)*(j+1))
a(ieq,ivr+2)  = sgn*2._DP*muC(i+iupp)*(j-1)*arad(i)**(j-2)/dfloat(j)
a(ieq,ivr+4)  = sgn*muC(i+iupp)*(j**2-1)/arad(i)**(j+1)/dfloat(j*(2*j-1))
a(ieq,ivr+5)  = sgn*2._DP*muC(i+iupp)*(j+2)/arad(i)**(j+3)/dfloat(j+1)

! Continuity of gravity potential
ieq  = 6*(i-1)+5
a(ieq,ivr+3)  = sgn*arad(i)**j
a(ieq,ivr+6)  = sgn*arad(i)**(-j-1)

! Discontinuity of its first derivative (due to surface density)
ieq  = 6*(i-1)+6
a(ieq,ivr+1)  = sgn*4._DP*pi*kappa*rho(i+iupp)*j*arad(i)**(j+1)/(2._DP*(2*j+3))
a(ieq,ivr+2)  = sgn*4._DP*pi*kappa*rho(i+iupp)*arad(i)**(j-1)
a(ieq,ivr+3)  = sgn*arad(i)**(j-1)*(2*j+1)
a(ieq,ivr+4)  = sgn*4._DP*pi*kappa*rho(i+iupp)*(j+1)/arad(i)**j/(2._DP*(2*j-1))
a(ieq,ivr+5)  = sgn*4._DP*pi*kappa*rho(i+iupp)/arad(i)**(j+2)

END SUBROUTINE fill_matrix

!========================================
! Fill matrix A - innermost layer only
!========================================
SUBROUTINE fill_matrix1(j)
USE mtype
USE mconst
IMPLICIT NONE

INTEGER :: j

grav = kappa*4._DP*pi*arad(1)*rho(1)/3.
ksi  = grav/arad(1)

! Continuity of radial displacement
a(1,1)  = -j*arad(1)**(j+1)/(2._DP*(2*j+3))
a(1,2)  = -arad(1)**(j-1)
a(1,4)  = j*arad(1)**(j+1)/(2._DP*(2*j+3))
a(1,5)  = arad(1)**(j-1)
a(1,7)  = (j+1)/arad(1)**j/(2._DP*(2*j-1))
a(1,8)  = 1._DP/arad(1)**(j+2)

! Continuity of tangential displacement
a(2,1)  = -(j+3)*arad(1)**(j+1)/(2._DP*(2*j+3)*(j+1))
a(2,2)  = -arad(1)**(j-1)/dfloat(j)
a(2,4)  = (j+3)*arad(1)**(j+1)/(2._DP*(2*j+3)*(j+1))
a(2,5)  = arad(1)**(j-1)/dfloat(j)
a(2,7)  = dfloat(2-j)/arad(1)**j/(2._DP*j*(2*j-1))
a(2,8)  = -1._DP/arad(1)**(j+2)/dfloat(j+1)

! Continuity of radial traction
a(3,1)  = rho(1)*ksi*j*arad(1)**(j+2) + muC(1)*2._DP*(j**2-j-3)*arad(1)**j
a(3,1)  = -a(3,1)/(2._DP*(2*j+3))
a(3,2)  = rho(1)*ksi*arad(1)**j + 2._DP*muC(1)*(j-1)*arad(1)**(j-2)
a(3,2)  = -a(3,2)
a(3,3)  = -rho(1)*arad(1)**j
a(3,4)  = rho(2)*ksi*j*arad(1)**(j+2) + muC(2)*2._DP*(j**2-j-3)*arad(1)**j
a(3,4)  = a(3,4)/(2._DP*(2*j+3))
a(3,5)  = rho(2)*ksi*arad(1)**j + 2._DP*muC(2)*(j-1)*arad(1)**(j-2)
a(3,6)  = rho(2)*arad(1)**j
a(3,7)  = rho(2)*ksi*(j+1)/arad(1)**(j-1) - muC(2)*2._DP*(j**2+3*j-1)/arad(1)**(j+1)
a(3,7)  = a(3,7)/(2._DP*(2*j-1))
a(3,8)  = rho(2)*ksi/arad(1)**(j+1) - 2._DP*muC(2)*(j+2)/arad(1)**(j+3)
a(3,9)  = rho(2)/arad(1)**(j+1)

! Continuity of tangential traction
a(4,1)  = -muC(1)*j*(j+2)*arad(1)**j/dfloat((2*j+3)*(j+1))
a(4,2)  = -2._DP*muC(1)*(j-1)*arad(1)**(j-2)/dfloat(j)
a(4,4)  = muC(2)*j*(j+2)*arad(1)**j/dfloat((2*j+3)*(j+1))
a(4,5)  = 2._DP*muC(2)*(j-1)*arad(1)**(j-2)/dfloat(j)
a(4,7)  = muC(2)*(j**2-1)/arad(1)**(j+1)/dfloat(j*(2*j-1))
a(4,8)  = 2._DP*muC(2)*(j+2)/arad(1)**(j+3)/dfloat(j+1)

! Continuity of gravity potential
a(5,3)  = -arad(1)**j
a(5,6)  = arad(1)**j
a(5,9)  = arad(1)**(-j-1)

! Discontinuity of its first derivative (due to surface density)
a(6,1)  = -4._DP*pi*kappa*rho(1)*j*arad(1)**(j+1)/(2._DP*(2*j+3))
a(6,2)  = -4._DP*pi*kappa*rho(1)*arad(1)**(j-1)
a(6,3)  = -arad(1)**(j-1)*(2*j+1)
a(6,4)  = 4._DP*pi*kappa*rho(2)*j*arad(1)**(j+1)/(2._DP*(2*j+3))
a(6,5)  = 4._DP*pi*kappa*rho(2)*arad(1)**(j-1)
a(6,6)  = arad(1)**(j-1)*(2*j+1)
a(6,7)  = 4._DP*pi*kappa*rho(2)*(j+1)/arad(1)**j/(2._DP*(2*j-1))
a(6,8)  = 4._DP*pi*kappa*rho(2)/arad(1)**(j+2)

END SUBROUTINE fill_matrix1

!========================================
! Fill matrix A - uppermost layer only
!========================================
SUBROUTINE fill_matrixN(j)
USE mtype
USE mconst
IMPLICIT NONE

INTEGER :: j,ivr,ieq

ivr  = 6*nl-9
grav = (arad(nl-1)/arad(nl))**2*grav &
&     + kappa*4._DP*pi*rho(nl)*(arad(nl)**3-arad(nl-1)**3)/3./arad(nl)**2
ksi  = grav/arad(nl)

! Continuity of radial traction
ieq  = 6*nl-5
a(ieq,ivr+1)  = rho(nl)*ksi*j*arad(nl)**(j+2) + muC(nl)*2._DP*(j**2-j-3)*arad(nl)**j
a(ieq,ivr+1)  = a(ieq,ivr+1)/(2._DP*(2*j+3))
a(ieq,ivr+2)  = rho(nl)*ksi*arad(nl)**j + 2._DP*muC(nl)*(j-1)*arad(nl)**(j-2)
a(ieq,ivr+3)  = rho(nl)*arad(nl)**j
a(ieq,ivr+4)  = rho(nl)*ksi*(j+1)/arad(nl)**(j-1) - muC(nl)*2._DP*(j**2+3*j-1)/arad(nl)**(j+1)
a(ieq,ivr+4)  = a(ieq,ivr+4)/(2._DP*(2*j-1))
a(ieq,ivr+5)  = rho(nl)*ksi/arad(nl)**(j+1) - 2._DP*muC(nl)*(j+2)/arad(nl)**(j+3)
a(ieq,ivr+6)  = rho(nl)/arad(nl)**(j+1)

! Continuity of tangential traction
ieq  = 6*nl-4
a(ieq,ivr+1)  = muC(nl)*j*(j+2)*arad(nl)**j/dfloat((2*j+3)*(j+1))
a(ieq,ivr+2)  = 2._DP*muC(nl)*(j-1)*arad(nl)**(j-2)/dfloat(j)
a(ieq,ivr+4)  = muC(nl)*(j**2-1)/arad(nl)**(j+1)/dfloat(j*(2*j-1))
a(ieq,ivr+5)  = 2._DP*muC(nl)*(j+2)/arad(nl)**(j+3)/dfloat(j+1)

! Discontinuity of its first derivative (due to surface density)
ieq  = 6*nl-3
a(ieq,ivr+1)  = 4._DP*pi*kappa*rho(nl)*j*arad(nl)**(j+1)/(2._DP*(2*j+3))
a(ieq,ivr+2)  = 4._DP*pi*kappa*rho(nl)*arad(nl)**(j-1)
a(ieq,ivr+3)  = arad(nl)**(j-1)*(2*j+1)
a(ieq,ivr+4)  = 4._DP*pi*kappa*rho(nl)*(j+1)/arad(nl)**j/(2._DP*(2*j-1))
a(ieq,ivr+5)  = 4._DP*pi*kappa*rho(nl)/arad(nl)**(j+2)

END SUBROUTINE fill_matrixN

!========================================
! Fill matrix A - homogeneous model
!========================================
SUBROUTINE fill_matrix_hom(j)
USE mtype
USE mconst
IMPLICIT NONE

INTEGER :: j

grav = kappa*4._DP*pi*arad(1)*rho(1)/3.
ksi  = grav/arad(1)

! Continuity of radial traction
a(1,1)  = rho(1)*ksi*j*arad(1)**(j+2) + muC(1)*2._DP*(j**2-j-3)*arad(1)**j
a(1,1)  = a(1,1)/(2._DP*(2*j+3))
a(1,2)  = rho(1)*ksi*arad(1)**j + 2._DP*muC(1)*(j-1)*arad(1)**(j-2)
a(1,3)  = rho(1)*arad(1)**j

! Continuity of tangential traction
a(2,1)  = muC(1)*j*(j+2)*arad(1)**j/dfloat((2*j+3)*(j+1))
a(2,2)  = 2._DP*muC(1)*(j-1)*arad(1)**(j-2)/dfloat(j)

! Discontinuity of its first derivative (due to surface density)
a(3,1)  = 4._DP*pi*kappa*rho(1)*j*arad(1)**(j+1)/(2._DP*(2*j+3))
a(3,2)  = 4._DP*pi*kappa*rho(1)*arad(1)**(j-1)
a(3,3)  = arad(1)**(j-1)*(2*j+1)

END SUBROUTINE fill_matrix_hom

!========================================

END MODULE mrheo