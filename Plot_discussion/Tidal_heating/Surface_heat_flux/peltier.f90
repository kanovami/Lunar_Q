SUBROUTINE peltier(j,freq,k2,muC)
USE mconst
IMPLICIT NONE

!==========================================================================
! VARIABLES
!==========================================================================
INTEGER               :: ivr,ieq,i
INTEGER               :: j ! SH degree
REAL(DP)              :: freq,ksi,grav
COMPLEX(DPC)          :: k2,JC,trn,ttn,dru
COMPLEX(DPC)          :: muC(nl)
!==========================================================================

do i=1,nl-1
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
if (abs(eta(nl)).le.epsilon(eta(nl))) then
  muC(nl)=dcmplx(mu(nl))
else
  call compliance(JC,abs(freq),nl)
  if (freq.gt.0._DP) then
    muC(i)  = 1._DP/JC
  else
    muC(i)  = dconjg(1._DP/JC)
  end if
end if

a=0._DP
b=0._DP

!--------------------------------------------------------------------------
! FILLING THE MATRIX
!--------------------------------------------------------------------------
! i=1
!-----
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

! other boundaries
!------------------
do i=2,nl-1
  ivr  = 6*(i-1)-3
  grav = (arad(i-1)/arad(i))**2*grav &
 &     + kappa*4._DP*pi*rho(i)*(arad(i)**3-arad(i-1)**3)/3./arad(i)**2
  ksi  = grav/arad(i)

  ! Continuity of radial displacement
  ieq  = 6*(i-1)+1
  a(ieq,ivr+1)  = -j*arad(i)**(j+1)/(2._DP*(2*j+3))
  a(ieq,ivr+2)  = -arad(i)**(j-1)
  a(ieq,ivr+4)  = -(j+1)/arad(i)**j/(2._DP*(2*j-1))
  a(ieq,ivr+5)  = -1._DP/arad(i)**(j+2)
  a(ieq,ivr+7)  = j*arad(i)**(j+1)/(2._DP*(2*j+3))
  a(ieq,ivr+8)  = arad(i)**(j-1)
  a(ieq,ivr+10) = (j+1)/arad(i)**j/(2._DP*(2*j-1))
  a(ieq,ivr+11) = 1._DP/arad(i)**(j+2)

  ! Continuity of tangential displacement
  ieq  = 6*(i-1)+2
  a(ieq,ivr+1)  = -(j+3)*arad(i)**(j+1)/(2._DP*(2*j+3)*(j+1))
  a(ieq,ivr+2)  = -arad(i)**(j-1)/dfloat(j)
  a(ieq,ivr+4)  = -dfloat(2-j)/arad(i)**j/(2._DP*j*(2*j-1))
  a(ieq,ivr+5)  = 1._DP/arad(i)**(j+2)/dfloat(j+1)
  a(ieq,ivr+7)  = (j+3)*arad(i)**(j+1)/(2._DP*(2*j+3)*(j+1))
  a(ieq,ivr+8)  = arad(i)**(j-1)/dfloat(j)
  a(ieq,ivr+10) = dfloat(2-j)/arad(i)**j/(2._DP*j*(2*j-1))
  a(ieq,ivr+11) = -1._DP/arad(i)**(j+2)/dfloat(j+1)

  ! Continuity of radial traction
  ieq  = 6*(i-1)+3
  a(ieq,ivr+1)  = rho(i)*ksi*j*arad(i)**(j+2) + muC(i)*2._DP*(j**2-j-3)*arad(i)**j
  a(ieq,ivr+1)  = -a(ieq,ivr+1)/(2._DP*(2*j+3))
  a(ieq,ivr+2)  = rho(i)*ksi*arad(i)**j + 2._DP*muC(i)*(j-1)*arad(i)**(j-2)
  a(ieq,ivr+2)  = -a(ieq,ivr+2)
  a(ieq,ivr+3)  = -rho(i)*arad(i)**j
  a(ieq,ivr+4)  = rho(i)*ksi*(j+1)/arad(i)**(j-1) - muC(i)*2._DP*(j**2+3*j-1)/arad(i)**(j+1)
  a(ieq,ivr+4)  = -a(ieq,ivr+4)/(2._DP*(2*j-1))
  a(ieq,ivr+5)  = rho(i)*ksi/arad(i)**(j+1) - 2._DP*muC(i)*(j+2)/arad(i)**(j+3)
  a(ieq,ivr+5)  = -a(ieq,ivr+5)
  a(ieq,ivr+6)  = -rho(i)/arad(i)**(j+1)
  a(ieq,ivr+7)  = rho(i+1)*ksi*j*arad(i)**(j+2) + muC(i+1)*2._DP*(j**2-j-3)*arad(i)**j
  a(ieq,ivr+7)  = a(ieq,ivr+7)/(2._DP*(2*j+3))
  a(ieq,ivr+8)  = rho(i+1)*ksi*arad(i)**j + 2._DP*muC(i+1)*(j-1)*arad(i)**(j-2)
  a(ieq,ivr+9)  = rho(i+1)*arad(i)**j
  a(ieq,ivr+10) = rho(i+1)*ksi*(j+1)/arad(i)**(j-1) - muC(i+1)*2._DP*(j**2+3*j-1)/arad(i)**(j+1)
  a(ieq,ivr+10) = a(ieq,ivr+10)/(2._DP*(2*j-1))
  a(ieq,ivr+11) = rho(i+1)*ksi/arad(i)**(j+1) - 2._DP*muC(i+1)*(j+2)/arad(i)**(j+3)
  a(ieq,ivr+12) = rho(i+1)/arad(i)**(j+1)

  ! Continuity of tangential traction
  ieq  = 6*(i-1)+4
  a(ieq,ivr+1)  = -muC(i)*j*(j+2)*arad(i)**j/dfloat((2*j+3)*(j+1))
  a(ieq,ivr+2)  = -2._DP*muC(i)*(j-1)*arad(i)**(j-2)/dfloat(j)
  a(ieq,ivr+4)  = -muC(i)*(j**2-1)/arad(i)**(j+1)/dfloat(j*(2*j-1))
  a(ieq,ivr+5)  = -2._DP*muC(i)*(j+2)/arad(i)**(j+3)/dfloat(j+1)
  a(ieq,ivr+7)  = muC(i+1)*j*(j+2)*arad(i)**j/dfloat((2*j+3)*(j+1))
  a(ieq,ivr+8)  = 2._DP*muC(i+1)*(j-1)*arad(i)**(j-2)/dfloat(j)
  a(ieq,ivr+10) = muC(i+1)*(j**2-1)/arad(i)**(j+1)/dfloat(j*(2*j-1))
  a(ieq,ivr+11) = 2._DP*muC(i+1)*(j+2)/arad(i)**(j+3)/dfloat(j+1)

  ! Continuity of gravity potential
  ieq  = 6*(i-1)+5
  a(ieq,ivr+3)  = -arad(i)**j
  a(ieq,ivr+6)  = -arad(i)**(-j-1)
  a(ieq,ivr+9)  = arad(i)**j
  a(ieq,ivr+12) = arad(i)**(-j-1)

  ! Discontinuity of its first derivative (due to surface density)
  ieq  = 6*(i-1)+6
  a(ieq,ivr+1)  = -4._DP*pi*kappa*rho(i)*j*arad(i)**(j+1)/(2._DP*(2*j+3))
  a(ieq,ivr+2)  = -4._DP*pi*kappa*rho(i)*arad(i)**(j-1)
  a(ieq,ivr+3)  = -arad(i)**(j-1)*(2*j+1)
  a(ieq,ivr+4)  = -4._DP*pi*kappa*rho(i)*(j+1)/arad(i)**j/(2._DP*(2*j-1))
  a(ieq,ivr+5)  = -4._DP*pi*kappa*rho(i)/arad(i)**(j+2)
  a(ieq,ivr+7)  = 4._DP*pi*kappa*rho(i+1)*j*arad(i)**(j+1)/(2._DP*(2*j+3))
  a(ieq,ivr+8)  = 4._DP*pi*kappa*rho(i+1)*arad(i)**(j-1)
  a(ieq,ivr+9)  = arad(i)**(j-1)*(2*j+1)
  a(ieq,ivr+10) = 4._DP*pi*kappa*rho(i+1)*(j+1)/arad(i)**j/(2._DP*(2*j-1))
  a(ieq,ivr+11) = 4._DP*pi*kappa*rho(i+1)/arad(i)**(j+2)

end do

! i=nl
!------
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

!--------------------------------------------------------------------------

! Fill RHS
!----------
b = 0._DP
b(neq) = -(2*j+1)/arad(nl)

! LU decomposition
!------------------
call ludcmp(a,neq,np,indx,d)
call lubksb(a,neq,np,indx,b)

! Love number k_j
!-----------------
k2 = -(b(neq-3)*arad(nl)**j+b(neq)*arad(nl)**(-j-1))-dcmplx(1._DP)

CONTAINS
!======================================
! Compute compliance
!======================================
SUBROUTINE compliance(Jc,freq,ilayer)
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
  Jc = dcmplx(1.d10,0._DP)
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
    imj=-1._DP/(freq*taum)-(freq*zeta(ilayer)*taum)**(-alpha(ilayer))*sina*exp(gammln(alpha(ilayer)+1._DP)) - cdef(ilayer)**2/mu(ilayer)*eta_kv_rel(ilayer)*eta(ilayer)*freq/denom
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

END SUBROUTINE peltier
