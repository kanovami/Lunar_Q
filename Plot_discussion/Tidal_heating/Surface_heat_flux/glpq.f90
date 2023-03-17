SUBROUTINE glpq(l,p,q,ecc,hans)
USE mconst
IMPLICIT NONE
!---------------------------
! Variables
!---------------------------
INTEGER   :: tau,rh,sigma,k,n,j,l,p,q,igl
INTEGER   :: crst,rsmax
REAL(DP)  :: ct,hans,dhans,ecc
REAL(DP)  :: ver,rstar
REAL(DP),ALLOCATABLE :: x(:),w(:)

n=-(l+1)
k=l-2*p
j=l-2*p+q

! Computation of Hansen coefficients
!------------------------------------
hans=0._DP
DO sigma=0,nhans
  rh=j-k+sigma
  IF (rh.lt.0) CYCLE
  dhans=newcomb(n,k,rh,sigma)*ecc**(rh+sigma)
  hans=hans+dhans
  IF (abs(hans).lt.epsilon(hans)) THEN
    if ((sigma.gt.0).and.(abs(dhans).lt.epsilon(dhans))) EXIT
  END IF
  IF (abs(dhans/hans).lt.1d-10) EXIT
END DO

CONTAINS

!---------------------------
! Combinatorial numbers
!---------------------------
FUNCTION comb(as,bs)
USE mconst
IMPLICIT NONE
REAL(DP) :: as
REAL(DP) :: comb
INTEGER  :: bs
IF (bs.eq.0) THEN
  comb=1._DP
ELSE
  comb=pochhammer(as,bs)/factorial(bs)
END IF
END FUNCTION comb
!---------------------------

!---------------------------
! Compute factorial
!---------------------------
RECURSIVE FUNCTION factorial(n)  RESULT(fact)
USE mconst
IMPLICIT NONE
INTEGER :: fact
INTEGER, INTENT(IN) :: n
IF (n.eq.0) THEN
  fact=1
ELSE
  fact=n*factorial(n-1)
END IF
END FUNCTION factorial
!---------------------------

!---------------------------
! Falling factorial, real r
!---------------------------
RECURSIVE FUNCTION pochhammer(r,k)  RESULT(poch)
USE mconst
IMPLICIT NONE
INTEGER  :: k
REAL(DP) :: poch,r
IF (k.eq.0) THEN
  poch=1._DP
ELSE
  poch=(r-1._DP*(k-1))*pochhammer(r,k-1)
END IF
END FUNCTION pochhammer
!---------------------------

!------------------------------
! Falling factorial, integer r
!------------------------------
RECURSIVE FUNCTION ipochhammer(n,k)  RESULT(ipoch)
USE mconst
IMPLICIT NONE
INTEGER(8) :: ipoch
INTEGER    :: k,n
IF (k.eq.1) THEN
  ipoch=1
ELSE
  ipoch=(n-(k-1))*ipochhammer(n,k-1)
END IF
END FUNCTION ipochhammer
!---------------------------

!---------------------------
! J polynomials for s=0
!---------------------------
RECURSIVE FUNCTION Jnull(n,k,r)  RESULT(Jnul)
USE mconst
IMPLICIT NONE
INTEGER(8)          :: Jnul
INTEGER             :: n,k,r
IF (r.eq.0) THEN
  Jnul=1
ELSEIF (r.eq.1) THEN
  Jnul=2*k-n
ELSEIF (r.lt.0) THEN
  Jnul=0
ELSE
  Jnul=(2*k-n)*Jnull(n,k+1,r-1) + (r-1)*(k-n)*Jnull(n,k+2,r-2)
END IF
END FUNCTION Jnull
!---------------------------

!---------------------------
! J polynomials for s>0
!---------------------------
RECURSIVE FUNCTION Jpositive(n,k,r,s)  RESULT(Jpos)
USE mconst
IMPLICIT NONE
INTEGER(8)          :: Jpos
INTEGER, INTENT(IN) :: n,k,r,s
REAL(DP)   :: ct
INTEGER    :: tau
INTEGER(8) :: crst,sumt
IF (s.eq.0) THEN
  Jpos=Jnull(n,k,r)
ELSEIF (s.lt.0) THEN
  Jpos=0
ELSE
  sumt=0
  DO tau=2,min(r,s)
    ct   = (-1._DP)**tau*comb(1.5_DP,tau)*2._DP**(2.*tau-1.)
    crst = ipochhammer(r,tau)*ipochhammer(s,tau)*int(ct)
    sumt=sumt+crst*Jpositive(n,k,r-tau,s-tau)
    IF (Jpositive(n,k,r-tau,s-tau).eq.0) EXIT
  END DO
  Jpos=-(2*k+n)*Jpositive(n,k-1,r,s-1) - (s-1)*(k+n)*Jpositive(n,k-2,r,s-2) &
       -r*(r-5*s+4+4*k+n)*Jpositive(n,k,r-1,s-1)+r*(r-s+k)*sumt
END IF
END FUNCTION Jpositive
!---------------------------

!-------------------------------------
! Computation of Newcomb coefficients
!-------------------------------------
FUNCTION newcomb(n,k,r,s)
USE mconst
IMPLICIT NONE
INTEGER   :: n,k,r,s
REAL(DP)  :: newcomb
IF ((r.lt.0).or.(s.lt.0)) THEN
  newcomb=0._DP
  RETURN
END IF

IF (r.lt.s) THEN
  newcomb=dfloat(Jpositive(n,-k,s,r))/2._DP**(r+s)/dfloat(factorial(r))/dfloat(factorial(s))
ELSE
  newcomb=dfloat(Jpositive(n,k,r,s))/2._DP**(r+s)/dfloat(factorial(r))/dfloat(factorial(s))
END IF
END FUNCTION newcomb
!-------------------------------------

END SUBROUTINE glpq
