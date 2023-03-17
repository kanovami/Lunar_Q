MODULE mconst
IMPLICIT NONE

INTEGER, PARAMETER  :: DP  = KIND(1.0D0)
INTEGER, PARAMETER  :: DPC = KIND((1.0D0,1.0D0))
INTEGER, PARAMETER  :: ifile=11,ifileP=12
INTEGER, PARAMETER  :: ngl=100  ! # of abscissae for Gauss-Legendre quadrature
! nlmax..maximum # of layers (for allocation)
! nhead..# of rows in the head of rheo.in (please do not change)
INTEGER, PARAMETER  :: nlmax=110,nhead=20
INTEGER, PARAMETER  :: nhans=10 ! # of Hansen coefficients (qmax=(nhans-1)/2)
! tidal potential expanded to the qmax-th order in orbital eccentricity
REAL(DP), PARAMETER :: kappa=6.6732d-11,pi=4._DP*datan(1._DP)
REAL(DP), PARAMETER :: msun=1.9891d30,rsun=6.957d8,mearth=5.97237d24
REAL(DP), PARAMETER :: day=86400._DP,rok=365.25_DP*day,AUm=1.496d11
COMPLEX(DPC),PARAMETER     :: imag=(0._DP,1._DP)
CHARACTER(LEN=*),PARAMETER :: filein='rheo_sc.in'

! Tidal heating calculation
INTEGER, PARAMETER  :: nph=72,nth=36
INTEGER, PARAMETER  :: ntime=100
INTEGER, PARAMETER  :: nr=100 ! # of layers for radially dependent heating computation

! Star and orbit
REAL(DP),PARAMETER  :: mstar=1._DP*mearth
REAL(DP),PARAMETER  :: astar=384399d3
REAL(DP),PARAMETER  :: estar=0.0549_DP
REAL(DP),PARAMETER  :: res=1._DP

! Degree of SH decomposition
!INTEGER, PARAMETER  :: jdeg=2

! Numerical Recipes
INTEGER                  :: neq,np
REAL(DP)                 :: d
INTEGER,ALLOCATABLE      :: indx(:)
COMPLEX(DPC),ALLOCATABLE :: a(:,:),b(:)

! From rheo.in
INTEGER                  :: iost,nl
INTEGER,ALLOCATABLE      :: ire(:)
REAL(DP),ALLOCATABLE     :: rho(:),arad(:),eta(:),mu(:),alpha(:),zeta(:)
REAL(DP),ALLOCATABLE     :: cdef(:),eta_kv_rel(:)
CHARACTER                :: adum

END MODULE mconst
