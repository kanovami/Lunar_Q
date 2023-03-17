!==========================================================================
! Code calculating the surface tidal heat flux pattern (M. Walterova, 2019)
!==========================================================================
! - orbital parameters and resolution can be changed in 'mconst.f90'
! - planet parameters (structure, densities, rheology) can be changed
!   in the input file (parameter 'filein' in 'mconst.f90')
! - please folow the structure of the existing input files
!==========================================================================
PROGRAM main
USE mconst
IMPLICIT NONE
!-------------------------------------
INTEGER  :: i,il
INTEGER  :: iph,ith,indtp
INTEGER  :: ir,it,ip
INTEGER, ALLOCATABLE  :: q(:),p(:)
REAL(DP) :: mpl,rpl,norb
REAL(DP) :: freq,phi,thet,rad
REAL(DP) :: f220,f201,b22,b20,dheat
REAL(DP) :: fact
REAL(DP) :: hflux,heat_aux
REAL(DP),ALLOCATABLE :: g20(:),g21(:),pot220sin(:),pot220cos(:),pot201(:),heatr_total(:)
REAL(DP),ALLOCATABLE :: x(:),w(:)
COMPLEX(DPC)         :: klove
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
! INITIALIZATION AND ALLOCATIONS
!--------------------------------------------------------------------------
ALLOCATE(rho(nlmax),arad(0:nlmax),eta(nlmax),mu(nlmax),ire(nlmax))
ALLOCATE(q(0:nhans),g20(0:nhans),g21(0:nhans),p(0:nhans))
ALLOCATE(pot220sin(nth*nph),pot220cos(nth*nph),pot201(nth*nph))
ALLOCATE(x(ngl_rad),w(ngl_rad))

arad(0) = 0._DP

! Read planet's parameters
!--------------------------
nl=0
open(ifileP,file=filein)
do i=1,nlmax+nhead
  read(ifileP,'(a)',IOSTAT=iost,ADVANCE='no') adum
  if (iost.lt.0) exit
  if (adum.eq.'#') then
    read(ifileP,*)
  else
    nl=nl+1
    read(ifileP,'(i3)',ADVANCE='no') ire(nl)
    if (ire(nl).eq.0) then
      read(ifileP,*) rho(nl),arad(nl),eta(nl),mu(nl)
    else if (ire(nl).eq.1) then
      if (.not.ALLOCATED(alpha)) then
        ALLOCATE(alpha(nlmax))
        alpha=0._DP
      end if
      if (.not.ALLOCATED(zeta)) then
        ALLOCATE(zeta(nlmax))
        zeta =0._DP
      end if
      read(ifileP,*) rho(nl),arad(nl),eta(nl),mu(nl),alpha(nl),zeta(nl)
    else if (ire(nl).eq.2) then ! Sundberg-Cooper
      if (.not.ALLOCATED(alpha)) then
        ALLOCATE(alpha(nlmax))
        alpha=0._DP
      end if
      if (.not.ALLOCATED(zeta)) then
        ALLOCATE(zeta(nlmax))
        zeta =0._DP
      end if
      if (.not.ALLOCATED(eta_kv_rel)) then
        ALLOCATE(eta_kv_rel(nlmax))
        eta_kv_rel=0._DP
      end if
      if (.not.ALLOCATED(cdef)) then
        ALLOCATE(cdef(nlmax))
        cdef =0._DP
      end if
      read(ifileP,*) rho(nl),arad(nl),eta(nl),mu(nl),alpha(nl),zeta(nl),eta_kv_rel(nl),cdef(nl)
      !print *,rho(nl),arad(nl),eta(nl),mu(nl),alpha(nl),zeta(nl),eta_kv_rel(nl),cdef(nl)
    else
      print *,'Unknown rheology in rheo.in'
    end if
  end if
end do
close(ifileP)

! Allocate arrays for k2 calculation
neq = 6*nl-3
np  = neq
ALLOCATE(a(np,np),b(neq),indx(neq))

! Planet radius and mass
rpl = arad(nl)
mpl = 0._DP
do i=1,nl-1
  mpl = mpl + arad(i)**3*(rho(i)-rho(i+1))
end do
mpl = mpl + arad(nl)**3*rho(nl)
mpl = 4._DP/3.*pi*mpl

! Mean motion
norb=sqrt(kappa*(mstar+mpl)/astar**3)

! Planetographic factors for tidal potential
fact=0.5_DP*kappa*mstar*rpl**2
DO iph=1,nph
  phi = 2._DP*pi*(iph-1)/dfloat(nph)
  DO ith=1,nth
    thet = pi*(1._DP*ith-0.5_DP)/dfloat(nth)
    indtp=(iph-1)*nth+ith
    pot220sin(indtp)=0.5_DP*fact*sin(thet)**2*sin(2._DP*phi)
    pot220cos(indtp)=0.5_DP*fact*sin(thet)**2*cos(2._DP*phi)
    pot201(indtp)=fact*(3._DP*cos(thet)**2-1._DP)
  END DO
END DO

! Gauss-Legendre weights for integration over radius
call gauleg(0._DP,rpl,x,w,ngl_rad)

!===================================================
! Calculate and output the heat fluxes (in mW/m**2)
!===================================================
OPEN(ifile,file='Moon_tidflux_MELT_aux.dat')

heat_aux=0._DP

do it=1,nth
  thet = pi*(1._DP*it-0.5_DP)/dfloat(nth)
  do ip=1,nph
    phi=(ip-1)*2._DP*pi/dfloat(nph)
    hflux=0._DP
    do ir=1,ngl_rad
      rad=x(ir)
      do il=1,nl
        if ((rad.ge.arad(il-1)).and.(rad.le.arad(il))) exit
      end do
      call unit_heating(dheat,rad,thet,phi,il)
      hflux=hflux+w(ir)*dheat*(rad/rpl)**2
    end do
    write(ifile, '(3g20.10)') phi,thet,hflux*1d3
    heat_aux=heat_aux+hflux
  end do
end do
heat_aux=heat_aux/nth/nph*4*pi*rpl**2
print *,heat_aux

CLOSE(ifile)
!===================================================

CONTAINS

!==============================================
! Kaula's F_lmp and G_lpq coefficients for i=0
!==============================================
SUBROUTINE kaula_coeffs(estar)
REAL(DP)  :: estar

! DEGREE 2, ORDER 0,2
b20 = kappa*mstar
b22 = 1._DP/12.*kappa*mstar

! Inclination functions F_201 and F_220
!---------------------------------------
f201 =-0.5_DP
f220 = 3._DP

! Eccentricity functions G_21q and G_20q
DO i=0,nhans
  q(i)=i-nhans/2
  call glpq(2,1,q(i),estar,g21(i))
  call glpq(2,0,q(i),estar,g20(i))
END DO

END SUBROUTINE kaula_coeffs
!==============================================

!==============================================
! Tidal heating in a unit volume (following Tobie, 2003)
!==============================================
SUBROUTINE unit_heating(dheat,rad,thet,phi,il)
USE mconst
IMPLICIT NONE
INTEGER   :: il,ieq,itime,j,ig
REAL(DP)  :: rad,thet,phi,arg,arg1,arg2,arg3,delta,delta_m
REAL(DP)  :: freq(2*nhans+2),fact,dheat,grav,ksi
REAL(DP)  :: lag(2*nhans+2) 
COMPLEX(DPC)  :: c1,c2,c3,c4,c5,c6
COMPLEX(DPC)  :: sum1,sum2,sum3,sum4,sum5,sum6
COMPLEX(DPC)  :: un,vn,trn,ttn,dundr
COMPLEX(DPC)  :: u,du_theta,ddu_theta,du_phi,ddu_phi,ddu_phth
COMPLEX(DPC),DIMENSION(2*nhans+2) :: stress1,stress2,stress3,stress5,stress6,stress4
COMPLEX(DPC),DIMENSION(2*nhans+2) :: strain1,strain2,strain3,strain5,strain6,strain4
COMPLEX(DPC),ALLOCATABLE :: muC(:)

j=2

norb=sqrt(kappa*(mstar+mpl)/astar**3)
call kaula_coeffs(estar)

if (il.eq.1) then
  grav = 0._DP
else if (il.eq.2) then
  grav = kappa*4._DP*pi*arad(1)*rho(1)/3.
else
  grav = kappa*4._DP*pi*arad(1)*rho(1)/3.
  do ig=2,il-1
    grav = (arad(ig-1)/arad(ig))**2*grav &
   &     + kappa*4._DP*pi*rho(ig)*(arad(ig)**3-arad(ig-1)**3)/3./arad(ig)**2
  end do
end if

if (rad.lt.epsilon(rad)) then
  grav = 0._DP
else
  grav = (arad(il-1)/rad)**2*grav &
  &    + kappa*4._DP*pi*rho(il)*(rad**3-arad(il-1)**3)/3./rad**2
end if

IF (.not.ALLOCATED(muC)) ALLOCATE(muC(nl))

fact=0.5_DP*kappa*mstar*rpl**2/astar**3

DO i=0,nhans ! iterate over frequencies
  q(i)=i-nhans/2

  ! 220q
  freq(i+1)   = (1._DP*(2+q(i))-2._DP*res)*norb
  call peltier(2,freq(i+1),klove,muC)
  lag(i+1) = atan(aimag(muc(il))/real(muc(il)))

  ! Tidal potential and its derivatives
  u         = 0.5_DP*fact*f220*g20(i)*dsin(thet)**2
  u         = 0.5_DP*u
  du_theta  = fact*f220*g20(i)*dcos(thet)*dsin(thet)
  du_theta  = 0.5_DP*du_theta
  du_phi    = fact*f220*g20(i)*dsin(thet)**2
  du_phi    = 0.5_DP*du_phi
  ddu_theta = fact*f220*g20(i)*(dcos(thet)**2-dsin(thet)**2)
  ddu_theta = 0.5_DP*ddu_theta
  ddu_phi   = -2._DP*fact*f220*g20(i)*dsin(thet)**2
  ddu_phi   = 0.5_DP*ddu_phi
  ddu_phth  = 2._DP*fact*f220*g20(i)*dsin(thet)*dcos(thet)
  ddu_phth  = 0.5_DP*ddu_phth

  if (il.gt.1) then
    ieq = 6*(il-1)-3
    c1 = b(ieq+1)
    c2 = b(ieq+2)
    c3 = b(ieq+3)
    c4 = b(ieq+4)
    c5 = b(ieq+5)
    c6 = b(ieq+6)
  else
    c1 = b(1)
    c2 = b(2)
    c3 = b(3)
    c4 = dcmplx(0.)
    c5 = dcmplx(0.)
    c6 = dcmplx(0.)
  end if

  ksi = grav/rad

  ! Displacements and tractions
  un = c1*j*rad**(j+1)/(2._DP*(2*j+3)) + c2*rad**(j-1) + c4*(j+1)/rad**j/(2._DP*(2*j-1)) + c5/rad**(j+2)
  vn = c1*(j+3)*rad**(j+1)/(2._DP*(2*j+3)*(j+1)) + c2*rad**(j-1)/j + c4*dfloat(2-j)/rad**j/(2._DP*j*(2*j-1)) - c5/rad**(j+2)/dfloat(j+1)
  dundr = c1*j*(j+1)*rad**j/2/(2*j+3) + c2*(j-1)*rad**(j-2) - c4*j*(j+1)/2/(2*j-1)/rad**(j+1) - c5*(j+2)/rad**(j+3)
  trn = c1*(rho(il)*ksi*j*rad**(j+2)+muC(il)*2*(j**2-j-3)*rad**j)/2/(2*j+3) &
      + c2*(rho(il)*ksi*rad**j+2*muC(il)*(j-1)*rad**(j-2)) + c3*rho(il)*rad**j &
      + c4*(rho(il)*ksi*(j+1)*rad**(-j+1)-muC(il)*2*(j**2+3*j-1)*rad**(-j-1))/2/(2*j-1) &
      + c5*(rho(il)*ksi/rad**(j+1)-2*muC(il)*(j+2)/rad**(j+3)) + c6*rho(il)/rad**(j+1)
  ttn = c1*muC(il)*j*(j+2)*rad**j/(2*j+3)/(j+1) + c2*2*muC(il)*(j-1)*rad**(j-2)/j &
      + c4*muC(il)*(j-1)*(j+1)/j/(2*j-1)/rad**(j+1) + c5*2*muC(il)*(j+2)/(j+1)/rad**(j+3)

  strain1(i+1) = dundr*u
  strain2(i+1) = 1._DP/rad*(vn*ddu_theta+un*u)
  strain3(i+1) = 1._DP/rad*(vn*ddu_phi/dsin(thet)**2+un*u+vn*dcos(thet)/dsin(thet)*du_theta)
  strain4(i+1) = 1._DP/rad*vn*(ddu_phth/dsin(thet)-dcos(thet)/dsin(thet)**2*du_phi)
  strain5(i+1) = 1._DP*ttn*du_phi/dsin(thet)/muC(il)/2.
  strain6(i+1) = 1._DP*ttn*du_theta/muC(il)/2.

  stress1(i+1) = 2*abs(muC(il))*strain1(i+1) !trn*u 
  stress2(i+1) = 2*abs(muC(il))*strain2(i+1) !(trn-2*muC(il)*dundr)*u+2*muC(il)*strain2(i)
  stress3(i+1) = 2*abs(muC(il))*strain3(i+1) !(trn-2*muC(il)*dundr)*u+2*muC(il)*strain3(i)
  stress4(i+1) = 2*abs(muC(il))*strain4(i+1)
  stress5(i+1) = 2*abs(muC(il))*strain5(i+1)
  stress6(i+1) = 2*abs(muC(il))*strain6(i+1)

  ! 201q
  freq(i+nhans+2)   = q(i)*norb
  call peltier(2,freq(i+nhans+2),klove,muC)
  lag(i+nhans+2) = atan(aimag(muc(il))/real(muc(il)))

  ! Tidal potential and its derivatives
  u         = fact*f201*g21(i)*(3*dcos(thet)**2-1._DP)
  u         = 0.5_DP*u
  du_theta  = -6*fact*f201*g21(i)*dcos(thet)*dsin(thet)
  du_theta  = 0.5_DP*du_theta
  du_phi    = 0._DP
  ddu_theta = -6*fact*f201*g21(i)*(dcos(thet)**2-dsin(thet)**2)
  ddu_theta = 0.5_DP*ddu_theta
  ddu_phi   = 0._DP
  ddu_phth  = 0._DP

  if (il.gt.1) then
    ieq = 6*(il-1)-3
    c1 = b(ieq+1)
    c2 = b(ieq+2)
    c3 = b(ieq+3)
    c4 = b(ieq+4)
    c5 = b(ieq+5)
    c6 = b(ieq+6)
  else
    c1 = b(1)
    c2 = b(2)
    c3 = b(3)
    c4 = dcmplx(0.)
    c5 = dcmplx(0.)
    c6 = dcmplx(0.)
  end if

  ksi = grav/rad

  ! Displacements and tractions
  un = c1*j*rad**(j+1)/(2._DP*(2*j+3)) + c2*rad**(j-1) + c4*(j+1)/rad**j/(2._DP*(2*j-1)) + c5/rad**(j+2)
  vn = c1*(j+3)*rad**(j+1)/(2._DP*(2*j+3)*(j+1)) + c2*rad**(j-1)/j + c4*dfloat(2-j)/rad**j/(2._DP*j*(2*j-1)) - c5/rad**(j+2)/dfloat(j+1)
  dundr = c1*j*(j+1)*rad**j/2/(2*j+3) + c2*(j-1)*rad**(j-2) - c4*j*(j+1)/2/(2*j-1)/rad**(j+1) - c5*(j+2)/rad**(j+3)
  trn = c1*(rho(il)*ksi*j*rad**(j+2)+muC(il)*2*(j**2-j-3)*rad**j)/2/(2*j+3) &
      + c2*(rho(il)*ksi*rad**j+2*muC(il)*(j-1)*rad**(j-2)) + c3*rho(il)*rad**j &
      + c4*(rho(il)*ksi*(j+1)*rad**(-j+1)-muC(il)*2*(j**2+3*j-1)*rad**(-j-1))/2/(2*j-1) &
      + c5*(rho(il)*ksi/rad**(j+1)-2*muC(il)*(j+2)/rad**(j+3)) + c6*rho(il)/rad**(j+1)
  ttn = c1*muC(il)*j*(j+2)*rad**j/(2*j+3)/(j+1) + c2*2*muC(il)*(j-1)*rad**(j-2)/j &
      + c4*muC(il)*(j-1)*(j+1)/j/(2*j-1)/rad**(j+1) + c5*2*muC(il)*(j+2)/(j+1)/rad**(j+3)

  strain1(i+nhans+2) = dundr*u
  strain2(i+nhans+2) = 1._DP/rad*(vn*ddu_theta+un*u)
  strain3(i+nhans+2) = 1._DP/rad*(vn*dcos(thet)/dsin(thet)*du_theta+un*u)
  strain4(i+nhans+2) = 0._DP
  strain5(i+nhans+2) = 0._DP
  strain6(i+nhans+2) = 1._DP*ttn*du_theta/muC(il)/2.

  stress1(i+nhans+2) = 2*abs(muC(il))*strain1(i+nhans+2) !trn*u 
  stress2(i+nhans+2) = 2*abs(muC(il))*strain2(i+nhans+2) !(trn-2*muC(il)*dundr)*u+2*muC(il)*strain2(i+nhans+1)
  stress3(i+nhans+2) = 2*abs(muC(il))*strain3(i+nhans+2) !(trn-2*muC(il)*dundr)*u+2*muC(il)*strain3(i+nhans+1)
  stress4(i+nhans+2) = 0._DP
  stress5(i+nhans+2) = 0._DP
  stress6(i+nhans+2) = 2*abs(muC(il))*strain6(i+nhans+2)

END DO

sum1 = dcmplx(0)
sum2 = dcmplx(0)
sum3 = dcmplx(0)
sum4 = dcmplx(0)
sum5 = dcmplx(0)
sum6 = dcmplx(0)

! Unit heating integrated over mean anomaly (-pi,pi)
DO i=0,nhans
  q(i)=i-nhans/2
  DO j=0,nhans
    p(j)=j-nhans/2

    if (i.eq.j) then
      delta=1._DP
    else
      delta=0._DP
    end if

    if (q(i).eq.(-p(j))) then
      delta_m=1._DP
    else
      delta_m=0._DP
    end if

    ! 220q
    arg = 4*(1._DP-res)+q(i)+p(j)

    IF ((abs(freq(i+1)).gt.epsilon(freq(i+1))).and.(abs(freq(j+1)).gt.epsilon(freq(j+1)))) THEN
      sum1 = sum1 + stress1(i+1)*conjg(strain1(j+1))*freq(j+1)*dsin(4*phi+lag(j+1))*delta_sin(arg) + &
        &  stress1(i+1)*conjg(strain1(j+1))*freq(j+1)*delta*dsin(lag(j+1))
      sum2 = sum2 + stress2(i+1)*conjg(strain2(j+1))*freq(j+1)*dsin(4*phi+lag(j+1))*delta_sin(arg) + &
        &  stress2(i+1)*conjg(strain2(j+1))*freq(j+1)*delta*dsin(lag(j+1))
      sum3 = sum3 + stress3(i+1)*conjg(strain3(j+1))*freq(j+1)*dsin(4*phi+lag(j+1))*delta_sin(arg) + &
        &  stress3(i+1)*conjg(strain3(j+1))*freq(j+1)*delta*dsin(lag(j+1))
      sum4 = sum4 - stress4(i+1)*conjg(strain4(j+1))*freq(j+1)*dsin(4*phi+lag(j+1))*delta_sin(arg) + &
        &  stress4(i+1)*conjg(strain4(j+1))*freq(j+1)*delta*dsin(lag(j+1))
      sum5 = sum5 - stress5(i+1)*conjg(strain5(j+1))*freq(j+1)*dsin(4*phi+lag(j+1))*delta_sin(arg) + &
        &  stress5(i+1)*conjg(strain5(j+1))*freq(j+1)*delta*dsin(lag(j+1))
      sum6 = sum6 + stress6(i+1)*conjg(strain6(j+1))*freq(j+1)*dsin(4*phi+lag(j+1))*delta_sin(arg) + &
        &  stress6(i+1)*conjg(strain6(j+1))*freq(j+1)*delta*dsin(lag(j+1))
    END IF

    ! 201q
    IF ((abs(freq(i+nhans+2)).gt.epsilon(freq(i+nhans+2))).and.(abs(freq(j+nhans+2)).gt.epsilon(freq(j+nhans+2)))) THEN
      sum1 = sum1 + stress1(i+nhans+2)*conjg(strain1(j+nhans+2))*freq(j+nhans+2)*(delta+delta_m)*dsin(lag(j+nhans+2))
      sum2 = sum2 + stress2(i+nhans+2)*conjg(strain2(j+nhans+2))*freq(j+nhans+2)*(delta+delta_m)*dsin(lag(j+nhans+2))
      sum3 = sum3 + stress3(i+nhans+2)*conjg(strain3(j+nhans+2))*freq(j+nhans+2)*(delta+delta_m)*dsin(lag(j+nhans+2))
      sum6 = sum6 + stress6(i+nhans+2)*conjg(strain6(j+nhans+2))*freq(j+nhans+2)*(delta+delta_m)*dsin(lag(j+nhans+2))
    END IF

    ! 220q x 201q
    arg1 = q(i)+p(j)+2*(1._DP-res)
    arg2 = -q(i)+p(j)-2*(1._DP-res)
    arg3 = -q(i)+p(j)+2*(1._DP-res)

    IF ((abs(freq(i+1)).gt.epsilon(freq(i+1))).and.(abs(freq(j+nhans+2)).gt.epsilon(freq(j+nhans+2)))) THEN
      sum1 = sum1 + stress1(i+1)*conjg(strain1(j+nhans+2))*freq(j+nhans+2)*(dsin(2*phi+lag(j+nhans+2))*delta_sin(arg1)-dsin(2*phi-lag(j+nhans+2))*delta_sin(arg2))
      sum2 = sum2 + stress2(i+1)*conjg(strain2(j+nhans+2))*freq(j+nhans+2)*(dsin(2*phi+lag(j+nhans+2))*delta_sin(arg1)-dsin(2*phi-lag(j+nhans+2))*delta_sin(arg2))
      sum3 = sum3 + stress3(i+1)*conjg(strain3(j+nhans+2))*freq(j+nhans+2)*(dsin(2*phi+lag(j+nhans+2))*delta_sin(arg1)-dsin(2*phi-lag(j+nhans+2))*delta_sin(arg2))
      sum6 = sum6 + stress6(i+1)*conjg(strain6(j+nhans+2))*freq(j+nhans+2)*(dsin(2*phi+lag(j+nhans+2))*delta_sin(arg1)-dsin(2*phi-lag(j+nhans+2))*delta_sin(arg2))
    END IF

    IF ((abs(freq(j+1)).gt.epsilon(freq(j+1))).and.(abs(freq(i+nhans+2)).gt.epsilon(freq(i+nhans+2)))) THEN
      sum1 = sum1 + dsin(2*phi+lag(j+1))*stress1(i+nhans+2)*conjg(strain1(j+1))*freq(j+1)*(delta_sin(arg1)+delta_sin(arg3))
      sum2 = sum2 + dsin(2*phi+lag(j+1))*stress2(i+nhans+2)*conjg(strain2(j+1))*freq(j+1)*(delta_sin(arg1)+delta_sin(arg3))
      sum3 = sum3 + dsin(2*phi+lag(j+1))*stress3(i+nhans+2)*conjg(strain3(j+1))*freq(j+1)*(delta_sin(arg1)+delta_sin(arg3))
      sum6 = sum6 + dsin(2*phi+lag(j+1))*stress6(i+nhans+2)*conjg(strain6(j+1))*freq(j+1)*(delta_sin(arg1)+delta_sin(arg3))
    END IF

  END DO
END DO

! Unit heating
dheat = 2*real(sum1 + sum2 + sum3 + 2*sum4 + 2*sum5 + 2*sum6)


END SUBROUTINE unit_heating
!==============================================

FUNCTION delta_sin(arg)
USE mconst
IMPLICIT NONE
REAL(DP)  :: arg
COMPLEX(DPC) :: delta_sin

IF (abs(arg).lt.EPSILON(arg)) THEN
  delta_sin = 1._DP
ELSE
  delta_sin = exp(imag*pi*arg)*dsin(arg*pi)/arg/pi
  if (abs(delta_sin).lt.1.d-10) delta_sin=0._DP
END IF

END FUNCTION delta_sin

END PROGRAM
