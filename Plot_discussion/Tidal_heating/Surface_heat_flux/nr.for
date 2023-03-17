C --- Gauss-Legendre abscissae and weights
      SUBROUTINE gauleg(x1,x2,x,w,n)
      INTEGER n
      DOUBLE PRECISION x1,x2,x(n),w(n)
      DOUBLE PRECISION EPS
      PARAMETER (EPS=1.d-14)
      INTEGER i,j,m
      DOUBLE PRECISION p1,p2,p3,pp,xl,xm,z,z1
      m=(n+1)/2
      xm=0.5d0*(x2+x1)
      xl=0.5d0*(x2-x1)
      do 12 i=1,m
        z=cos(3.141592654d0*(i-.25d0)/(n+.5d0))
1       continue
          p1=1.d0
          p2=0.d0
          do 11 j=1,n
            p3=p2
            p2=p1
            p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
11        continue
          pp=n*(z*p1-p2)/(z*z-1.d0)
          z1=z
          z=z1-p1/pp
        if(abs(z-z1).gt.EPS)goto 1
        x(i)=xm-xl*z
        x(n+1-i)=xm+xl*z
        w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
        w(n+1-i)=w(i)
12    continue
      return
      END


      SUBROUTINE ludcmp(a,n,np,indx,d)
      INTEGER n,np,indx(n),NMAX
      COMPLEX*16 a(np,np),sum,dumc
      DOUBLE PRECISION d,TINY
      PARAMETER (NMAX=500,TINY=1.0d-20)
      INTEGER i,imax,j,k
      DOUBLE PRECISION aamax,dum,vv(NMAX)

      !print *,a(1,:)

      d=1.d0
      do 12 i=1,n
        aamax=0.d0
        do 11 j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11      continue
        !print *,i,a(i,:),aamax
        if (aamax.eq.0.d0) pause 'singular matrix in ludcmp'
        vv(i)=1.d0/aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.d0
        do 16 i=j,n
          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dumc=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dumc
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(abs(a(j,j)).eq.0.d0)a(j,j)=TINY
        if(j.ne.n)then
          dumc=1.d0/a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dumc
18        continue
        endif
19    continue
      return
      END


      SUBROUTINE lubksb(a,n,np,indx,b)
      INTEGER n,np,indx(n)
      COMPLEX*16 a(np,np),b(n)
      INTEGER i,ii,j,ll
      COMPLEX*16 sum
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.d0) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14    continue
      return
      END
