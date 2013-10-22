      subroutine bcgstab(n,itmax,crit,x,ubar,b,p,r,s,t,v,suba)
      implicit real*8(a-h,o-z)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Code supplied by Laurette Tuckerman, July 2011.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  solves Ax=b by the biconjugate gradient stabilized algorithm (van der Vorst)
c  requires a subroutine suba(x,y) which has x as input and y=Ax as output
c  n is the size of the system
c  itmax is the maximum number of iterations allowed
c  crit is the convergence criterion
c  x is the last solution vector produced
c  ubar is the best solution vector (smallest residual) produced
c  b is the right-hand-side vector
c  p,r,s,t,v are work arrays
c
      dimension b(n),x(n),p(n),r(n),s(n),t(n),v(n),ubar(n)
      EXTERNAL suba
      write(6,*)
      write(6,*) '  bcgstab'
      write(6,*)      
      write(6,*) ' intermediate output after each iteration'
      write(6,*) ' iteration    test'
      rho_old=1.
      alpha=1.
      omega=1.
      do j=1,n
         x(j)=0.
         p(j)=0.
         r(j)=b(j)
         v(j)=0.
      end do
      call dot(n,b,r,rho)
      call dot(n,b,b,bb)
      iter=0
      call dot(n,r,r,test)
      test=sqrt(test/bb)
      write(6,100) iter,test
      do iter=1,itmax
         beta=(rho/rho_old)*(alpha/omega)
         do j=1,n
            p(j)=r(j)+beta*(p(j)-omega*v(j))
         end do
         call suba(1.0, p, 0.0, v)
	 call dot(n,b,v,alpha)
         write(6,*) ' got to call dot'
         alpha=rho/alpha
         do j=1,n
            s(j)=r(j)-alpha*v(j)
         end do
         write(6,*) ' S=RLOOP'
         call suba(1.0, s, 0.0, t)
	 call dot(n,t,s,omeganum)
	 call dot(n,t,t,omegaden)
         omega=omeganum/omegaden
         rho_old=rho
	 call dot(n,b,t,rho)
         rho=-omega*rho
         do j=1,n
            x(j)=x(j)+alpha*p(j)+omega*s(j)
         end do
         write(6,*) ' x(j)=x(j) loop'
         do j=1,n
            r(j)=s(j)-omega*t(j)
         end do
         call dot(n,r,r,test)
         write(6,*) ' call dot'
         write(6,*) ' test = ',test,' bb = ',bb
	 test=sqrt(test/bb)
         write(6,100) iter,test
         write(6,*) ' write iter,test'
         if (test.lt.wfac.or.iter.eq.1) then
            wfac=test
	    jwfac=iter
	    do j=1,n
	       ubar(j)=x(j)
            end do
            write(6,*) ' inside if'
         endif
         write(6,*) ' before test crit test'
         write(6,*) ' test = ',test,' crit = ',crit
	 if (test.lt.crit) return
         write(6,*) ' after test crit test'
      end do
      write(6,*)
      return
 100  format(i4,2x,1pg12.4)
      end
      subroutine dot(n,x,y,s)
      implicit real*8(a-h,o-z)
      dimension x(n),y(n)
      s=0.
      do i=1,n
         s=s+x(i)*y(i)
      end do
      s=s/n
      return
      end
