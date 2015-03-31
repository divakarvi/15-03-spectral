c     This file contains the Fast Legendre Transform code for double
c     precision.  There are 3 external entries, exactly analogous to
c     an FFT.  They are plini, plf, and plb, which perform the initialization,
c     forward computation, and backward computation.  The calling sequence
c     is as follows:
c         plini(n, wsave)
c             n is the size of the input vector to be transformed.  n must
c             be a power of 2 greater than or equal to 64.  wsave is a
c             double-precision array of length at least 199.375*n-2295, which
c             is used to store values computed in the initialization.
c         plf(n, v, wsave)
c             n is as above.  v is a double-precision array of length at
c             least n.  The n entries of v represent a function f
c             and contain the coefficients of f's expansion in 
c             the Legendre polynomials of degree 0 through n-1.  wsave
c             is as it was initialized in plini.  Upon return, v contains
c             the coefficients of the expansion in Chebyshev polynomials,
c             degrees 0 through n-1.
c         plb(n, v, wsave)
c             n and wsave are as above.  v contains the coefficients of
c             the expansion of a function f in Chebyshev polynomials
c             of degree 0 through n-1.  Upon return, v contains coefficients
c             of f's expansion in the Legendre polynomials of degree 0
c             through n-1.
c
      subroutine plini(n, wsave)
      implicit real *8 (a-h,o-z)
      real *8 wsave(*)
      integer *4 i(10)
      call plhat(n, i)
      call plinj(n, i(1), i(2), wsave(i(3)), wsave(i(4)),
     +     wsave(i(5)), wsave(i(6)), wsave(i(7)), wsave(i(8)))
      return
      end

      subroutine plf(n, v, wsave)
      real *8 v(*), wsave(*)
      integer *4 i(10)
      call plhat(n, i)
      call plfb(v, n, i(1), i(2),
     +     wsave(i(3)), wsave(i(4)), wsave(i(5)), wsave(i(6)),
     +     wsave(i(9)), wsave(i(10)))
      v(1)=v(1)/2
      return
      end

      subroutine plb(n, v, wsave)
      real *8 v(*), wsave(*)
      integer *4 i(10)
      call plhat(n, i)
      v(1)=v(1)*2
      call plfb(v, n, i(1), i(2),
     +     wsave(i(3)), wsave(i(4)), wsave(i(7)), wsave(i(8)),
     +     wsave(i(9)), wsave(i(10)))
      return
      end


      subroutine plhat(n, i)
      integer *4 i(*)
      nc=18
      n0=64
      i(1)=nc
      i(2)=n0
      i(3)=1
      i(4)=i(3)+2*nc*nc
      i(5)=i(4)+n0*nc
      i(6)=i(5)+3*n/n0*nc*nc
      i(7)=i(6)+n/2+(3*n-2*n0)*n0/4
      i(8)=i(7)+3*n/n0*nc*nc
      i(9)=i(8)+n/2+(3*n-2*n0)*n0/4
      i(10)=i(9)+2*n*nc
ccc      print *, 'In plhat, n, space:', n, i(10)+2*n*nc
      return
      end


      subroutine plinj(n, nc, n0, ucoef, fralp,
     +     axyfar, axyner, bxyfar, bxyner)
      external plaxy, plaxyd, plbxy, plbxyd
      real *8 ucoef(*), fralp(*)
      real *8 axyfar(*), axyner(*), bxyfar(*), bxyner(*)
      real *8 cnode(30), dcnode(60), unode(128)
      nm1=n-1
      call plnode(cnode, dcnode, unode, n0, nc)
      call plcoef(ucoef, cnode, dcnode, nc, 2*nc)
      call plcoef(fralp, cnode, unode, nc, n0)
      call plfar(axyfar, nm1, n0, cnode, nc, plaxy)
      call plner(axyner, nm1, n0, plaxy, plaxyd)
      call plfar(bxyfar, nm1, n0, cnode, nc, plbxy)
      call plner(bxyner, nm1, n0, plbxy, plbxyd)
      return
      end


      subroutine plfb(v, n, nc, n0, upcoef, fralph, ffar, fnear,
     +                atree, btree)
      real *8 v(*), upcoef(*), fralph(*), ffar(*), fnear(*)
      real *8 atree(*), btree(*)
      nm1=n-1
      np1=n+1
      call plup(v, atree, upcoef, fralph, nm1, n0, nc)
      call placr(atree, btree, ffar, nm1, n0, nc)
      call plnei(v, fnear, nm1, n0, nc)
      call pldn(btree, v, upcoef, fralph, nm1, n0, nc)
      return
      end


c     Build atree from the bottom up.
c
c     Input:
c       alpha(0:n)        the coefficients of Pi(cos x)
c       fralph(nc,n0)     the interpolation coefficients from alphas
c                         to finest grain Chebyshev nodes
c       upcoef(nc,2*nc)   the interpolation coefficients from finer to coarser
c                         Chebyshev nodes
c       n0                number of alphas per finest subinterval
c       nc                number of Chebyshev nodes per subinterval
c       n                 last alpha index (number of alphas is n+1)
c     Output:
c       atree(nc,2,*)     tree of interpolation coefficients from the alphas
c
      subroutine plup(alpha, atree, upcoef, fralph, n, n0, nc)
      integer n, n0, nc, nintv, koff, i, j, k, iptr0, iptr1, lev
      real *8 alpha(0:n), atree(nc,2,*), upcoef(nc,2*nc)
      real *8 fralph(nc, n0), sum1, sum2
c
c     Build first level of atree from alpha
c     Treat even and odd-indexed alphas separately
c
      nintv=(n+1)/n0
      koff=-1
      do 1600 i=1, nintv
         do 1400 j=1, nc
            sum1=0
            sum2=0
            do 1000 k=1, n0, 2
               sum1=sum1+fralph(j,k)*alpha(k+koff)
               sum2=sum2+fralph(j,k+1)*alpha(k+1+koff)
 1000       continue
            atree(j,1,i)=sum1
            atree(j,2,i)=sum2
 1400    continue
         koff=koff+n0
 1600 continue
c
c     Build rest of atree, each level from the one below
c
      iptr0=1
      iptr1=nintv+1
      do 2400 lev=2, 30
         nintv=nintv/2
         if (nintv .lt. 4)  goto 2600
         do 2200 i=1, nintv
            do 2000 j=1, nc
               atree(j,1,iptr1)=0
               atree(j,2,iptr1)=0
 2000       continue
            call plmatv(upcoef, atree(1,1,iptr0), atree(1,1,iptr1), nc)
            call plmatv(upcoef(1,1+nc), atree(1,1,iptr0+1),
     +           atree(1,1,iptr1), nc)
            call plmatv(upcoef, atree(1,2,iptr0), atree(1,2,iptr1), nc)
            call plmatv(upcoef(1,1+nc), atree(1,2,iptr0+1),
     +           atree(1,2,iptr1), nc)
            iptr0=iptr0+2
            iptr1=iptr1+1
 2200    continue
 2400 continue
 2600 continue
      return
      end


c     Add matrix-vector product to another vector.
c
      subroutine plmatv(xmat, vcin, vcout, n)
      integer n, i, j
      real *8 xmat(n,n), vcin(n), vcout(n), s
      do 1200 i=1, n
         s=0
         do 1000 j=1, n
            s=s+xmat(i,j)*vcin(j)
 1000    continue
         vcout(i)=vcout(i)+s
 1200 continue
      return
      end


c     Do the little matrix-vector multiplies for each square of the matrix axy.
c
c     Input:
c       atree(nc,2,*)     tree of interpolation coefficients from the alphas
c       axyfar(nc,nc,*)   the interpolated matrix squares 
c       nc                number of Chebyshev nodes per subinterval
c       n0                number of alphas per finest subinterval
c       n                 last alpha index (number of alphas is n+1)
c     Output:
c       btree(nc,2,*)     tree of interpolation coefficients for the betas
c
      subroutine placr(atree, btree, axyfar, n, n0, nc)
      integer n, n0, nc, iptr1, iptr2, nintv, lev, i, ilim
      real *8 atree(*), btree(*), axyfar(nc,nc,*)
      logical left
      ilim=4*nc*(n+1)/n0
      do 1000 i=1, ilim
         btree(i)=0
 1000 continue
c
c     For each level of atree, produce the same level of btree by
c     doing miniature matrix-vector multiplies.
c       iptr1          interval number of interval in 'btree'
c       iptr1+2 or 3   interval number of corresp. interval in 'atree'
c       iptr2          box number of square in 'axyfar'
c
      iptr1=1
      iptr2=1
      nintv=(n+1)/n0-2
      do 1400 lev=2, 30
         if (nintv .lt. 2)  goto 1600
         left= .true.
         do 1200 i=1, nintv
            call            plsacr(atree,btree,axyfar,nc,iptr1,iptr2,2)
            if (left)  call plsacr(atree,btree,axyfar,nc,iptr1,iptr2,3)
            iptr1=iptr1+1
            left=.not. left
 1200    continue
         iptr1=iptr1+2
         nintv=nintv/2-1
 1400 continue
 1600 continue
      return
      end


c     Perform matrix-vector multiplies for one square block of axy.
c     
      subroutine plsacr(atree,btree,axyfar,nc,iptr1,iptr2,ioff)
      integer nc, iptr1, iptr2, ioff, iptr3, j, k
      real *8 atree(nc,2,*), btree(nc,2,*), axyfar(nc,nc,*), sum1, sum2
      iptr3=iptr1+ioff
      do 1200 j=1, nc
         sum1=0
         sum2=0
         do 1000 k=1, nc
            sum1=sum1 + axyfar(j,k,iptr2)*atree(k,1,iptr3)
            sum2=sum2 + axyfar(j,k,iptr2)*atree(k,2,iptr3)
 1000    continue
         btree(j,1,iptr1)=btree(j,1,iptr1)+sum1
         btree(j,2,iptr1)=btree(j,2,iptr1)+sum2
 1200 continue
      iptr2=iptr2+1
      return
      end


c     Reduce btree from the top down.
c
c     Input:
c       btree(nc,2,*)     tree of interpolation coefficients to the betas
c       tobeta(n0,nc)     the interpolation coefficients from finest grain
c                         Chebyshev nodes to beta
c       upcoef(nc,2*nc)   the interpolation coefficients from finer to coarser
c                         Chebyshev nodes
c       n0                number of alphas per finest subinterval
c       nc                number of Chebyshev nodes per subinterval
c       n                 last alpha index (number of alphas is n+1)
c     Output:
c       beta(0:n)         the coefficients of cos ix
c
      subroutine pldn(btree, beta, upcoef, fralph, n, n0, nc)
      integer n, n0, nc, iptr0, iptr1, nintv, lev, i, j, k, joff
      real *8 btree(nc,2,*), beta(0:n), upcoef(nc,2*nc)
      real *8 fralph(nc, n0), sum1, sum2, sum3, sum4
c
c     Trickle the information at higher intervals in btree down to lower.
c
      iptr1=2*(n+1)/n0-4
      iptr0=iptr1-4
      nintv=2
      do 1600 lev=2, 30
         if (iptr0 .lt. 2)  goto 1800
         nintv=nintv*2
         do 1400 i=1, nintv
            do 1200 j=1, nc
               sum1=0
               sum2=0
               sum3=0
               sum4=0
               do 1000 k=1, nc
                  sum1=sum1 + upcoef(k,j   )*btree(k,1,iptr1)
                  sum2=sum2 + upcoef(k,j   )*btree(k,2,iptr1)
                  sum3=sum3 + upcoef(k,j+nc)*btree(k,1,iptr1)
                  sum4=sum4 + upcoef(k,j+nc)*btree(k,2,iptr1)
 1000          continue
               btree(j,1,iptr0-1)=btree(j,1,iptr0-1)+sum1
               btree(j,2,iptr0-1)=btree(j,2,iptr0-1)+sum2
               btree(j,1,iptr0  )=btree(j,1,iptr0  )+sum3
               btree(j,2,iptr0  )=btree(j,2,iptr0  )+sum4
 1200       continue
            iptr0=iptr0-2
            iptr1=iptr1-1
 1400    continue
 1600 continue
 1800 continue
c
c     Build betas from first level of btree.
c     Treat even and odd-indexed betas separately.
c
      nintv=(n+1)/n0
      joff=-1
      do 2400 i=1, nintv
         do 2200 j=1, n0, 2
            sum1=0
            sum2=0
            do 2000 k=1, nc
               sum1=sum1+fralph(k,j  )*btree(k,1,i)
               sum2=sum2+fralph(k,j+1)*btree(k,2,i)
 2000       continue
            beta(j+joff  )=beta(j+joff  )+sum1
            beta(j+joff+1)=beta(j+joff+1)+sum2
 2200    continue
         joff=joff+n0
 2400 continue
      return
      end


c     Handle  near  interactions directly
c
c     Input:
c       alpha(0:n)        coefficients of  Pi(cos x)
c       beta(0:n)         coefficients of  cos ix
c       axyner(nc,nc,*)   the matrix triangles and squares near the diagonal
c       nc                number of Chebyshev nodes per subinterval
c       n0                number of alphas per finest subinterval
c       n                 last alpha index (number of alphas is n+1)
c     Output:
c       alpha(0:n)        these get the effect of  near  interactions
c
      subroutine plnei(alpha, axyner, n, n0, nc)
      integer n, n0, nc, jhigh, i, j, k
      real *8 alpha(0:n), axyner(*), sum
      k=1
      do 1200 i=0, n
         jhigh=min0(n, n0*(2+i/n0)-1)
         sum=0
         do 1000 j=i, jhigh, 2
            sum=sum + axyner(k)*alpha(j)
            k=k+1
 1000    continue
         alpha(i)=sum
 1200 continue
      return
      end
      


c     Precompute all values of 'fun' corresponding to "far" interactions.
c
c     Input:
c       n            highest alpha index (n+1 is number of alpha values)
c       n0           number of alpha values in smallest interval
c       cnode        Chebyshev nodes for interval [-1,1]
c       nc           number of Chebyshev nodes
c       fun          function name: calling sequence is fun(x,y)
c     Output:
c       v(nc,nc,*)   resulting function values
c
      subroutine plfar(v, n, n0, cnode, nc, fun)
      integer n, n0, nc, intv, nintv, iwidth, i, j, ncmax
      parameter (ncmax=30)
      real *8 v(*), cnode(nc), x(ncmax), y(ncmax), fun
      intv=1      
      iwidth=n0
      nintv=(n+1)/(2*n0)-1
      do 1400 i=1, 30
         if (nintv .lt. 1)  goto 1600
         do 1000 j=1, nc
            x(j)=iwidth*(cnode(j)+1)/2
            y(j)=x(j)+2*iwidth
 1000    continue
         do 1200 j=1, nintv
            call plfsqr(v, x, y, nc, intv, fun)
            call plshft(y, nc, iwidth)
            call plfsqr(v, x, y, nc, intv, fun)
            call plshft(x, nc, iwidth)
            call plfsqr(v, x, y, nc, intv, fun)
            call plshft(x, nc, iwidth)
            call plshft(y, nc, iwidth)
 1200    continue
         iwidth=iwidth*2
         nintv=(nintv-1)/2
 1400 continue
 1600 continue
      return
      end


c     Shift the array x of n elements by i.
c
      subroutine plshft(x, n, i)
      real *8 x(n)
      integer n, i, j
      do 1000 j=1, n
         x(j)=x(j)+i
 1000 continue
      return
      end


c     Record values of 'fun' on a square
c
      subroutine plfsqr(v, x, y, n, intv, fun)
      real *8 v(n,n,*), x(n), y(n), fun
      integer n, intv, i, j
      do 1200 i=1, n
         do 1000 j=1, n
            v(i,j,intv)=fun(x(i),y(j))
 1000    continue
 1200 continue
      intv=intv+1
      return
      end


c     Precompute all values of 'fun' corresponding to "near" interactions.
c
c     Input:
c       n            highest alpha index (n+1 is number of alpha values)
c       n0           number of alpha values in smallest interval
c       fun          function name: calling sequence is fun(x,y)
c       fundg        function on diagonal: calling sequence is fundg(x)
c     Output:
c       v(n0,n0,*)   resulting function values
c
      subroutine plner(v, n, n0, fun, fundg)
      integer n, n0, i, j, k, jhigh
      real *8 v(*), fun, fundg, x, y
      k=1
      do 1200 i=0, n
         x=i
         v(k)=fundg(x)
         k=k+1
         jhigh=min0(n, n0*(2+i/n0)-1)
         do 1000 j=i+2, jhigh, 2
            y=j
            v(k)=fun(x,y)
            k=k+1
 1000    continue
 1200 continue
      return
      end



c     Calculate Chebyshev, double Chebyshev, and uniform-spaced nodes
c
c     Input:
c       n0            number of uniform-spaced nodes
c       nc            number of Chebyshev nodes
c     Output:
c       cnode(nc)     Chebyshev nodes on [-1,1]
c       dcnode(2*nc)  double Chebyshev nodes (i.e. nodes on [-1,0] and [0,1])
c       unode(n0)     uniform-spaced nodes on [-1,1]  (-1,-1+2/n0,...,1-2/n0)
c
      subroutine plnode(cnode, dcnode, unode, n0, nc)
      integer n0, nc, i
      real *8 cnode(nc), dcnode(2*nc), unode(n0), one, pi2
      one=1
      pi2=atan(one)*2
      do 1000 i=1, nc
         cnode(i)=cos((2*(nc-i)+1)*pi2/nc)
         dcnode(i)=(cnode(i)-1)/2
         dcnode(i+nc)=(cnode(i)+1)/2
 1000 continue
      do 1200 i=1, n0
         unode(i)=2*(i-one)/n0-1
 1200 continue
      return
      end


c     Calculate polynomial interpolation coefficients for a function
c     evaluated at specified nodes to obtain its value at specified points.
c
c     Input:
c       cnode(nc)    the nodes at which the function is available
c       point(np)    the points for which the function is to be interpolated
c       nc
c       np
c     Output:                               nc
c       coef         so that f(point(i)) ~ SUM coef(i,j)*f(cnode(j))
c                                          j=1
c
      subroutine plcoef(coef, cnode, point, nc, np)
      integer nc, np, ncmax, i, j, k
      parameter (ncmax=30)
      real *8 coef(nc,np), cnode(nc), point(np), denom(ncmax), d,p,eps
      eps=1
      eps=eps/10
      eps=eps**(13+nc)
c
c     Construct denominators of Lagrange interpolation formula
c
      do 1200 i=1, nc
         d=1
         do 1000 j=1, nc
            if (i .eq. j)  goto 1000
            d=d*(cnode(i)-cnode(j))
 1000    continue
         denom(i)=d
 1200 continue
c
c     Calculate coefficients
c
      do 1900 i=1, np
         p=point(i)
         d=1
         do 1400 j=1, nc
            d=d*(p-cnode(j))
 1400    continue
         if (d .lt. eps) goto 1600
         do 1500 j=1, nc
            coef(j,i)=d/((p-cnode(j))*denom(j))
 1500    continue
 1600    continue
c
c     p nearly coincides with one of the nodes.  React accordingly.
c
         do 1800 j=1, nc
            d=1
            do 1700 k=1, nc
               if (j .eq. k)  goto 1700
               d=d*(p-cnode(k))
 1700       continue
            coef(j,i)=d/denom(j)
 1800    continue
 1900 continue
      return
      end



c     Axy, for forward transform.
c
      function plaxy(x, y)
      real *8 plaxy, x, y, plfff, twobpi, half
      data twobpi/.6366197723675814d0/,half/.5d0/
      plaxy = twobpi * plfff((y+x)*half) * plfff((y-x)*half)
      return
      end

      function plaxyd(x)
      real *8 x, plaxyd, plfff, sqtp2i
      data sqtp2i/.11283791670955126d1/
      plaxyd=sqtp2i*plfff(x)
      return
      end



c     Bxy, for inverse transform.
c
      function plbxy(x, y)
      implicit real *8 (a-h,o-z)
      real *8 plbxy, x, y, plfff, half
      data half/.5d0/
      plbxy=-y*(x+half) / ((y+x+1)*(y-x))
     +      * plfff((y+x-1)*half)*plfff((y-x-2)*half)
      return
      end

      function plbxyd(x)
      real *8 x, plbxyd, plfff, sqtpi2
      data sqtpi2/.88622692545275800d0/
      plbxyd=sqtpi2/plfff(x)
      return
      end



c     Double precision function to compute gamma(x+.5)/gamma(x+1) for
c     x in the set {0, .5, 1, 1.5, 2,...,14.5} or the interval [15,infinity).
c
      function plfff(x)
      integer index
      real *8 x, y, plfff
      real *8 a0,a1,a2,a3,a4,a5, b0,b1,b2,b3,b4,b5
      real *8 c0,c1,c2,c3,c4,c5, d0,d1,d2,d3,d4,d5, fn(30)
      data a0,a1,a2,a3,a4,a5/
     +      0.99999999999999999d+00,-0.12499999999996888d+00,
     +      0.78124999819509772d-02, 0.48828163023526451d-02,
     +     -0.64122689844951054d-03,-0.15070098356496836d-02/
      data b0,b1,b2,b3,b4,b5/
     +      0.99999999999974725d+00,-0.12499999994706490d+00,
     +      0.78124954632722315d-02, 0.48830152125039076d-02,
     +     -0.64579205161159155d-03,-0.14628616278637035d-02/
      data c0,c1,c2,c3,c4,c5/
     +      0.99999999999298844d+00,-0.12499999914033463d+00,
     +      0.78124565447111342d-02, 0.48839648427432678d-02,
     +     -0.65752249058053233d-03,-0.14041419931494052d-02/
      data d0,d1,d2,d3,d4,d5/
     +      0.99999999996378690d+00,-0.12499999657282661d+00,
     +      0.78123659464717666d-02, 0.48855685911602214d-02,
     +     -0.67176366234107532d-03,-0.13533949520771154d-02/
      data fn/0.17724538509055159d+01,0.11283791670955126d+01,
     +        0.88622692545275794d+00,0.75225277806367508d+00,
     +        0.66467019408956851d+00,0.60180222245094006d+00,
     +        0.55389182840797380d+00,0.51583047638652002d+00,
     +        0.48465534985697706d+00,0.45851597901023999d+00,
     +        0.43618981487127934d+00,0.41683270819112728d+00,
     +        0.39984066363200604d+00,0.38476865371488672d+00,
     +        0.37128061622971992d+00,0.35911741013389425d+00,
     +        0.34807557771536241d+00,0.33799285659660638d+00,
     +        0.32873804562006448d+00,0.32020375888099550d+00,
     +        0.31230114333906123d+00,0.30495596083904336d+00,
     +        0.29810563682364938d+00,0.29169700601995452d+00,
     +        0.28568456862266400d+00,0.28002912577915634d+00,
     +        0.27469670059871537d+00,0.26965767667622464d+00,
     +        0.26488610414876124d+00,0.26035913610118244d+00/
      if (x .lt. .01)  then
         y=100
      else
         y=1/x
      end if
c
c     The computation is divided into 34 cases, corresponding
c     to the 30 elements of the set {0, .5, 1, 1.5,...,14.5} and 4 subintervals
c     of the interval [1/.067,infinity).  The subintervals, expressed in
c     y=1/x, are (0,.02], (.02,.04], (.04,.058], (.058,.067].
c     (They were chosen to maintain double precision using polynomials of
c     degree 5.)
c
      if (y .le. .02) then
         plfff=(((((a5*y+a4)*y+a3)*y+a2)*y+a1)*y+a0)*sqrt(y)
      else if (y .le. .04) then
         plfff=(((((b5*y+b4)*y+b3)*y+b2)*y+b1)*y+b0)*sqrt(y)
      else if (y .le. .058) then
         plfff=(((((c5*y+c4)*y+c3)*y+c2)*y+c1)*y+c0)*sqrt(y)
      else if (y .le. .067) then
         plfff=(((((d5*y+d4)*y+d3)*y+d2)*y+d1)*y+d0)*sqrt(y)
      else
         index=x*2+1.01
         plfff=fn(index)
      end if
      return
      end
