      subroutine izamaxsub(n,x,incx,iamax)
c
      external izamax
      integer  izamax,iamax
      integer n,incx
      double complex x(*)
c
      iamax=izamax(n,x,incx)
      return
      end