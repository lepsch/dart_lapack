      subroutine dznrm2sub(n,x,incx,nrm2)
c
      external dznrm2
      double           dznrm2,nrm2;
      int     n,incx
      double complex x(*)
c
      nrm2=dznrm2(n,x,incx)
      return
      end
