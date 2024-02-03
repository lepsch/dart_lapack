      subroutine dnrm2sub(n,x,incx,nrm2)
c
      external dnrm2
      double precision dnrm2,nrm2
      int     n,incx
      double precision x(*)
c
      nrm2=dnrm2(n,x,incx)
      return
      end
