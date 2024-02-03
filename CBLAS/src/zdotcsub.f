      subroutine zdotcsub(n,x,incx,y,incy,dotc)
c
      external zdotc
      double complex zdotc,dotc
      int     n,incx,incy
      double complex x(*),y(*)
c
      dotc=zdotc(n,x,incx,y,incy)
      return
      end
