      subroutine cdotcsub(n,x,incx,y,incy,dotc)
c
      // external cdotc
      complex cdotc,dotc
      int     n,incx,incy;
      complex x(*),y(*)
c
      dotc=cdotc(n,x,incx,y,incy)
      return
      end
