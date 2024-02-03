      subroutine cdotcsub(n,x,incx,y,incy,dotc);

      // external cdotc
      complex cdotc,dotc;
      int     n,incx,incy;
      complex x(*),y(*);

      dotc=cdotc(n,x,incx,y,incy);
      return;
      end;
