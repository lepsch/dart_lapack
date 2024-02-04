      void cdotcsub(n,x,incx,y,incy,dotc) {

      // external cdotc
      Complex cdotc,dotc;
      int     n,incx,incy;
      Complex x(*),y(*);

      dotc=cdotc(n,x,incx,y,incy);
      return;
      end;
