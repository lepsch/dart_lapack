      void zdotcsub(n,x,incx,y,incy,dotc) {

      // external zdotc
      Complex zdotc,dotc;
      int     n,incx,incy;
      Complex x(*),y(*);

      dotc=zdotc(n,x,incx,y,incy);
      return;
      end;
