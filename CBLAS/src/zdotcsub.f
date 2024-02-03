      void zdotcsub(n,x,incx,y,incy,dotc) {

      // external zdotc
      double complex zdotc,dotc;
      int     n,incx,incy;
      double complex x(*),y(*);

      dotc=zdotc(n,x,incx,y,incy);
      return;
      end;
