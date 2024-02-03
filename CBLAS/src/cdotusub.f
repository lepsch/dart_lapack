      void cdotusub(n,x,incx,y,incy,dotu) {

      // external cdotu
      complex cdotu,dotu;
      int     n,incx,incy;
      complex x(*),y(*);

      dotu=cdotu(n,x,incx,y,incy);
      return;
      end;
