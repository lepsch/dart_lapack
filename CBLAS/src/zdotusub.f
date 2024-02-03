      subroutine zdotusub(n,x,incx,y,incy,dotu)

      // external zdotu
      double complex zdotu,dotu
      int     n,incx,incy;
      double complex x(*),y(*)

      dotu=zdotu(n,x,incx,y,incy)
      return
      end
