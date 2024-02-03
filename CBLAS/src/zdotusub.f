      subroutine zdotusub(n,x,incx,y,incy,dotu)
c
      // external zdotu
      double complex zdotu,dotu
      int     n,incx,incy;
      double complex x(*),y(*)
c
      dotu=zdotu(n,x,incx,y,incy)
      return
      end
