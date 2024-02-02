      subroutine cdotusub(n,x,incx,y,incy,dotu)
c
      external cdotu
      complex cdotu,dotu
      integer n,incx,incy
      complex x(*),y(*)
c
      dotu=cdotu(n,x,incx,y,incy)
      return
      end