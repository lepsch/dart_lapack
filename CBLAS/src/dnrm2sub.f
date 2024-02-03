      void dnrm2sub(n,x,incx,nrm2) {

      // external dnrm2
      double           dnrm2,nrm2;
      int     n,incx;
      double           x(*);

      nrm2=dnrm2(n,x,incx);
      return;
      end;
