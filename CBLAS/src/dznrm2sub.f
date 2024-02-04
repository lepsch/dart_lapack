      void dznrm2sub(n,x,incx,nrm2) {

      // external dznrm2
      double           dznrm2,nrm2;
      int     n,incx;
      Complex x(*);

      nrm2=dznrm2(n,x,incx);
      return;
      end;
