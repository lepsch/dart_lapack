      void scnrm2sub(n,x,incx,nrm2) {

      // external scnrm2
      real scnrm2,nrm2;
      int     n,incx;
      complex x(*);

      nrm2=scnrm2(n,x,incx);
      return;
      end;
