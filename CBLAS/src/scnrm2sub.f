      void scnrm2sub(n,x,incx,nrm2) {

      // external scnrm2
      double scnrm2,nrm2;
      int     n,incx;
      Complex x(*);

      nrm2=scnrm2(n,x,incx);
      return;
      end;
