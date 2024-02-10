      void zerrab(infoc.NUNIT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                infoc.NUNIT;
      // ..

      int                NMAX;
      const              NMAX = 4 ;
      int                I, INFO, ITER, J;
      int                IP( NMAX );
      Complex         A( NMAX, NMAX ), AF( NMAX, NMAX ), B( NMAX ), C( NMAX ), R( NMAX ), R1( NMAX ), R2( NMAX ), W( 2*NMAX ), X( NMAX );
      Complex         WORK(1);
      Complex            SWORK(1);
      double             RWORK(1);
      // ..
      // .. External Functions ..
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHKXER, ZCGESV
      // ..
      // .. Scalars in Common ..
      bool               infoc.LERR, infoc.OK;
      String            srnamc.SRNAMT;
      int                infoc.INFOT, NOUT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / infoc.INFOT, NOUT, infoc.OK, infoc.LERR
      // COMMON / SRNAMC /srnamc.SRNAMT
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE

      NOUT = infoc.NUNIT;
      WRITE( NOUT, FMT = * );

      // Set the variables to innocuous values.

      for (J = 1; J <= NMAX; J++) { // 20
         for (I = 1; I <= NMAX; I++) { // 10
            A[I][J] = 1.0 / (I+J).toDouble();
            AF[I][J] = 1.0 / (I+J).toDouble();
         } // 10
         B[J] = 0.0;
         R1[J] = 0.0;
         R2[J] = 0.0;
         W[J] = 0.0;
         X[J] = 0.0;
         C[J] = 0.0;
         R[J] = 0.0;
         IP[J] = J;
      } // 20
      infoc.OK = true;

     srnamc.SRNAMT = 'ZCGESV';
      infoc.INFOT = 1;
      zcgesv(-1,0,A,1,IP,B,1,X,1,WORK,SWORK,RWORK,ITER,INFO);
      chkxer('ZCGESV', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      zcgesv(0,-1,A,1,IP,B,1,X,1,WORK,SWORK,RWORK,ITER,INFO);
      chkxer('ZCGESV', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 4;
      zcgesv(2,1,A,1,IP,B,2,X,2,WORK,SWORK,RWORK,ITER,INFO);
      chkxer('ZCGESV', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 7;
      zcgesv(2,1,A,2,IP,B,1,X,2,WORK,SWORK,RWORK,ITER,INFO);
      chkxer('ZCGESV', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 9;
      zcgesv(2,1,A,2,IP,B,2,X,1,WORK,SWORK,RWORK,ITER,INFO);
      chkxer('ZCGESV', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

      // Print a summary line.

      if ( infoc.OK ) {
         WRITE( NOUT, FMT = 9999 )'ZCGESV';
      } else {
         WRITE( NOUT, FMT = 9998 )'ZCGESV';
      }

 9999 FORMAT(' ${.a6} drivers passed the tests of the error exits' );
 9998 FORMAT( ' *** ${.a6} drivers failed the tests of the error exits ***' );

      }
