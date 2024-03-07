      void zerrqr(PATH, infoc.NUNIT ) {

// -- LAPACK test routine (--
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             PATH;
      int                infoc.NUNIT;
      // ..

      int                NMAX;
      const              NMAX = 2 ;
      int                I, INFO, J;
      Complex         A( NMAX, NMAX ), AF( NMAX, NMAX ), B( NMAX ), W( NMAX ), X( NMAX );
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, ZGEQR2, ZGEQR2P, ZGEQRF, ZGEQRFP, ZUNG2R, ZUNGQR, ZUNM2R, ZUNMQR
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
      // INTRINSIC DBLE, DCMPLX

      NOUT = infoc.NUNIT;
      WRITE( NOUT, FMT = * );

      // Set the variables to innocuous values.

      for (J = 1; J <= NMAX; J++) { // 20
         for (I = 1; I <= NMAX; I++) { // 10
            A[I][J] = DCMPLX( 1.0 / DBLE( I+J ), -1.0 / DBLE( I+J ) )             AF( I, J ) = DCMPLX( 1.0 / DBLE( I+J ), -1.0 / (I+J).toDouble() );
         } // 10
         B[J] = 0.0;
         W[J] = 0.0;
         X[J] = 0.0;
      } // 20
      infoc.OK.value = true;

      // Error exits for QR factorization

      // ZGEQRF

     srnamc.SRNAMT = 'ZGEQRF';
      infoc.INFOT = 1;
      zgeqrf(-1, 0, A, 1, B, W, 1, INFO );
      chkxer('ZGEQRF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      zgeqrf(0, -1, A, 1, B, W, 1, INFO );
      chkxer('ZGEQRF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 4;
      zgeqrf(2, 1, A, 1, B, W, 1, INFO );
      chkxer('ZGEQRF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 7;
      zgeqrf(1, 2, A, 1, B, W, 1, INFO );
      chkxer('ZGEQRF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

      // ZGEQRFP

     srnamc.SRNAMT = 'ZGEQRFP';
      infoc.INFOT = 1;
      zgeqrfp(-1, 0, A, 1, B, W, 1, INFO );
      chkxer('ZGEQRFP', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      zgeqrfp(0, -1, A, 1, B, W, 1, INFO );
      chkxer('ZGEQRFP', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 4;
      zgeqrfp(2, 1, A, 1, B, W, 1, INFO );
      chkxer('ZGEQRFP', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 7;
      zgeqrfp(1, 2, A, 1, B, W, 1, INFO );
      chkxer('ZGEQRFP', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

      // ZGEQR2

     srnamc.SRNAMT = 'ZGEQR2';
      infoc.INFOT = 1;
      zgeqr2(-1, 0, A, 1, B, W, INFO );
      chkxer('ZGEQR2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      zgeqr2(0, -1, A, 1, B, W, INFO );
      chkxer('ZGEQR2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 4;
      zgeqr2(2, 1, A, 1, B, W, INFO );
      chkxer('ZGEQR2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

      // ZGEQR2P

     srnamc.SRNAMT = 'ZGEQR2P';
      infoc.INFOT = 1;
      zgeqr2p(-1, 0, A, 1, B, W, INFO );
      chkxer('ZGEQR2P', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      zgeqr2p(0, -1, A, 1, B, W, INFO );
      chkxer('ZGEQR2P', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 4;
      zgeqr2p(2, 1, A, 1, B, W, INFO );
      chkxer('ZGEQR2P', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

      // ZUNGQR

     srnamc.SRNAMT = 'ZUNGQR';
      infoc.INFOT = 1;
      zungqr(-1, 0, 0, A, 1, X, W, 1, INFO );
      chkxer('ZUNGQR', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      zungqr(0, -1, 0, A, 1, X, W, 1, INFO );
      chkxer('ZUNGQR', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      zungqr(1, 2, 0, A, 1, X, W, 2, INFO );
      chkxer('ZUNGQR', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 3;
      zungqr(0, 0, -1, A, 1, X, W, 1, INFO );
      chkxer('ZUNGQR', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 3;
      zungqr(1, 1, 2, A, 1, X, W, 1, INFO );
      chkxer('ZUNGQR', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      zungqr(2, 2, 0, A, 1, X, W, 2, INFO );
      chkxer('ZUNGQR', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 8;
      zungqr(2, 2, 0, A, 2, X, W, 1, INFO );
      chkxer('ZUNGQR', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

      // ZUNG2R

     srnamc.SRNAMT = 'ZUNG2R';
      infoc.INFOT = 1;
      zung2r(-1, 0, 0, A, 1, X, W, INFO );
      chkxer('ZUNG2R', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      zung2r(0, -1, 0, A, 1, X, W, INFO );
      chkxer('ZUNG2R', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      zung2r(1, 2, 0, A, 1, X, W, INFO );
      chkxer('ZUNG2R', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 3;
      zung2r(0, 0, -1, A, 1, X, W, INFO );
      chkxer('ZUNG2R', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 3;
      zung2r(2, 1, 2, A, 2, X, W, INFO );
      chkxer('ZUNG2R', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      zung2r(2, 1, 0, A, 1, X, W, INFO );
      chkxer('ZUNG2R', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

      // ZUNMQR

     srnamc.SRNAMT = 'ZUNMQR';
      infoc.INFOT = 1;
      zunmqr('/', 'N', 0, 0, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('ZUNMQR', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      zunmqr('L', '/', 0, 0, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('ZUNMQR', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 3;
      zunmqr('L', 'N', -1, 0, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('ZUNMQR', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 4;
      zunmqr('L', 'N', 0, -1, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('ZUNMQR', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      zunmqr('L', 'N', 0, 0, -1, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('ZUNMQR', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      zunmqr('L', 'N', 0, 1, 1, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('ZUNMQR', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      zunmqr('R', 'N', 1, 0, 1, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('ZUNMQR', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 7;
      zunmqr('L', 'N', 2, 1, 0, A, 1, X, AF, 2, W, 1, INFO );
      chkxer('ZUNMQR', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 7;
      zunmqr('R', 'N', 1, 2, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('ZUNMQR', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 10;
      zunmqr('L', 'N', 2, 1, 0, A, 2, X, AF, 1, W, 1, INFO );
      chkxer('ZUNMQR', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 12;
      zunmqr('L', 'N', 1, 2, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('ZUNMQR', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 12;
      zunmqr('R', 'N', 2, 1, 0, A, 1, X, AF, 2, W, 1, INFO );
      chkxer('ZUNMQR', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

      // ZUNM2R

     srnamc.SRNAMT = 'ZUNM2R';
      infoc.INFOT = 1;
      zunm2r('/', 'N', 0, 0, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('ZUNM2R', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      zunm2r('L', '/', 0, 0, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('ZUNM2R', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 3;
      zunm2r('L', 'N', -1, 0, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('ZUNM2R', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 4;
      zunm2r('L', 'N', 0, -1, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('ZUNM2R', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      zunm2r('L', 'N', 0, 0, -1, A, 1, X, AF, 1, W, INFO );
      chkxer('ZUNM2R', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      zunm2r('L', 'N', 0, 1, 1, A, 1, X, AF, 1, W, INFO );
      chkxer('ZUNM2R', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      zunm2r('R', 'N', 1, 0, 1, A, 1, X, AF, 1, W, INFO );
      chkxer('ZUNM2R', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 7;
      zunm2r('L', 'N', 2, 1, 0, A, 1, X, AF, 2, W, INFO );
      chkxer('ZUNM2R', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 7;
      zunm2r('R', 'N', 1, 2, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('ZUNM2R', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 10;
      zunm2r('L', 'N', 2, 1, 0, A, 2, X, AF, 1, W, INFO );
      chkxer('ZUNM2R', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

      // Print a summary line.

      alaesm(PATH, infoc.OK, NOUT );

      }
