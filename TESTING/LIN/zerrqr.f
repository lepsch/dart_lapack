      void zerrqr(PATH, NUNIT ) {

// -- LAPACK test routine (--
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             PATH;
      int                NUNIT;
      // ..

// =====================================================================

      // .. Parameters ..
      int                NMAX;
      const              NMAX = 2 ;
      // ..
      // .. Local Scalars ..
      int                I, INFO, J;
      // ..
      // .. Local Arrays ..
      COMPLEX*16         A( NMAX, NMAX ), AF( NMAX, NMAX ), B( NMAX ), W( NMAX ), X( NMAX );
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, ZGEQR2, ZGEQR2P, ZGEQRF, ZGEQRFP, ZUNG2R, ZUNGQR, ZUNM2R, ZUNMQR
      // ..
      // .. Scalars in Common ..
      bool               LERR, OK;
      String             SRNAMT;
      int                INFOT, NOUT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / INFOT, NOUT, OK, LERR
      // COMMON / SRNAMC / SRNAMT
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DCMPLX
      // ..
      // .. Executable Statements ..

      NOUT = NUNIT;
      WRITE( NOUT, FMT = * );

      // Set the variables to innocuous values.

      for (J = 1; J <= NMAX; J++) { // 20
         for (I = 1; I <= NMAX; I++) { // 10
            A( I, J ) = DCMPLX( 1.0 / DBLE( I+J ), -1.0 / DBLE( I+J ) )             AF( I, J ) = DCMPLX( 1.0 / DBLE( I+J ), -1.0 / DBLE( I+J ) );
         } // 10
         B( J ) = 0.0;
         W( J ) = 0.0;
         X( J ) = 0.0;
      } // 20
      OK = true;

      // Error exits for QR factorization

      // ZGEQRF

      SRNAMT = 'ZGEQRF';
      INFOT = 1;
      zgeqrf(-1, 0, A, 1, B, W, 1, INFO );
      chkxer('ZGEQRF', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      zgeqrf(0, -1, A, 1, B, W, 1, INFO );
      chkxer('ZGEQRF', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      zgeqrf(2, 1, A, 1, B, W, 1, INFO );
      chkxer('ZGEQRF', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      zgeqrf(1, 2, A, 1, B, W, 1, INFO );
      chkxer('ZGEQRF', INFOT, NOUT, LERR, OK );

      // ZGEQRFP

      SRNAMT = 'ZGEQRFP';
      INFOT = 1;
      zgeqrfp(-1, 0, A, 1, B, W, 1, INFO );
      chkxer('ZGEQRFP', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      zgeqrfp(0, -1, A, 1, B, W, 1, INFO );
      chkxer('ZGEQRFP', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      zgeqrfp(2, 1, A, 1, B, W, 1, INFO );
      chkxer('ZGEQRFP', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      zgeqrfp(1, 2, A, 1, B, W, 1, INFO );
      chkxer('ZGEQRFP', INFOT, NOUT, LERR, OK );

      // ZGEQR2

      SRNAMT = 'ZGEQR2';
      INFOT = 1;
      zgeqr2(-1, 0, A, 1, B, W, INFO );
      chkxer('ZGEQR2', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      zgeqr2(0, -1, A, 1, B, W, INFO );
      chkxer('ZGEQR2', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      zgeqr2(2, 1, A, 1, B, W, INFO );
      chkxer('ZGEQR2', INFOT, NOUT, LERR, OK );

      // ZGEQR2P

      SRNAMT = 'ZGEQR2P';
      INFOT = 1;
      zgeqr2p(-1, 0, A, 1, B, W, INFO );
      chkxer('ZGEQR2P', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      zgeqr2p(0, -1, A, 1, B, W, INFO );
      chkxer('ZGEQR2P', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      zgeqr2p(2, 1, A, 1, B, W, INFO );
      chkxer('ZGEQR2P', INFOT, NOUT, LERR, OK );

      // ZUNGQR

      SRNAMT = 'ZUNGQR';
      INFOT = 1;
      zungqr(-1, 0, 0, A, 1, X, W, 1, INFO );
      chkxer('ZUNGQR', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      zungqr(0, -1, 0, A, 1, X, W, 1, INFO );
      chkxer('ZUNGQR', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      zungqr(1, 2, 0, A, 1, X, W, 2, INFO );
      chkxer('ZUNGQR', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      zungqr(0, 0, -1, A, 1, X, W, 1, INFO );
      chkxer('ZUNGQR', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      zungqr(1, 1, 2, A, 1, X, W, 1, INFO );
      chkxer('ZUNGQR', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      zungqr(2, 2, 0, A, 1, X, W, 2, INFO );
      chkxer('ZUNGQR', INFOT, NOUT, LERR, OK );
      INFOT = 8;
      zungqr(2, 2, 0, A, 2, X, W, 1, INFO );
      chkxer('ZUNGQR', INFOT, NOUT, LERR, OK );

      // ZUNG2R

      SRNAMT = 'ZUNG2R';
      INFOT = 1;
      zung2r(-1, 0, 0, A, 1, X, W, INFO );
      chkxer('ZUNG2R', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      zung2r(0, -1, 0, A, 1, X, W, INFO );
      chkxer('ZUNG2R', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      zung2r(1, 2, 0, A, 1, X, W, INFO );
      chkxer('ZUNG2R', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      zung2r(0, 0, -1, A, 1, X, W, INFO );
      chkxer('ZUNG2R', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      zung2r(2, 1, 2, A, 2, X, W, INFO );
      chkxer('ZUNG2R', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      zung2r(2, 1, 0, A, 1, X, W, INFO );
      chkxer('ZUNG2R', INFOT, NOUT, LERR, OK );

      // ZUNMQR

      SRNAMT = 'ZUNMQR';
      INFOT = 1;
      zunmqr('/', 'N', 0, 0, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('ZUNMQR', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      zunmqr('L', '/', 0, 0, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('ZUNMQR', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      zunmqr('L', 'N', -1, 0, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('ZUNMQR', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      zunmqr('L', 'N', 0, -1, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('ZUNMQR', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      zunmqr('L', 'N', 0, 0, -1, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('ZUNMQR', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      zunmqr('L', 'N', 0, 1, 1, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('ZUNMQR', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      zunmqr('R', 'N', 1, 0, 1, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('ZUNMQR', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      zunmqr('L', 'N', 2, 1, 0, A, 1, X, AF, 2, W, 1, INFO );
      chkxer('ZUNMQR', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      zunmqr('R', 'N', 1, 2, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('ZUNMQR', INFOT, NOUT, LERR, OK );
      INFOT = 10;
      zunmqr('L', 'N', 2, 1, 0, A, 2, X, AF, 1, W, 1, INFO );
      chkxer('ZUNMQR', INFOT, NOUT, LERR, OK );
      INFOT = 12;
      zunmqr('L', 'N', 1, 2, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('ZUNMQR', INFOT, NOUT, LERR, OK );
      INFOT = 12;
      zunmqr('R', 'N', 2, 1, 0, A, 1, X, AF, 2, W, 1, INFO );
      chkxer('ZUNMQR', INFOT, NOUT, LERR, OK );

      // ZUNM2R

      SRNAMT = 'ZUNM2R';
      INFOT = 1;
      zunm2r('/', 'N', 0, 0, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('ZUNM2R', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      zunm2r('L', '/', 0, 0, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('ZUNM2R', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      zunm2r('L', 'N', -1, 0, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('ZUNM2R', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      zunm2r('L', 'N', 0, -1, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('ZUNM2R', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      zunm2r('L', 'N', 0, 0, -1, A, 1, X, AF, 1, W, INFO );
      chkxer('ZUNM2R', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      zunm2r('L', 'N', 0, 1, 1, A, 1, X, AF, 1, W, INFO );
      chkxer('ZUNM2R', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      zunm2r('R', 'N', 1, 0, 1, A, 1, X, AF, 1, W, INFO );
      chkxer('ZUNM2R', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      zunm2r('L', 'N', 2, 1, 0, A, 1, X, AF, 2, W, INFO );
      chkxer('ZUNM2R', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      zunm2r('R', 'N', 1, 2, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('ZUNM2R', INFOT, NOUT, LERR, OK );
      INFOT = 10;
      zunm2r('L', 'N', 2, 1, 0, A, 2, X, AF, 1, W, INFO );
      chkxer('ZUNM2R', INFOT, NOUT, LERR, OK );

      // Print a summary line.

      alaesm(PATH, OK, NOUT );

      return;

      // End of ZERRQR

      }
