      void cerrqr(PATH, NUNIT ) {

// -- LAPACK test routine --
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
      COMPLEX            A( NMAX, NMAX ), AF( NMAX, NMAX ), B( NMAX ), W( NMAX ), X( NMAX );
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CGEQR2, CGEQR2P, CGEQRF, CGEQRFP, CHKXER, CUNG2R, CUNGQR, CUNM2R, CUNMQR
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
      // INTRINSIC CMPLX, REAL
      // ..
      // .. Executable Statements ..

      NOUT = NUNIT;
      WRITE( NOUT, FMT = * );

      // Set the variables to innocuous values.

      for (J = 1; J <= NMAX; J++) { // 20
         for (I = 1; I <= NMAX; I++) { // 10
            A[I, J] = CMPLX( 1. / REAL( I+J ), -1. / REAL( I+J ) );
            AF[I, J] = CMPLX( 1. / REAL( I+J ), -1. / REAL( I+J ) );
         } // 10
         B[J] = 0.;
         W[J] = 0.;
         X[J] = 0.;
      } // 20
      OK = true;

      // Error exits for QR factorization

      // CGEQRF

      SRNAMT = 'CGEQRF';
      INFOT = 1;
      cgeqrf(-1, 0, A, 1, B, W, 1, INFO );
      chkxer('CGEQRF', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      cgeqrf(0, -1, A, 1, B, W, 1, INFO );
      chkxer('CGEQRF', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      cgeqrf(2, 1, A, 1, B, W, 1, INFO );
      chkxer('CGEQRF', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      cgeqrf(1, 2, A, 1, B, W, 1, INFO );
      chkxer('CGEQRF', INFOT, NOUT, LERR, OK );

      // CGEQRFP

      SRNAMT = 'CGEQRFP';
      INFOT = 1;
      cgeqrfp(-1, 0, A, 1, B, W, 1, INFO );
      chkxer('CGEQRFP', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      cgeqrfp(0, -1, A, 1, B, W, 1, INFO );
      chkxer('CGEQRFP', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      cgeqrfp(2, 1, A, 1, B, W, 1, INFO );
      chkxer('CGEQRFP', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      cgeqrfp(1, 2, A, 1, B, W, 1, INFO );
      chkxer('CGEQRFP', INFOT, NOUT, LERR, OK );

      // CGEQR2

      SRNAMT = 'CGEQR2';
      INFOT = 1;
      cgeqr2(-1, 0, A, 1, B, W, INFO );
      chkxer('CGEQR2', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      cgeqr2(0, -1, A, 1, B, W, INFO );
      chkxer('CGEQR2', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      cgeqr2(2, 1, A, 1, B, W, INFO );
      chkxer('CGEQR2', INFOT, NOUT, LERR, OK );

      // CGEQR2P

      SRNAMT = 'CGEQR2P';
      INFOT = 1;
      cgeqr2p(-1, 0, A, 1, B, W, INFO );
      chkxer('CGEQR2P', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      cgeqr2p(0, -1, A, 1, B, W, INFO );
      chkxer('CGEQR2P', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      cgeqr2p(2, 1, A, 1, B, W, INFO );
      chkxer('CGEQR2P', INFOT, NOUT, LERR, OK );

      // CUNGQR

      SRNAMT = 'CUNGQR';
      INFOT = 1;
      cungqr(-1, 0, 0, A, 1, X, W, 1, INFO );
      chkxer('CUNGQR', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      cungqr(0, -1, 0, A, 1, X, W, 1, INFO );
      chkxer('CUNGQR', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      cungqr(1, 2, 0, A, 1, X, W, 2, INFO );
      chkxer('CUNGQR', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      cungqr(0, 0, -1, A, 1, X, W, 1, INFO );
      chkxer('CUNGQR', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      cungqr(1, 1, 2, A, 1, X, W, 1, INFO );
      chkxer('CUNGQR', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      cungqr(2, 2, 0, A, 1, X, W, 2, INFO );
      chkxer('CUNGQR', INFOT, NOUT, LERR, OK );
      INFOT = 8;
      cungqr(2, 2, 0, A, 2, X, W, 1, INFO );
      chkxer('CUNGQR', INFOT, NOUT, LERR, OK );

      // CUNG2R

      SRNAMT = 'CUNG2R';
      INFOT = 1;
      cung2r(-1, 0, 0, A, 1, X, W, INFO );
      chkxer('CUNG2R', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      cung2r(0, -1, 0, A, 1, X, W, INFO );
      chkxer('CUNG2R', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      cung2r(1, 2, 0, A, 1, X, W, INFO );
      chkxer('CUNG2R', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      cung2r(0, 0, -1, A, 1, X, W, INFO );
      chkxer('CUNG2R', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      cung2r(2, 1, 2, A, 2, X, W, INFO );
      chkxer('CUNG2R', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      cung2r(2, 1, 0, A, 1, X, W, INFO );
      chkxer('CUNG2R', INFOT, NOUT, LERR, OK );

      // CUNMQR

      SRNAMT = 'CUNMQR';
      INFOT = 1;
      cunmqr('/', 'N', 0, 0, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('CUNMQR', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      cunmqr('L', '/', 0, 0, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('CUNMQR', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      cunmqr('L', 'N', -1, 0, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('CUNMQR', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      cunmqr('L', 'N', 0, -1, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('CUNMQR', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      cunmqr('L', 'N', 0, 0, -1, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('CUNMQR', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      cunmqr('L', 'N', 0, 1, 1, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('CUNMQR', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      cunmqr('R', 'N', 1, 0, 1, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('CUNMQR', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      cunmqr('L', 'N', 2, 1, 0, A, 1, X, AF, 2, W, 1, INFO );
      chkxer('CUNMQR', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      cunmqr('R', 'N', 1, 2, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('CUNMQR', INFOT, NOUT, LERR, OK );
      INFOT = 10;
      cunmqr('L', 'N', 2, 1, 0, A, 2, X, AF, 1, W, 1, INFO );
      chkxer('CUNMQR', INFOT, NOUT, LERR, OK );
      INFOT = 12;
      cunmqr('L', 'N', 1, 2, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('CUNMQR', INFOT, NOUT, LERR, OK );
      INFOT = 12;
      cunmqr('R', 'N', 2, 1, 0, A, 1, X, AF, 2, W, 1, INFO );
      chkxer('CUNMQR', INFOT, NOUT, LERR, OK );

      // CUNM2R

      SRNAMT = 'CUNM2R';
      INFOT = 1;
      cunm2r('/', 'N', 0, 0, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('CUNM2R', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      cunm2r('L', '/', 0, 0, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('CUNM2R', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      cunm2r('L', 'N', -1, 0, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('CUNM2R', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      cunm2r('L', 'N', 0, -1, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('CUNM2R', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      cunm2r('L', 'N', 0, 0, -1, A, 1, X, AF, 1, W, INFO );
      chkxer('CUNM2R', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      cunm2r('L', 'N', 0, 1, 1, A, 1, X, AF, 1, W, INFO );
      chkxer('CUNM2R', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      cunm2r('R', 'N', 1, 0, 1, A, 1, X, AF, 1, W, INFO );
      chkxer('CUNM2R', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      cunm2r('L', 'N', 2, 1, 0, A, 1, X, AF, 2, W, INFO );
      chkxer('CUNM2R', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      cunm2r('R', 'N', 1, 2, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('CUNM2R', INFOT, NOUT, LERR, OK );
      INFOT = 10;
      cunm2r('L', 'N', 2, 1, 0, A, 2, X, AF, 1, W, INFO );
      chkxer('CUNM2R', INFOT, NOUT, LERR, OK );

      // Print a summary line.

      alaesm(PATH, OK, NOUT );

      return;
      }
