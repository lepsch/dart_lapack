      void derrqrtp(PATH, NUNIT ) {
      // IMPLICIT NONE

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
      double             A( NMAX, NMAX ), T( NMAX, NMAX ), W( NMAX ), B( NMAX, NMAX ), C( NMAX, NMAX );
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, DTPQRT2, DTPQRT, DTPMQRT
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
      // INTRINSIC DBLE
      // ..
      // .. Executable Statements ..

      NOUT = NUNIT;
      WRITE( NOUT, FMT = * );

      // Set the variables to innocuous values.

      for (J = 1; J <= NMAX; J++) {
         for (I = 1; I <= NMAX; I++) {
            A[I, J] = 1.0 / DBLE( I+J );
            C[I, J] = 1.0 / DBLE( I+J );
            T[I, J] = 1.0 / DBLE( I+J );
         }
         W[J] = 0.0;
      }
      OK = true;

      // Error exits for TPQRT factorization

      // DTPQRT

      SRNAMT = 'DTPQRT';
      INFOT = 1;
      dtpqrt(-1, 1, 0, 1, A, 1, B, 1, T, 1, W, INFO );
      chkxer('DTPQRT', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      dtpqrt(1, -1, 0, 1, A, 1, B, 1, T, 1, W, INFO );
      chkxer('DTPQRT', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      dtpqrt(0, 1, -1, 1, A, 1, B, 1, T, 1, W, INFO );
      chkxer('DTPQRT', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      dtpqrt(0, 1, 1, 1, A, 1, B, 1, T, 1, W, INFO );
      chkxer('DTPQRT', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      dtpqrt(0, 1, 0, 0, A, 1, B, 1, T, 1, W, INFO );
      chkxer('DTPQRT', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      dtpqrt(0, 1, 0, 2, A, 1, B, 1, T, 1, W, INFO );
      chkxer('DTPQRT', INFOT, NOUT, LERR, OK );
      INFOT = 6;
      dtpqrt(1, 2, 0, 2, A, 1, B, 1, T, 1, W, INFO );
      chkxer('DTPQRT', INFOT, NOUT, LERR, OK );
      INFOT = 8;
      dtpqrt(2, 1, 0, 1, A, 1, B, 1, T, 1, W, INFO );
      chkxer('DTPQRT', INFOT, NOUT, LERR, OK );
      INFOT = 10;
      dtpqrt(2, 2, 1, 2, A, 2, B, 2, T, 1, W, INFO );
      chkxer('DTPQRT', INFOT, NOUT, LERR, OK );

      // DTPQRT2

      SRNAMT = 'DTPQRT2';
      INFOT = 1;
      dtpqrt2(-1, 0, 0, A, 1, B, 1, T, 1, INFO );
      chkxer('DTPQRT2', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      dtpqrt2(0, -1, 0, A, 1, B, 1, T, 1, INFO );
      chkxer('DTPQRT2', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      dtpqrt2(0, 0, -1, A, 1, B, 1, T, 1, INFO );
      chkxer('DTPQRT2', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      dtpqrt2(2, 2, 0, A, 1, B, 2, T, 2, INFO );
      chkxer('DTPQRT2', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      dtpqrt2(2, 2, 0, A, 2, B, 1, T, 2, INFO );
      chkxer('DTPQRT2', INFOT, NOUT, LERR, OK );
      INFOT = 9;
      dtpqrt2(2, 2, 0, A, 2, B, 2, T, 1, INFO );
      chkxer('DTPQRT2', INFOT, NOUT, LERR, OK );

      // DTPMQRT

      SRNAMT = 'DTPMQRT';
      INFOT = 1;
      dtpmqrt('/', 'N', 0, 0, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('DTPMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      dtpmqrt('L', '/', 0, 0, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('DTPMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      dtpmqrt('L', 'N', -1, 0, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('DTPMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      dtpmqrt('L', 'N', 0, -1, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('DTPMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      dtpmqrt('L', 'N', 0, 0, -1, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      INFOT = 6;
      dtpmqrt('L', 'N', 0, 0, 0, -1, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('DTPMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      dtpmqrt('L', 'N', 0, 0, 0, 0, 0, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('DTPMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 9;
      dtpmqrt('R', 'N', 1, 2, 1, 1, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('DTPMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 9;
      dtpmqrt('L', 'N', 2, 1, 1, 1, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('DTPMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 11;
      dtpmqrt('R', 'N', 1, 1, 1, 1, 1, A, 1, T, 0, B, 1, C, 1, W, INFO );
      chkxer('DTPMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 13;
      dtpmqrt('L', 'N', 1, 1, 1, 1, 1, A, 1, T, 1, B, 0, C, 1, W, INFO );
      chkxer('DTPMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 15;
      dtpmqrt('L', 'N', 1, 1, 1, 1, 1, A, 1, T, 1, B, 1, C, 0, W, INFO );
      chkxer('DTPMQRT', INFOT, NOUT, LERR, OK );

      // Print a summary line.

      alaesm(PATH, OK, NOUT );

      return;
      }
