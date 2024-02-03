      void zerrqrt(PATH, NUNIT ) {
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
      COMPLEX*16         A( NMAX, NMAX ), T( NMAX, NMAX ), W( NMAX ), C( NMAX, NMAX );
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, ZGEQRT2, ZGEQRT3, ZGEQRT, ZGEMQRT
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

      for (J = 1; J <= NMAX; J++) {
         for (I = 1; I <= NMAX; I++) {
            A( I, J ) = 1.0 / DCMPLX( DBLE( I+J ), 0.0 );
            C( I, J ) = 1.0 / DCMPLX( DBLE( I+J ), 0.0 );
            T( I, J ) = 1.0 / DCMPLX( DBLE( I+J ), 0.0 );
         }
         W( J ) = 0.0;
      }
      OK = true;

      // Error exits for QRT factorization

      // ZGEQRT

      SRNAMT = 'ZGEQRT';
      INFOT = 1;
      zgeqrt(-1, 0, 1, A, 1, T, 1, W, INFO );
      chkxer('ZGEQRT', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      zgeqrt(0, -1, 1, A, 1, T, 1, W, INFO );
      chkxer('ZGEQRT', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      zgeqrt(0, 0, 0, A, 1, T, 1, W, INFO );
      chkxer('ZGEQRT', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      zgeqrt(2, 1, 1, A, 1, T, 1, W, INFO );
      chkxer('ZGEQRT', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      zgeqrt(2, 2, 2, A, 2, T, 1, W, INFO );
      chkxer('ZGEQRT', INFOT, NOUT, LERR, OK );

      // ZGEQRT2

      SRNAMT = 'ZGEQRT2';
      INFOT = 1;
      zgeqrt2(-1, 0, A, 1, T, 1, INFO );
      chkxer('ZGEQRT2', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      zgeqrt2(0, -1, A, 1, T, 1, INFO );
      chkxer('ZGEQRT2', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      zgeqrt2(2, 1, A, 1, T, 1, INFO );
      chkxer('ZGEQRT2', INFOT, NOUT, LERR, OK );
      INFOT = 6;
      zgeqrt2(2, 2, A, 2, T, 1, INFO );
      chkxer('ZGEQRT2', INFOT, NOUT, LERR, OK );

      // ZGEQRT3

      SRNAMT = 'ZGEQRT3';
      INFOT = 1;
      zgeqrt3(-1, 0, A, 1, T, 1, INFO );
      chkxer('ZGEQRT3', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      zgeqrt3(0, -1, A, 1, T, 1, INFO );
      chkxer('ZGEQRT3', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      zgeqrt3(2, 1, A, 1, T, 1, INFO );
      chkxer('ZGEQRT3', INFOT, NOUT, LERR, OK );
      INFOT = 6;
      zgeqrt3(2, 2, A, 2, T, 1, INFO );
      chkxer('ZGEQRT3', INFOT, NOUT, LERR, OK );

      // ZGEMQRT

      SRNAMT = 'ZGEMQRT';
      INFOT = 1;
      zgemqrt('/', 'N', 0, 0, 0, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('ZGEMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      zgemqrt('L', '/', 0, 0, 0, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('ZGEMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      zgemqrt('L', 'N', -1, 0, 0, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('ZGEMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      zgemqrt('L', 'N', 0, -1, 0, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('ZGEMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      zgemqrt('L', 'N', 0, 0, -1, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('ZGEMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      zgemqrt('R', 'N', 0, 0, -1, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('ZGEMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 6;
      zgemqrt('L', 'N', 0, 0, 0, 0, A, 1, T, 1, C, 1, W, INFO );
      chkxer('ZGEMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 8;
      zgemqrt('R', 'N', 1, 2, 1, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('ZGEMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 8;
      zgemqrt('L', 'N', 2, 1, 1, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('ZGEMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 10;
      zgemqrt('R', 'N', 1, 1, 1, 1, A, 1, T, 0, C, 1, W, INFO );
      chkxer('ZGEMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 12;
      zgemqrt('L', 'N', 1, 1, 1, 1, A, 1, T, 1, C, 0, W, INFO );
      chkxer('ZGEMQRT', INFOT, NOUT, LERR, OK );

      // Print a summary line.

      alaesm(PATH, OK, NOUT );

      return;

      // End of ZERRQRT

      }
