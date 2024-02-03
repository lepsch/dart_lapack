      SUBROUTINE ZERRQRTP( PATH, NUNIT );
      // IMPLICIT NONE

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             PATH;
      int                NUNIT;
      // ..

*  =====================================================================

      // .. Parameters ..
      int                NMAX;
      const              NMAX = 2 ;
      // ..
      // .. Local Scalars ..
      int                I, INFO, J;
      // ..
      // .. Local Arrays ..
      COMPLEX*16         A( NMAX, NMAX ), T( NMAX, NMAX ), W( NMAX ), B( NMAX, NMAX ), C( NMAX, NMAX );
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, ZTPQRT2, ZTPQRT, ZTPMQRT
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
            A( I, J ) = 1.0 / DCMPLX(DBLE( I+J ),0.0);
            C( I, J ) = 1.0 / DCMPLX(DBLE( I+J ),0.0);
            T( I, J ) = 1.0 / DCMPLX(DBLE( I+J ),0.0);
         }
         W( J ) = DCMPLX(0.0,0.0);
      }
      OK = true;

      // Error exits for TPQRT factorization

      // ZTPQRT

      SRNAMT = 'ZTPQRT';
      INFOT = 1;
      ztpqrt(-1, 1, 0, 1, A, 1, B, 1, T, 1, W, INFO );
      chkxer('ZTPQRT', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      ztpqrt(1, -1, 0, 1, A, 1, B, 1, T, 1, W, INFO );
      chkxer('ZTPQRT', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      ztpqrt(0, 1, -1, 1, A, 1, B, 1, T, 1, W, INFO );
      chkxer('ZTPQRT', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      ztpqrt(0, 1, 1, 1, A, 1, B, 1, T, 1, W, INFO );
      chkxer('ZTPQRT', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      ztpqrt(0, 1, 0, 0, A, 1, B, 1, T, 1, W, INFO );
      chkxer('ZTPQRT', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      ztpqrt(0, 1, 0, 2, A, 1, B, 1, T, 1, W, INFO );
      chkxer('ZTPQRT', INFOT, NOUT, LERR, OK );
      INFOT = 6;
      ztpqrt(1, 2, 0, 2, A, 1, B, 1, T, 1, W, INFO );
      chkxer('ZTPQRT', INFOT, NOUT, LERR, OK );
      INFOT = 8;
      ztpqrt(2, 1, 0, 1, A, 1, B, 1, T, 1, W, INFO );
      chkxer('ZTPQRT', INFOT, NOUT, LERR, OK );
      INFOT = 10;
      ztpqrt(2, 2, 1, 2, A, 2, B, 2, T, 1, W, INFO );
      chkxer('ZTPQRT', INFOT, NOUT, LERR, OK );

      // ZTPQRT2

      SRNAMT = 'ZTPQRT2';
      INFOT = 1;
      ztpqrt2(-1, 0, 0, A, 1, B, 1, T, 1, INFO );
      chkxer('ZTPQRT2', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      ztpqrt2(0, -1, 0, A, 1, B, 1, T, 1, INFO );
      chkxer('ZTPQRT2', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      ztpqrt2(0, 0, -1, A, 1, B, 1, T, 1, INFO );
      chkxer('ZTPQRT2', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      ztpqrt2(2, 2, 0, A, 1, B, 2, T, 2, INFO );
      chkxer('ZTPQRT2', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      ztpqrt2(2, 2, 0, A, 2, B, 1, T, 2, INFO );
      chkxer('ZTPQRT2', INFOT, NOUT, LERR, OK );
      INFOT = 9;
      ztpqrt2(2, 2, 0, A, 2, B, 2, T, 1, INFO );
      chkxer('ZTPQRT2', INFOT, NOUT, LERR, OK );

      // ZTPMQRT

      SRNAMT = 'ZTPMQRT';
      INFOT = 1;
      ztpmqrt('/', 'N', 0, 0, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('ZTPMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      ztpmqrt('L', '/', 0, 0, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('ZTPMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      ztpmqrt('L', 'N', -1, 0, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('ZTPMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      ztpmqrt('L', 'N', 0, -1, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('ZTPMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      ztpmqrt('L', 'N', 0, 0, -1, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      INFOT = 6;
      ztpmqrt('L', 'N', 0, 0, 0, -1, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('ZTPMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      ztpmqrt('L', 'N', 0, 0, 0, 0, 0, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('ZTPMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 9;
      ztpmqrt('R', 'N', 1, 2, 1, 1, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('ZTPMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 9;
      ztpmqrt('L', 'N', 2, 1, 1, 1, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('ZTPMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 11;
      ztpmqrt('R', 'N', 1, 1, 1, 1, 1, A, 1, T, 0, B, 1, C, 1, W, INFO );
      chkxer('ZTPMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 13;
      ztpmqrt('L', 'N', 1, 1, 1, 1, 1, A, 1, T, 1, B, 0, C, 1, W, INFO );
      chkxer('ZTPMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 15;
      ztpmqrt('L', 'N', 1, 1, 1, 1, 1, A, 1, T, 1, B, 1, C, 0, W, INFO );
      chkxer('ZTPMQRT', INFOT, NOUT, LERR, OK );

      // Print a summary line.

      alaesm(PATH, OK, NOUT );

      RETURN;

      // End of ZERRQRTP

      }
