      SUBROUTINE CERRQRTP( PATH, NUNIT );
      IMPLICIT NONE;

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
      COMPLEX            A( NMAX, NMAX ), T( NMAX, NMAX ), W( NMAX ), B( NMAX, NMAX ), C( NMAX, NMAX );
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, CTPQRT2, CTPQRT, CTPMQRT
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
      // INTRINSIC FLOAT, CMPLX
      // ..
      // .. Executable Statements ..

      NOUT = NUNIT;
      WRITE( NOUT, FMT = * );

      // Set the variables to innocuous values.

      for (J = 1; J <= NMAX; J++) {
         for (I = 1; I <= NMAX; I++) {
            A( I, J ) = 1.0 / CMPLX(FLOAT( I+J ),0.0);
            C( I, J ) = 1.0 / CMPLX(FLOAT( I+J ),0.0);
            T( I, J ) = 1.0 / CMPLX(FLOAT( I+J ),0.0);
         }
         W( J ) = CMPLX(0.0,0.0);
      }
      OK = true;

      // Error exits for TPQRT factorization

      // CTPQRT

      SRNAMT = 'CTPQRT';
      INFOT = 1;
      ctpqrt(-1, 1, 0, 1, A, 1, B, 1, T, 1, W, INFO );
      chkxer('CTPQRT', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      ctpqrt(1, -1, 0, 1, A, 1, B, 1, T, 1, W, INFO );
      chkxer('CTPQRT', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      ctpqrt(0, 1, -1, 1, A, 1, B, 1, T, 1, W, INFO );
      chkxer('CTPQRT', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      ctpqrt(0, 1, 1, 1, A, 1, B, 1, T, 1, W, INFO );
      chkxer('CTPQRT', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      ctpqrt(0, 1, 0, 0, A, 1, B, 1, T, 1, W, INFO );
      chkxer('CTPQRT', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      ctpqrt(0, 1, 0, 2, A, 1, B, 1, T, 1, W, INFO );
      chkxer('CTPQRT', INFOT, NOUT, LERR, OK );
      INFOT = 6;
      ctpqrt(1, 2, 0, 2, A, 1, B, 1, T, 1, W, INFO );
      chkxer('CTPQRT', INFOT, NOUT, LERR, OK );
      INFOT = 8;
      ctpqrt(2, 1, 0, 1, A, 1, B, 1, T, 1, W, INFO );
      chkxer('CTPQRT', INFOT, NOUT, LERR, OK );
      INFOT = 10;
      ctpqrt(2, 2, 1, 2, A, 2, B, 2, T, 1, W, INFO );
      chkxer('CTPQRT', INFOT, NOUT, LERR, OK );

      // CTPQRT2

      SRNAMT = 'CTPQRT2';
      INFOT = 1;
      ctpqrt2(-1, 0, 0, A, 1, B, 1, T, 1, INFO );
      chkxer('CTPQRT2', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      ctpqrt2(0, -1, 0, A, 1, B, 1, T, 1, INFO );
      chkxer('CTPQRT2', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      ctpqrt2(0, 0, -1, A, 1, B, 1, T, 1, INFO );
      chkxer('CTPQRT2', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      ctpqrt2(2, 2, 0, A, 1, B, 2, T, 2, INFO );
      chkxer('CTPQRT2', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      ctpqrt2(2, 2, 0, A, 2, B, 1, T, 2, INFO );
      chkxer('CTPQRT2', INFOT, NOUT, LERR, OK );
      INFOT = 9;
      ctpqrt2(2, 2, 0, A, 2, B, 2, T, 1, INFO );
      chkxer('CTPQRT2', INFOT, NOUT, LERR, OK );

      // CTPMQRT

      SRNAMT = 'CTPMQRT';
      INFOT = 1;
      ctpmqrt('/', 'N', 0, 0, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('CTPMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      ctpmqrt('L', '/', 0, 0, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('CTPMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      ctpmqrt('L', 'N', -1, 0, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('CTPMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      ctpmqrt('L', 'N', 0, -1, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('CTPMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      ctpmqrt('L', 'N', 0, 0, -1, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      INFOT = 6;
      ctpmqrt('L', 'N', 0, 0, 0, -1, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('CTPMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      ctpmqrt('L', 'N', 0, 0, 0, 0, 0, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('CTPMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 9;
      ctpmqrt('R', 'N', 1, 2, 1, 1, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('CTPMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 9;
      ctpmqrt('L', 'N', 2, 1, 1, 1, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('CTPMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 11;
      ctpmqrt('R', 'N', 1, 1, 1, 1, 1, A, 1, T, 0, B, 1, C, 1, W, INFO );
      chkxer('CTPMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 13;
      ctpmqrt('L', 'N', 1, 1, 1, 1, 1, A, 1, T, 1, B, 0, C, 1, W, INFO );
      chkxer('CTPMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 15;
      ctpmqrt('L', 'N', 1, 1, 1, 1, 1, A, 1, T, 1, B, 1, C, 0, W, INFO );
      chkxer('CTPMQRT', INFOT, NOUT, LERR, OK );

      // Print a summary line.

      alaesm(PATH, OK, NOUT );

      RETURN;

      // End of CERRQRTP

      }
