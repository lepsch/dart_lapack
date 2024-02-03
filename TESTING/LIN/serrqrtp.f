      SUBROUTINE SERRQRTP( PATH, NUNIT );
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
      REAL               A( NMAX, NMAX ), T( NMAX, NMAX ), W( NMAX ), B( NMAX, NMAX ), C( NMAX, NMAX );
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, STPQRT2, STPQRT, STPMQRT
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
      // INTRINSIC FLOAT
      // ..
      // .. Executable Statements ..

      NOUT = NUNIT;
      WRITE( NOUT, FMT = * );

      // Set the variables to innocuous values.

      for (J = 1; J <= NMAX; J++) {
         for (I = 1; I <= NMAX; I++) {
            A( I, J ) = 1.0 / FLOAT( I+J );
            C( I, J ) = 1.0 / FLOAT( I+J );
            T( I, J ) = 1.0 / FLOAT( I+J );
         }
         W( J ) = 0.0;
      }
      OK = true;

      // Error exits for TPQRT factorization

      // STPQRT

      SRNAMT = 'STPQRT';
      INFOT = 1;
      stpqrt(-1, 1, 0, 1, A, 1, B, 1, T, 1, W, INFO );
      chkxer('STPQRT', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      stpqrt(1, -1, 0, 1, A, 1, B, 1, T, 1, W, INFO );
      chkxer('STPQRT', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      stpqrt(0, 1, -1, 1, A, 1, B, 1, T, 1, W, INFO );
      chkxer('STPQRT', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      stpqrt(0, 1, 1, 1, A, 1, B, 1, T, 1, W, INFO );
      chkxer('STPQRT', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      stpqrt(0, 1, 0, 0, A, 1, B, 1, T, 1, W, INFO );
      chkxer('STPQRT', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      stpqrt(0, 1, 0, 2, A, 1, B, 1, T, 1, W, INFO );
      chkxer('STPQRT', INFOT, NOUT, LERR, OK );
      INFOT = 6;
      stpqrt(1, 2, 0, 2, A, 1, B, 1, T, 1, W, INFO );
      chkxer('STPQRT', INFOT, NOUT, LERR, OK );
      INFOT = 8;
      stpqrt(2, 1, 0, 1, A, 1, B, 1, T, 1, W, INFO );
      chkxer('STPQRT', INFOT, NOUT, LERR, OK );
      INFOT = 10;
      stpqrt(2, 2, 1, 2, A, 2, B, 2, T, 1, W, INFO );
      chkxer('STPQRT', INFOT, NOUT, LERR, OK );

      // STPQRT2

      SRNAMT = 'STPQRT2';
      INFOT = 1;
      stpqrt2(-1, 0, 0, A, 1, B, 1, T, 1, INFO );
      chkxer('STPQRT2', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      stpqrt2(0, -1, 0, A, 1, B, 1, T, 1, INFO );
      chkxer('STPQRT2', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      stpqrt2(0, 0, -1, A, 1, B, 1, T, 1, INFO );
      chkxer('STPQRT2', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      stpqrt2(2, 2, 0, A, 1, B, 2, T, 2, INFO );
      chkxer('STPQRT2', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      stpqrt2(2, 2, 0, A, 2, B, 1, T, 2, INFO );
      chkxer('STPQRT2', INFOT, NOUT, LERR, OK );
      INFOT = 9;
      stpqrt2(2, 2, 0, A, 2, B, 2, T, 1, INFO );
      chkxer('STPQRT2', INFOT, NOUT, LERR, OK );

      // STPMQRT

      SRNAMT = 'STPMQRT';
      INFOT = 1;
      stpmqrt('/', 'N', 0, 0, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('STPMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      stpmqrt('L', '/', 0, 0, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('STPMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      stpmqrt('L', 'N', -1, 0, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('STPMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      stpmqrt('L', 'N', 0, -1, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('STPMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      stpmqrt('L', 'N', 0, 0, -1, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      INFOT = 6;
      stpmqrt('L', 'N', 0, 0, 0, -1, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('STPMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      stpmqrt('L', 'N', 0, 0, 0, 0, 0, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('STPMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 9;
      stpmqrt('R', 'N', 1, 2, 1, 1, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('STPMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 9;
      stpmqrt('L', 'N', 2, 1, 1, 1, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('STPMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 11;
      stpmqrt('R', 'N', 1, 1, 1, 1, 1, A, 1, T, 0, B, 1, C, 1, W, INFO );
      chkxer('STPMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 13;
      stpmqrt('L', 'N', 1, 1, 1, 1, 1, A, 1, T, 1, B, 0, C, 1, W, INFO );
      chkxer('STPMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 15;
      stpmqrt('L', 'N', 1, 1, 1, 1, 1, A, 1, T, 1, B, 1, C, 0, W, INFO );
      chkxer('STPMQRT', INFOT, NOUT, LERR, OK );

      // Print a summary line.

      alaesm(PATH, OK, NOUT );

      RETURN;

      // End of SERRQRTP

      }
