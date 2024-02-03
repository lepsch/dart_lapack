      SUBROUTINE DERRQRT( PATH, NUNIT )
      IMPLICIT NONE

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
      double             A( NMAX, NMAX ), T( NMAX, NMAX ), W( NMAX ), C( NMAX, NMAX );
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, DGEQRT2, DGEQRT3, DGEQRT, DGEMQRT
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

      NOUT = NUNIT
      WRITE( NOUT, FMT = * )

      // Set the variables to innocuous values.

      for (J = 1; J <= NMAX; J++) {
         for (I = 1; I <= NMAX; I++) {
            A( I, J ) = 1.D0 / DBLE( I+J )
            C( I, J ) = 1.D0 / DBLE( I+J )
            T( I, J ) = 1.D0 / DBLE( I+J )
         }
         W( J ) = 0.D0
      }
      OK = .TRUE.

      // Error exits for QRT factorization

      // DGEQRT

      SRNAMT = 'DGEQRT'
      INFOT = 1
      dgeqrt(-1, 0, 1, A, 1, T, 1, W, INFO );
      chkxer('DGEQRT', INFOT, NOUT, LERR, OK );
      INFOT = 2
      dgeqrt(0, -1, 1, A, 1, T, 1, W, INFO );
      chkxer('DGEQRT', INFOT, NOUT, LERR, OK );
      INFOT = 3
      dgeqrt(0, 0, 0, A, 1, T, 1, W, INFO );
      chkxer('DGEQRT', INFOT, NOUT, LERR, OK );
      INFOT = 5
      dgeqrt(2, 1, 1, A, 1, T, 1, W, INFO );
      chkxer('DGEQRT', INFOT, NOUT, LERR, OK );
      INFOT = 7
      dgeqrt(2, 2, 2, A, 2, T, 1, W, INFO );
      chkxer('DGEQRT', INFOT, NOUT, LERR, OK );

      // DGEQRT2

      SRNAMT = 'DGEQRT2'
      INFOT = 1
      dgeqrt2(-1, 0, A, 1, T, 1, INFO );
      chkxer('DGEQRT2', INFOT, NOUT, LERR, OK );
      INFOT = 2
      dgeqrt2(0, -1, A, 1, T, 1, INFO );
      chkxer('DGEQRT2', INFOT, NOUT, LERR, OK );
      INFOT = 4
      dgeqrt2(2, 1, A, 1, T, 1, INFO );
      chkxer('DGEQRT2', INFOT, NOUT, LERR, OK );
      INFOT = 6
      dgeqrt2(2, 2, A, 2, T, 1, INFO );
      chkxer('DGEQRT2', INFOT, NOUT, LERR, OK );

      // DGEQRT3

      SRNAMT = 'DGEQRT3'
      INFOT = 1
      dgeqrt3(-1, 0, A, 1, T, 1, INFO );
      chkxer('DGEQRT3', INFOT, NOUT, LERR, OK );
      INFOT = 2
      dgeqrt3(0, -1, A, 1, T, 1, INFO );
      chkxer('DGEQRT3', INFOT, NOUT, LERR, OK );
      INFOT = 4
      dgeqrt3(2, 1, A, 1, T, 1, INFO );
      chkxer('DGEQRT3', INFOT, NOUT, LERR, OK );
      INFOT = 6
      dgeqrt3(2, 2, A, 2, T, 1, INFO );
      chkxer('DGEQRT3', INFOT, NOUT, LERR, OK );

      // DGEMQRT

      SRNAMT = 'DGEMQRT'
      INFOT = 1
      dgemqrt('/', 'N', 0, 0, 0, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('DGEMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 2
      dgemqrt('L', '/', 0, 0, 0, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('DGEMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 3
      dgemqrt('L', 'N', -1, 0, 0, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('DGEMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 4
      dgemqrt('L', 'N', 0, -1, 0, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('DGEMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 5
      dgemqrt('L', 'N', 0, 0, -1, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('DGEMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 5
      dgemqrt('R', 'N', 0, 0, -1, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('DGEMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 6
      dgemqrt('L', 'N', 0, 0, 0, 0, A, 1, T, 1, C, 1, W, INFO );
      chkxer('DGEMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 8
      dgemqrt('R', 'N', 1, 2, 1, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('DGEMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 8
      dgemqrt('L', 'N', 2, 1, 1, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('DGEMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 10
      dgemqrt('R', 'N', 1, 1, 1, 1, A, 1, T, 0, C, 1, W, INFO );
      chkxer('DGEMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 12
      dgemqrt('L', 'N', 1, 1, 1, 1, A, 1, T, 1, C, 0, W, INFO );
      chkxer('DGEMQRT', INFOT, NOUT, LERR, OK );

      // Print a summary line.

      alaesm(PATH, OK, NOUT );

      RETURN

      // End of DERRQRT

      }
