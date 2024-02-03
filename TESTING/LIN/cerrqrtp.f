      SUBROUTINE CERRQRTP( PATH, NUNIT )
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
      COMPLEX            A( NMAX, NMAX ), T( NMAX, NMAX ), W( NMAX ), B( NMAX, NMAX ), C( NMAX, NMAX )
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
      COMMON             / INFOC / INFOT, NOUT, OK, LERR
      COMMON             / SRNAMC / SRNAMT
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC FLOAT, CMPLX
      // ..
      // .. Executable Statements ..

      NOUT = NUNIT
      WRITE( NOUT, FMT = * )

      // Set the variables to innocuous values.

      DO J = 1, NMAX
         DO I = 1, NMAX
            A( I, J ) = 1.0 / CMPLX(FLOAT( I+J ),0.0)
            C( I, J ) = 1.0 / CMPLX(FLOAT( I+J ),0.0)
            T( I, J ) = 1.0 / CMPLX(FLOAT( I+J ),0.0)
         END DO
         W( J ) = CMPLX(0.0,0.0)
      END DO
      OK = .TRUE.

      // Error exits for TPQRT factorization

      // CTPQRT

      SRNAMT = 'CTPQRT'
      INFOT = 1
      CALL CTPQRT( -1, 1, 0, 1, A, 1, B, 1, T, 1, W, INFO )
      CALL CHKXER( 'CTPQRT', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL CTPQRT( 1, -1, 0, 1, A, 1, B, 1, T, 1, W, INFO )
      CALL CHKXER( 'CTPQRT', INFOT, NOUT, LERR, OK )
      INFOT = 3
      CALL CTPQRT( 0, 1, -1, 1, A, 1, B, 1, T, 1, W, INFO )
      CALL CHKXER( 'CTPQRT', INFOT, NOUT, LERR, OK )
      INFOT = 3
      CALL CTPQRT( 0, 1, 1, 1, A, 1, B, 1, T, 1, W, INFO )
      CALL CHKXER( 'CTPQRT', INFOT, NOUT, LERR, OK )
      INFOT = 4
      CALL CTPQRT( 0, 1, 0, 0, A, 1, B, 1, T, 1, W, INFO )
      CALL CHKXER( 'CTPQRT', INFOT, NOUT, LERR, OK )
      INFOT = 4
      CALL CTPQRT( 0, 1, 0, 2, A, 1, B, 1, T, 1, W, INFO )
      CALL CHKXER( 'CTPQRT', INFOT, NOUT, LERR, OK )
      INFOT = 6
      CALL CTPQRT( 1, 2, 0, 2, A, 1, B, 1, T, 1, W, INFO )
      CALL CHKXER( 'CTPQRT', INFOT, NOUT, LERR, OK )
      INFOT = 8
      CALL CTPQRT( 2, 1, 0, 1, A, 1, B, 1, T, 1, W, INFO )
      CALL CHKXER( 'CTPQRT', INFOT, NOUT, LERR, OK )
      INFOT = 10
      CALL CTPQRT( 2, 2, 1, 2, A, 2, B, 2, T, 1, W, INFO )
      CALL CHKXER( 'CTPQRT', INFOT, NOUT, LERR, OK )

      // CTPQRT2

      SRNAMT = 'CTPQRT2'
      INFOT = 1
      CALL CTPQRT2( -1, 0, 0, A, 1, B, 1, T, 1, INFO )
      CALL CHKXER( 'CTPQRT2', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL CTPQRT2( 0, -1, 0, A, 1, B, 1, T, 1, INFO )
      CALL CHKXER( 'CTPQRT2', INFOT, NOUT, LERR, OK )
      INFOT = 3
      CALL CTPQRT2( 0, 0, -1, A, 1, B, 1, T, 1, INFO )
      CALL CHKXER( 'CTPQRT2', INFOT, NOUT, LERR, OK )
      INFOT = 5
      CALL CTPQRT2( 2, 2, 0, A, 1, B, 2, T, 2, INFO )
      CALL CHKXER( 'CTPQRT2', INFOT, NOUT, LERR, OK )
      INFOT = 7
      CALL CTPQRT2( 2, 2, 0, A, 2, B, 1, T, 2, INFO )
      CALL CHKXER( 'CTPQRT2', INFOT, NOUT, LERR, OK )
      INFOT = 9
      CALL CTPQRT2( 2, 2, 0, A, 2, B, 2, T, 1, INFO )
      CALL CHKXER( 'CTPQRT2', INFOT, NOUT, LERR, OK )

      // CTPMQRT

      SRNAMT = 'CTPMQRT'
      INFOT = 1
      CALL CTPMQRT( '/', 'N', 0, 0, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO )
      CALL CHKXER( 'CTPMQRT', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL CTPMQRT( 'L', '/', 0, 0, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO )
      CALL CHKXER( 'CTPMQRT', INFOT, NOUT, LERR, OK )
      INFOT = 3
      CALL CTPMQRT( 'L', 'N', -1, 0, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO )
      CALL CHKXER( 'CTPMQRT', INFOT, NOUT, LERR, OK )
      INFOT = 4
      CALL CTPMQRT( 'L', 'N', 0, -1, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO )
      CALL CHKXER( 'CTPMQRT', INFOT, NOUT, LERR, OK )
      INFOT = 5
      CALL CTPMQRT( 'L', 'N', 0, 0, -1, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO )
      INFOT = 6
      CALL CTPMQRT( 'L', 'N', 0, 0, 0, -1, 1, A, 1, T, 1, B, 1, C, 1, W, INFO )
      CALL CHKXER( 'CTPMQRT', INFOT, NOUT, LERR, OK )
      INFOT = 7
      CALL CTPMQRT( 'L', 'N', 0, 0, 0, 0, 0, A, 1, T, 1, B, 1, C, 1, W, INFO )
      CALL CHKXER( 'CTPMQRT', INFOT, NOUT, LERR, OK )
      INFOT = 9
      CALL CTPMQRT( 'R', 'N', 1, 2, 1, 1, 1, A, 1, T, 1, B, 1, C, 1, W, INFO )
      CALL CHKXER( 'CTPMQRT', INFOT, NOUT, LERR, OK )
      INFOT = 9
      CALL CTPMQRT( 'L', 'N', 2, 1, 1, 1, 1, A, 1, T, 1, B, 1, C, 1, W, INFO )
      CALL CHKXER( 'CTPMQRT', INFOT, NOUT, LERR, OK )
      INFOT = 11
      CALL CTPMQRT( 'R', 'N', 1, 1, 1, 1, 1, A, 1, T, 0, B, 1, C, 1, W, INFO )
      CALL CHKXER( 'CTPMQRT', INFOT, NOUT, LERR, OK )
      INFOT = 13
      CALL CTPMQRT( 'L', 'N', 1, 1, 1, 1, 1, A, 1, T, 1, B, 0, C, 1, W, INFO )
      CALL CHKXER( 'CTPMQRT', INFOT, NOUT, LERR, OK )
      INFOT = 15
      CALL CTPMQRT( 'L', 'N', 1, 1, 1, 1, 1, A, 1, T, 1, B, 1, C, 0, W, INFO )
      CALL CHKXER( 'CTPMQRT', INFOT, NOUT, LERR, OK )

      // Print a summary line.

      CALL ALAESM( PATH, OK, NOUT )

      RETURN

      // End of CERRQRTP

      }
