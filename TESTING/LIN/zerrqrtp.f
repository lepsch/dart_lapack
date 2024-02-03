      SUBROUTINE ZERRQRTP( PATH, NUNIT )
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
      PARAMETER          ( NMAX = 2 )
      // ..
      // .. Local Scalars ..
      int                I, INFO, J;
      // ..
      // .. Local Arrays ..
      COMPLEX*16         A( NMAX, NMAX ), T( NMAX, NMAX ), W( NMAX ), B( NMAX, NMAX ), C( NMAX, NMAX )
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
      COMMON             / INFOC / INFOT, NOUT, OK, LERR
      COMMON             / SRNAMC / SRNAMT
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DCMPLX
      // ..
      // .. Executable Statements ..

      NOUT = NUNIT
      WRITE( NOUT, FMT = * )

      // Set the variables to innocuous values.

      DO J = 1, NMAX
         DO I = 1, NMAX
            A( I, J ) = 1.D0 / DCMPLX(DBLE( I+J ),0.D0)
            C( I, J ) = 1.D0 / DCMPLX(DBLE( I+J ),0.D0)
            T( I, J ) = 1.D0 / DCMPLX(DBLE( I+J ),0.D0)
         END DO
         W( J ) = DCMPLX(0.D0,0.D0)
      END DO
      OK = .TRUE.

      // Error exits for TPQRT factorization

      // ZTPQRT

      SRNAMT = 'ZTPQRT'
      INFOT = 1
      CALL ZTPQRT( -1, 1, 0, 1, A, 1, B, 1, T, 1, W, INFO )
      CALL CHKXER( 'ZTPQRT', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL ZTPQRT( 1, -1, 0, 1, A, 1, B, 1, T, 1, W, INFO )
      CALL CHKXER( 'ZTPQRT', INFOT, NOUT, LERR, OK )
      INFOT = 3
      CALL ZTPQRT( 0, 1, -1, 1, A, 1, B, 1, T, 1, W, INFO )
      CALL CHKXER( 'ZTPQRT', INFOT, NOUT, LERR, OK )
      INFOT = 3
      CALL ZTPQRT( 0, 1, 1, 1, A, 1, B, 1, T, 1, W, INFO )
      CALL CHKXER( 'ZTPQRT', INFOT, NOUT, LERR, OK )
      INFOT = 4
      CALL ZTPQRT( 0, 1, 0, 0, A, 1, B, 1, T, 1, W, INFO )
      CALL CHKXER( 'ZTPQRT', INFOT, NOUT, LERR, OK )
      INFOT = 4
      CALL ZTPQRT( 0, 1, 0, 2, A, 1, B, 1, T, 1, W, INFO )
      CALL CHKXER( 'ZTPQRT', INFOT, NOUT, LERR, OK )
      INFOT = 6
      CALL ZTPQRT( 1, 2, 0, 2, A, 1, B, 1, T, 1, W, INFO )
      CALL CHKXER( 'ZTPQRT', INFOT, NOUT, LERR, OK )
      INFOT = 8
      CALL ZTPQRT( 2, 1, 0, 1, A, 1, B, 1, T, 1, W, INFO )
      CALL CHKXER( 'ZTPQRT', INFOT, NOUT, LERR, OK )
      INFOT = 10
      CALL ZTPQRT( 2, 2, 1, 2, A, 2, B, 2, T, 1, W, INFO )
      CALL CHKXER( 'ZTPQRT', INFOT, NOUT, LERR, OK )

      // ZTPQRT2

      SRNAMT = 'ZTPQRT2'
      INFOT = 1
      CALL ZTPQRT2( -1, 0, 0, A, 1, B, 1, T, 1, INFO )
      CALL CHKXER( 'ZTPQRT2', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL ZTPQRT2( 0, -1, 0, A, 1, B, 1, T, 1, INFO )
      CALL CHKXER( 'ZTPQRT2', INFOT, NOUT, LERR, OK )
      INFOT = 3
      CALL ZTPQRT2( 0, 0, -1, A, 1, B, 1, T, 1, INFO )
      CALL CHKXER( 'ZTPQRT2', INFOT, NOUT, LERR, OK )
      INFOT = 5
      CALL ZTPQRT2( 2, 2, 0, A, 1, B, 2, T, 2, INFO )
      CALL CHKXER( 'ZTPQRT2', INFOT, NOUT, LERR, OK )
      INFOT = 7
      CALL ZTPQRT2( 2, 2, 0, A, 2, B, 1, T, 2, INFO )
      CALL CHKXER( 'ZTPQRT2', INFOT, NOUT, LERR, OK )
      INFOT = 9
      CALL ZTPQRT2( 2, 2, 0, A, 2, B, 2, T, 1, INFO )
      CALL CHKXER( 'ZTPQRT2', INFOT, NOUT, LERR, OK )

      // ZTPMQRT

      SRNAMT = 'ZTPMQRT'
      INFOT = 1
      CALL ZTPMQRT( '/', 'N', 0, 0, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO )
      CALL CHKXER( 'ZTPMQRT', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL ZTPMQRT( 'L', '/', 0, 0, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO )
      CALL CHKXER( 'ZTPMQRT', INFOT, NOUT, LERR, OK )
      INFOT = 3
      CALL ZTPMQRT( 'L', 'N', -1, 0, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO )
      CALL CHKXER( 'ZTPMQRT', INFOT, NOUT, LERR, OK )
      INFOT = 4
      CALL ZTPMQRT( 'L', 'N', 0, -1, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO )
      CALL CHKXER( 'ZTPMQRT', INFOT, NOUT, LERR, OK )
      INFOT = 5
      CALL ZTPMQRT( 'L', 'N', 0, 0, -1, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO )
      INFOT = 6
      CALL ZTPMQRT( 'L', 'N', 0, 0, 0, -1, 1, A, 1, T, 1, B, 1, C, 1, W, INFO )
      CALL CHKXER( 'ZTPMQRT', INFOT, NOUT, LERR, OK )
      INFOT = 7
      CALL ZTPMQRT( 'L', 'N', 0, 0, 0, 0, 0, A, 1, T, 1, B, 1, C, 1, W, INFO )
      CALL CHKXER( 'ZTPMQRT', INFOT, NOUT, LERR, OK )
      INFOT = 9
      CALL ZTPMQRT( 'R', 'N', 1, 2, 1, 1, 1, A, 1, T, 1, B, 1, C, 1, W, INFO )
      CALL CHKXER( 'ZTPMQRT', INFOT, NOUT, LERR, OK )
      INFOT = 9
      CALL ZTPMQRT( 'L', 'N', 2, 1, 1, 1, 1, A, 1, T, 1, B, 1, C, 1, W, INFO )
      CALL CHKXER( 'ZTPMQRT', INFOT, NOUT, LERR, OK )
      INFOT = 11
      CALL ZTPMQRT( 'R', 'N', 1, 1, 1, 1, 1, A, 1, T, 0, B, 1, C, 1, W, INFO )
      CALL CHKXER( 'ZTPMQRT', INFOT, NOUT, LERR, OK )
      INFOT = 13
      CALL ZTPMQRT( 'L', 'N', 1, 1, 1, 1, 1, A, 1, T, 1, B, 0, C, 1, W, INFO )
      CALL CHKXER( 'ZTPMQRT', INFOT, NOUT, LERR, OK )
      INFOT = 15
      CALL ZTPMQRT( 'L', 'N', 1, 1, 1, 1, 1, A, 1, T, 1, B, 1, C, 0, W, INFO )
      CALL CHKXER( 'ZTPMQRT', INFOT, NOUT, LERR, OK )

      // Print a summary line.

      CALL ALAESM( PATH, OK, NOUT )

      RETURN

      // End of ZERRQRTP

      END
