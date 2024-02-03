      SUBROUTINE ZERRTSQR( PATH, NUNIT )
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
      int                I, INFO, J, MB, NB;
      // ..
      // .. Local Arrays ..
      COMPLEX*16         A( NMAX, NMAX ), T( NMAX, NMAX ), W( NMAX ), C( NMAX, NMAX ), TAU(NMAX)
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, ZGEQR, ZGEMQR, ZGELQ, ZGEMLQ
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
      // INTRINSIC DBLE
      // ..
      // .. Executable Statements ..

      NOUT = NUNIT
      WRITE( NOUT, FMT = * )

      // Set the variables to innocuous values.

      DO J = 1, NMAX
         DO I = 1, NMAX
            A( I, J ) = 1.D0 / DBLE( I+J )
            C( I, J ) = 1.D0 / DBLE( I+J )
            T( I, J ) = 1.D0 / DBLE( I+J )
         END DO
         W( J ) = 0.D0
      END DO
      OK = .TRUE.

      // Error exits for TS factorization

      // ZGEQR

      SRNAMT = 'ZGEQR'
      INFOT = 1
      CALL ZGEQR( -1, 0, A, 1, TAU, 1, W, 1, INFO )
      CALL CHKXER( 'ZGEQR', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL ZGEQR( 0, -1, A, 1, TAU, 1, W, 1, INFO )
      CALL CHKXER( 'ZGEQR', INFOT, NOUT, LERR, OK )
      INFOT = 4
      CALL ZGEQR( 1, 1, A, 0, TAU, 1, W, 1, INFO )
      CALL CHKXER( 'ZGEQR', INFOT, NOUT, LERR, OK )
      INFOT = 6
      CALL ZGEQR( 3, 2, A, 3, TAU, 1, W, 1, INFO )
      CALL CHKXER( 'ZGEQR', INFOT, NOUT, LERR, OK )
      INFOT = 8
      CALL ZGEQR( 3, 2, A, 3, TAU, 8, W, 0, INFO )
      CALL CHKXER( 'ZGEQR', INFOT, NOUT, LERR, OK )

      // ZLATSQR

      MB = 1
      NB = 1
      SRNAMT = 'ZLATSQR'
      INFOT = 1
      CALL ZLATSQR( -1, 0, MB, NB, A, 1, TAU, 1, W, 1, INFO )
      CALL CHKXER( 'ZLATSQR', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL ZLATSQR( 1, 2, MB, NB, A, 1, TAU, 1, W, 1, INFO )
      CALL CHKXER( 'ZLATSQR', INFOT, NOUT, LERR, OK )
      CALL ZLATSQR( 0, -1, MB, NB, A, 1, TAU, 1, W, 1, INFO )
      CALL CHKXER( 'ZLATSQR', INFOT, NOUT, LERR, OK )
      INFOT = 3
      CALL ZLATSQR( 2, 1, -1, NB, A, 2, TAU, 1, W, 1, INFO )
      CALL CHKXER( 'ZLATSQR', INFOT, NOUT, LERR, OK )
      INFOT = 4
      CALL ZLATSQR( 2, 1, MB, 2, A, 2, TAU, 1, W, 1, INFO )
      CALL CHKXER( 'ZLATSQR', INFOT, NOUT, LERR, OK )
      INFOT = 6
      CALL ZLATSQR( 2, 1, MB, NB, A, 1, TAU, 1, W, 1, INFO )
      CALL CHKXER( 'ZLATSQR', INFOT, NOUT, LERR, OK )
      INFOT = 8
      CALL ZLATSQR( 2, 1, MB, NB, A, 2, TAU, 0, W, 1, INFO )
      CALL CHKXER( 'ZLATSQR', INFOT, NOUT, LERR, OK )
      INFOT = 10
      CALL ZLATSQR( 2, 1, MB, NB, A, 2, TAU, 2, W, 0, INFO )
      CALL CHKXER( 'ZLATSQR', INFOT, NOUT, LERR, OK )

      // ZGEMQR

      TAU(1)=1
      TAU(2)=1
      SRNAMT = 'ZGEMQR'
      NB=1
      INFOT = 1
      CALL ZGEMQR( '/', 'N', 0, 0, 0, A, 1, TAU, 1, C, 1, W, 1,INFO)
      CALL CHKXER( 'ZGEMQR', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL ZGEMQR( 'L', '/', 0, 0, 0, A, 1, TAU, 1, C, 1, W, 1,INFO)
      CALL CHKXER( 'ZGEMQR', INFOT, NOUT, LERR, OK )
      INFOT = 3
      CALL ZGEMQR( 'L', 'N', -1, 0, 0, A, 1, TAU, 1, C, 1, W,1,INFO)
      CALL CHKXER( 'ZGEMQR', INFOT, NOUT, LERR, OK )
      INFOT = 4
      CALL ZGEMQR( 'L', 'N', 0, -1, 0, A, 1, TAU, 1, C, 1, W,1,INFO)
      CALL CHKXER( 'ZGEMQR', INFOT, NOUT, LERR, OK )
      INFOT = 5
      CALL ZGEMQR( 'L', 'N', 0, 0, -1, A, 1, TAU, 1, C, 1, W,1,INFO)
      CALL CHKXER( 'ZGEMQR', INFOT, NOUT, LERR, OK )
      INFOT = 5
      CALL ZGEMQR( 'R', 'N', 0, 0, -1, A, 1, TAU, 1, C, 1, W,1,INFO)
      CALL CHKXER( 'ZGEMQR', INFOT, NOUT, LERR, OK )
      INFOT = 7
      CALL ZGEMQR( 'L', 'N', 2, 1, 0, A, 0, TAU, 1, C, 1, W, 1,INFO)
      CALL CHKXER( 'ZGEMQR', INFOT, NOUT, LERR, OK )
      INFOT = 9
      CALL ZGEMQR( 'R', 'N', 2, 2, 1, A, 2, TAU, 0, C, 1, W, 1,INFO)
      CALL CHKXER( 'ZGEMQR', INFOT, NOUT, LERR, OK )
      INFOT = 9
      CALL ZGEMQR( 'L', 'N', 2, 2, 1, A, 2, TAU, 0, C, 1, W, 1,INFO)
      CALL CHKXER( 'ZGEMQR', INFOT, NOUT, LERR, OK )
      INFOT = 11
      CALL ZGEMQR( 'L', 'N', 2, 1, 1, A, 2, TAU, 6, C, 0, W, 1,INFO)
      CALL CHKXER( 'ZGEMQR', INFOT, NOUT, LERR, OK )
      INFOT = 13
      CALL ZGEMQR( 'L', 'N', 2, 2, 1, A, 2, TAU, 6, C, 2, W, 0,INFO)
      CALL CHKXER( 'ZGEMQR', INFOT, NOUT, LERR, OK )

      // ZGELQ

      SRNAMT = 'ZGELQ'
      INFOT = 1
      CALL ZGELQ( -1, 0, A, 1, TAU, 1, W, 1, INFO )
      CALL CHKXER( 'ZGELQ', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL ZGELQ( 0, -1, A, 1, TAU, 1, W, 1, INFO )
      CALL CHKXER( 'ZGELQ', INFOT, NOUT, LERR, OK )
      INFOT = 4
      CALL ZGELQ( 1, 1, A, 0, TAU, 1, W, 1, INFO )
      CALL CHKXER( 'ZGELQ', INFOT, NOUT, LERR, OK )
      INFOT = 6
      CALL ZGELQ( 2, 3, A, 3, TAU, 1, W, 1, INFO )
      CALL CHKXER( 'ZGELQ', INFOT, NOUT, LERR, OK )
      INFOT = 8
      CALL ZGELQ( 2, 3, A, 3, TAU, 8, W, 0, INFO )
      CALL CHKXER( 'ZGELQ', INFOT, NOUT, LERR, OK )

      // ZLASWLQ

      MB = 1
      NB = 1
      SRNAMT = 'ZLASWLQ'
      INFOT = 1
      CALL ZLASWLQ( -1, 0, MB, NB, A, 1, TAU, 1, W, 1, INFO )
      CALL CHKXER( 'ZLASWLQ', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL ZLASWLQ( 2, 1, MB, NB, A, 1, TAU, 1, W, 1, INFO )
      CALL CHKXER( 'ZLASWLQ', INFOT, NOUT, LERR, OK )
      CALL ZLASWLQ( 0, -1, MB, NB, A, 1, TAU, 1, W, 1, INFO )
      CALL CHKXER( 'ZLASWLQ', INFOT, NOUT, LERR, OK )
      INFOT = 3
      CALL ZLASWLQ( 1, 2, -1, NB, A, 1, TAU, 1, W, 1, INFO )
      CALL CHKXER( 'ZLASWLQ', INFOT, NOUT, LERR, OK )
      CALL ZLASWLQ( 1, 1, 2, NB, A, 1, TAU, 1, W, 1, INFO )
      CALL CHKXER( 'ZLASWLQ', INFOT, NOUT, LERR, OK )
      INFOT = 4
      CALL ZLASWLQ( 1, 2, MB, -1, A, 1, TAU, 1, W, 1, INFO )
      CALL CHKXER( 'ZLASWLQ', INFOT, NOUT, LERR, OK )
      INFOT = 6
      CALL ZLASWLQ( 1, 2, MB, NB, A, 0, TAU, 1, W, 1, INFO )
      CALL CHKXER( 'ZLASWLQ', INFOT, NOUT, LERR, OK )
      INFOT = 8
      CALL ZLASWLQ( 1, 2, MB, NB, A, 1, TAU, 0, W, 1, INFO )
      CALL CHKXER( 'ZLASWLQ', INFOT, NOUT, LERR, OK )
      INFOT = 10
      CALL ZLASWLQ( 1, 2, MB, NB, A, 1, TAU, 1, W, 0, INFO )
      CALL CHKXER( 'ZLASWLQ', INFOT, NOUT, LERR, OK )

      // ZGEMLQ

      TAU(1)=1
      TAU(2)=1
      SRNAMT = 'ZGEMLQ'
      NB=1
      INFOT = 1
      CALL ZGEMLQ( '/', 'N', 0, 0, 0, A, 1, TAU, 1, C, 1, W, 1,INFO)
      CALL CHKXER( 'ZGEMLQ', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL ZGEMLQ( 'L', '/', 0, 0, 0, A, 1, TAU, 1, C, 1, W, 1,INFO)
      CALL CHKXER( 'ZGEMLQ', INFOT, NOUT, LERR, OK )
      INFOT = 3
      CALL ZGEMLQ( 'L', 'N', -1, 0, 0, A, 1, TAU, 1, C, 1, W,1,INFO)
      CALL CHKXER( 'ZGEMLQ', INFOT, NOUT, LERR, OK )
      INFOT = 4
      CALL ZGEMLQ( 'L', 'N', 0, -1, 0, A, 1, TAU, 1, C, 1, W,1,INFO)
      CALL CHKXER( 'ZGEMLQ', INFOT, NOUT, LERR, OK )
      INFOT = 5
      CALL ZGEMLQ( 'L', 'N', 0, 0, -1, A, 1, TAU, 1, C, 1, W,1,INFO)
      CALL CHKXER( 'ZGEMLQ', INFOT, NOUT, LERR, OK )
      INFOT = 5
      CALL ZGEMLQ( 'R', 'N', 0, 0, -1, A, 1, TAU, 1, C, 1, W,1,INFO)
      CALL CHKXER( 'ZGEMLQ', INFOT, NOUT, LERR, OK )
      INFOT = 7
      CALL ZGEMLQ( 'L', 'N', 1, 2, 0, A, 0, TAU, 1, C, 1, W, 1,INFO)
      CALL CHKXER( 'ZGEMLQ', INFOT, NOUT, LERR, OK )
      INFOT = 9
      CALL ZGEMLQ( 'R', 'N', 2, 2, 1, A, 1, TAU, 0, C, 1, W, 1,INFO)
      CALL CHKXER( 'ZGEMLQ', INFOT, NOUT, LERR, OK )
      INFOT = 9
      CALL ZGEMLQ( 'L', 'N', 2, 2, 1, A, 1, TAU, 0, C, 1, W, 1,INFO)
      CALL CHKXER( 'ZGEMLQ', INFOT, NOUT, LERR, OK )
      INFOT = 11
      CALL ZGEMLQ( 'L', 'N', 1, 2, 1, A, 1, TAU, 6, C, 0, W, 1,INFO)
      CALL CHKXER( 'ZGEMLQ', INFOT, NOUT, LERR, OK )
      INFOT = 13
      CALL ZGEMLQ( 'L', 'N', 2, 2, 1, A, 2, TAU, 6, C, 2, W, 0,INFO)
      CALL CHKXER( 'ZGEMLQ', INFOT, NOUT, LERR, OK )

      // Print a summary line.

      CALL ALAESM( PATH, OK, NOUT )

      RETURN

      // End of ZERRTSQR

      }
