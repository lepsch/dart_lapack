      SUBROUTINE ZERRQP( PATH, NUNIT )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER*3        PATH
      int                NUNIT
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      int                NMAX
      PARAMETER          ( NMAX = 3 )
*     ..
*     .. Local Scalars ..
      CHARACTER*2        C2
      int                INFO, LW
*     ..
*     .. Local Arrays ..
      int                IP( NMAX )
      DOUBLE PRECISION   RW( 2*NMAX )
      COMPLEX*16         A( NMAX, NMAX ), TAU( NMAX ), W( 2*NMAX+3*NMAX )
*     ..
*     .. External Functions ..
      LOGICAL            LSAMEN
      EXTERNAL           LSAMEN
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALAESM, CHKXER, ZGEQP3
*     ..
*     .. Scalars in Common ..
      LOGICAL            LERR, OK
      CHARACTER*32       SRNAMT
      int                INFOT, NOUT
*     ..
*     .. Common blocks ..
      COMMON             / INFOC / INFOT, NOUT, OK, LERR
      COMMON             / SRNAMC / SRNAMT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DCMPLX
*     ..
*     .. Executable Statements ..
*
      NOUT = NUNIT
      C2 = PATH( 2: 3 )
      LW = NMAX + 1
      A( 1, 1 ) = DCMPLX( 1.0D+0, -1.0D+0 )
      A( 1, 2 ) = DCMPLX( 2.0D+0, -2.0D+0 )
      A( 2, 2 ) = DCMPLX( 3.0D+0, -3.0D+0 )
      A( 2, 1 ) = DCMPLX( 4.0D+0, -4.0D+0 )
      OK = .TRUE.
      WRITE( NOUT, FMT = * )
*
*     Test error exits for QR factorization with pivoting
*
      IF( LSAMEN( 2, C2, 'QP' ) ) THEN
*
*        ZGEQP3
*
         SRNAMT = 'ZGEQP3'
         INFOT = 1
         CALL ZGEQP3( -1, 0, A, 1, IP, TAU, W, LW, RW, INFO )
         CALL CHKXER( 'ZGEQP3', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL ZGEQP3( 1, -1, A, 1, IP, TAU, W, LW, RW, INFO )
         CALL CHKXER( 'ZGEQP3', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL ZGEQP3( 2, 3, A, 1, IP, TAU, W, LW, RW, INFO )
         CALL CHKXER( 'ZGEQP3', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL ZGEQP3( 2, 2, A, 2, IP, TAU, W, LW-10, RW, INFO )
         CALL CHKXER( 'ZGEQP3', INFOT, NOUT, LERR, OK )
      END IF
*
*     Print a summary line.
*
      CALL ALAESM( PATH, OK, NOUT )
*
      RETURN
*
*     End of ZERRQP
*
      END
