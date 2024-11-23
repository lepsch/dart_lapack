      SUBROUTINE ZERRLS( PATH, NUNIT )
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
      PARAMETER          ( NMAX = 2 )
*     ..
*     .. Local Scalars ..
      CHARACTER*2        C2
      int                INFO, IRNK
      DOUBLE PRECISION   RCOND
*     ..
*     .. Local Arrays ..
      int                IP( NMAX )
      DOUBLE PRECISION   RW( NMAX ), S( NMAX )
      COMPLEX*16         A( NMAX, NMAX ), B( NMAX, NMAX ), W( NMAX )
*     ..
*     .. External Functions ..
      LOGICAL            LSAMEN
      EXTERNAL           LSAMEN
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALAESM, CHKXER, ZGELS, ZGELSD, ZGELSS, ZGELST, ZGELSY, ZGETSLS
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
*     .. Executable Statements ..
*
      NOUT = NUNIT
      C2 = PATH( 2: 3 )
      A( 1, 1 ) = ( 1.0D+0, 0.0D+0 )
      A( 1, 2 ) = ( 2.0D+0, 0.0D+0 )
      A( 2, 2 ) = ( 3.0D+0, 0.0D+0 )
      A( 2, 1 ) = ( 4.0D+0, 0.0D+0 )
      OK = .TRUE.
      WRITE( NOUT, FMT = * )
*
*     Test error exits for the least squares driver routines.
*
      IF( LSAMEN( 2, C2, 'LS' ) ) THEN
*
*        ZGELS
*
         SRNAMT = 'ZGELS '
         INFOT = 1
         CALL ZGELS( '/', 0, 0, 0, A, 1, B, 1, W, 1, INFO )
         CALL CHKXER( 'ZGELS ', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL ZGELS( 'N', -1, 0, 0, A, 1, B, 1, W, 1, INFO )
         CALL CHKXER( 'ZGELS ', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL ZGELS( 'N', 0, -1, 0, A, 1, B, 1, W, 1, INFO )
         CALL CHKXER( 'ZGELS ', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL ZGELS( 'N', 0, 0, -1, A, 1, B, 1, W, 1, INFO )
         CALL CHKXER( 'ZGELS ', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL ZGELS( 'N', 2, 0, 0, A, 1, B, 2, W, 2, INFO )
         CALL CHKXER( 'ZGELS ', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL ZGELS( 'N', 2, 0, 0, A, 2, B, 1, W, 2, INFO )
         CALL CHKXER( 'ZGELS ', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL ZGELS( 'N', 0, 2, 0, A, 1, B, 1, W, 2, INFO )
         CALL CHKXER( 'ZGELS', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL ZGELS( 'N', 1, 1, 0, A, 1, B, 1, W, 1, INFO )
         CALL CHKXER( 'ZGELS ', INFOT, NOUT, LERR, OK )
*
*        ZGELST
*
         SRNAMT = 'ZGELST'
         INFOT = 1
         CALL ZGELST( '/', 0, 0, 0, A, 1, B, 1, W, 1, INFO )
         CALL CHKXER( 'ZGELST', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL ZGELST( 'N', -1, 0, 0, A, 1, B, 1, W, 1, INFO )
         CALL CHKXER( 'ZGELST', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL ZGELST( 'N', 0, -1, 0, A, 1, B, 1, W, 1, INFO )
         CALL CHKXER( 'ZGELST', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL ZGELST( 'N', 0, 0, -1, A, 1, B, 1, W, 1, INFO )
         CALL CHKXER( 'ZGELST', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL ZGELST( 'N', 2, 0, 0, A, 1, B, 2, W, 2, INFO )
         CALL CHKXER( 'ZGELST', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL ZGELST( 'N', 2, 0, 0, A, 2, B, 1, W, 2, INFO )
         CALL CHKXER( 'ZGELST', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL ZGELST( 'N', 0, 2, 0, A, 1, B, 1, W, 2, INFO )
         CALL CHKXER( 'ZGELST', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL ZGELST( 'N', 1, 1, 0, A, 1, B, 1, W, 1, INFO )
         CALL CHKXER( 'ZGELST', INFOT, NOUT, LERR, OK )
*
*        ZGETSLS
*
         SRNAMT = 'ZGETSLS'
         INFOT = 1
         CALL ZGETSLS( '/', 0, 0, 0, A, 1, B, 1, W, 1, INFO )
         CALL CHKXER( 'ZGETSLS', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL ZGETSLS( 'N', -1, 0, 0, A, 1, B, 1, W, 1, INFO )
         CALL CHKXER( 'ZGETSLS', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL ZGETSLS( 'N', 0, -1, 0, A, 1, B, 1, W, 1, INFO )
         CALL CHKXER( 'ZGETSLS', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL ZGETSLS( 'N', 0, 0, -1, A, 1, B, 1, W, 1, INFO )
         CALL CHKXER( 'ZGETSLS', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL ZGETSLS( 'N', 2, 0, 0, A, 1, B, 2, W, 2, INFO )
         CALL CHKXER( 'ZGETSLS', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL ZGETSLS( 'N', 2, 0, 0, A, 2, B, 1, W, 2, INFO )
         CALL CHKXER( 'ZGETSLS', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL ZGETSLS( 'N', 0, 2, 0, A, 1, B, 1, W, 2, INFO )
         CALL CHKXER( 'ZGETSLS', INFOT, NOUT, LERR, OK )
*
*        ZGELSS
*
         SRNAMT = 'ZGELSS'
         INFOT = 1
         CALL ZGELSS( -1, 0, 0, A, 1, B, 1, S, RCOND, IRNK, W, 1, RW, INFO )
         CALL CHKXER( 'ZGELSS', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL ZGELSS( 0, -1, 0, A, 1, B, 1, S, RCOND, IRNK, W, 1, RW, INFO )
         CALL CHKXER( 'ZGELSS', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL ZGELSS( 0, 0, -1, A, 1, B, 1, S, RCOND, IRNK, W, 1, RW, INFO )
         CALL CHKXER( 'ZGELSS', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL ZGELSS( 2, 0, 0, A, 1, B, 2, S, RCOND, IRNK, W, 2, RW, INFO )
         CALL CHKXER( 'ZGELSS', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL ZGELSS( 2, 0, 0, A, 2, B, 1, S, RCOND, IRNK, W, 2, RW, INFO )
         CALL CHKXER( 'ZGELSS', INFOT, NOUT, LERR, OK )
*
*        ZGELSY
*
         SRNAMT = 'ZGELSY'
         INFOT = 1
         CALL ZGELSY( -1, 0, 0, A, 1, B, 1, IP, RCOND, IRNK, W, 10, RW, INFO )
         CALL CHKXER( 'ZGELSY', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL ZGELSY( 0, -1, 0, A, 1, B, 1, IP, RCOND, IRNK, W, 10, RW, INFO )
         CALL CHKXER( 'ZGELSY', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL ZGELSY( 0, 0, -1, A, 1, B, 1, IP, RCOND, IRNK, W, 10, RW, INFO )
         CALL CHKXER( 'ZGELSY', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL ZGELSY( 2, 0, 0, A, 1, B, 2, IP, RCOND, IRNK, W, 10, RW, INFO )
         CALL CHKXER( 'ZGELSY', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL ZGELSY( 2, 0, 0, A, 2, B, 1, IP, RCOND, IRNK, W, 10, RW, INFO )
         CALL CHKXER( 'ZGELSY', INFOT, NOUT, LERR, OK )
         INFOT = 12
         CALL ZGELSY( 0, 3, 0, A, 1, B, 3, IP, RCOND, IRNK, W, 1, RW, INFO )
         CALL CHKXER( 'ZGELSY', INFOT, NOUT, LERR, OK )
*
*        ZGELSD
*
         SRNAMT = 'ZGELSD'
         INFOT = 1
         CALL ZGELSD( -1, 0, 0, A, 1, B, 1, S, RCOND, IRNK, W, 10, RW, IP, INFO )
         CALL CHKXER( 'ZGELSD', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL ZGELSD( 0, -1, 0, A, 1, B, 1, S, RCOND, IRNK, W, 10, RW, IP, INFO )
         CALL CHKXER( 'ZGELSD', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL ZGELSD( 0, 0, -1, A, 1, B, 1, S, RCOND, IRNK, W, 10, RW, IP, INFO )
         CALL CHKXER( 'ZGELSD', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL ZGELSD( 2, 0, 0, A, 1, B, 2, S, RCOND, IRNK, W, 10, RW, IP, INFO )
         CALL CHKXER( 'ZGELSD', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL ZGELSD( 2, 0, 0, A, 2, B, 1, S, RCOND, IRNK, W, 10, RW, IP, INFO )
         CALL CHKXER( 'ZGELSD', INFOT, NOUT, LERR, OK )
         INFOT = 12
         CALL ZGELSD( 2, 2, 1, A, 2, B, 2, S, RCOND, IRNK, W, 1, RW, IP, INFO )
         CALL CHKXER( 'ZGELSD', INFOT, NOUT, LERR, OK )
      END IF
*
*     Print a summary line.
*
      CALL ALAESM( PATH, OK, NOUT )
*
      RETURN
*
*     End of ZERRLS
*
      END
