      SUBROUTINE CERRLS( PATH, NUNIT )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             PATH;
      int                NUNIT;
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      int                NMAX;
      PARAMETER          ( NMAX = 2 )
*     ..
*     .. Local Scalars ..
      String             C2;
      int                INFO, IRNK;
      REAL               RCOND
*     ..
*     .. Local Arrays ..
      int                IP( NMAX );
      REAL               RW( NMAX ), S( NMAX )
      COMPLEX            A( NMAX, NMAX ), B( NMAX, NMAX ), W( NMAX )
*     ..
*     .. External Functions ..
      bool               LSAMEN;
      // EXTERNAL LSAMEN
*     ..
*     .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, CGELS, CGELSD, CGELSS, CGELST, CGELSY, CGETSLS
*     ..
*     .. Scalars in Common ..
      bool               LERR, OK;
      String             SRNAMT;
      int                INFOT, NOUT;
*     ..
*     .. Common blocks ..
      COMMON             / INFOC / INFOT, NOUT, OK, LERR
      COMMON             / SRNAMC / SRNAMT
*     ..
*     .. Executable Statements ..
*
      NOUT = NUNIT
      C2 = PATH( 2: 3 )
      A( 1, 1 ) = ( 1.0E+0, 0.0E+0 )
      A( 1, 2 ) = ( 2.0E+0, 0.0E+0 )
      A( 2, 2 ) = ( 3.0E+0, 0.0E+0 )
      A( 2, 1 ) = ( 4.0E+0, 0.0E+0 )
      OK = .TRUE.
      WRITE( NOUT, FMT = * )
*
*     Test error exits for the least squares driver routines.
*
      IF( LSAMEN( 2, C2, 'LS' ) ) THEN
*
*        CGELS
*
         SRNAMT = 'CGELS '
         INFOT = 1
         CALL CGELS( '/', 0, 0, 0, A, 1, B, 1, W, 1, INFO )
         CALL CHKXER( 'CGELS ', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL CGELS( 'N', -1, 0, 0, A, 1, B, 1, W, 1, INFO )
         CALL CHKXER( 'CGELS ', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL CGELS( 'N', 0, -1, 0, A, 1, B, 1, W, 1, INFO )
         CALL CHKXER( 'CGELS ', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL CGELS( 'N', 0, 0, -1, A, 1, B, 1, W, 1, INFO )
         CALL CHKXER( 'CGELS ', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL CGELS( 'N', 2, 0, 0, A, 1, B, 2, W, 2, INFO )
         CALL CHKXER( 'CGELS ', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL CGELS( 'N', 2, 0, 0, A, 2, B, 1, W, 2, INFO )
         CALL CHKXER( 'CGELS ', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL CGELS( 'N', 0, 2, 0, A, 1, B, 1, W, 2, INFO )
         CALL CHKXER( 'CGELS', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL CGELS( 'N', 1, 1, 0, A, 1, B, 1, W, 1, INFO )
         CALL CHKXER( 'CGELS ', INFOT, NOUT, LERR, OK )
*
*        CGELST
*
         SRNAMT = 'CGELST'
         INFOT = 1
         CALL CGELST( '/', 0, 0, 0, A, 1, B, 1, W, 1, INFO )
         CALL CHKXER( 'CGELST', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL CGELST( 'N', -1, 0, 0, A, 1, B, 1, W, 1, INFO )
         CALL CHKXER( 'CGELST', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL CGELST( 'N', 0, -1, 0, A, 1, B, 1, W, 1, INFO )
         CALL CHKXER( 'CGELST', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL CGELST( 'N', 0, 0, -1, A, 1, B, 1, W, 1, INFO )
         CALL CHKXER( 'CGELST', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL CGELST( 'N', 2, 0, 0, A, 1, B, 2, W, 2, INFO )
         CALL CHKXER( 'CGELST', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL CGELST( 'N', 2, 0, 0, A, 2, B, 1, W, 2, INFO )
         CALL CHKXER( 'CGELST', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL CGELST( 'N', 0, 2, 0, A, 1, B, 1, W, 2, INFO )
         CALL CHKXER( 'CGELST', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL CGELST( 'N', 1, 1, 0, A, 1, B, 1, W, 1, INFO )
         CALL CHKXER( 'CGELST', INFOT, NOUT, LERR, OK )
*
*        CGETSLS
*
         SRNAMT = 'CGETSLS'
         INFOT = 1
         CALL CGETSLS( '/', 0, 0, 0, A, 1, B, 1, W, 1, INFO )
         CALL CHKXER( 'CGETSLS', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL CGETSLS( 'N', -1, 0, 0, A, 1, B, 1, W, 1, INFO )
         CALL CHKXER( 'CGETSLS', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL CGETSLS( 'N', 0, -1, 0, A, 1, B, 1, W, 1, INFO )
         CALL CHKXER( 'CGETSLS', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL CGETSLS( 'N', 0, 0, -1, A, 1, B, 1, W, 1, INFO )
         CALL CHKXER( 'CGETSLS', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL CGETSLS( 'N', 2, 0, 0, A, 1, B, 2, W, 2, INFO )
         CALL CHKXER( 'CGETSLS', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL CGETSLS( 'N', 2, 0, 0, A, 2, B, 1, W, 2, INFO )
         CALL CHKXER( 'CGETSLS', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL CGETSLS( 'N', 0, 2, 0, A, 1, B, 1, W, 2, INFO )
         CALL CHKXER( 'CGETSLS', INFOT, NOUT, LERR, OK )
*
*        CGELSS
*
         SRNAMT = 'CGELSS'
         INFOT = 1
         CALL CGELSS( -1, 0, 0, A, 1, B, 1, S, RCOND, IRNK, W, 1, RW, INFO )
         CALL CHKXER( 'CGELSS', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL CGELSS( 0, -1, 0, A, 1, B, 1, S, RCOND, IRNK, W, 1, RW, INFO )
         CALL CHKXER( 'CGELSS', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL CGELSS( 0, 0, -1, A, 1, B, 1, S, RCOND, IRNK, W, 1, RW, INFO )
         CALL CHKXER( 'CGELSS', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL CGELSS( 2, 0, 0, A, 1, B, 2, S, RCOND, IRNK, W, 2, RW, INFO )
         CALL CHKXER( 'CGELSS', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL CGELSS( 2, 0, 0, A, 2, B, 1, S, RCOND, IRNK, W, 2, RW, INFO )
         CALL CHKXER( 'CGELSS', INFOT, NOUT, LERR, OK )
*
*        CGELSY
*
         SRNAMT = 'CGELSY'
         INFOT = 1
         CALL CGELSY( -1, 0, 0, A, 1, B, 1, IP, RCOND, IRNK, W, 10, RW, INFO )
         CALL CHKXER( 'CGELSY', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL CGELSY( 0, -1, 0, A, 1, B, 1, IP, RCOND, IRNK, W, 10, RW, INFO )
         CALL CHKXER( 'CGELSY', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL CGELSY( 0, 0, -1, A, 1, B, 1, IP, RCOND, IRNK, W, 10, RW, INFO )
         CALL CHKXER( 'CGELSY', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL CGELSY( 2, 0, 0, A, 1, B, 2, IP, RCOND, IRNK, W, 10, RW, INFO )
         CALL CHKXER( 'CGELSY', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL CGELSY( 2, 0, 0, A, 2, B, 1, IP, RCOND, IRNK, W, 10, RW, INFO )
         CALL CHKXER( 'CGELSY', INFOT, NOUT, LERR, OK )
         INFOT = 12
         CALL CGELSY( 0, 3, 0, A, 1, B, 3, IP, RCOND, IRNK, W, 1, RW, INFO )
         CALL CHKXER( 'CGELSY', INFOT, NOUT, LERR, OK )
*
*        CGELSD
*
         SRNAMT = 'CGELSD'
         INFOT = 1
         CALL CGELSD( -1, 0, 0, A, 1, B, 1, S, RCOND, IRNK, W, 10, RW, IP, INFO )
         CALL CHKXER( 'CGELSD', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL CGELSD( 0, -1, 0, A, 1, B, 1, S, RCOND, IRNK, W, 10, RW, IP, INFO )
         CALL CHKXER( 'CGELSD', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL CGELSD( 0, 0, -1, A, 1, B, 1, S, RCOND, IRNK, W, 10, RW, IP, INFO )
         CALL CHKXER( 'CGELSD', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL CGELSD( 2, 0, 0, A, 1, B, 2, S, RCOND, IRNK, W, 10, RW, IP, INFO )
         CALL CHKXER( 'CGELSD', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL CGELSD( 2, 0, 0, A, 2, B, 1, S, RCOND, IRNK, W, 10, RW, IP, INFO )
         CALL CHKXER( 'CGELSD', INFOT, NOUT, LERR, OK )
         INFOT = 12
         CALL CGELSD( 2, 2, 1, A, 2, B, 2, S, RCOND, IRNK, W, 1, RW, IP, INFO )
         CALL CHKXER( 'CGELSD', INFOT, NOUT, LERR, OK )
      END IF
*
*     Print a summary line.
*
      CALL ALAESM( PATH, OK, NOUT )
*
      RETURN
*
*     End of CERRLS
*
      END
