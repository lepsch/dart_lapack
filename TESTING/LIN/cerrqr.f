      SUBROUTINE CERRQR( PATH, NUNIT )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             PATH;
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
      int                I, INFO, J
*     ..
*     .. Local Arrays ..
      COMPLEX            A( NMAX, NMAX ), AF( NMAX, NMAX ), B( NMAX ), W( NMAX ), X( NMAX )
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALAESM, CGEQR2, CGEQR2P, CGEQRF, CGEQRFP, CHKXER, CUNG2R, CUNGQR, CUNM2R, CUNMQR
*     ..
*     .. Scalars in Common ..
      bool               LERR, OK;
      String             SRNAMT;
      int                INFOT, NOUT
*     ..
*     .. Common blocks ..
      COMMON             / INFOC / INFOT, NOUT, OK, LERR
      COMMON             / SRNAMC / SRNAMT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CMPLX, REAL
*     ..
*     .. Executable Statements ..
*
      NOUT = NUNIT
      WRITE( NOUT, FMT = * )
*
*     Set the variables to innocuous values.
*
      DO 20 J = 1, NMAX
         DO 10 I = 1, NMAX
            A( I, J ) = CMPLX( 1. / REAL( I+J ), -1. / REAL( I+J ) )
            AF( I, J ) = CMPLX( 1. / REAL( I+J ), -1. / REAL( I+J ) )
   10    CONTINUE
         B( J ) = 0.
         W( J ) = 0.
         X( J ) = 0.
   20 CONTINUE
      OK = .TRUE.
*
*     Error exits for QR factorization
*
*     CGEQRF
*
      SRNAMT = 'CGEQRF'
      INFOT = 1
      CALL CGEQRF( -1, 0, A, 1, B, W, 1, INFO )
      CALL CHKXER( 'CGEQRF', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL CGEQRF( 0, -1, A, 1, B, W, 1, INFO )
      CALL CHKXER( 'CGEQRF', INFOT, NOUT, LERR, OK )
      INFOT = 4
      CALL CGEQRF( 2, 1, A, 1, B, W, 1, INFO )
      CALL CHKXER( 'CGEQRF', INFOT, NOUT, LERR, OK )
      INFOT = 7
      CALL CGEQRF( 1, 2, A, 1, B, W, 1, INFO )
      CALL CHKXER( 'CGEQRF', INFOT, NOUT, LERR, OK )
*
*     CGEQRFP
*
      SRNAMT = 'CGEQRFP'
      INFOT = 1
      CALL CGEQRFP( -1, 0, A, 1, B, W, 1, INFO )
      CALL CHKXER( 'CGEQRFP', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL CGEQRFP( 0, -1, A, 1, B, W, 1, INFO )
      CALL CHKXER( 'CGEQRFP', INFOT, NOUT, LERR, OK )
      INFOT = 4
      CALL CGEQRFP( 2, 1, A, 1, B, W, 1, INFO )
      CALL CHKXER( 'CGEQRFP', INFOT, NOUT, LERR, OK )
      INFOT = 7
      CALL CGEQRFP( 1, 2, A, 1, B, W, 1, INFO )
      CALL CHKXER( 'CGEQRFP', INFOT, NOUT, LERR, OK )
*
*     CGEQR2
*
      SRNAMT = 'CGEQR2'
      INFOT = 1
      CALL CGEQR2( -1, 0, A, 1, B, W, INFO )
      CALL CHKXER( 'CGEQR2', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL CGEQR2( 0, -1, A, 1, B, W, INFO )
      CALL CHKXER( 'CGEQR2', INFOT, NOUT, LERR, OK )
      INFOT = 4
      CALL CGEQR2( 2, 1, A, 1, B, W, INFO )
      CALL CHKXER( 'CGEQR2', INFOT, NOUT, LERR, OK )
*
*     CGEQR2P
*
      SRNAMT = 'CGEQR2P'
      INFOT = 1
      CALL CGEQR2P( -1, 0, A, 1, B, W, INFO )
      CALL CHKXER( 'CGEQR2P', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL CGEQR2P( 0, -1, A, 1, B, W, INFO )
      CALL CHKXER( 'CGEQR2P', INFOT, NOUT, LERR, OK )
      INFOT = 4
      CALL CGEQR2P( 2, 1, A, 1, B, W, INFO )
      CALL CHKXER( 'CGEQR2P', INFOT, NOUT, LERR, OK )
*
*     CUNGQR
*
      SRNAMT = 'CUNGQR'
      INFOT = 1
      CALL CUNGQR( -1, 0, 0, A, 1, X, W, 1, INFO )
      CALL CHKXER( 'CUNGQR', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL CUNGQR( 0, -1, 0, A, 1, X, W, 1, INFO )
      CALL CHKXER( 'CUNGQR', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL CUNGQR( 1, 2, 0, A, 1, X, W, 2, INFO )
      CALL CHKXER( 'CUNGQR', INFOT, NOUT, LERR, OK )
      INFOT = 3
      CALL CUNGQR( 0, 0, -1, A, 1, X, W, 1, INFO )
      CALL CHKXER( 'CUNGQR', INFOT, NOUT, LERR, OK )
      INFOT = 3
      CALL CUNGQR( 1, 1, 2, A, 1, X, W, 1, INFO )
      CALL CHKXER( 'CUNGQR', INFOT, NOUT, LERR, OK )
      INFOT = 5
      CALL CUNGQR( 2, 2, 0, A, 1, X, W, 2, INFO )
      CALL CHKXER( 'CUNGQR', INFOT, NOUT, LERR, OK )
      INFOT = 8
      CALL CUNGQR( 2, 2, 0, A, 2, X, W, 1, INFO )
      CALL CHKXER( 'CUNGQR', INFOT, NOUT, LERR, OK )
*
*     CUNG2R
*
      SRNAMT = 'CUNG2R'
      INFOT = 1
      CALL CUNG2R( -1, 0, 0, A, 1, X, W, INFO )
      CALL CHKXER( 'CUNG2R', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL CUNG2R( 0, -1, 0, A, 1, X, W, INFO )
      CALL CHKXER( 'CUNG2R', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL CUNG2R( 1, 2, 0, A, 1, X, W, INFO )
      CALL CHKXER( 'CUNG2R', INFOT, NOUT, LERR, OK )
      INFOT = 3
      CALL CUNG2R( 0, 0, -1, A, 1, X, W, INFO )
      CALL CHKXER( 'CUNG2R', INFOT, NOUT, LERR, OK )
      INFOT = 3
      CALL CUNG2R( 2, 1, 2, A, 2, X, W, INFO )
      CALL CHKXER( 'CUNG2R', INFOT, NOUT, LERR, OK )
      INFOT = 5
      CALL CUNG2R( 2, 1, 0, A, 1, X, W, INFO )
      CALL CHKXER( 'CUNG2R', INFOT, NOUT, LERR, OK )
*
*     CUNMQR
*
      SRNAMT = 'CUNMQR'
      INFOT = 1
      CALL CUNMQR( '/', 'N', 0, 0, 0, A, 1, X, AF, 1, W, 1, INFO )
      CALL CHKXER( 'CUNMQR', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL CUNMQR( 'L', '/', 0, 0, 0, A, 1, X, AF, 1, W, 1, INFO )
      CALL CHKXER( 'CUNMQR', INFOT, NOUT, LERR, OK )
      INFOT = 3
      CALL CUNMQR( 'L', 'N', -1, 0, 0, A, 1, X, AF, 1, W, 1, INFO )
      CALL CHKXER( 'CUNMQR', INFOT, NOUT, LERR, OK )
      INFOT = 4
      CALL CUNMQR( 'L', 'N', 0, -1, 0, A, 1, X, AF, 1, W, 1, INFO )
      CALL CHKXER( 'CUNMQR', INFOT, NOUT, LERR, OK )
      INFOT = 5
      CALL CUNMQR( 'L', 'N', 0, 0, -1, A, 1, X, AF, 1, W, 1, INFO )
      CALL CHKXER( 'CUNMQR', INFOT, NOUT, LERR, OK )
      INFOT = 5
      CALL CUNMQR( 'L', 'N', 0, 1, 1, A, 1, X, AF, 1, W, 1, INFO )
      CALL CHKXER( 'CUNMQR', INFOT, NOUT, LERR, OK )
      INFOT = 5
      CALL CUNMQR( 'R', 'N', 1, 0, 1, A, 1, X, AF, 1, W, 1, INFO )
      CALL CHKXER( 'CUNMQR', INFOT, NOUT, LERR, OK )
      INFOT = 7
      CALL CUNMQR( 'L', 'N', 2, 1, 0, A, 1, X, AF, 2, W, 1, INFO )
      CALL CHKXER( 'CUNMQR', INFOT, NOUT, LERR, OK )
      INFOT = 7
      CALL CUNMQR( 'R', 'N', 1, 2, 0, A, 1, X, AF, 1, W, 1, INFO )
      CALL CHKXER( 'CUNMQR', INFOT, NOUT, LERR, OK )
      INFOT = 10
      CALL CUNMQR( 'L', 'N', 2, 1, 0, A, 2, X, AF, 1, W, 1, INFO )
      CALL CHKXER( 'CUNMQR', INFOT, NOUT, LERR, OK )
      INFOT = 12
      CALL CUNMQR( 'L', 'N', 1, 2, 0, A, 1, X, AF, 1, W, 1, INFO )
      CALL CHKXER( 'CUNMQR', INFOT, NOUT, LERR, OK )
      INFOT = 12
      CALL CUNMQR( 'R', 'N', 2, 1, 0, A, 1, X, AF, 2, W, 1, INFO )
      CALL CHKXER( 'CUNMQR', INFOT, NOUT, LERR, OK )
*
*     CUNM2R
*
      SRNAMT = 'CUNM2R'
      INFOT = 1
      CALL CUNM2R( '/', 'N', 0, 0, 0, A, 1, X, AF, 1, W, INFO )
      CALL CHKXER( 'CUNM2R', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL CUNM2R( 'L', '/', 0, 0, 0, A, 1, X, AF, 1, W, INFO )
      CALL CHKXER( 'CUNM2R', INFOT, NOUT, LERR, OK )
      INFOT = 3
      CALL CUNM2R( 'L', 'N', -1, 0, 0, A, 1, X, AF, 1, W, INFO )
      CALL CHKXER( 'CUNM2R', INFOT, NOUT, LERR, OK )
      INFOT = 4
      CALL CUNM2R( 'L', 'N', 0, -1, 0, A, 1, X, AF, 1, W, INFO )
      CALL CHKXER( 'CUNM2R', INFOT, NOUT, LERR, OK )
      INFOT = 5
      CALL CUNM2R( 'L', 'N', 0, 0, -1, A, 1, X, AF, 1, W, INFO )
      CALL CHKXER( 'CUNM2R', INFOT, NOUT, LERR, OK )
      INFOT = 5
      CALL CUNM2R( 'L', 'N', 0, 1, 1, A, 1, X, AF, 1, W, INFO )
      CALL CHKXER( 'CUNM2R', INFOT, NOUT, LERR, OK )
      INFOT = 5
      CALL CUNM2R( 'R', 'N', 1, 0, 1, A, 1, X, AF, 1, W, INFO )
      CALL CHKXER( 'CUNM2R', INFOT, NOUT, LERR, OK )
      INFOT = 7
      CALL CUNM2R( 'L', 'N', 2, 1, 0, A, 1, X, AF, 2, W, INFO )
      CALL CHKXER( 'CUNM2R', INFOT, NOUT, LERR, OK )
      INFOT = 7
      CALL CUNM2R( 'R', 'N', 1, 2, 0, A, 1, X, AF, 1, W, INFO )
      CALL CHKXER( 'CUNM2R', INFOT, NOUT, LERR, OK )
      INFOT = 10
      CALL CUNM2R( 'L', 'N', 2, 1, 0, A, 2, X, AF, 1, W, INFO )
      CALL CHKXER( 'CUNM2R', INFOT, NOUT, LERR, OK )
*
*     Print a summary line.
*
      CALL ALAESM( PATH, OK, NOUT )
*
      RETURN
*
*     End of CERRQR
*
      END
