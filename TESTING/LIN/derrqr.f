      SUBROUTINE DERRQR( PATH, NUNIT )

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
      double             A( NMAX, NMAX ), AF( NMAX, NMAX ), B( NMAX ), W( NMAX ), X( NMAX );
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, DGEQR2, DGEQR2P, DGEQRF, DGEQRFP, DORG2R, DORGQR, DORM2R, DORMQR
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

      DO 20 J = 1, NMAX
         DO 10 I = 1, NMAX
            A( I, J ) = 1.D0 / DBLE( I+J )
            AF( I, J ) = 1.D0 / DBLE( I+J )
   10    CONTINUE
         B( J ) = 0.D0
         W( J ) = 0.D0
         X( J ) = 0.D0
   20 CONTINUE
      OK = .TRUE.

      // Error exits for QR factorization

      // DGEQRF

      SRNAMT = 'DGEQRF'
      INFOT = 1
      dgeqrf(-1, 0, A, 1, B, W, 1, INFO );
      chkxer('DGEQRF', INFOT, NOUT, LERR, OK );
      INFOT = 2
      dgeqrf(0, -1, A, 1, B, W, 1, INFO );
      chkxer('DGEQRF', INFOT, NOUT, LERR, OK );
      INFOT = 4
      dgeqrf(2, 1, A, 1, B, W, 1, INFO );
      chkxer('DGEQRF', INFOT, NOUT, LERR, OK );
      INFOT = 7
      dgeqrf(1, 2, A, 1, B, W, 1, INFO );
      chkxer('DGEQRF', INFOT, NOUT, LERR, OK );

      // DGEQRFP

      SRNAMT = 'DGEQRFP'
      INFOT = 1
      dgeqrfp(-1, 0, A, 1, B, W, 1, INFO );
      chkxer('DGEQRFP', INFOT, NOUT, LERR, OK );
      INFOT = 2
      dgeqrfp(0, -1, A, 1, B, W, 1, INFO );
      chkxer('DGEQRFP', INFOT, NOUT, LERR, OK );
      INFOT = 4
      dgeqrfp(2, 1, A, 1, B, W, 1, INFO );
      chkxer('DGEQRFP', INFOT, NOUT, LERR, OK );
      INFOT = 7
      dgeqrfp(1, 2, A, 1, B, W, 1, INFO );
      chkxer('DGEQRFP', INFOT, NOUT, LERR, OK );

      // DGEQR2

      SRNAMT = 'DGEQR2'
      INFOT = 1
      dgeqr2(-1, 0, A, 1, B, W, INFO );
      chkxer('DGEQR2', INFOT, NOUT, LERR, OK );
      INFOT = 2
      dgeqr2(0, -1, A, 1, B, W, INFO );
      chkxer('DGEQR2', INFOT, NOUT, LERR, OK );
      INFOT = 4
      dgeqr2(2, 1, A, 1, B, W, INFO );
      chkxer('DGEQR2', INFOT, NOUT, LERR, OK );

      // DGEQR2P

      SRNAMT = 'DGEQR2P'
      INFOT = 1
      dgeqr2p(-1, 0, A, 1, B, W, INFO );
      chkxer('DGEQR2P', INFOT, NOUT, LERR, OK );
      INFOT = 2
      dgeqr2p(0, -1, A, 1, B, W, INFO );
      chkxer('DGEQR2P', INFOT, NOUT, LERR, OK );
      INFOT = 4
      dgeqr2p(2, 1, A, 1, B, W, INFO );
      chkxer('DGEQR2P', INFOT, NOUT, LERR, OK );

      // DORGQR

      SRNAMT = 'DORGQR'
      INFOT = 1
      dorgqr(-1, 0, 0, A, 1, X, W, 1, INFO );
      chkxer('DORGQR', INFOT, NOUT, LERR, OK );
      INFOT = 2
      dorgqr(0, -1, 0, A, 1, X, W, 1, INFO );
      chkxer('DORGQR', INFOT, NOUT, LERR, OK );
      INFOT = 2
      dorgqr(1, 2, 0, A, 1, X, W, 2, INFO );
      chkxer('DORGQR', INFOT, NOUT, LERR, OK );
      INFOT = 3
      dorgqr(0, 0, -1, A, 1, X, W, 1, INFO );
      chkxer('DORGQR', INFOT, NOUT, LERR, OK );
      INFOT = 3
      dorgqr(1, 1, 2, A, 1, X, W, 1, INFO );
      chkxer('DORGQR', INFOT, NOUT, LERR, OK );
      INFOT = 5
      dorgqr(2, 2, 0, A, 1, X, W, 2, INFO );
      chkxer('DORGQR', INFOT, NOUT, LERR, OK );
      INFOT = 8
      dorgqr(2, 2, 0, A, 2, X, W, 1, INFO );
      chkxer('DORGQR', INFOT, NOUT, LERR, OK );

      // DORG2R

      SRNAMT = 'DORG2R'
      INFOT = 1
      dorg2r(-1, 0, 0, A, 1, X, W, INFO );
      chkxer('DORG2R', INFOT, NOUT, LERR, OK );
      INFOT = 2
      dorg2r(0, -1, 0, A, 1, X, W, INFO );
      chkxer('DORG2R', INFOT, NOUT, LERR, OK );
      INFOT = 2
      dorg2r(1, 2, 0, A, 1, X, W, INFO );
      chkxer('DORG2R', INFOT, NOUT, LERR, OK );
      INFOT = 3
      dorg2r(0, 0, -1, A, 1, X, W, INFO );
      chkxer('DORG2R', INFOT, NOUT, LERR, OK );
      INFOT = 3
      dorg2r(2, 1, 2, A, 2, X, W, INFO );
      chkxer('DORG2R', INFOT, NOUT, LERR, OK );
      INFOT = 5
      dorg2r(2, 1, 0, A, 1, X, W, INFO );
      chkxer('DORG2R', INFOT, NOUT, LERR, OK );

      // DORMQR

      SRNAMT = 'DORMQR'
      INFOT = 1
      dormqr('/', 'N', 0, 0, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('DORMQR', INFOT, NOUT, LERR, OK );
      INFOT = 2
      dormqr('L', '/', 0, 0, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('DORMQR', INFOT, NOUT, LERR, OK );
      INFOT = 3
      dormqr('L', 'N', -1, 0, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('DORMQR', INFOT, NOUT, LERR, OK );
      INFOT = 4
      dormqr('L', 'N', 0, -1, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('DORMQR', INFOT, NOUT, LERR, OK );
      INFOT = 5
      dormqr('L', 'N', 0, 0, -1, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('DORMQR', INFOT, NOUT, LERR, OK );
      INFOT = 5
      dormqr('L', 'N', 0, 1, 1, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('DORMQR', INFOT, NOUT, LERR, OK );
      INFOT = 5
      dormqr('R', 'N', 1, 0, 1, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('DORMQR', INFOT, NOUT, LERR, OK );
      INFOT = 7
      dormqr('L', 'N', 2, 1, 0, A, 1, X, AF, 2, W, 1, INFO );
      chkxer('DORMQR', INFOT, NOUT, LERR, OK );
      INFOT = 7
      dormqr('R', 'N', 1, 2, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('DORMQR', INFOT, NOUT, LERR, OK );
      INFOT = 10
      dormqr('L', 'N', 2, 1, 0, A, 2, X, AF, 1, W, 1, INFO );
      chkxer('DORMQR', INFOT, NOUT, LERR, OK );
      INFOT = 12
      dormqr('L', 'N', 1, 2, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('DORMQR', INFOT, NOUT, LERR, OK );
      INFOT = 12
      dormqr('R', 'N', 2, 1, 0, A, 1, X, AF, 2, W, 1, INFO );
      chkxer('DORMQR', INFOT, NOUT, LERR, OK );

      // DORM2R

      SRNAMT = 'DORM2R'
      INFOT = 1
      dorm2r('/', 'N', 0, 0, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('DORM2R', INFOT, NOUT, LERR, OK );
      INFOT = 2
      dorm2r('L', '/', 0, 0, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('DORM2R', INFOT, NOUT, LERR, OK );
      INFOT = 3
      dorm2r('L', 'N', -1, 0, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('DORM2R', INFOT, NOUT, LERR, OK );
      INFOT = 4
      dorm2r('L', 'N', 0, -1, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('DORM2R', INFOT, NOUT, LERR, OK );
      INFOT = 5
      dorm2r('L', 'N', 0, 0, -1, A, 1, X, AF, 1, W, INFO );
      chkxer('DORM2R', INFOT, NOUT, LERR, OK );
      INFOT = 5
      dorm2r('L', 'N', 0, 1, 1, A, 1, X, AF, 1, W, INFO );
      chkxer('DORM2R', INFOT, NOUT, LERR, OK );
      INFOT = 5
      dorm2r('R', 'N', 1, 0, 1, A, 1, X, AF, 1, W, INFO );
      chkxer('DORM2R', INFOT, NOUT, LERR, OK );
      INFOT = 7
      dorm2r('L', 'N', 2, 1, 0, A, 1, X, AF, 2, W, INFO );
      chkxer('DORM2R', INFOT, NOUT, LERR, OK );
      INFOT = 7
      dorm2r('R', 'N', 1, 2, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('DORM2R', INFOT, NOUT, LERR, OK );
      INFOT = 10
      dorm2r('L', 'N', 2, 1, 0, A, 2, X, AF, 1, W, INFO );
      chkxer('DORM2R', INFOT, NOUT, LERR, OK );

      // Print a summary line.

      alaesm(PATH, OK, NOUT );

      RETURN

      // End of DERRQR

      }
