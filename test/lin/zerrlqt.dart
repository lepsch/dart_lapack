      void zerrlqt(PATH, NUNIT ) {
      // IMPLICIT NONE

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             PATH;
      int                NUNIT;
      // ..

// =====================================================================

      // .. Parameters ..
      int                NMAX;
      const              NMAX = 2 ;
      // ..
      // .. Local Scalars ..
      int                I, INFO, J;
      // ..
      // .. Local Arrays ..
      Complex   A( NMAX, NMAX ), T( NMAX, NMAX ), W( NMAX ), C( NMAX, NMAX );
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, ZGELQT3, ZGELQT, ZGEMLQT
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
      // INTRINSIC DBLE, DCMPLX
      // ..
      // .. Executable Statements ..

      NOUT = NUNIT;
      WRITE( NOUT, FMT = * );

      // Set the variables to innocuous values.

      for (J = 1; J <= NMAX; J++) {
         for (I = 1; I <= NMAX; I++) {
            A[I, J] = 1.0 / DCMPLX( DBLE( I+J ), 0.0 );
            C[I, J] = 1.0 / DCMPLX( DBLE( I+J ), 0.0 );
            T[I, J] = 1.0 / DCMPLX( DBLE( I+J ), 0.0 );
         }
         W[J] = 0.0;
      }
      OK = true;

      // Error exits for LQT factorization

      // ZGELQT

      SRNAMT = 'ZGELQT';
      INFOT = 1;
      zgelqt(-1, 0, 1, A, 1, T, 1, W, INFO );
      chkxer('ZGELQT', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      zgelqt(0, -1, 1, A, 1, T, 1, W, INFO );
      chkxer('ZGELQT', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      zgelqt(0, 0, 0, A, 1, T, 1, W, INFO );
      chkxer('ZGELQT', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      zgelqt(2, 1, 1, A, 1, T, 1, W, INFO );
      chkxer('ZGELQT', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      zgelqt(2, 2, 2, A, 2, T, 1, W, INFO );
      chkxer('ZGELQT', INFOT, NOUT, LERR, OK );

      // ZGELQT3

      SRNAMT = 'ZGELQT3';
      INFOT = 1;
      zgelqt3(-1, 0, A, 1, T, 1, INFO );
      chkxer('ZGELQT3', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      zgelqt3(0, -1, A, 1, T, 1, INFO );
      chkxer('ZGELQT3', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      zgelqt3(2, 2, A, 1, T, 1, INFO );
      chkxer('ZGELQT3', INFOT, NOUT, LERR, OK );
      INFOT = 6;
      zgelqt3(2, 2, A, 2, T, 1, INFO );
      chkxer('ZGELQT3', INFOT, NOUT, LERR, OK );

      // ZGEMLQT

      SRNAMT = 'ZGEMLQT';
      INFOT = 1;
      zgemlqt('/', 'N', 0, 0, 0, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('ZGEMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      zgemlqt('L', '/', 0, 0, 0, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('ZGEMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      zgemlqt('L', 'N', -1, 0, 0, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('ZGEMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      zgemlqt('L', 'N', 0, -1, 0, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('ZGEMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      zgemlqt('L', 'N', 0, 0, -1, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('ZGEMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      zgemlqt('R', 'N', 0, 0, -1, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('ZGEMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 6;
      zgemlqt('L', 'N', 0, 0, 0, 0, A, 1, T, 1, C, 1, W, INFO );
      chkxer('ZGEMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 8;
      zgemlqt('R', 'N', 2, 2, 2, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('ZGEMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 8;
      zgemlqt('L', 'N', 2, 2, 2, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('ZGEMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 10;
      zgemlqt('R', 'N', 1, 1, 1, 1, A, 1, T, 0, C, 1, W, INFO );
      chkxer('ZGEMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 12;
      zgemlqt('L', 'N', 1, 1, 1, 1, A, 1, T, 1, C, 0, W, INFO );
      chkxer('ZGEMLQT', INFOT, NOUT, LERR, OK );

      // Print a summary line.

      alaesm(PATH, OK, NOUT );

      return;
      }
