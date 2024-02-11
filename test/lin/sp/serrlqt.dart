      void serrlqt(final int PATH, final int NUNIT,) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             PATH;
      int                NUNIT;
      // ..

      int                NMAX;
      const              NMAX = 2 ;
      int                I, INFO, J;
      double               A( NMAX, NMAX ), T( NMAX, NMAX ), W( NMAX ), C( NMAX, NMAX );
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, SGELQT3, SGELQT, SGEMLQT
      // ..
      // .. Scalars in Common ..
      bool               LERR, OK;
      String            srnamc.SRNAMT;
      int                INFOT, NOUT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / INFOT, NOUT, OK, LERR
      // COMMON / SRNAMC /srnamc.SRNAMT
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC REAL

      NOUT = NUNIT;
      WRITE( NOUT, FMT = * );

      // Set the variables to innocuous values.

      for (J = 1; J <= NMAX; J++) {
         for (I = 1; I <= NMAX; I++) {
            A[I][J] = 1. / REAL( I+J );
            C[I][J] = 1. / REAL( I+J );
            T[I][J] = 1. / REAL( I+J );
         }
         W[J] = 0.;
      }
      OK = true;

      // Error exits for LQT factorization

      // SGELQT

     srnamc.SRNAMT = 'SGELQT';
      INFOT = 1;
      sgelqt(-1, 0, 1, A, 1, T, 1, W, INFO );
      chkxer('SGELQT', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      sgelqt(0, -1, 1, A, 1, T, 1, W, INFO );
      chkxer('SGELQT', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      sgelqt(0, 0, 0, A, 1, T, 1, W, INFO );
      chkxer('SGELQT', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      sgelqt(2, 1, 1, A, 1, T, 1, W, INFO );
      chkxer('SGELQT', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      sgelqt(2, 2, 2, A, 2, T, 1, W, INFO );
      chkxer('SGELQT', INFOT, NOUT, LERR, OK );

      // SGELQT3

     srnamc.SRNAMT = 'SGELQT3';
      INFOT = 1;
      sgelqt3(-1, 0, A, 1, T, 1, INFO );
      chkxer('SGELQT3', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      sgelqt3(0, -1, A, 1, T, 1, INFO );
      chkxer('SGELQT3', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      sgelqt3(2, 2, A, 1, T, 1, INFO );
      chkxer('SGELQT3', INFOT, NOUT, LERR, OK );
      INFOT = 6;
      sgelqt3(2, 2, A, 2, T, 1, INFO );
      chkxer('SGELQT3', INFOT, NOUT, LERR, OK );

      // SGEMLQT

     srnamc.SRNAMT = 'SGEMLQT';
      INFOT = 1;
      sgemlqt('/', 'N', 0, 0, 0, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('SGEMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      sgemlqt('L', '/', 0, 0, 0, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('SGEMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      sgemlqt('L', 'N', -1, 0, 0, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('SGEMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      sgemlqt('L', 'N', 0, -1, 0, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('SGEMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      sgemlqt('L', 'N', 0, 0, -1, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('SGEMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      sgemlqt('R', 'N', 0, 0, -1, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('SGEMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 6;
      sgemlqt('L', 'N', 0, 0, 0, 0, A, 1, T, 1, C, 1, W, INFO );
      chkxer('SGEMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 8;
      sgemlqt('R', 'N', 2, 2, 2, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('SGEMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 8;
      sgemlqt('L', 'N', 2, 2, 2, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('SGEMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 10;
      sgemlqt('R', 'N', 1, 1, 1, 1, A, 1, T, 0, C, 1, W, INFO );
      chkxer('SGEMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 12;
      sgemlqt('L', 'N', 1, 1, 1, 1, A, 1, T, 1, C, 0, W, INFO );
      chkxer('SGEMLQT', INFOT, NOUT, LERR, OK );

      // Print a summary line.

      alaesm(PATH, OK, NOUT );

      }
