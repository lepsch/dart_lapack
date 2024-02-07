      void cerrlqtp(PATH, NUNIT ) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             PATH;
      int                NUNIT;
      // ..

      int                NMAX;
      const              NMAX = 2 ;
      int                I, INFO, J;
      Complex            A( NMAX, NMAX ), T( NMAX, NMAX ), W( NMAX ), B( NMAX, NMAX ), C( NMAX, NMAX );
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, CTPLQT2, CTPLQT, CTPMLQT
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
      // INTRINSIC REAL, CMPLX

      NOUT = NUNIT;
      WRITE( NOUT, FMT = * );

      // Set the variables to innocuous values.

      for (J = 1; J <= NMAX; J++) {
         for (I = 1; I <= NMAX; I++) {
            A[I][J] = 1.0 / CMPLX( REAL( I+J ), 0.0 );
            C[I][J] = 1.0 / CMPLX( REAL( I+J ), 0.0 );
            T[I][J] = 1.0 / CMPLX( REAL( I+J ), 0.0 );
         }
         W[J] = 0.0;
      }
      OK = true;

      // Error exits for TPLQT factorization

      // CTPLQT

     srnamc.SRNAMT = 'CTPLQT';
      INFOT = 1;
      ctplqt(-1, 1, 0, 1, A, 1, B, 1, T, 1, W, INFO );
      chkxer('CTPLQT', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      ctplqt(1, -1, 0, 1, A, 1, B, 1, T, 1, W, INFO );
      chkxer('CTPLQT', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      ctplqt(0, 1, -1, 1, A, 1, B, 1, T, 1, W, INFO );
      chkxer('CTPLQT', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      ctplqt(0, 1, 1, 1, A, 1, B, 1, T, 1, W, INFO );
      chkxer('CTPLQT', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      ctplqt(0, 1, 0, 0, A, 1, B, 1, T, 1, W, INFO );
      chkxer('CTPLQT', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      ctplqt(1, 1, 0, 2, A, 1, B, 1, T, 1, W, INFO );
      chkxer('CTPLQT', INFOT, NOUT, LERR, OK );
      INFOT = 6;
      ctplqt(2, 1, 0, 2, A, 1, B, 1, T, 1, W, INFO );
      chkxer('CTPLQT', INFOT, NOUT, LERR, OK );
      INFOT = 8;
      ctplqt(2, 1, 0, 1, A, 2, B, 1, T, 1, W, INFO );
      chkxer('CTPLQT', INFOT, NOUT, LERR, OK );
      INFOT = 10;
      ctplqt(2, 2, 1, 2, A, 2, B, 2, T, 1, W, INFO );
      chkxer('CTPLQT', INFOT, NOUT, LERR, OK );

      // CTPLQT2

     srnamc.SRNAMT = 'CTPLQT2';
      INFOT = 1;
      ctplqt2(-1, 0, 0, A, 1, B, 1, T, 1, INFO );
      chkxer('CTPLQT2', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      ctplqt2(0, -1, 0, A, 1, B, 1, T, 1, INFO );
      chkxer('CTPLQT2', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      ctplqt2(0, 0, -1, A, 1, B, 1, T, 1, INFO );
      chkxer('CTPLQT2', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      ctplqt2(2, 2, 0, A, 1, B, 2, T, 2, INFO );
      chkxer('CTPLQT2', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      ctplqt2(2, 2, 0, A, 2, B, 1, T, 2, INFO );
      chkxer('CTPLQT2', INFOT, NOUT, LERR, OK );
      INFOT = 9;
      ctplqt2(2, 2, 0, A, 2, B, 2, T, 1, INFO );
      chkxer('CTPLQT2', INFOT, NOUT, LERR, OK );

      // CTPMLQT

     srnamc.SRNAMT = 'CTPMLQT';
      INFOT = 1;
      ctpmlqt('/', 'N', 0, 0, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('CTPMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      ctpmlqt('L', '/', 0, 0, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('CTPMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      ctpmlqt('L', 'N', -1, 0, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('CTPMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      ctpmlqt('L', 'N', 0, -1, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('CTPMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      ctpmlqt('L', 'N', 0, 0, -1, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      INFOT = 6;
      ctpmlqt('L', 'N', 0, 0, 0, -1, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('CTPMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      ctpmlqt('L', 'N', 0, 0, 0, 0, 0, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('CTPMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 9;
      ctpmlqt('R', 'N', 2, 2, 2, 1, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('CTPMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 11;
      ctpmlqt('R', 'N', 1, 1, 1, 1, 1, A, 1, T, 0, B, 1, C, 1, W, INFO );
      chkxer('CTPMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 13;
      ctpmlqt('L', 'N', 1, 1, 1, 1, 1, A, 1, T, 1, B, 0, C, 1, W, INFO );
      chkxer('CTPMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 15;
      ctpmlqt('L', 'N', 1, 1, 1, 1, 1, A, 1, T, 1, B, 1, C, 0, W, INFO );
      chkxer('CTPMLQT', INFOT, NOUT, LERR, OK );

      // Print a summary line.

      alaesm(PATH, OK, NOUT );

      return;
      }
