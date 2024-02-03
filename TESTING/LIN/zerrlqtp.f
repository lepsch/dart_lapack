      void zerrlqtp(PATH, NUNIT ) {
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
      Complex         A( NMAX, NMAX ), T( NMAX, NMAX ), W( NMAX ), B( NMAX, NMAX ), C( NMAX, NMAX );
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, ZTPLQT2, ZTPLQT, ZTPMLQT
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
            A( I, J ) = 1.0 / DCMPLX( DBLE( I+J ), 0.0 );
            C( I, J ) = 1.0 / DCMPLX( DBLE( I+J ), 0.0 );
            T( I, J ) = 1.0 / DCMPLX( DBLE( I+J ), 0.0 );
         }
         W( J ) = 0.0;
      }
      OK = true;

      // Error exits for TPLQT factorization

      // ZTPLQT

      SRNAMT = 'ZTPLQT';
      INFOT = 1;
      ztplqt(-1, 1, 0, 1, A, 1, B, 1, T, 1, W, INFO );
      chkxer('ZTPLQT', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      ztplqt(1, -1, 0, 1, A, 1, B, 1, T, 1, W, INFO );
      chkxer('ZTPLQT', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      ztplqt(0, 1, -1, 1, A, 1, B, 1, T, 1, W, INFO );
      chkxer('ZTPLQT', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      ztplqt(0, 1, 1, 1, A, 1, B, 1, T, 1, W, INFO );
      chkxer('ZTPLQT', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      ztplqt(0, 1, 0, 0, A, 1, B, 1, T, 1, W, INFO );
      chkxer('ZTPLQT', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      ztplqt(1, 1, 0, 2, A, 1, B, 1, T, 1, W, INFO );
      chkxer('ZTPLQT', INFOT, NOUT, LERR, OK );
      INFOT = 6;
      ztplqt(2, 1, 0, 2, A, 1, B, 1, T, 1, W, INFO );
      chkxer('ZTPLQT', INFOT, NOUT, LERR, OK );
      INFOT = 8;
      ztplqt(2, 1, 0, 1, A, 2, B, 1, T, 1, W, INFO );
      chkxer('ZTPLQT', INFOT, NOUT, LERR, OK );
      INFOT = 10;
      ztplqt(2, 2, 1, 2, A, 2, B, 2, T, 1, W, INFO );
      chkxer('ZTPLQT', INFOT, NOUT, LERR, OK );

      // ZTPLQT2

      SRNAMT = 'ZTPLQT2';
      INFOT = 1;
      ztplqt2(-1, 0, 0, A, 1, B, 1, T, 1, INFO );
      chkxer('ZTPLQT2', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      ztplqt2(0, -1, 0, A, 1, B, 1, T, 1, INFO );
      chkxer('ZTPLQT2', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      ztplqt2(0, 0, -1, A, 1, B, 1, T, 1, INFO );
      chkxer('ZTPLQT2', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      ztplqt2(2, 2, 0, A, 1, B, 2, T, 2, INFO );
      chkxer('ZTPLQT2', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      ztplqt2(2, 2, 0, A, 2, B, 1, T, 2, INFO );
      chkxer('ZTPLQT2', INFOT, NOUT, LERR, OK );
      INFOT = 9;
      ztplqt2(2, 2, 0, A, 2, B, 2, T, 1, INFO );
      chkxer('ZTPLQT2', INFOT, NOUT, LERR, OK );

      // ZTPMLQT

      SRNAMT = 'ZTPMLQT';
      INFOT = 1;
      ztpmlqt('/', 'N', 0, 0, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('ZTPMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      ztpmlqt('L', '/', 0, 0, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('ZTPMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      ztpmlqt('L', 'N', -1, 0, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('ZTPMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      ztpmlqt('L', 'N', 0, -1, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('ZTPMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      ztpmlqt('L', 'N', 0, 0, -1, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      INFOT = 6;
      ztpmlqt('L', 'N', 0, 0, 0, -1, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('ZTPMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      ztpmlqt('L', 'N', 0, 0, 0, 0, 0, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('ZTPMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 9;
      ztpmlqt('R', 'N', 2, 2, 2, 1, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('ZTPMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 11;
      ztpmlqt('R', 'N', 1, 1, 1, 1, 1, A, 1, T, 0, B, 1, C, 1, W, INFO );
      chkxer('ZTPMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 13;
      ztpmlqt('L', 'N', 1, 1, 1, 1, 1, A, 1, T, 1, B, 0, C, 1, W, INFO );
      chkxer('ZTPMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 15;
      ztpmlqt('L', 'N', 1, 1, 1, 1, 1, A, 1, T, 1, B, 1, C, 0, W, INFO );
      chkxer('ZTPMLQT', INFOT, NOUT, LERR, OK );

      // Print a summary line.

      alaesm(PATH, OK, NOUT );

      return;
      }
