      void zerrlqtp(PATH, infoc.NUNIT ) {
      // IMPLICIT NONE

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             PATH;
      int                infoc.NUNIT;
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
      bool               infoc.LERR, infoc.OK;
      String            srnamc.SRNAMT;
      int                infoc.INFOT, NOUT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / infoc.INFOT, NOUT, infoc.OK, infoc.LERR
      // COMMON / SRNAMC /srnamc.SRNAMT
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DCMPLX
      // ..
      // .. Executable Statements ..

      NOUT = infoc.NUNIT;
      WRITE( NOUT, FMT = * );

      // Set the variables to innocuous values.

      for (J = 1; J <= NMAX; J++) {
         for (I = 1; I <= NMAX; I++) {
            A[I][J] = 1.0 / DCMPLX( (I+J).toDouble(), 0.0 );
            C[I][J] = 1.0 / DCMPLX( (I+J).toDouble(), 0.0 );
            T[I][J] = 1.0 / DCMPLX( (I+J).toDouble(), 0.0 );
         }
         W[J] = 0.0;
      }
      infoc.OK = true;

      // Error exits for TPLQT factorization

      // ZTPLQT

     srnamc.SRNAMT = 'ZTPLQT';
      infoc.INFOT = 1;
      ztplqt(-1, 1, 0, 1, A, 1, B, 1, T, 1, W, INFO );
      chkxer('ZTPLQT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      ztplqt(1, -1, 0, 1, A, 1, B, 1, T, 1, W, INFO );
      chkxer('ZTPLQT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 3;
      ztplqt(0, 1, -1, 1, A, 1, B, 1, T, 1, W, INFO );
      chkxer('ZTPLQT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 3;
      ztplqt(0, 1, 1, 1, A, 1, B, 1, T, 1, W, INFO );
      chkxer('ZTPLQT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 4;
      ztplqt(0, 1, 0, 0, A, 1, B, 1, T, 1, W, INFO );
      chkxer('ZTPLQT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 4;
      ztplqt(1, 1, 0, 2, A, 1, B, 1, T, 1, W, INFO );
      chkxer('ZTPLQT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 6;
      ztplqt(2, 1, 0, 2, A, 1, B, 1, T, 1, W, INFO );
      chkxer('ZTPLQT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 8;
      ztplqt(2, 1, 0, 1, A, 2, B, 1, T, 1, W, INFO );
      chkxer('ZTPLQT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 10;
      ztplqt(2, 2, 1, 2, A, 2, B, 2, T, 1, W, INFO );
      chkxer('ZTPLQT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

      // ZTPLQT2

     srnamc.SRNAMT = 'ZTPLQT2';
      infoc.INFOT = 1;
      ztplqt2(-1, 0, 0, A, 1, B, 1, T, 1, INFO );
      chkxer('ZTPLQT2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      ztplqt2(0, -1, 0, A, 1, B, 1, T, 1, INFO );
      chkxer('ZTPLQT2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 3;
      ztplqt2(0, 0, -1, A, 1, B, 1, T, 1, INFO );
      chkxer('ZTPLQT2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      ztplqt2(2, 2, 0, A, 1, B, 2, T, 2, INFO );
      chkxer('ZTPLQT2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 7;
      ztplqt2(2, 2, 0, A, 2, B, 1, T, 2, INFO );
      chkxer('ZTPLQT2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 9;
      ztplqt2(2, 2, 0, A, 2, B, 2, T, 1, INFO );
      chkxer('ZTPLQT2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

      // ZTPMLQT

     srnamc.SRNAMT = 'ZTPMLQT';
      infoc.INFOT = 1;
      ztpmlqt('/', 'N', 0, 0, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('ZTPMLQT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      ztpmlqt('L', '/', 0, 0, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('ZTPMLQT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 3;
      ztpmlqt('L', 'N', -1, 0, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('ZTPMLQT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 4;
      ztpmlqt('L', 'N', 0, -1, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('ZTPMLQT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      ztpmlqt('L', 'N', 0, 0, -1, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      infoc.INFOT = 6;
      ztpmlqt('L', 'N', 0, 0, 0, -1, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('ZTPMLQT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 7;
      ztpmlqt('L', 'N', 0, 0, 0, 0, 0, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('ZTPMLQT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 9;
      ztpmlqt('R', 'N', 2, 2, 2, 1, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('ZTPMLQT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 11;
      ztpmlqt('R', 'N', 1, 1, 1, 1, 1, A, 1, T, 0, B, 1, C, 1, W, INFO );
      chkxer('ZTPMLQT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 13;
      ztpmlqt('L', 'N', 1, 1, 1, 1, 1, A, 1, T, 1, B, 0, C, 1, W, INFO );
      chkxer('ZTPMLQT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 15;
      ztpmlqt('L', 'N', 1, 1, 1, 1, 1, A, 1, T, 1, B, 1, C, 0, W, INFO );
      chkxer('ZTPMLQT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

      // Print a summary line.

      alaesm(PATH, infoc.OK, NOUT );

      return;
      }
