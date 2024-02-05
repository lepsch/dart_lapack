      void zerrtsqr(PATH, infoc.NUNIT ) {
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
      int                I, INFO, J, MB, NB;
      // ..
      // .. Local Arrays ..
      Complex         A( NMAX, NMAX ), T( NMAX, NMAX ), W( NMAX ), C( NMAX, NMAX ), TAU(NMAX);
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, ZGEQR, ZGEMQR, ZGELQ, ZGEMLQ
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
      // INTRINSIC DBLE
      // ..
      // .. Executable Statements ..

      NOUT = infoc.NUNIT;
      WRITE( NOUT, FMT = * );

      // Set the variables to innocuous values.

      for (J = 1; J <= NMAX; J++) {
         for (I = 1; I <= NMAX; I++) {
            A[I, J] = 1.0 / (I+J).toDouble();
            C[I, J] = 1.0 / (I+J).toDouble();
            T[I, J] = 1.0 / (I+J).toDouble();
         }
         W[J] = 0.0;
      }
      infoc.OK = true;

      // Error exits for TS factorization

      // ZGEQR

     srnamc.SRNAMT = 'ZGEQR';
      infoc.INFOT = 1;
      zgeqr(-1, 0, A, 1, TAU, 1, W, 1, INFO );
      chkxer('ZGEQR', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      zgeqr(0, -1, A, 1, TAU, 1, W, 1, INFO );
      chkxer('ZGEQR', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 4;
      zgeqr(1, 1, A, 0, TAU, 1, W, 1, INFO );
      chkxer('ZGEQR', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 6;
      zgeqr(3, 2, A, 3, TAU, 1, W, 1, INFO );
      chkxer('ZGEQR', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 8;
      zgeqr(3, 2, A, 3, TAU, 8, W, 0, INFO );
      chkxer('ZGEQR', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

      // ZLATSQR

      MB = 1;
      NB = 1;
     srnamc.SRNAMT = 'ZLATSQR';
      infoc.INFOT = 1;
      zlatsqr(-1, 0, MB, NB, A, 1, TAU, 1, W, 1, INFO );
      chkxer('ZLATSQR', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      zlatsqr(1, 2, MB, NB, A, 1, TAU, 1, W, 1, INFO );
      chkxer('ZLATSQR', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      zlatsqr(0, -1, MB, NB, A, 1, TAU, 1, W, 1, INFO );
      chkxer('ZLATSQR', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 3;
      zlatsqr(2, 1, -1, NB, A, 2, TAU, 1, W, 1, INFO );
      chkxer('ZLATSQR', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 4;
      zlatsqr(2, 1, MB, 2, A, 2, TAU, 1, W, 1, INFO );
      chkxer('ZLATSQR', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 6;
      zlatsqr(2, 1, MB, NB, A, 1, TAU, 1, W, 1, INFO );
      chkxer('ZLATSQR', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 8;
      zlatsqr(2, 1, MB, NB, A, 2, TAU, 0, W, 1, INFO );
      chkxer('ZLATSQR', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 10;
      zlatsqr(2, 1, MB, NB, A, 2, TAU, 2, W, 0, INFO );
      chkxer('ZLATSQR', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

      // ZGEMQR

      TAU(1)=1;
      TAU(2)=1;
     srnamc.SRNAMT = 'ZGEMQR';
      NB=1;
      infoc.INFOT = 1;
      zgemqr('/', 'N', 0, 0, 0, A, 1, TAU, 1, C, 1, W, 1,INFO);
      chkxer('ZGEMQR', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      zgemqr('L', '/', 0, 0, 0, A, 1, TAU, 1, C, 1, W, 1,INFO);
      chkxer('ZGEMQR', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 3;
      zgemqr('L', 'N', -1, 0, 0, A, 1, TAU, 1, C, 1, W,1,INFO);
      chkxer('ZGEMQR', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 4;
      zgemqr('L', 'N', 0, -1, 0, A, 1, TAU, 1, C, 1, W,1,INFO);
      chkxer('ZGEMQR', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      zgemqr('L', 'N', 0, 0, -1, A, 1, TAU, 1, C, 1, W,1,INFO);
      chkxer('ZGEMQR', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      zgemqr('R', 'N', 0, 0, -1, A, 1, TAU, 1, C, 1, W,1,INFO);
      chkxer('ZGEMQR', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 7;
      zgemqr('L', 'N', 2, 1, 0, A, 0, TAU, 1, C, 1, W, 1,INFO);
      chkxer('ZGEMQR', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 9;
      zgemqr('R', 'N', 2, 2, 1, A, 2, TAU, 0, C, 1, W, 1,INFO);
      chkxer('ZGEMQR', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 9;
      zgemqr('L', 'N', 2, 2, 1, A, 2, TAU, 0, C, 1, W, 1,INFO);
      chkxer('ZGEMQR', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 11;
      zgemqr('L', 'N', 2, 1, 1, A, 2, TAU, 6, C, 0, W, 1,INFO);
      chkxer('ZGEMQR', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 13;
      zgemqr('L', 'N', 2, 2, 1, A, 2, TAU, 6, C, 2, W, 0,INFO);
      chkxer('ZGEMQR', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

      // ZGELQ

     srnamc.SRNAMT = 'ZGELQ';
      infoc.INFOT = 1;
      zgelq(-1, 0, A, 1, TAU, 1, W, 1, INFO );
      chkxer('ZGELQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      zgelq(0, -1, A, 1, TAU, 1, W, 1, INFO );
      chkxer('ZGELQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 4;
      zgelq(1, 1, A, 0, TAU, 1, W, 1, INFO );
      chkxer('ZGELQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 6;
      zgelq(2, 3, A, 3, TAU, 1, W, 1, INFO );
      chkxer('ZGELQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 8;
      zgelq(2, 3, A, 3, TAU, 8, W, 0, INFO );
      chkxer('ZGELQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

      // ZLASWLQ

      MB = 1;
      NB = 1;
     srnamc.SRNAMT = 'ZLASWLQ';
      infoc.INFOT = 1;
      zlaswlq(-1, 0, MB, NB, A, 1, TAU, 1, W, 1, INFO );
      chkxer('ZLASWLQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      zlaswlq(2, 1, MB, NB, A, 1, TAU, 1, W, 1, INFO );
      chkxer('ZLASWLQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      zlaswlq(0, -1, MB, NB, A, 1, TAU, 1, W, 1, INFO );
      chkxer('ZLASWLQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 3;
      zlaswlq(1, 2, -1, NB, A, 1, TAU, 1, W, 1, INFO );
      chkxer('ZLASWLQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      zlaswlq(1, 1, 2, NB, A, 1, TAU, 1, W, 1, INFO );
      chkxer('ZLASWLQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 4;
      zlaswlq(1, 2, MB, -1, A, 1, TAU, 1, W, 1, INFO );
      chkxer('ZLASWLQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 6;
      zlaswlq(1, 2, MB, NB, A, 0, TAU, 1, W, 1, INFO );
      chkxer('ZLASWLQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 8;
      zlaswlq(1, 2, MB, NB, A, 1, TAU, 0, W, 1, INFO );
      chkxer('ZLASWLQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 10;
      zlaswlq(1, 2, MB, NB, A, 1, TAU, 1, W, 0, INFO );
      chkxer('ZLASWLQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

      // ZGEMLQ

      TAU(1)=1;
      TAU(2)=1;
     srnamc.SRNAMT = 'ZGEMLQ';
      NB=1;
      infoc.INFOT = 1;
      zgemlq('/', 'N', 0, 0, 0, A, 1, TAU, 1, C, 1, W, 1,INFO);
      chkxer('ZGEMLQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      zgemlq('L', '/', 0, 0, 0, A, 1, TAU, 1, C, 1, W, 1,INFO);
      chkxer('ZGEMLQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 3;
      zgemlq('L', 'N', -1, 0, 0, A, 1, TAU, 1, C, 1, W,1,INFO);
      chkxer('ZGEMLQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 4;
      zgemlq('L', 'N', 0, -1, 0, A, 1, TAU, 1, C, 1, W,1,INFO);
      chkxer('ZGEMLQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      zgemlq('L', 'N', 0, 0, -1, A, 1, TAU, 1, C, 1, W,1,INFO);
      chkxer('ZGEMLQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      zgemlq('R', 'N', 0, 0, -1, A, 1, TAU, 1, C, 1, W,1,INFO);
      chkxer('ZGEMLQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 7;
      zgemlq('L', 'N', 1, 2, 0, A, 0, TAU, 1, C, 1, W, 1,INFO);
      chkxer('ZGEMLQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 9;
      zgemlq('R', 'N', 2, 2, 1, A, 1, TAU, 0, C, 1, W, 1,INFO);
      chkxer('ZGEMLQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 9;
      zgemlq('L', 'N', 2, 2, 1, A, 1, TAU, 0, C, 1, W, 1,INFO);
      chkxer('ZGEMLQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 11;
      zgemlq('L', 'N', 1, 2, 1, A, 1, TAU, 6, C, 0, W, 1,INFO);
      chkxer('ZGEMLQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 13;
      zgemlq('L', 'N', 2, 2, 1, A, 2, TAU, 6, C, 2, W, 0,INFO);
      chkxer('ZGEMLQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

      // Print a summary line.

      alaesm(PATH, infoc.OK, NOUT );

      return;
      }
