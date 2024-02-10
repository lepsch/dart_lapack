      void derrtsqr(final int PATH, final int NUNIT) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             PATH;
      int                NUNIT;
      // ..

      int                NMAX;
      const              NMAX = 2 ;
      int                I, INFO, J, MB, NB;
      double             A( NMAX, NMAX ), T( NMAX, NMAX ), W( NMAX ), C( NMAX, NMAX ), TAU(NMAX*2);
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, DGEQR, DGEMQR, DGELQ, DGEMLQ
      // ..
      // .. Scalars in Common ..
      // bool               infoc.LERR, infoc.OK;
      // String            srnamc.SRNAMT;
      // int                infoc.INFOT, infoc.NOUT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / infoc.INFOT, infoc.NOUT, infoc.OK, infoc.LERR
      // COMMON / SRNAMC /srnamc.SRNAMT
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE

      infoc.NOUT = NUNIT;
      WRITE( infoc.NOUT, FMT = * );

      // Set the variables to innocuous values.

      for (J = 1; J <= NMAX; J++) {
         for (I = 1; I <= NMAX; I++) {
            A[I][J] = 1.0 / (I+J).toDouble();
            C[I][J] = 1.0 / (I+J).toDouble();
            T[I][J] = 1.0 / (I+J).toDouble();
         }
         W[J] = 0.0;
      }
      infoc.OK = true;

      // Error exits for TS factorization

      // DGEQR

     srnamc.SRNAMT = 'DGEQR';
      infoc.INFOT = 1;
      dgeqr(-1, 0, A, 1, TAU, 1, W, 1, INFO );
      chkxer('DGEQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      dgeqr(0, -1, A, 1, TAU, 1, W, 1, INFO );
      chkxer('DGEQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 4;
      dgeqr(1, 1, A, 0, TAU, 1, W, 1, INFO );
      chkxer('DGEQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 6;
      dgeqr(3, 2, A, 3, TAU, 1, W, 1, INFO );
      chkxer('DGEQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 8;
      dgeqr(3, 2, A, 3, TAU, 7, W, 0, INFO );
      chkxer('DGEQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

      // DLATSQR

      MB = 1;
      NB = 1;
     srnamc.SRNAMT = 'DLATSQR';
      infoc.INFOT = 1;
      dlatsqr(-1, 0, MB, NB, A, 1, TAU, 1, W, 1, INFO );
      chkxer('DLATSQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      dlatsqr(1, 2, MB, NB, A, 1, TAU, 1, W, 1, INFO );
      chkxer('DLATSQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      dlatsqr(0, -1, MB, NB, A, 1, TAU, 1, W, 1, INFO );
      chkxer('DLATSQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 3;
      dlatsqr(2, 1, -1, NB, A, 2, TAU, 1, W, 1, INFO );
      chkxer('DLATSQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 4;
      dlatsqr(2, 1, MB, 2, A, 2, TAU, 1, W, 1, INFO );
      chkxer('DLATSQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 6;
      dlatsqr(2, 1, MB, NB, A, 1, TAU, 1, W, 1, INFO );
      chkxer('DLATSQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 8;
      dlatsqr(2, 1, MB, NB, A, 2, TAU, 0, W, 1, INFO );
      chkxer('DLATSQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 10;
      dlatsqr(2, 1, MB, NB, A, 2, TAU, 2, W, 0, INFO );
      chkxer('DLATSQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

      // DGEMQR

      TAU(1)=1;
      TAU(2)=1;
      TAU(3)=1;
      TAU(4)=1;
     srnamc.SRNAMT = 'DGEMQR';
      NB=1;
      infoc.INFOT = 1;
      dgemqr('/', 'N', 0, 0, 0, A, 1, TAU, 1, C, 1, W, 1,INFO);
      chkxer('DGEMQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      dgemqr('L', '/', 0, 0, 0, A, 1, TAU, 1, C, 1, W, 1,INFO);
      chkxer('DGEMQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 3;
      dgemqr('L', 'N', -1, 0, 0, A, 1, TAU, 1, C, 1, W,1,INFO);
      chkxer('DGEMQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 4;
      dgemqr('L', 'N', 0, -1, 0, A, 1, TAU, 1, C, 1, W,1,INFO);
      chkxer('DGEMQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      dgemqr('L', 'N', 0, 0, -1, A, 1, TAU, 1, C, 1, W,1,INFO);
      chkxer('DGEMQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      dgemqr('R', 'N', 0, 0, -1, A, 1, TAU, 1, C, 1, W,1,INFO);
      chkxer('DGEMQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 7;
      dgemqr('L', 'N', 2, 1, 0, A, 0, TAU, 1, C, 1, W, 1,INFO);
      chkxer('DGEMQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 9;
      dgemqr('R', 'N', 2, 2, 1, A, 2, TAU, 0, C, 1, W, 1,INFO);
      chkxer('DGEMQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 9;
      dgemqr('L', 'N', 2, 2, 1, A, 2, TAU, 0, C, 1, W, 1,INFO);
      chkxer('DGEMQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 11;
      dgemqr('L', 'N', 2, 1, 1, A, 2, TAU, 6, C, 0, W, 1,INFO);
      chkxer('DGEMQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 13;
      dgemqr('L', 'N', 2, 2, 1, A, 2, TAU, 6, C, 2, W, 0,INFO);
      chkxer('DGEMQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

      // DGELQ

     srnamc.SRNAMT = 'DGELQ';
      infoc.INFOT = 1;
      dgelq(-1, 0, A, 1, TAU, 1, W, 1, INFO );
      chkxer('DGELQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      dgelq(0, -1, A, 1, TAU, 1, W, 1, INFO );
      chkxer('DGELQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 4;
      dgelq(1, 1, A, 0, TAU, 1, W, 1, INFO );
      chkxer('DGELQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 6;
      dgelq(2, 3, A, 3, TAU, 1, W, 1, INFO );
      chkxer('DGELQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 8;
      dgelq(2, 3, A, 3, TAU, 7, W, 0, INFO );
      chkxer('DGELQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

      // DLASWLQ

      MB = 1;
      NB = 1;
     srnamc.SRNAMT = 'DLASWLQ';
      infoc.INFOT = 1;
      dlaswlq(-1, 0, MB, NB, A, 1, TAU, 1, W, 1, INFO );
      chkxer('DLASWLQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      dlaswlq(2, 1, MB, NB, A, 1, TAU, 1, W, 1, INFO );
      chkxer('DLASWLQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      dlaswlq(0, -1, MB, NB, A, 1, TAU, 1, W, 1, INFO );
      chkxer('DLASWLQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 3;
      dlaswlq(1, 2, -1, NB, A, 1, TAU, 1, W, 1, INFO );
      chkxer('DLASWLQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      dlaswlq(1, 1, 2, NB, A, 1, TAU, 1, W, 1, INFO );
      chkxer('DLASWLQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 4;
      dlaswlq(1, 2, MB, -1, A, 1, TAU, 1, W, 1, INFO );
      chkxer('DLASWLQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 6;
      dlaswlq(1, 2, MB, NB, A, 0, TAU, 1, W, 1, INFO );
      chkxer('DLASWLQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 8;
      dlaswlq(1, 2, MB, NB, A, 1, TAU, 0, W, 1, INFO );
      chkxer('DLASWLQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 10;
      dlaswlq(1, 2, MB, NB, A, 1, TAU, 1, W, 0, INFO );
      chkxer('DLASWLQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

      // DGEMLQ

      TAU(1)=1;
      TAU(2)=1;
     srnamc.SRNAMT = 'DGEMLQ';
      NB=1;
      infoc.INFOT = 1;
      dgemlq('/', 'N', 0, 0, 0, A, 1, TAU, 1, C, 1, W, 1,INFO);
      chkxer('DGEMLQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      dgemlq('L', '/', 0, 0, 0, A, 1, TAU, 1, C, 1, W, 1,INFO);
      chkxer('DGEMLQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 3;
      dgemlq('L', 'N', -1, 0, 0, A, 1, TAU, 1, C, 1, W,1,INFO);
      chkxer('DGEMLQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 4;
      dgemlq('L', 'N', 0, -1, 0, A, 1, TAU, 1, C, 1, W,1,INFO);
      chkxer('DGEMLQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      dgemlq('L', 'N', 0, 0, -1, A, 1, TAU, 1, C, 1, W,1,INFO);
      chkxer('DGEMLQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      dgemlq('R', 'N', 0, 0, -1, A, 1, TAU, 1, C, 1, W,1,INFO);
      chkxer('DGEMLQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 7;
      dgemlq('L', 'N', 1, 2, 0, A, 0, TAU, 1, C, 1, W, 1,INFO);
      chkxer('DGEMLQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 9;
      dgemlq('R', 'N', 2, 2, 1, A, 1, TAU, 0, C, 1, W, 1,INFO);
      chkxer('DGEMLQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 9;
      dgemlq('L', 'N', 2, 2, 1, A, 1, TAU, 0, C, 1, W, 1,INFO);
      chkxer('DGEMLQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 11;
      dgemlq('L', 'N', 1, 2, 1, A, 1, TAU, 6, C, 0, W, 1,INFO);
      chkxer('DGEMLQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 13;
      dgemlq('L', 'N', 2, 2, 1, A, 2, TAU, 6, C, 2, W, 0,INFO);
      chkxer('DGEMLQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

      // Print a summary line.

      alaesm(PATH, infoc.OK, infoc.NOUT );

      }
