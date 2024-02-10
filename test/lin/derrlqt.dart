import 'common.dart';

      void derrlqt(PATH, final int NUNIT) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             PATH;
      int                NUNIT;
      // ..

      int                NMAX;
      const              NMAX = 2 ;
      int                I, INFO, J;
      double             A( NMAX, NMAX ), T( NMAX, NMAX ), W( NMAX ), C( NMAX, NMAX );
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, DGELQT3, DGELQT, DGEMLQT
      // ..
      // .. Scalars in Common ..
      // bool               infoc.LERR, infoc.OK;
      // String             srnamc.SRNAMT;
      // int                infoc.INFOT, infoc.NOUT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / infoc.INFOT, infoc.NOUT, infoc.OK, infoc.LERR
      // COMMON / SRNAMC / srnamc.SRNAMT
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

      // Error exits for LQT factorization

      // DGELQT

      srnamc.SRNAMT = 'DGELQT';
      infoc.INFOT = 1;
      dgelqt(-1, 0, 1, A, 1, T, 1, W, INFO );
      chkxer('DGELQT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      dgelqt(0, -1, 1, A, 1, T, 1, W, INFO );
      chkxer('DGELQT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 3;
      dgelqt(0, 0, 0, A, 1, T, 1, W, INFO );
      chkxer('DGELQT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      dgelqt(2, 1, 1, A, 1, T, 1, W, INFO );
      chkxer('DGELQT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 7;
      dgelqt(2, 2, 2, A, 2, T, 1, W, INFO );
      chkxer('DGELQT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

      // DGELQT3

      srnamc.SRNAMT = 'DGELQT3';
      infoc.INFOT = 1;
      dgelqt3(-1, 0, A, 1, T, 1, INFO );
      chkxer('DGELQT3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      dgelqt3(0, -1, A, 1, T, 1, INFO );
      chkxer('DGELQT3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 4;
      dgelqt3(2, 2, A, 1, T, 1, INFO );
      chkxer('DGELQT3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 6;
      dgelqt3(2, 2, A, 2, T, 1, INFO );
      chkxer('DGELQT3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

      // DGEMLQT

      srnamc.SRNAMT = 'DGEMLQT';
      infoc.INFOT = 1;
      dgemlqt('/', 'N', 0, 0, 0, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('DGEMLQT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      dgemlqt('L', '/', 0, 0, 0, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('DGEMLQT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 3;
      dgemlqt('L', 'N', -1, 0, 0, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('DGEMLQT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 4;
      dgemlqt('L', 'N', 0, -1, 0, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('DGEMLQT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      dgemlqt('L', 'N', 0, 0, -1, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('DGEMLQT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      dgemlqt('R', 'N', 0, 0, -1, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('DGEMLQT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 6;
      dgemlqt('L', 'N', 0, 0, 0, 0, A, 1, T, 1, C, 1, W, INFO );
      chkxer('DGEMLQT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 8;
      dgemlqt('R', 'N', 2, 2, 2, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('DGEMLQT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 8;
      dgemlqt('L', 'N', 2, 2, 2, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('DGEMLQT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 10;
      dgemlqt('R', 'N', 1, 1, 1, 1, A, 1, T, 0, C, 1, W, INFO );
      chkxer('DGEMLQT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 12;
      dgemlqt('L', 'N', 1, 1, 1, 1, A, 1, T, 1, C, 0, W, INFO );
      chkxer('DGEMLQT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

      // Print a summary line.

      alaesm(PATH, infoc.OK, infoc.NOUT );

      }
