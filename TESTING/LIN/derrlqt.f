      SUBROUTINE DERRLQT( PATH, NUNIT );
      // IMPLICIT NONE

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
      double             A( NMAX, NMAX ), T( NMAX, NMAX ), W( NMAX ), C( NMAX, NMAX );
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, DGELQT3, DGELQT, DGEMLQT
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
      // INTRINSIC DBLE
      // ..
      // .. Executable Statements ..

      NOUT = NUNIT;
      WRITE( NOUT, FMT = * );

      // Set the variables to innocuous values.

      for (J = 1; J <= NMAX; J++) {
         for (I = 1; I <= NMAX; I++) {
            A( I, J ) = 1.0 / DBLE( I+J );
            C( I, J ) = 1.0 / DBLE( I+J );
            T( I, J ) = 1.0 / DBLE( I+J );
         }
         W( J ) = 0.0;
      }
      OK = true;

      // Error exits for LQT factorization

      // DGELQT

      SRNAMT = 'DGELQT';
      INFOT = 1;
      dgelqt(-1, 0, 1, A, 1, T, 1, W, INFO );
      chkxer('DGELQT', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      dgelqt(0, -1, 1, A, 1, T, 1, W, INFO );
      chkxer('DGELQT', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      dgelqt(0, 0, 0, A, 1, T, 1, W, INFO );
      chkxer('DGELQT', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      dgelqt(2, 1, 1, A, 1, T, 1, W, INFO );
      chkxer('DGELQT', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      dgelqt(2, 2, 2, A, 2, T, 1, W, INFO );
      chkxer('DGELQT', INFOT, NOUT, LERR, OK );

      // DGELQT3

      SRNAMT = 'DGELQT3';
      INFOT = 1;
      dgelqt3(-1, 0, A, 1, T, 1, INFO );
      chkxer('DGELQT3', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      dgelqt3(0, -1, A, 1, T, 1, INFO );
      chkxer('DGELQT3', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      dgelqt3(2, 2, A, 1, T, 1, INFO );
      chkxer('DGELQT3', INFOT, NOUT, LERR, OK );
      INFOT = 6;
      dgelqt3(2, 2, A, 2, T, 1, INFO );
      chkxer('DGELQT3', INFOT, NOUT, LERR, OK );

      // DGEMLQT

      SRNAMT = 'DGEMLQT';
      INFOT = 1;
      dgemlqt('/', 'N', 0, 0, 0, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('DGEMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      dgemlqt('L', '/', 0, 0, 0, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('DGEMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      dgemlqt('L', 'N', -1, 0, 0, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('DGEMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      dgemlqt('L', 'N', 0, -1, 0, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('DGEMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      dgemlqt('L', 'N', 0, 0, -1, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('DGEMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      dgemlqt('R', 'N', 0, 0, -1, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('DGEMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 6;
      dgemlqt('L', 'N', 0, 0, 0, 0, A, 1, T, 1, C, 1, W, INFO );
      chkxer('DGEMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 8;
      dgemlqt('R', 'N', 2, 2, 2, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('DGEMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 8;
      dgemlqt('L', 'N', 2, 2, 2, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('DGEMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 10;
      dgemlqt('R', 'N', 1, 1, 1, 1, A, 1, T, 0, C, 1, W, INFO );
      chkxer('DGEMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 12;
      dgemlqt('L', 'N', 1, 1, 1, 1, A, 1, T, 1, C, 0, W, INFO );
      chkxer('DGEMLQT', INFOT, NOUT, LERR, OK );

      // Print a summary line.

      alaesm(PATH, OK, NOUT );

      return;

      // End of DERRLQT

      }
