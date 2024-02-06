      void serrqrt(PATH, NUNIT ) {
      // IMPLICIT NONE

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
      // EXTERNAL ALAESM, CHKXER, SGEQRT2, SGEQRT3, SGEQRT, SGEMQRT
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
      // INTRINSIC FLOAT

      NOUT = NUNIT;
      WRITE( NOUT, FMT = * );

      // Set the variables to innocuous values.

      for (J = 1; J <= NMAX; J++) {
         for (I = 1; I <= NMAX; I++) {
            A[I][J] = 1.0 / FLOAT( I+J );
            C[I][J] = 1.0 / FLOAT( I+J );
            T[I][J] = 1.0 / FLOAT( I+J );
         }
         W[J] = 0.0;
      }
      OK = true;

      // Error exits for QRT factorization

      // SGEQRT

     srnamc.SRNAMT = 'SGEQRT';
      INFOT = 1;
      sgeqrt(-1, 0, 1, A, 1, T, 1, W, INFO );
      chkxer('SGEQRT', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      sgeqrt(0, -1, 1, A, 1, T, 1, W, INFO );
      chkxer('SGEQRT', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      sgeqrt(0, 0, 0, A, 1, T, 1, W, INFO );
      chkxer('SGEQRT', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      sgeqrt(2, 1, 1, A, 1, T, 1, W, INFO );
      chkxer('SGEQRT', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      sgeqrt(2, 2, 2, A, 2, T, 1, W, INFO );
      chkxer('SGEQRT', INFOT, NOUT, LERR, OK );

      // SGEQRT2

     srnamc.SRNAMT = 'SGEQRT2';
      INFOT = 1;
      sgeqrt2(-1, 0, A, 1, T, 1, INFO );
      chkxer('SGEQRT2', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      sgeqrt2(0, -1, A, 1, T, 1, INFO );
      chkxer('SGEQRT2', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      sgeqrt2(2, 1, A, 1, T, 1, INFO );
      chkxer('SGEQRT2', INFOT, NOUT, LERR, OK );
      INFOT = 6;
      sgeqrt2(2, 2, A, 2, T, 1, INFO );
      chkxer('SGEQRT2', INFOT, NOUT, LERR, OK );

      // SGEQRT3

     srnamc.SRNAMT = 'SGEQRT3';
      INFOT = 1;
      sgeqrt3(-1, 0, A, 1, T, 1, INFO );
      chkxer('SGEQRT3', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      sgeqrt3(0, -1, A, 1, T, 1, INFO );
      chkxer('SGEQRT3', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      sgeqrt3(2, 1, A, 1, T, 1, INFO );
      chkxer('SGEQRT3', INFOT, NOUT, LERR, OK );
      INFOT = 6;
      sgeqrt3(2, 2, A, 2, T, 1, INFO );
      chkxer('SGEQRT3', INFOT, NOUT, LERR, OK );

      // SGEMQRT

     srnamc.SRNAMT = 'SGEMQRT';
      INFOT = 1;
      sgemqrt('/', 'N', 0, 0, 0, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('SGEMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      sgemqrt('L', '/', 0, 0, 0, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('SGEMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      sgemqrt('L', 'N', -1, 0, 0, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('SGEMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      sgemqrt('L', 'N', 0, -1, 0, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('SGEMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      sgemqrt('L', 'N', 0, 0, -1, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('SGEMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      sgemqrt('R', 'N', 0, 0, -1, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('SGEMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 6;
      sgemqrt('L', 'N', 0, 0, 0, 0, A, 1, T, 1, C, 1, W, INFO );
      chkxer('SGEMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 8;
      sgemqrt('R', 'N', 1, 2, 1, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('SGEMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 8;
      sgemqrt('L', 'N', 2, 1, 1, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('SGEMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 10;
      sgemqrt('R', 'N', 1, 1, 1, 1, A, 1, T, 0, C, 1, W, INFO );
      chkxer('SGEMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 12;
      sgemqrt('L', 'N', 1, 1, 1, 1, A, 1, T, 1, C, 0, W, INFO );
      chkxer('SGEMQRT', INFOT, NOUT, LERR, OK );

      // Print a summary line.

      alaesm(PATH, OK, NOUT );

      return;
      }
