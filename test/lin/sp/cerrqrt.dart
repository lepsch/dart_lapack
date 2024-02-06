      void cerrqrt(PATH, NUNIT ) {
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
      Complex            A( NMAX, NMAX ), T( NMAX, NMAX ), W( NMAX ), C( NMAX, NMAX );
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, CGEQRT2, CGEQRT3, CGEQRT, CGEMQRT
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
      // INTRINSIC FLOAT, CMPLX

      NOUT = NUNIT;
      WRITE( NOUT, FMT = * );

      // Set the variables to innocuous values.

      for (J = 1; J <= NMAX; J++) {
         for (I = 1; I <= NMAX; I++) {
            A[I][J] = 1.0 / CMPLX( FLOAT(I+J), 0.0 );
            C[I][J] = 1.0 / CMPLX( FLOAT(I+J), 0.0 );
            T[I][J] = 1.0 / CMPLX( FLOAT(I+J), 0.0 );
         }
         W[J] = 0.0;
      }
      OK = true;

      // Error exits for QRT factorization

      // CGEQRT

     srnamc.SRNAMT = 'CGEQRT';
      INFOT = 1;
      cgeqrt(-1, 0, 1, A, 1, T, 1, W, INFO );
      chkxer('CGEQRT', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      cgeqrt(0, -1, 1, A, 1, T, 1, W, INFO );
      chkxer('CGEQRT', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      cgeqrt(0, 0, 0, A, 1, T, 1, W, INFO );
      chkxer('CGEQRT', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      cgeqrt(2, 1, 1, A, 1, T, 1, W, INFO );
      chkxer('CGEQRT', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      cgeqrt(2, 2, 2, A, 2, T, 1, W, INFO );
      chkxer('CGEQRT', INFOT, NOUT, LERR, OK );

      // CGEQRT2

     srnamc.SRNAMT = 'CGEQRT2';
      INFOT = 1;
      cgeqrt2(-1, 0, A, 1, T, 1, INFO );
      chkxer('CGEQRT2', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      cgeqrt2(0, -1, A, 1, T, 1, INFO );
      chkxer('CGEQRT2', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      cgeqrt2(2, 1, A, 1, T, 1, INFO );
      chkxer('CGEQRT2', INFOT, NOUT, LERR, OK );
      INFOT = 6;
      cgeqrt2(2, 2, A, 2, T, 1, INFO );
      chkxer('CGEQRT2', INFOT, NOUT, LERR, OK );

      // CGEQRT3

     srnamc.SRNAMT = 'CGEQRT3';
      INFOT = 1;
      cgeqrt3(-1, 0, A, 1, T, 1, INFO );
      chkxer('CGEQRT3', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      cgeqrt3(0, -1, A, 1, T, 1, INFO );
      chkxer('CGEQRT3', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      cgeqrt3(2, 1, A, 1, T, 1, INFO );
      chkxer('CGEQRT3', INFOT, NOUT, LERR, OK );
      INFOT = 6;
      cgeqrt3(2, 2, A, 2, T, 1, INFO );
      chkxer('CGEQRT3', INFOT, NOUT, LERR, OK );

      // CGEMQRT

     srnamc.SRNAMT = 'CGEMQRT';
      INFOT = 1;
      cgemqrt('/', 'N', 0, 0, 0, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('CGEMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      cgemqrt('L', '/', 0, 0, 0, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('CGEMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      cgemqrt('L', 'N', -1, 0, 0, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('CGEMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      cgemqrt('L', 'N', 0, -1, 0, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('CGEMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      cgemqrt('L', 'N', 0, 0, -1, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('CGEMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      cgemqrt('R', 'N', 0, 0, -1, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('CGEMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 6;
      cgemqrt('L', 'N', 0, 0, 0, 0, A, 1, T, 1, C, 1, W, INFO );
      chkxer('CGEMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 8;
      cgemqrt('R', 'N', 1, 2, 1, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('CGEMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 8;
      cgemqrt('L', 'N', 2, 1, 1, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('CGEMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 10;
      cgemqrt('R', 'N', 1, 1, 1, 1, A, 1, T, 0, C, 1, W, INFO );
      chkxer('CGEMQRT', INFOT, NOUT, LERR, OK );
      INFOT = 12;
      cgemqrt('L', 'N', 1, 1, 1, 1, A, 1, T, 1, C, 0, W, INFO );
      chkxer('CGEMQRT', INFOT, NOUT, LERR, OK );

      // Print a summary line.

      alaesm(PATH, OK, NOUT );

      return;
      }
