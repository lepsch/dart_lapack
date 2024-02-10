      void serrqr(final int PATH, final int NUNIT) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             PATH;
      int                NUNIT;
      // ..

      int                NMAX;
      const              NMAX = 2 ;
      int                I, INFO, J;
      double               A( NMAX, NMAX ), AF( NMAX, NMAX ), B( NMAX ), W( NMAX ), X( NMAX );
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, SGEQR2, SGEQR2P, SGEQRF, SGEQRFP, SORG2R, SORGQR, SORM2R, SORMQR
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

      for (J = 1; J <= NMAX; J++) { // 20
         for (I = 1; I <= NMAX; I++) { // 10
            A[I][J] = 1. / REAL( I+J );
            AF[I][J] = 1. / REAL( I+J );
         } // 10
         B[J] = 0.;
         W[J] = 0.;
         X[J] = 0.;
      } // 20
      OK = true;

      // Error exits for QR factorization

      // SGEQRF

     srnamc.SRNAMT = 'SGEQRF';
      INFOT = 1;
      sgeqrf(-1, 0, A, 1, B, W, 1, INFO );
      chkxer('SGEQRF', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      sgeqrf(0, -1, A, 1, B, W, 1, INFO );
      chkxer('SGEQRF', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      sgeqrf(2, 1, A, 1, B, W, 1, INFO );
      chkxer('SGEQRF', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      sgeqrf(1, 2, A, 1, B, W, 1, INFO );
      chkxer('SGEQRF', INFOT, NOUT, LERR, OK );

      // SGEQRFP

     srnamc.SRNAMT = 'SGEQRFP';
      INFOT = 1;
      sgeqrfp(-1, 0, A, 1, B, W, 1, INFO );
      chkxer('SGEQRFP', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      sgeqrfp(0, -1, A, 1, B, W, 1, INFO );
      chkxer('SGEQRFP', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      sgeqrfp(2, 1, A, 1, B, W, 1, INFO );
      chkxer('SGEQRFP', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      sgeqrfp(1, 2, A, 1, B, W, 1, INFO );
      chkxer('SGEQRFP', INFOT, NOUT, LERR, OK );

      // SGEQR2

     srnamc.SRNAMT = 'SGEQR2';
      INFOT = 1;
      sgeqr2(-1, 0, A, 1, B, W, INFO );
      chkxer('SGEQR2', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      sgeqr2(0, -1, A, 1, B, W, INFO );
      chkxer('SGEQR2', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      sgeqr2(2, 1, A, 1, B, W, INFO );
      chkxer('SGEQR2', INFOT, NOUT, LERR, OK );

      // SGEQR2P

     srnamc.SRNAMT = 'SGEQR2P';
      INFOT = 1;
      sgeqr2p(-1, 0, A, 1, B, W, INFO );
      chkxer('SGEQR2P', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      sgeqr2p(0, -1, A, 1, B, W, INFO );
      chkxer('SGEQR2P', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      sgeqr2p(2, 1, A, 1, B, W, INFO );
      chkxer('SGEQR2P', INFOT, NOUT, LERR, OK );

      // SORGQR

     srnamc.SRNAMT = 'SORGQR';
      INFOT = 1;
      sorgqr(-1, 0, 0, A, 1, X, W, 1, INFO );
      chkxer('SORGQR', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      sorgqr(0, -1, 0, A, 1, X, W, 1, INFO );
      chkxer('SORGQR', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      sorgqr(1, 2, 0, A, 1, X, W, 2, INFO );
      chkxer('SORGQR', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      sorgqr(0, 0, -1, A, 1, X, W, 1, INFO );
      chkxer('SORGQR', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      sorgqr(1, 1, 2, A, 1, X, W, 1, INFO );
      chkxer('SORGQR', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      sorgqr(2, 2, 0, A, 1, X, W, 2, INFO );
      chkxer('SORGQR', INFOT, NOUT, LERR, OK );
      INFOT = 8;
      sorgqr(2, 2, 0, A, 2, X, W, 1, INFO );
      chkxer('SORGQR', INFOT, NOUT, LERR, OK );

      // SORG2R

     srnamc.SRNAMT = 'SORG2R';
      INFOT = 1;
      sorg2r(-1, 0, 0, A, 1, X, W, INFO );
      chkxer('SORG2R', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      sorg2r(0, -1, 0, A, 1, X, W, INFO );
      chkxer('SORG2R', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      sorg2r(1, 2, 0, A, 1, X, W, INFO );
      chkxer('SORG2R', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      sorg2r(0, 0, -1, A, 1, X, W, INFO );
      chkxer('SORG2R', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      sorg2r(2, 1, 2, A, 2, X, W, INFO );
      chkxer('SORG2R', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      sorg2r(2, 1, 0, A, 1, X, W, INFO );
      chkxer('SORG2R', INFOT, NOUT, LERR, OK );

      // SORMQR

     srnamc.SRNAMT = 'SORMQR';
      INFOT = 1;
      sormqr('/', 'N', 0, 0, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('SORMQR', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      sormqr('L', '/', 0, 0, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('SORMQR', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      sormqr('L', 'N', -1, 0, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('SORMQR', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      sormqr('L', 'N', 0, -1, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('SORMQR', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      sormqr('L', 'N', 0, 0, -1, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('SORMQR', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      sormqr('L', 'N', 0, 1, 1, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('SORMQR', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      sormqr('R', 'N', 1, 0, 1, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('SORMQR', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      sormqr('L', 'N', 2, 1, 0, A, 1, X, AF, 2, W, 1, INFO );
      chkxer('SORMQR', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      sormqr('R', 'N', 1, 2, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('SORMQR', INFOT, NOUT, LERR, OK );
      INFOT = 10;
      sormqr('L', 'N', 2, 1, 0, A, 2, X, AF, 1, W, 1, INFO );
      chkxer('SORMQR', INFOT, NOUT, LERR, OK );
      INFOT = 12;
      sormqr('L', 'N', 1, 2, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('SORMQR', INFOT, NOUT, LERR, OK );
      INFOT = 12;
      sormqr('R', 'N', 2, 1, 0, A, 1, X, AF, 2, W, 1, INFO );
      chkxer('SORMQR', INFOT, NOUT, LERR, OK );

      // SORM2R

     srnamc.SRNAMT = 'SORM2R';
      INFOT = 1;
      sorm2r('/', 'N', 0, 0, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('SORM2R', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      sorm2r('L', '/', 0, 0, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('SORM2R', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      sorm2r('L', 'N', -1, 0, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('SORM2R', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      sorm2r('L', 'N', 0, -1, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('SORM2R', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      sorm2r('L', 'N', 0, 0, -1, A, 1, X, AF, 1, W, INFO );
      chkxer('SORM2R', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      sorm2r('L', 'N', 0, 1, 1, A, 1, X, AF, 1, W, INFO );
      chkxer('SORM2R', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      sorm2r('R', 'N', 1, 0, 1, A, 1, X, AF, 1, W, INFO );
      chkxer('SORM2R', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      sorm2r('L', 'N', 2, 1, 0, A, 1, X, AF, 2, W, INFO );
      chkxer('SORM2R', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      sorm2r('R', 'N', 1, 2, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('SORM2R', INFOT, NOUT, LERR, OK );
      INFOT = 10;
      sorm2r('L', 'N', 2, 1, 0, A, 2, X, AF, 1, W, INFO );
      chkxer('SORM2R', INFOT, NOUT, LERR, OK );

      // Print a summary line.

      alaesm(PATH, OK, NOUT );

      }
