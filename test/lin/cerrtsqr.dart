      void cerrtsqr(PATH, NUNIT ) {
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
      int                I, INFO, J, MB, NB;
      // ..
      // .. Local Arrays ..
      COMPLEX            A( NMAX, NMAX ), T( NMAX, NMAX ), W( NMAX ), C( NMAX, NMAX ), TAU(NMAX);
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, CGEQR, CGEMQR, CGELQ, CGEMLQ
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
      // INTRINSIC REAL
      // ..
      // .. Executable Statements ..

      NOUT = NUNIT;
      WRITE( NOUT, FMT = * );

      // Set the variables to innocuous values.

      for (J = 1; J <= NMAX; J++) {
         for (I = 1; I <= NMAX; I++) {
            A( I, J ) = 1.0 / CMPLX( REAL( I+J ), 0.0 );
            C( I, J ) = 1.0 / CMPLX( REAL( I+J ), 0.0 );
            T( I, J ) = 1.0 / CMPLX( REAL( I+J ), 0.0 );
         }
         W( J ) = 0.0;
      }
      OK = true;

      // Error exits for TS factorization

      // CGEQR

      SRNAMT = 'CGEQR';
      INFOT = 1;
      cgeqr(-1, 0, A, 1, TAU, 1, W, 1, INFO );
      chkxer('CGEQR', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      cgeqr(0, -1, A, 1, TAU, 1, W, 1, INFO );
      chkxer('CGEQR', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      cgeqr(1, 1, A, 0, TAU, 1, W, 1, INFO );
      chkxer('CGEQR', INFOT, NOUT, LERR, OK );
      INFOT = 6;
      cgeqr(3, 2, A, 3, TAU, 1, W, 1, INFO );
      chkxer('CGEQR', INFOT, NOUT, LERR, OK );
      INFOT = 8;
      cgeqr(3, 2, A, 3, TAU, 8, W, 0, INFO );
      chkxer('CGEQR', INFOT, NOUT, LERR, OK );

      // CLATSQR

      MB = 1;
      NB = 1;
      SRNAMT = 'CLATSQR';
      INFOT = 1;
      clatsqr(-1, 0, MB, NB, A, 1, TAU, 1, W, 1, INFO );
      chkxer('CLATSQR', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      clatsqr(1, 2, MB, NB, A, 1, TAU, 1, W, 1, INFO );
      chkxer('CLATSQR', INFOT, NOUT, LERR, OK );
      clatsqr(0, -1, MB, NB, A, 1, TAU, 1, W, 1, INFO );
      chkxer('CLATSQR', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      clatsqr(2, 1, -1, NB, A, 2, TAU, 1, W, 1, INFO );
      chkxer('CLATSQR', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      clatsqr(2, 1, MB, 2, A, 2, TAU, 1, W, 1, INFO );
      chkxer('CLATSQR', INFOT, NOUT, LERR, OK );
      INFOT = 6;
      clatsqr(2, 1, MB, NB, A, 1, TAU, 1, W, 1, INFO );
      chkxer('CLATSQR', INFOT, NOUT, LERR, OK );
      INFOT = 8;
      clatsqr(2, 1, MB, NB, A, 2, TAU, 0, W, 1, INFO );
      chkxer('CLATSQR', INFOT, NOUT, LERR, OK );
      INFOT = 10;
      clatsqr(2, 1, MB, NB, A, 2, TAU, 2, W, 0, INFO );
      chkxer('CLATSQR', INFOT, NOUT, LERR, OK );

      // CGEMQR

      TAU(1)=1;
      TAU(2)=1;
      SRNAMT = 'CGEMQR';
      NB=1;
      INFOT = 1;
      cgemqr('/', 'N', 0, 0, 0, A, 1, TAU, 1, C, 1, W, 1,INFO);
      chkxer('CGEMQR', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      cgemqr('L', '/', 0, 0, 0, A, 1, TAU, 1, C, 1, W, 1,INFO);
      chkxer('CGEMQR', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      cgemqr('L', 'N', -1, 0, 0, A, 1, TAU, 1, C, 1, W,1,INFO);
      chkxer('CGEMQR', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      cgemqr('L', 'N', 0, -1, 0, A, 1, TAU, 1, C, 1, W,1,INFO);
      chkxer('CGEMQR', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      cgemqr('L', 'N', 0, 0, -1, A, 1, TAU, 1, C, 1, W,1,INFO);
      chkxer('CGEMQR', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      cgemqr('R', 'N', 0, 0, -1, A, 1, TAU, 1, C, 1, W,1,INFO);
      chkxer('CGEMQR', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      cgemqr('L', 'N', 2, 1, 0, A, 0, TAU, 1, C, 1, W, 1,INFO);
      chkxer('CGEMQR', INFOT, NOUT, LERR, OK );
      INFOT = 9;
      cgemqr('R', 'N', 2, 2, 1, A, 2, TAU, 0, C, 1, W, 1,INFO);
      chkxer('CGEMQR', INFOT, NOUT, LERR, OK );
      INFOT = 9;
      cgemqr('L', 'N', 2, 2, 1, A, 2, TAU, 0, C, 1, W, 1,INFO);
      chkxer('CGEMQR', INFOT, NOUT, LERR, OK );
      INFOT = 11;
      cgemqr('L', 'N', 2, 1, 1, A, 2, TAU, 6, C, 0, W, 1,INFO);
      chkxer('CGEMQR', INFOT, NOUT, LERR, OK );
      INFOT = 13;
      cgemqr('L', 'N', 2, 2, 1, A, 2, TAU, 6, C, 2, W, 0,INFO);
      chkxer('CGEMQR', INFOT, NOUT, LERR, OK );

      // CGELQ

      SRNAMT = 'CGELQ';
      INFOT = 1;
      cgelq(-1, 0, A, 1, TAU, 1, W, 1, INFO );
      chkxer('CGELQ', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      cgelq(0, -1, A, 1, TAU, 1, W, 1, INFO );
      chkxer('CGELQ', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      cgelq(1, 1, A, 0, TAU, 1, W, 1, INFO );
      chkxer('CGELQ', INFOT, NOUT, LERR, OK );
      INFOT = 6;
      cgelq(2, 3, A, 3, TAU, 1, W, 1, INFO );
      chkxer('CGELQ', INFOT, NOUT, LERR, OK );
      INFOT = 8;
      cgelq(2, 3, A, 3, TAU, 8, W, 0, INFO );
      chkxer('CGELQ', INFOT, NOUT, LERR, OK );

      // CLASWLQ

      MB = 1;
      NB = 1;
      SRNAMT = 'CLASWLQ';
      INFOT = 1;
      claswlq(-1, 0, MB, NB, A, 1, TAU, 1, W, 1, INFO );
      chkxer('CLASWLQ', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      claswlq(2, 1, MB, NB, A, 1, TAU, 1, W, 1, INFO );
      chkxer('CLASWLQ', INFOT, NOUT, LERR, OK );
      claswlq(0, -1, MB, NB, A, 1, TAU, 1, W, 1, INFO );
      chkxer('CLASWLQ', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      claswlq(1, 2, -1, NB, A, 1, TAU, 1, W, 1, INFO );
      chkxer('CLASWLQ', INFOT, NOUT, LERR, OK );
      claswlq(1, 1, 2, NB, A, 1, TAU, 1, W, 1, INFO );
      chkxer('CLASWLQ', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      claswlq(1, 2, MB, -1, A, 1, TAU, 1, W, 1, INFO );
      chkxer('CLASWLQ', INFOT, NOUT, LERR, OK );
      INFOT = 6;
      claswlq(1, 2, MB, NB, A, 0, TAU, 1, W, 1, INFO );
      chkxer('CLASWLQ', INFOT, NOUT, LERR, OK );
      INFOT = 8;
      claswlq(1, 2, MB, NB, A, 1, TAU, 0, W, 1, INFO );
      chkxer('CLASWLQ', INFOT, NOUT, LERR, OK );
      INFOT = 10;
      claswlq(1, 2, MB, NB, A, 1, TAU, 1, W, 0, INFO );
      chkxer('CLASWLQ', INFOT, NOUT, LERR, OK );

      // CGEMLQ

      TAU(1)=1;
      TAU(2)=1;
      SRNAMT = 'CGEMLQ';
      NB=1;
      INFOT = 1;
      cgemlq('/', 'N', 0, 0, 0, A, 1, TAU, 1, C, 1, W, 1,INFO);
      chkxer('CGEMLQ', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      cgemlq('L', '/', 0, 0, 0, A, 1, TAU, 1, C, 1, W, 1,INFO);
      chkxer('CGEMLQ', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      cgemlq('L', 'N', -1, 0, 0, A, 1, TAU, 1, C, 1, W,1,INFO);
      chkxer('CGEMLQ', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      cgemlq('L', 'N', 0, -1, 0, A, 1, TAU, 1, C, 1, W,1,INFO);
      chkxer('CGEMLQ', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      cgemlq('L', 'N', 0, 0, -1, A, 1, TAU, 1, C, 1, W,1,INFO);
      chkxer('CGEMLQ', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      cgemlq('R', 'N', 0, 0, -1, A, 1, TAU, 1, C, 1, W,1,INFO);
      chkxer('CGEMLQ', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      cgemlq('L', 'N', 1, 2, 0, A, 0, TAU, 1, C, 1, W, 1,INFO);
      chkxer('CGEMLQ', INFOT, NOUT, LERR, OK );
      INFOT = 9;
      cgemlq('R', 'N', 2, 2, 1, A, 1, TAU, 0, C, 1, W, 1,INFO);
      chkxer('CGEMLQ', INFOT, NOUT, LERR, OK );
      INFOT = 9;
      cgemlq('L', 'N', 2, 2, 1, A, 1, TAU, 0, C, 1, W, 1,INFO);
      chkxer('CGEMLQ', INFOT, NOUT, LERR, OK );
      INFOT = 11;
      cgemlq('L', 'N', 1, 2, 1, A, 1, TAU, 6, C, 0, W, 1,INFO);
      chkxer('CGEMLQ', INFOT, NOUT, LERR, OK );
      INFOT = 13;
      cgemlq('L', 'N', 2, 2, 1, A, 2, TAU, 6, C, 2, W, 0,INFO);
      chkxer('CGEMLQ', INFOT, NOUT, LERR, OK );

      // Print a summary line.

      alaesm(PATH, OK, NOUT );

      return;
      }
