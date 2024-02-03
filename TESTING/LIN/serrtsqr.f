      SUBROUTINE SERRTSQR( PATH, NUNIT );
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
      REAL               A( NMAX, NMAX ), T( NMAX, NMAX ), W( NMAX ), C( NMAX, NMAX ), TAU(NMAX*2);
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, SGEQR, SGEMQR, SGELQ, SGEMLQ
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
            A( I, J ) = 1. / REAL( I+J );
            C( I, J ) = 1. / REAL( I+J );
            T( I, J ) = 1. / REAL( I+J );
         }
         W( J ) = 0.;
      }
      OK = true;

      // Error exits for TS factorization

      // SGEQR

      SRNAMT = 'SGEQR';
      INFOT = 1;
      sgeqr(-1, 0, A, 1, TAU, 1, W, 1, INFO );
      chkxer('SGEQR', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      sgeqr(0, -1, A, 1, TAU, 1, W, 1, INFO );
      chkxer('SGEQR', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      sgeqr(1, 1, A, 0, TAU, 1, W, 1, INFO );
      chkxer('SGEQR', INFOT, NOUT, LERR, OK );
      INFOT = 6;
      sgeqr(3, 2, A, 3, TAU, 1, W, 1, INFO );
      chkxer('SGEQR', INFOT, NOUT, LERR, OK );
      INFOT = 8;
      sgeqr(3, 2, A, 3, TAU, 7, W, 0, INFO );
      chkxer('SGEQR', INFOT, NOUT, LERR, OK );

      // SLATSQR

      MB = 1;
      NB = 1;
      SRNAMT = 'SLATSQR';
      INFOT = 1;
      slatsqr(-1, 0, MB, NB, A, 1, TAU, 1, W, 1, INFO );
      chkxer('SLATSQR', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      slatsqr(1, 2, MB, NB, A, 1, TAU, 1, W, 1, INFO );
      chkxer('SLATSQR', INFOT, NOUT, LERR, OK );
      slatsqr(0, -1, MB, NB, A, 1, TAU, 1, W, 1, INFO );
      chkxer('SLATSQR', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      slatsqr(2, 1, -1, NB, A, 2, TAU, 1, W, 1, INFO );
      chkxer('SLATSQR', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      slatsqr(2, 1, MB, 2, A, 2, TAU, 1, W, 1, INFO );
      chkxer('SLATSQR', INFOT, NOUT, LERR, OK );
      INFOT = 6;
      slatsqr(2, 1, MB, NB, A, 1, TAU, 1, W, 1, INFO );
      chkxer('SLATSQR', INFOT, NOUT, LERR, OK );
      INFOT = 8;
      slatsqr(2, 1, MB, NB, A, 2, TAU, 0, W, 1, INFO );
      chkxer('SLATSQR', INFOT, NOUT, LERR, OK );
      INFOT = 10;
      slatsqr(2, 1, MB, NB, A, 2, TAU, 2, W, 0, INFO );
      chkxer('SLATSQR', INFOT, NOUT, LERR, OK );

      // SGEMQR

      TAU(1)=1;
      TAU(2)=1;
      TAU(3)=1;
      TAU(4)=1;
      SRNAMT = 'SGEMQR';
      NB=1;
      INFOT = 1;
      sgemqr('/', 'N', 0, 0, 0, A, 1, TAU, 1, C, 1, W, 1,INFO);
      chkxer('SGEMQR', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      sgemqr('L', '/', 0, 0, 0, A, 1, TAU, 1, C, 1, W, 1,INFO);
      chkxer('SGEMQR', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      sgemqr('L', 'N', -1, 0, 0, A, 1, TAU, 1, C, 1, W,1,INFO);
      chkxer('SGEMQR', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      sgemqr('L', 'N', 0, -1, 0, A, 1, TAU, 1, C, 1, W,1,INFO);
      chkxer('SGEMQR', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      sgemqr('L', 'N', 0, 0, -1, A, 1, TAU, 1, C, 1, W,1,INFO);
      chkxer('SGEMQR', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      sgemqr('R', 'N', 0, 0, -1, A, 1, TAU, 1, C, 1, W,1,INFO);
      chkxer('SGEMQR', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      sgemqr('L', 'N', 2, 1, 0, A, 0, TAU, 1, C, 1, W, 1,INFO);
      chkxer('SGEMQR', INFOT, NOUT, LERR, OK );
      INFOT = 9;
      sgemqr('R', 'N', 2, 2, 1, A, 2, TAU, 0, C, 1, W, 1,INFO);
      chkxer('SGEMQR', INFOT, NOUT, LERR, OK );
      INFOT = 9;
      sgemqr('L', 'N', 2, 2, 1, A, 2, TAU, 0, C, 1, W, 1,INFO);
      chkxer('SGEMQR', INFOT, NOUT, LERR, OK );
      INFOT = 11;
      sgemqr('L', 'N', 2, 1, 1, A, 2, TAU, 6, C, 0, W, 1,INFO);
      chkxer('SGEMQR', INFOT, NOUT, LERR, OK );
      INFOT = 13;
      sgemqr('L', 'N', 2, 2, 1, A, 2, TAU, 6, C, 2, W, 0,INFO);
      chkxer('SGEMQR', INFOT, NOUT, LERR, OK );

      // SGELQ

      SRNAMT = 'SGELQ';
      INFOT = 1;
      sgelq(-1, 0, A, 1, TAU, 1, W, 1, INFO );
      chkxer('SGELQ', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      sgelq(0, -1, A, 1, TAU, 1, W, 1, INFO );
      chkxer('SGELQ', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      sgelq(1, 1, A, 0, TAU, 1, W, 1, INFO );
      chkxer('SGELQ', INFOT, NOUT, LERR, OK );
      INFOT = 6;
      sgelq(2, 3, A, 3, TAU, 1, W, 1, INFO );
      chkxer('SGELQ', INFOT, NOUT, LERR, OK );
      INFOT = 8;
      sgelq(2, 3, A, 3, TAU, 7, W, 0, INFO );
      chkxer('SGELQ', INFOT, NOUT, LERR, OK );

      // SLASWLQ

      MB = 1;
      NB = 1;
      SRNAMT = 'SLASWLQ';
      INFOT = 1;
      slaswlq(-1, 0, MB, NB, A, 1, TAU, 1, W, 1, INFO );
      chkxer('SLASWLQ', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      slaswlq(2, 1, MB, NB, A, 1, TAU, 1, W, 1, INFO );
      chkxer('SLASWLQ', INFOT, NOUT, LERR, OK );
      slaswlq(0, -1, MB, NB, A, 1, TAU, 1, W, 1, INFO );
      chkxer('SLASWLQ', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      slaswlq(1, 2, -1, NB, A, 1, TAU, 1, W, 1, INFO );
      chkxer('SLASWLQ', INFOT, NOUT, LERR, OK );
      slaswlq(1, 1, 2, NB, A, 1, TAU, 1, W, 1, INFO );
      chkxer('SLASWLQ', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      slaswlq(1, 2, MB, -1, A, 1, TAU, 1, W, 1, INFO );
      chkxer('SLASWLQ', INFOT, NOUT, LERR, OK );
      INFOT = 6;
      slaswlq(1, 2, MB, NB, A, 0, TAU, 1, W, 1, INFO );
      chkxer('SLASWLQ', INFOT, NOUT, LERR, OK );
      INFOT = 8;
      slaswlq(1, 2, MB, NB, A, 1, TAU, 0, W, 1, INFO );
      chkxer('SLASWLQ', INFOT, NOUT, LERR, OK );
      INFOT = 10;
      slaswlq(1, 2, MB, NB, A, 1, TAU, 1, W, 0, INFO );
      chkxer('SLASWLQ', INFOT, NOUT, LERR, OK );

      // SGEMLQ

      TAU(1)=1;
      TAU(2)=1;
      SRNAMT = 'SGEMLQ';
      NB=1;
      INFOT = 1;
      sgemlq('/', 'N', 0, 0, 0, A, 1, TAU, 1, C, 1, W, 1,INFO);
      chkxer('SGEMLQ', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      sgemlq('L', '/', 0, 0, 0, A, 1, TAU, 1, C, 1, W, 1,INFO);
      chkxer('SGEMLQ', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      sgemlq('L', 'N', -1, 0, 0, A, 1, TAU, 1, C, 1, W,1,INFO);
      chkxer('SGEMLQ', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      sgemlq('L', 'N', 0, -1, 0, A, 1, TAU, 1, C, 1, W,1,INFO);
      chkxer('SGEMLQ', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      sgemlq('L', 'N', 0, 0, -1, A, 1, TAU, 1, C, 1, W,1,INFO);
      chkxer('SGEMLQ', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      sgemlq('R', 'N', 0, 0, -1, A, 1, TAU, 1, C, 1, W,1,INFO);
      chkxer('SGEMLQ', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      sgemlq('L', 'N', 1, 2, 0, A, 0, TAU, 1, C, 1, W, 1,INFO);
      chkxer('SGEMLQ', INFOT, NOUT, LERR, OK );
      INFOT = 9;
      sgemlq('R', 'N', 2, 2, 1, A, 1, TAU, 0, C, 1, W, 1,INFO);
      chkxer('SGEMLQ', INFOT, NOUT, LERR, OK );
      INFOT = 9;
      sgemlq('L', 'N', 2, 2, 1, A, 1, TAU, 0, C, 1, W, 1,INFO);
      chkxer('SGEMLQ', INFOT, NOUT, LERR, OK );
      INFOT = 11;
      sgemlq('L', 'N', 1, 2, 1, A, 1, TAU, 6, C, 0, W, 1,INFO);
      chkxer('SGEMLQ', INFOT, NOUT, LERR, OK );
      INFOT = 13;
      sgemlq('L', 'N', 2, 2, 1, A, 2, TAU, 6, C, 2, W, 0,INFO);
      chkxer('SGEMLQ', INFOT, NOUT, LERR, OK );

      // Print a summary line.

      alaesm(PATH, OK, NOUT );

      return;

      // End of SERRTSQR

      }
