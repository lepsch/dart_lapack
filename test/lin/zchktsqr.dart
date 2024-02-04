      void zchktsqr(THRESH, TSTERR, NM, MVAL, NN, NVAL, NNB, NBVAL, NOUT ) {
      // IMPLICIT NONE

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               TSTERR;
      int                NM, NN, NNB, NOUT;
      double             THRESH;
      // ..
      // .. Array Arguments ..
      int                MVAL( * ), NBVAL( * ), NVAL( * );
      // ..

// =====================================================================

      // .. Parameters ..
      int                NTESTS;
      const              NTESTS = 6 ;
      // ..
      // .. Local Scalars ..
      String             PATH;
      int                I, J, K, T, M, N, NB, NFAIL, NERRS, NRUN, INB, MINMN, MB, IMB;

      // .. Local Arrays ..
      double             RESULT( NTESTS );
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAERH, ALAHD, ALASUM, ZERRTSQR, ZTSQR01, XLAENV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Scalars in Common ..
      bool               LERR, OK;
      String             SRNAMT;
      int                INFOT, NUNIT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / INFOT, NUNIT, OK, LERR
      // COMMON / SRNAMC / SRNAMT
      // ..
      // .. Executable Statements ..

      // Initialize constants

      PATH[1: 1] = 'Z';
      PATH[2: 3] = 'TS';
      NRUN = 0;
      NFAIL = 0;
      NERRS = 0;

      // Test the error exits

      xlaenv(1, 0 );
      xlaenv(2, 0 );
      if (TSTERR) zerrtsqr( PATH, NOUT );
      INFOT = 0;

      // Do for each value of M in MVAL.

      for (I = 1; I <= NM; I++) {
         M = MVAL( I );

         // Do for each value of N in NVAL.

         for (J = 1; J <= NN; J++) {
            N = NVAL( J );
              if (min(M,N) != 0) {
              for (INB = 1; INB <= NNB; INB++) {
                MB = NBVAL( INB );
                  xlaenv(1, MB );
                  for (IMB = 1; IMB <= NNB; IMB++) {
                    NB = NBVAL( IMB );
                    xlaenv(2, NB );

                  // Test ZGEQR and ZGEMQR

                    ztsqr01('TS', M, N, MB, NB, RESULT );

                  // Print information about the tests that did not
                  // pass the threshold.

                    for (T = 1; T <= NTESTS; T++) {
                      if ( RESULT( T ) >= THRESH ) {
                        if (NFAIL == 0 && NERRS == 0) alahd( NOUT, PATH );
                        WRITE( NOUT, FMT = 9999 )M, N, MB, NB, T, RESULT( T );
                        NFAIL = NFAIL + 1;
                      }
                    }
                    NRUN = NRUN + NTESTS;
                  }
              }
              }
         }
      }

      // Do for each value of M in MVAL.

      for (I = 1; I <= NM; I++) {
         M = MVAL( I );

         // Do for each value of N in NVAL.

         for (J = 1; J <= NN; J++) {
            N = NVAL( J );
              if (min(M,N) != 0) {
              for (INB = 1; INB <= NNB; INB++) {
                MB = NBVAL( INB );
                  xlaenv(1, MB );
                  for (IMB = 1; IMB <= NNB; IMB++) {
                    NB = NBVAL( IMB );
                    xlaenv(2, NB );

                  // Test ZGELQ and ZGEMLQ

                    ztsqr01('SW', M, N, MB, NB, RESULT );

                  // Print information about the tests that did not
                  // pass the threshold.

                    for (T = 1; T <= NTESTS; T++) {
                      if ( RESULT( T ) >= THRESH ) {
                        if (NFAIL == 0 && NERRS == 0) alahd( NOUT, PATH );
                           WRITE( NOUT, FMT = 9998 )M, N, MB, NB, T, RESULT( T );
                        NFAIL = NFAIL + 1;
                      }
                    }
                    NRUN = NRUN + NTESTS;
                  }
              }
              }
         }
      }

      // Print a summary of the results.

      alasum(PATH, NOUT, NFAIL, NRUN, NERRS );

 9999 FORMAT( 'TS: M=', I5, ', N=', I5, ', MB=', I5, ', NB=', I5,' test(', I2, ')=', G12.5 );
 9998 FORMAT( 'SW: M=', I5, ', N=', I5, ', MB=', I5, ', NB=', I5,' test(', I2, ')=', G12.5 );
      return;
      }
