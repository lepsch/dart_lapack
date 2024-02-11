import 'common.dart';

      void dchktsqr(final int THRESH, final int TSTERR, final int NM, final int MVAL, final int NN, final int NVAL, final int NNB, final int NBVAL, final int NOUT,) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      bool               TSTERR;
      int                NM, NN, NNB, NOUT;
      double             THRESH;
      int                MVAL( * ), NBVAL( * ), NVAL( * );
      // ..

      int                NTESTS;
      const              NTESTS = 6 ;
      String             PATH;
      int                I, J, K, T, M, N, NB, NFAIL, NERRS, NRUN, INB, MINMN, MB, IMB;

      // .. Local Arrays ..
      double             RESULT( NTESTS );
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAERH, ALAHD, ALASUM, DERRTSQR, DTSQR01, XLAENV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Scalars in Common ..
      // bool               infoc.LERR, infoc.OK;
      // String             srnamc.SRNAMT;
      // int                infoc.INFOT, infoc.NUNIT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / infoc.INFOT, infoc.NUNIT, infoc.OK, infoc.LERR
      // COMMON / SRNAMC / srnamc.SRNAMT

      // Initialize constants

      PATH[1: 1] = 'D';
      PATH[2: 3] = 'TS';
      NRUN = 0;
      NFAIL = 0;
      NERRS = 0;

      // Test the error exits

      xlaenv(1, 0 );
      xlaenv(2, 0 );
      if (TSTERR) derrtsqr( PATH, NOUT );
      infoc.INFOT = 0;

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

                  // Test DGEQR and DGEMQR

                    dtsqr01('TS', M, N, MB, NB, RESULT );

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

                  // Test DGEQR and DGEMQR

                    dtsqr01('SW', M, N, MB, NB, RESULT );

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

 9999 FORMAT( 'TS: M=${.i5}, N=${.i5}, MB=${.i5}, NB=', I5,' test(${.i2})=${.g12_5};
 9998 FORMAT( 'SW: M=${.i5}, N=${.i5}, MB=${.i5}, NB=', I5,' test(${.i2})=${.g12_5};
      }
