      void schkorhr_col(final int THRESH, final int TSTERR, final int NM, final int MVAL, final int NN, final int NVAL, final int NNB, final int NBVAL, final int NOUT,) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      bool               TSTERR;
      int                NM, NN, NNB, NOUT;
      double               THRESH;
      int                MVAL( * ), NBVAL( * ), NVAL( * );
      // ..

      int                NTESTS;
      const              NTESTS = 6 ;
      String   (LEN=3)   PATH;
      int                I, IMB1, INB1, INB2, J, T, M, N, MB1, NB1, NB2, NFAIL, NERRS, NRUN;

      // .. Local Arrays ..
      double               RESULT( NTESTS );
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAHD, ALASUM, SERRORHR_COL, SORHR_COL01, SORHR_COL02
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Scalars in Common ..
      bool               LERR, OK;
      String   (LEN=32) srnamc.SRNAMT;
      int                INFOT, NUNIT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / INFOT, NUNIT, OK, LERR
      // COMMON / SRNAMC /srnamc.SRNAMT

      // Initialize constants

      PATH[1: 1] = 'S';
      PATH[2: 3] = 'HH';
      NRUN = 0;
      NFAIL = 0;
      NERRS = 0;

      // Test the error exits

      if (TSTERR) serrorhr_col( PATH, NOUT );
      INFOT = 0;

      // Do for each value of M in MVAL.

      for (I = 1; I <= NM; I++) {
         M = MVAL( I );

         // Do for each value of N in NVAL.

         for (J = 1; J <= NN; J++) {
            N = NVAL( J );

            // Only for M >= N

            if ( min( M, N ) > 0 && M >= N ) {

               // Do for each possible value of MB1

               for (IMB1 = 1; IMB1 <= NNB; IMB1++) {
                  MB1 = NBVAL( IMB1 );

                  // Only for MB1 > N

                  if ( MB1 > N ) {

                     // Do for each possible value of NB1

                     for (INB1 = 1; INB1 <= NNB; INB1++) {
                        NB1 = NBVAL( INB1 );

                        // Do for each possible value of NB2

                        for (INB2 = 1; INB2 <= NNB; INB2++) {
                           NB2 = NBVAL( INB2 );

                           if ( NB1 > 0 && NB2 > 0 ) {

                              // Test SORHR_COL

                              sorhr_col01(M, N, MB1, NB1, NB2, RESULT );

                              // Print information about the tests that did
                              // not pass the threshold.

                              for (T = 1; T <= NTESTS; T++) {
                                 if ( RESULT( T ) >= THRESH ) {
                                    if (NFAIL == 0 && NERRS == 0) alahd( NOUT, PATH );
                                    WRITE( NOUT, FMT = 9999 ) M, N, MB1, NB1, NB2, T, RESULT( T );
                                    NFAIL = NFAIL + 1;
                                 }
                              }
                              NRUN = NRUN + NTESTS;
                           }
                        }
                     }
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

            // Only for M >= N

            if ( min( M, N ) > 0 && M >= N ) {

               // Do for each possible value of MB1

               for (IMB1 = 1; IMB1 <= NNB; IMB1++) {
                  MB1 = NBVAL( IMB1 );

                  // Only for MB1 > N

                  if ( MB1 > N ) {

                     // Do for each possible value of NB1

                     for (INB1 = 1; INB1 <= NNB; INB1++) {
                        NB1 = NBVAL( INB1 );

                        // Do for each possible value of NB2

                        for (INB2 = 1; INB2 <= NNB; INB2++) {
                           NB2 = NBVAL( INB2 );

                           if ( NB1 > 0 && NB2 > 0 ) {

                              // Test SORHR_COL

                              sorhr_col02(M, N, MB1, NB1, NB2, RESULT );

                              // Print information about the tests that did
                              // not pass the threshold.

                              for (T = 1; T <= NTESTS; T++) {
                                 if ( RESULT( T ) >= THRESH ) {
                                    if (NFAIL == 0 && NERRS == 0) alahd( NOUT, PATH );
                                    WRITE( NOUT, FMT = 9998 ) M, N, MB1, NB1, NB2, T, RESULT( T );
                                    NFAIL = NFAIL + 1;
                                 }
                              }
                              NRUN = NRUN + NTESTS;
                           }
                        }
                     }
                  }
                }
            }
         }
      }

      // Print a summary of the results.

      alasum(PATH, NOUT, NFAIL, NRUN, NERRS );

 9999 FORMAT( 'SORGTSQR and SORHR_COL: M=${.i5}, N=${.i5}, MB1=${.i5}, NB1=${.i5}, NB2=${.i5} test(${.i2})=${.g12_5};
 9998 FORMAT( 'SORGTSQR_ROW and SORHR_COL: M=${.i5}, N=${.i5}, MB1=${.i5}, NB1=${.i5}, NB2=${.i5} test(${.i2})=${.g12_5};
      }
