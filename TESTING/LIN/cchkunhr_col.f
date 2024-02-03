      SUBROUTINE CCHKUNHR_COL( THRESH, TSTERR, NM, MVAL, NN, NVAL, NNB, NBVAL, NOUT )
      IMPLICIT NONE

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               TSTERR;
      int                NM, NN, NNB, NOUT;
      REAL               THRESH
      // ..
      // .. Array Arguments ..
      int                MVAL( * ), NBVAL( * ), NVAL( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      int                NTESTS;
      const              NTESTS = 6 ;
      // ..
      // .. Local Scalars ..
      String   (LEN=3)   PATH;
      int                I, IMB1, INB1, INB2, J, T, M, N, MB1, NB1, NB2, NFAIL, NERRS, NRUN;

      // .. Local Arrays ..
      REAL               RESULT( NTESTS )
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAHD, ALASUM, CERRUNHR_COL, CUNHR_COL01, CUNHR_COL02
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Scalars in Common ..
      bool               LERR, OK;
      String   (LEN=32)  SRNAMT;
      int                INFOT, NUNIT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / INFOT, NUNIT, OK, LERR
      // COMMON / SRNAMC / SRNAMT
      // ..
      // .. Executable Statements ..

      // Initialize constants

      PATH( 1: 1 ) = 'C'
      PATH( 2: 3 ) = 'HH'
      NRUN = 0
      NFAIL = 0
      NERRS = 0

      // Test the error exits

      if (TSTERR) CALL CERRUNHR_COL( PATH, NOUT );
      INFOT = 0

      // Do for each value of M in MVAL.

      for (I = 1; I <= NM; I++) {
         M = MVAL( I )

         // Do for each value of N in NVAL.

         for (J = 1; J <= NN; J++) {
            N = NVAL( J )

            // Only for M >= N

            if ( MIN( M, N ).GT.0 && M.GE.N ) {

               // Do for each possible value of MB1

               for (IMB1 = 1; IMB1 <= NNB; IMB1++) {
                  MB1 = NBVAL( IMB1 )

                  // Only for MB1 > N

                  if ( MB1.GT.N ) {

                     // Do for each possible value of NB1

                     for (INB1 = 1; INB1 <= NNB; INB1++) {
                        NB1 = NBVAL( INB1 )

                        // Do for each possible value of NB2

                        for (INB2 = 1; INB2 <= NNB; INB2++) {
                           NB2 = NBVAL( INB2 )

                           if ( NB1.GT.0 && NB2.GT.0 ) {

                              // Test CUNHR_COL

                              cunhr_col01(M, N, MB1, NB1, NB2, RESULT );

                              // Print information about the tests that did
                              // not pass the threshold.

                              for (T = 1; T <= NTESTS; T++) {
                                 if ( RESULT( T ).GE.THRESH ) {
                                    if (NFAIL == 0 && NERRS == 0) CALL ALAHD( NOUT, PATH )                                     WRITE( NOUT, FMT = 9999 ) M, N, MB1, NB1, NB2, T, RESULT( T );
                                    NFAIL = NFAIL + 1
                                 }
                              }
                              NRUN = NRUN + NTESTS
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
         M = MVAL( I )

         // Do for each value of N in NVAL.

         for (J = 1; J <= NN; J++) {
            N = NVAL( J )

            // Only for M >= N

            if ( MIN( M, N ).GT.0 && M.GE.N ) {

               // Do for each possible value of MB1

               for (IMB1 = 1; IMB1 <= NNB; IMB1++) {
                  MB1 = NBVAL( IMB1 )

                  // Only for MB1 > N

                  if ( MB1.GT.N ) {

                     // Do for each possible value of NB1

                     for (INB1 = 1; INB1 <= NNB; INB1++) {
                        NB1 = NBVAL( INB1 )

                        // Do for each possible value of NB2

                        for (INB2 = 1; INB2 <= NNB; INB2++) {
                           NB2 = NBVAL( INB2 )

                           if ( NB1.GT.0 && NB2.GT.0 ) {

                              // Test CUNHR_COL

                              cunhr_col02(M, N, MB1, NB1, NB2, RESULT );

                              // Print information about the tests that did
                              // not pass the threshold.

                              for (T = 1; T <= NTESTS; T++) {
                                 if ( RESULT( T ).GE.THRESH ) {
                                    if (NFAIL == 0 && NERRS == 0) CALL ALAHD( NOUT, PATH )                                     WRITE( NOUT, FMT = 9998 ) M, N, MB1, NB1, NB2, T, RESULT( T );
                                    NFAIL = NFAIL + 1
                                 }
                              }
                              NRUN = NRUN + NTESTS
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

 9999 FORMAT( 'CUNGTSQR and CUNHR_COL: M=', I5, ', N=', I5, ', MB1=', I5, ', NB1=', I5, ', NB2=', I5, ' test(', I2, ')=', G12.5 )
 9998 FORMAT( 'CUNGTSQR_ROW and CUNHR_COL: M=', I5, ', N=', I5, ', MB1=', I5, ', NB1=', I5, ', NB2=', I5, ' test(', I2, ')=', G12.5 )
      RETURN

      // End of CCHKUNHR_COL

      }
