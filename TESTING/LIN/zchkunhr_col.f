      SUBROUTINE ZCHKUNHR_COL( THRESH, TSTERR, NM, MVAL, NN, NVAL, NNB, NBVAL, NOUT )
      IMPLICIT NONE

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               TSTERR;
      int                NM, NN, NNB, NOUT;
      double             THRESH;
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
      double             RESULT( NTESTS );
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAHD, ALASUM, ZERRUNHR_COL, ZUNHR_COL01, ZUNHR_COL02
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
      COMMON             / INFOC / INFOT, NUNIT, OK, LERR
      COMMON             / SRNAMC / SRNAMT
      // ..
      // .. Executable Statements ..

      // Initialize constants

      PATH( 1: 1 ) = 'Z'
      PATH( 2: 3 ) = 'HH'
      NRUN = 0
      NFAIL = 0
      NERRS = 0

      // Test the error exits

      IF( TSTERR ) CALL ZERRUNHR_COL( PATH, NOUT )
      INFOT = 0

      // Do for each value of M in MVAL.

      DO I = 1, NM
         M = MVAL( I )

         // Do for each value of N in NVAL.

         DO J = 1, NN
            N = NVAL( J )

            // Only for M >= N

            if ( MIN( M, N ).GT.0 .AND. M.GE.N ) {

               // Do for each possible value of MB1

               DO IMB1 = 1, NNB
                  MB1 = NBVAL( IMB1 )

                  // Only for MB1 > N

                  if ( MB1.GT.N ) {

                     // Do for each possible value of NB1

                     DO INB1 = 1, NNB
                        NB1 = NBVAL( INB1 )

                        // Do for each possible value of NB2

                        DO INB2 = 1, NNB
                           NB2 = NBVAL( INB2 )

                           if ( NB1.GT.0 .AND. NB2.GT.0 ) {

                              // Test ZUNHR_COL

                              zunhr_col01(M, N, MB1, NB1, NB2, RESULT );

                              // Print information about the tests that did
                              // not pass the threshold.

                              DO T = 1, NTESTS
                                 if ( RESULT( T ).GE.THRESH ) {
                                    IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALAHD( NOUT, PATH )                                     WRITE( NOUT, FMT = 9999 ) M, N, MB1, NB1, NB2, T, RESULT( T )
                                    NFAIL = NFAIL + 1
                                 }
                              END DO
                              NRUN = NRUN + NTESTS
                           }
                        END DO
                     END DO
                  }
                END DO
            }
         END DO
      END DO

      // Do for each value of M in MVAL.

      DO I = 1, NM
         M = MVAL( I )

         // Do for each value of N in NVAL.

         DO J = 1, NN
            N = NVAL( J )

            // Only for M >= N

            if ( MIN( M, N ).GT.0 .AND. M.GE.N ) {

               // Do for each possible value of MB1

               DO IMB1 = 1, NNB
                  MB1 = NBVAL( IMB1 )

                  // Only for MB1 > N

                  if ( MB1.GT.N ) {

                     // Do for each possible value of NB1

                     DO INB1 = 1, NNB
                        NB1 = NBVAL( INB1 )

                        // Do for each possible value of NB2

                        DO INB2 = 1, NNB
                           NB2 = NBVAL( INB2 )

                           if ( NB1.GT.0 .AND. NB2.GT.0 ) {

                              // Test ZUNHR_COL

                              zunhr_col02(M, N, MB1, NB1, NB2, RESULT );

                              // Print information about the tests that did
                              // not pass the threshold.

                              DO T = 1, NTESTS
                                 if ( RESULT( T ).GE.THRESH ) {
                                    IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALAHD( NOUT, PATH )                                     WRITE( NOUT, FMT = 9998 ) M, N, MB1, NB1, NB2, T, RESULT( T )
                                    NFAIL = NFAIL + 1
                                 }
                              END DO
                              NRUN = NRUN + NTESTS
                           }
                        END DO
                     END DO
                  }
                END DO
            }
         END DO
      END DO

      // Print a summary of the results.

      alasum(PATH, NOUT, NFAIL, NRUN, NERRS );

 9999 FORMAT( 'ZUNGTSQR and ZUNHR_COL: M=', I5, ', N=', I5, ', MB1=', I5, ', NB1=', I5, ', NB2=', I5, ' test(', I2, ')=', G12.5 )
 9998 FORMAT( 'ZUNGTSQR_ROW and ZUNHR_COL: M=', I5, ', N=', I5, ', MB1=', I5, ', NB1=', I5, ', NB2=', I5, ' test(', I2, ')=', G12.5 )
      RETURN

      // End of ZCHKUNHR_COL

      }
