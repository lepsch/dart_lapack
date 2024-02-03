      void dchklqtp(THRESH, TSTERR, NM, MVAL, NN, NVAL, NNB, NBVAL, NOUT ) {
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
      int                I, J, K, L, T, M, N, NB, NFAIL, NERRS, NRUN, MINMN;
      // ..
      // .. Local Arrays ..
      double             RESULT( NTESTS );
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAERH, ALAHD, ALASUM, DERRLQTP, DLQT05
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

      PATH( 1: 1 ) = 'D';
      PATH( 2: 3 ) = 'XQ';
      NRUN = 0;
      NFAIL = 0;
      NERRS = 0;

      // Test the error exits

      if (TSTERR) CALL DERRLQTP( PATH, NOUT );
      INFOT = 0;

      // Do for each value of M

      for (I = 1; I <= NM; I++) {
         M = MVAL( I );

         // Do for each value of N

         for (J = 1; J <= NN; J++) {
            N = NVAL( J );

            // Do for each value of L

            MINMN = min( M, N );
            DO L = 0, MINMN, max( MINMN, 1 );

               // Do for each possible value of NB

               for (K = 1; K <= NNB; K++) {
                  NB = NBVAL( K );

                  // Test DTPLQT and DTPMLQT

                  if ( (NB <= M) && (NB > 0) ) {
                     dlqt05(M, N, L, NB, RESULT );

                     // Print information about the tests that did not
                     // pass the threshold.

                     for (T = 1; T <= NTESTS; T++) {
                     if ( RESULT( T ) >= THRESH ) {
                        if (NFAIL == 0 && NERRS == 0) CALL ALAHD( NOUT, PATH );
                           WRITE( NOUT, FMT = 9999 )M, N, NB, L, T, RESULT( T );
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

 9999 FORMAT( ' M=', I5, ', N=', I5, ', NB=', I4,' L=', I4, ' test(', I2, ')=', G12.5 );
      return;

      // End of DCHKQRTP

      }
