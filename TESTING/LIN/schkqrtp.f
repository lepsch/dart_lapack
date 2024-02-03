      SUBROUTINE SCHKQRTP( THRESH, TSTERR, NM, MVAL, NN, NVAL, NNB, NBVAL, NOUT )
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
      String             PATH;
      int                I, J, K, T, L, M, N, NB, NFAIL, NERRS, NRUN, MINMN;
      // ..
      // .. Local Arrays ..
      REAL               RESULT( NTESTS )
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAERH, ALAHD, ALASUM, SERRQRTP, SQRT05
      // ..
      // .. Scalars in Common ..
      bool               LERR, OK;
      String             SRNAMT;
      int                INFOT, NUNIT;
      // ..
      // .. Common blocks ..
      COMMON             / INFOC / INFOT, NUNIT, OK, LERR
      COMMON             / SRNAMC / SRNAMT
      // ..
      // .. Executable Statements ..

      // Initialize constants

      PATH( 1: 1 ) = 'S'
      PATH( 2: 3 ) = 'QX'
      NRUN = 0
      NFAIL = 0
      NERRS = 0

      // Test the error exits

      IF( TSTERR ) CALL SERRQRTP( PATH, NOUT )
      INFOT = 0

      // Do for each value of M

      DO I = 1, NM
         M = MVAL( I )

         // Do for each value of N

         DO J = 1, NN
            N = NVAL( J )

            // Do for each value of L

            MINMN = MIN( M, N )
            DO L = 0, MINMN, MAX( MINMN, 1 )

               // Do for each possible value of NB

               DO K = 1, NNB
                  NB = NBVAL( K )

                     // Test STPQRT and STPMQRT

                  if ( (NB.LE.N).AND.(NB.GT.0) ) {
                     sqrt05(M, N, L, NB, RESULT );

                     // Print information about the tests that did not
                     // pass the threshold.

                     DO T = 1, NTESTS
                        if ( RESULT( T ).GE.THRESH ) {
                           IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALAHD( NOUT, PATH )                            WRITE( NOUT, FMT = 9999 )M, N, NB, L, T, RESULT( T )
                           NFAIL = NFAIL + 1
                        }
                     END DO
                     NRUN = NRUN + NTESTS
                  }
               END DO
            END DO
         END DO
      END DO

      // Print a summary of the results.

      alasum(PATH, NOUT, NFAIL, NRUN, NERRS );

 9999 FORMAT( ' M=', I5, ', N=', I5, ', NB=', I4, ', L=', I4, ' test(', I2, ')=', G12.5 )
      RETURN

      // End of SCHKQRTP

      }
