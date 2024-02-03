      SUBROUTINE ZCHKQRTP( THRESH, TSTERR, NM, MVAL, NN, NVAL, NNB, NBVAL, NOUT )
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
      PARAMETER          ( NTESTS = 6 )
      // ..
      // .. Local Scalars ..
      String             PATH;
      int                I, J, K, T, L, M, N, NB, NFAIL, NERRS, NRUN, MINMN;
      // ..
      // .. Local Arrays ..
      double             RESULT( NTESTS );
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAERH, ALAHD, ALASUM, ZERRQRTP
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

      PATH( 1: 1 ) = 'Z'
      PATH( 2: 3 ) = 'QX'
      NRUN = 0
      NFAIL = 0
      NERRS = 0

      // Test the error exits

      IF( TSTERR ) CALL ZERRQRTP( PATH, NOUT )
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

                  // Test ZTPQRT and ZTPMQRT

                  IF( (NB.LE.N).AND.(NB.GT.0) ) THEN
                     CALL ZQRT05( M, N, L, NB, RESULT )

                     // Print information about the tests that did not
                     // pass the threshold.

                     DO T = 1, NTESTS
                        IF( RESULT( T ).GE.THRESH ) THEN
                           IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALAHD( NOUT, PATH )                            WRITE( NOUT, FMT = 9999 )M, N, NB, T, RESULT( T )
                           NFAIL = NFAIL + 1
                        END IF
                     END DO
                     NRUN = NRUN + NTESTS
                  END IF
               END DO
            END DO
         END DO
      END DO

      // Print a summary of the results.

      CALL ALASUM( PATH, NOUT, NFAIL, NRUN, NERRS )

 9999 FORMAT( ' M=', I5, ', N=', I5, ', NB=', I4,
     $      ' test(', I2, ')=', G12.5 )
      RETURN

      // End of ZCHKQRTP

      END
