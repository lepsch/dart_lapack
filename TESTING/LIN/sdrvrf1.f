      SUBROUTINE SDRVRF1( NOUT, NN, NVAL, THRESH, A, LDA, ARF, WORK )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, NN, NOUT;
      REAL               THRESH
      // ..
      // .. Array Arguments ..
      int                NVAL( NN );
      REAL               A( LDA, * ), ARF( * ), WORK( * )
      // ..

*  =====================================================================
      // ..
      // .. Parameters ..
      REAL               ONE
      const              ONE = 1.0E+0 ;
      int                NTESTS;
      const              NTESTS = 1 ;
      // ..
      // .. Local Scalars ..
      String             UPLO, CFORM, NORM;
      int                I, IFORM, IIN, IIT, INFO, INORM, IUPLO, J, N, NERRS, NFAIL, NRUN;
      REAL               EPS, LARGE, NORMA, NORMARF, SMALL
      // ..
      // .. Local Arrays ..
      String             UPLOS( 2 ), FORMS( 2 ), NORMS( 4 );
      int                ISEED( 4 ), ISEEDY( 4 );
      REAL               RESULT( NTESTS )
      // ..
      // .. External Functions ..
      REAL               SLAMCH, SLANSY, SLANSF, SLARND
      // EXTERNAL SLAMCH, SLANSY, SLANSF, SLARND
      // ..
      // .. External Subroutines ..
      // EXTERNAL STRTTF
      // ..
      // .. Scalars in Common ..
      String             SRNAMT;
      // ..
      // .. Common blocks ..
      COMMON             / SRNAMC / SRNAMT
      // ..
      // .. Data statements ..
      DATA               ISEEDY / 1988, 1989, 1990, 1991 /
      DATA               UPLOS / 'U', 'L' /
      DATA               FORMS / 'N', 'T' /
      DATA               NORMS / 'M', '1', 'I', 'F' /
      // ..
      // .. Executable Statements ..

      // Initialize constants and the random number seed.

      NRUN = 0
      NFAIL = 0
      NERRS = 0
      INFO = 0
      DO 10 I = 1, 4
         ISEED( I ) = ISEEDY( I )
   10 CONTINUE

      EPS = SLAMCH( 'Precision' )
      SMALL = SLAMCH( 'Safe minimum' )
      LARGE = ONE / SMALL
      SMALL = SMALL * LDA * LDA
      LARGE = LARGE / LDA / LDA

      DO 130 IIN = 1, NN

         N = NVAL( IIN )

         DO 120 IIT = 1, 3
            // Nothing to do for N=0
            IF ( N .EQ. 0 ) EXIT

            // Quick Return if possible
            IF ( N .EQ. 0 ) EXIT

            // IIT = 1 : random matrix
            // IIT = 2 : random matrix scaled near underflow
            // IIT = 3 : random matrix scaled near overflow

            DO J = 1, N
               DO I = 1, N
                  A( I, J) = SLARND( 2, ISEED )
               END DO
            END DO

            if ( IIT.EQ.2 ) {
               DO J = 1, N
                  DO I = 1, N
                     A( I, J) = A( I, J ) * LARGE
                  END DO
               END DO
            }

            if ( IIT.EQ.3 ) {
               DO J = 1, N
                  DO I = 1, N
                     A( I, J) = A( I, J) * SMALL
                  END DO
               END DO
            }

            // Do first for UPLO = 'U', then for UPLO = 'L'

            DO 110 IUPLO = 1, 2

               UPLO = UPLOS( IUPLO )

               // Do first for CFORM = 'N', then for CFORM = 'C'

               DO 100 IFORM = 1, 2

                  CFORM = FORMS( IFORM )

                  SRNAMT = 'STRTTF'
                  strttf(CFORM, UPLO, N, A, LDA, ARF, INFO );

                  // Check error code from STRTTF

                  if ( INFO.NE.0 ) {
                     if ( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) {
                        WRITE( NOUT, * )
                        WRITE( NOUT, FMT = 9999 )
                     }
                     WRITE( NOUT, FMT = 9998 ) SRNAMT, UPLO, CFORM, N
                     NERRS = NERRS + 1
                     GO TO 100
                  }

                  DO 90 INORM = 1, 4

                     // Check all four norms: 'M', '1', 'I', 'F'

                     NORM = NORMS( INORM )
                     NORMARF = SLANSF( NORM, CFORM, UPLO, N, ARF, WORK )
                     NORMA = SLANSY( NORM, UPLO, N, A, LDA, WORK )

                     RESULT(1) = ( NORMA - NORMARF ) / NORMA / EPS
                     NRUN = NRUN + 1

                     if ( RESULT(1).GE.THRESH ) {
                        if ( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) {
                           WRITE( NOUT, * )
                           WRITE( NOUT, FMT = 9999 )
                        }
                        WRITE( NOUT, FMT = 9997 ) 'SLANSF', N, IIT, UPLO, CFORM, NORM, RESULT(1)
                        NFAIL = NFAIL + 1
                     }
   90             CONTINUE
  100          CONTINUE
  110       CONTINUE
  120    CONTINUE
  130 CONTINUE

      // Print a summary of the results.

      if ( NFAIL.EQ.0 ) {
         WRITE( NOUT, FMT = 9996 ) 'SLANSF', NRUN
      } else {
         WRITE( NOUT, FMT = 9995 ) 'SLANSF', NFAIL, NRUN
      }
      if ( NERRS.NE.0 ) {
         WRITE( NOUT, FMT = 9994 ) NERRS, 'SLANSF'
      }

 9999 FORMAT( 1X, ' *** Error(s) or Failure(s) while testing SLANSF ***')
 9998 FORMAT( 1X, '     Error in ',A6,' with UPLO=''',A1,''', FORM=''', A1,''', N=',I5)
 9997 FORMAT( 1X, '     Failure in ',A6,' N=',I5,' TYPE=',I5,' UPLO=''', A1, ''', FORM =''',A1,''', NORM=''',A1,''', test=',G12.5)
 9996 FORMAT( 1X, 'All tests for ',A6,' auxiliary routine passed the ', 'threshold ( ',I5,' tests run)')
 9995 FORMAT( 1X, A6, ' auxiliary routine: ',I5,' out of ',I5, ' tests failed to pass the threshold')
 9994 FORMAT( 26X, I5,' error message recorded (',A6,')')

      RETURN

      // End of SDRVRF1

      }
