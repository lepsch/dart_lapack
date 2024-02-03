      SUBROUTINE ZDRVRF1( NOUT, NN, NVAL, THRESH, A, LDA, ARF, WORK )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, NN, NOUT;
      double             THRESH;
      // ..
      // .. Array Arguments ..
      int                NVAL( NN );
      double             WORK( * );
      COMPLEX*16         A( LDA, * ), ARF( * )
      // ..

*  =====================================================================
      // ..
      // .. Parameters ..
      double             ONE;
      const              ONE = 1.0D+0 ;
      int                NTESTS;
      const              NTESTS = 1 ;
      // ..
      // .. Local Scalars ..
      String             UPLO, CFORM, NORM;
      int                I, IFORM, IIN, IIT, INFO, INORM, IUPLO, J, N, NERRS, NFAIL, NRUN;
      double             EPS, LARGE, NORMA, NORMARF, SMALL;
      // ..
      // .. Local Arrays ..
      String             UPLOS( 2 ), FORMS( 2 ), NORMS( 4 );
      int                ISEED( 4 ), ISEEDY( 4 );
      double             RESULT( NTESTS );
      // ..
      // .. External Functions ..
      COMPLEX*16         ZLARND
      double             DLAMCH, ZLANHE, ZLANHF;
      // EXTERNAL DLAMCH, ZLARND, ZLANHE, ZLANHF
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZTRTTF
      // ..
      // .. Scalars in Common ..
      String             SRNAMT;
      // ..
      // .. Common blocks ..
      // COMMON / SRNAMC / SRNAMT
      // ..
      // .. Data statements ..
      DATA               ISEEDY / 1988, 1989, 1990, 1991 /
      DATA               UPLOS / 'U', 'L' /
      DATA               FORMS / 'N', 'C' /
      DATA               NORMS / 'M', '1', 'I', 'F' /
      // ..
      // .. Executable Statements ..

      // Initialize constants and the random number seed.

      NRUN = 0
      NFAIL = 0
      NERRS = 0
      INFO = 0
      for (I = 1; I <= 4; I++) { // 10
         ISEED( I ) = ISEEDY( I )
      } // 10

      EPS = DLAMCH( 'Precision' )
      SMALL = DLAMCH( 'Safe minimum' )
      LARGE = ONE / SMALL
      SMALL = SMALL * LDA * LDA
      LARGE = LARGE / LDA / LDA

      for (IIN = 1; IIN <= NN; IIN++) { // 130

         N = NVAL( IIN )

         for (IIT = 1; IIT <= 3; IIT++) { // 120
            // Nothing to do for N=0
            if (N .EQ. 0) EXIT;

            // IIT = 1 : random matrix
            // IIT = 2 : random matrix scaled near underflow
            // IIT = 3 : random matrix scaled near overflow

            for (J = 1; J <= N; J++) {
               for (I = 1; I <= N; I++) {
                  A( I, J) = ZLARND( 4, ISEED )
               }
            }

            if ( IIT.EQ.2 ) {
               for (J = 1; J <= N; J++) {
                  for (I = 1; I <= N; I++) {
                     A( I, J) = A( I, J ) * LARGE
                  }
               }
            }

            if ( IIT.EQ.3 ) {
               for (J = 1; J <= N; J++) {
                  for (I = 1; I <= N; I++) {
                     A( I, J) = A( I, J) * SMALL
                  }
               }
            }

            // Do first for UPLO = 'U', then for UPLO = 'L'

            for (IUPLO = 1; IUPLO <= 2; IUPLO++) { // 110

               UPLO = UPLOS( IUPLO )

               // Do first for CFORM = 'N', then for CFORM = 'C'

               for (IFORM = 1; IFORM <= 2; IFORM++) { // 100

                  CFORM = FORMS( IFORM )

                  SRNAMT = 'ZTRTTF'
                  ztrttf(CFORM, UPLO, N, A, LDA, ARF, INFO );

                  // Check error code from ZTRTTF

                  if ( INFO.NE.0 ) {
                     if ( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) {
                        WRITE( NOUT, * )
                        WRITE( NOUT, FMT = 9999 )
                     }
                     WRITE( NOUT, FMT = 9998 ) SRNAMT, UPLO, CFORM, N
                     NERRS = NERRS + 1
                     GO TO 100
                  }

                  for (INORM = 1; INORM <= 4; INORM++) { // 90

                     // Check all four norms: 'M', '1', 'I', 'F'

                     NORM = NORMS( INORM )
                     NORMARF = ZLANHF( NORM, CFORM, UPLO, N, ARF, WORK )
                     NORMA = ZLANHE( NORM, UPLO, N, A, LDA, WORK )

                     RESULT(1) = ( NORMA - NORMARF ) / NORMA / EPS
                     NRUN = NRUN + 1

                     if ( RESULT(1).GE.THRESH ) {
                        if ( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) {
                           WRITE( NOUT, * )
                           WRITE( NOUT, FMT = 9999 )
                        }
                        WRITE( NOUT, FMT = 9997 ) 'ZLANHF', N, IIT, UPLO, CFORM, NORM, RESULT(1)
                        NFAIL = NFAIL + 1
                     }
                  } // 90
               } // 100
            } // 110
         } // 120
      } // 130

      // Print a summary of the results.

      if ( NFAIL.EQ.0 ) {
         WRITE( NOUT, FMT = 9996 ) 'ZLANHF', NRUN
      } else {
         WRITE( NOUT, FMT = 9995 ) 'ZLANHF', NFAIL, NRUN
      }
      if ( NERRS.NE.0 ) {
         WRITE( NOUT, FMT = 9994 ) NERRS, 'ZLANHF'
      }

 9999 FORMAT( 1X, ' *** Error(s) or Failure(s) while testing ZLANHF ***')
 9998 FORMAT( 1X, '     Error in ',A6,' with UPLO=''',A1,''', FORM=''', A1,''', N=',I5)
 9997 FORMAT( 1X, '     Failure in ',A6,' N=',I5,' TYPE=',I5,' UPLO=''', A1, ''', FORM =''',A1,''', NORM=''',A1,''', test=',G12.5)
 9996 FORMAT( 1X, 'All tests for ',A6,' auxiliary routine passed the ', 'threshold ( ',I5,' tests run)')
 9995 FORMAT( 1X, A6, ' auxiliary routine:',I5,' out of ',I5, ' tests failed to pass the threshold')
 9994 FORMAT( 26X, I5,' error message recorded (',A6,')')

      RETURN

      // End of ZDRVRF1

      }
