      void sdrvrf1(NOUT, NN, NVAL, THRESH, A, LDA, ARF, WORK ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, NN, NOUT;
      REAL               THRESH;
      // ..
      // .. Array Arguments ..
      int                NVAL( NN );
      REAL               A( LDA, * ), ARF( * ), WORK( * );
      // ..

// =====================================================================
      // ..
      // .. Parameters ..
      REAL               ONE;
      const              ONE = 1.0 ;
      int                NTESTS;
      const              NTESTS = 1 ;
      // ..
      // .. Local Scalars ..
      String             UPLO, CFORM, NORM;
      int                I, IFORM, IIN, IIT, INFO, INORM, IUPLO, J, N, NERRS, NFAIL, NRUN;
      REAL               EPS, LARGE, NORMA, NORMARF, SMALL;
      // ..
      // .. Local Arrays ..
      String             UPLOS( 2 ), FORMS( 2 ), NORMS( 4 );
      int                ISEED( 4 ), ISEEDY( 4 );
      REAL               RESULT( NTESTS );
      // ..
      // .. External Functions ..
      REAL               SLAMCH, SLANSY, SLANSF, SLARND;
      // EXTERNAL SLAMCH, SLANSY, SLANSF, SLARND
      // ..
      // .. External Subroutines ..
      // EXTERNAL STRTTF
      // ..
      // .. Scalars in Common ..
      String             SRNAMT;
      // ..
      // .. Common blocks ..
      // COMMON / SRNAMC / SRNAMT
      // ..
      // .. Data statements ..
      const ISEEDY = [ 1988, 1989, 1990, 1991 ];
      const UPLOS = [ 'U', 'L' ];
      const FORMS = [ 'N', 'T' ];
      const NORMS = [ 'M', '1', 'I', 'F' ];
      // ..
      // .. Executable Statements ..

      // Initialize constants and the random number seed.

      NRUN = 0;
      NFAIL = 0;
      NERRS = 0;
      INFO = 0;
      for (I = 1; I <= 4; I++) { // 10
         ISEED( I ) = ISEEDY( I );
      } // 10

      EPS = SLAMCH( 'Precision' );
      SMALL = SLAMCH( 'Safe minimum' );
      LARGE = ONE / SMALL;
      SMALL = SMALL * LDA * LDA;
      LARGE = LARGE / LDA / LDA;

      for (IIN = 1; IIN <= NN; IIN++) { // 130

         N = NVAL( IIN );

         for (IIT = 1; IIT <= 3; IIT++) { // 120
            // Nothing to do for N=0
            if (N == 0) EXIT;

            // Quick Return if possible
            if (N == 0) EXIT;

            // IIT = 1 : random matrix
            // IIT = 2 : random matrix scaled near underflow
            // IIT = 3 : random matrix scaled near overflow

            for (J = 1; J <= N; J++) {
               for (I = 1; I <= N; I++) {
                  A( I, J) = SLARND( 2, ISEED );
               }
            }

            if ( IIT == 2 ) {
               for (J = 1; J <= N; J++) {
                  for (I = 1; I <= N; I++) {
                     A( I, J) = A( I, J ) * LARGE;
                  }
               }
            }

            if ( IIT == 3 ) {
               for (J = 1; J <= N; J++) {
                  for (I = 1; I <= N; I++) {
                     A( I, J) = A( I, J) * SMALL;
                  }
               }
            }

            // Do first for UPLO = 'U', then for UPLO = 'L'

            for (IUPLO = 1; IUPLO <= 2; IUPLO++) { // 110

               UPLO = UPLOS( IUPLO );

               // Do first for CFORM = 'N', then for CFORM = 'C'

               for (IFORM = 1; IFORM <= 2; IFORM++) { // 100

                  CFORM = FORMS( IFORM );

                  SRNAMT = 'STRTTF';
                  strttf(CFORM, UPLO, N, A, LDA, ARF, INFO );

                  // Check error code from STRTTF

                  if ( INFO != 0 ) {
                     if ( NFAIL == 0 && NERRS == 0 ) {
                        WRITE( NOUT, * );
                        WRITE( NOUT, FMT = 9999 );
                     }
                     WRITE( NOUT, FMT = 9998 ) SRNAMT, UPLO, CFORM, N;
                     NERRS = NERRS + 1;
                     GO TO 100;
                  }

                  for (INORM = 1; INORM <= 4; INORM++) { // 90

                     // Check all four norms: 'M', '1', 'I', 'F'

                     NORM = NORMS( INORM );
                     NORMARF = SLANSF( NORM, CFORM, UPLO, N, ARF, WORK );
                     NORMA = SLANSY( NORM, UPLO, N, A, LDA, WORK );

                     RESULT(1) = ( NORMA - NORMARF ) / NORMA / EPS;
                     NRUN = NRUN + 1;

                     if ( RESULT(1) >= THRESH ) {
                        if ( NFAIL == 0 && NERRS == 0 ) {
                           WRITE( NOUT, * );
                           WRITE( NOUT, FMT = 9999 );
                        }
                        WRITE( NOUT, FMT = 9997 ) 'SLANSF', N, IIT, UPLO, CFORM, NORM, RESULT(1);
                        NFAIL = NFAIL + 1;
                     }
                  } // 90
               } // 100
            } // 110
         } // 120
      } // 130

      // Print a summary of the results.

      if ( NFAIL == 0 ) {
         WRITE( NOUT, FMT = 9996 ) 'SLANSF', NRUN;
      } else {
         WRITE( NOUT, FMT = 9995 ) 'SLANSF', NFAIL, NRUN;
      }
      if ( NERRS != 0 ) {
         WRITE( NOUT, FMT = 9994 ) NERRS, 'SLANSF';
      }

 9999 FORMAT( 1X, ' *** Error(s) or Failure(s) while testing SLANSF ***');
 9998 FORMAT( 1X, '     Error in ',A6,' with UPLO=''',A1,''', FORM=''', A1,''', N=',I5);
 9997 FORMAT( 1X, '     Failure in ',A6,' N=',I5,' TYPE=',I5,' UPLO=''', A1, ''', FORM =''',A1,''', NORM=''',A1,''', test=',G12.5);
 9996 FORMAT( 1X, 'All tests for ',A6,' auxiliary routine passed the ', 'threshold ( ',I5,' tests run)');
 9995 FORMAT( 1X, A6, ' auxiliary routine: ',I5,' out of ',I5, ' tests failed to pass the threshold');
 9994 FORMAT( 26X, I5,' error message recorded (',A6,')');

      return;
      }
