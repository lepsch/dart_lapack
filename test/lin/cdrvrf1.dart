      void cdrvrf1(final int NOUT, final int NN, final int NVAL, final int THRESH, final Matrix<double> A, final int LDA, final int ARF, final Array<double> WORK) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                LDA, NN, NOUT;
      double               THRESH;
      int                NVAL( NN );
      double               WORK( * );
      Complex            A( LDA, * ), ARF( * );
      // ..

// =====================================================================
      // ..
      // .. Parameters ..
      double               ONE;
      const              ONE = 1.0 ;
      int                NTESTS;
      const              NTESTS = 1 ;
      String             UPLO, CFORM, NORM;
      int                I, IFORM, IIN, IIT, INFO, INORM, IUPLO, J, N, NERRS, NFAIL, NRUN;
      double               EPS, LARGE, NORMA, NORMARF, SMALL;
      String             UPLOS( 2 ), FORMS( 2 ), NORMS( 4 );
      int                ISEED( 4 ), ISEEDY( 4 );
      double               RESULT( NTESTS );
      // ..
      // .. External Functions ..
      //- COMPLEX            CLARND;
      //- REAL               SLAMCH, CLANHE, CLANHF;
      // EXTERNAL SLAMCH, CLARND, CLANHE, CLANHF
      // ..
      // .. External Subroutines ..
      // EXTERNAL CTRTTF
      // ..
      // .. Scalars in Common ..
      String            srnamc.SRNAMT;
      // ..
      // .. Common blocks ..
      // COMMON / SRNAMC /srnamc.SRNAMT
      // ..
      // .. Data statements ..
      const ISEEDY = [ 1988, 1989, 1990, 1991 ];
      const UPLOS = [ 'U', 'L' ];
      const FORMS = [ 'N', 'C' ];
      const NORMS = [ 'M', '1', 'I', 'F' ];

      // Initialize constants and the random number seed.

      NRUN = 0;
      NFAIL = 0;
      NERRS = 0;
      INFO = 0;
      for (I = 1; I <= 4; I++) { // 10
         ISEED[I] = ISEEDY( I );
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
            if (N == 0) break;

            // IIT = 1 : random matrix
            // IIT = 2 : random matrix scaled near underflow
            // IIT = 3 : random matrix scaled near overflow

            for (J = 1; J <= N; J++) {
               for (I = 1; I <= N; I++) {
                  A[I][J] = CLARND( 4, ISEED );
               }
            }

            if ( IIT == 2 ) {
               for (J = 1; J <= N; J++) {
                  for (I = 1; I <= N; I++) {
                     A[I][J] = A( I, J ) * LARGE;
                  }
               }
            }

            if ( IIT == 3 ) {
               for (J = 1; J <= N; J++) {
                  for (I = 1; I <= N; I++) {
                     A[I][J] = A( I, J) * SMALL;
                  }
               }
            }

            // Do first for UPLO = 'U', then for UPLO = 'L'

            for (IUPLO = 1; IUPLO <= 2; IUPLO++) { // 110

               UPLO = UPLOS( IUPLO );

               // Do first for CFORM = 'N', then for CFORM = 'C'

               for (IFORM = 1; IFORM <= 2; IFORM++) { // 100

                  CFORM = FORMS( IFORM );

                 srnamc.SRNAMT = 'CTRTTF';
                  ctrttf(CFORM, UPLO, N, A, LDA, ARF, INFO );

                  // Check error code from CTRTTF

                  if ( INFO != 0 ) {
                     if ( NFAIL == 0 && NERRS == 0 ) {
                        WRITE( NOUT, * );
                        WRITE( NOUT, FMT = 9999 );
                     }
                     WRITE( NOUT, FMT = 9998 )srnamc.SRNAMT, UPLO, CFORM, N;
                     NERRS = NERRS + 1;
                     GO TO 100;
                  }

                  for (INORM = 1; INORM <= 4; INORM++) { // 90

                     // Check all four norms: 'M', '1', 'I', 'F'

                     NORM = NORMS( INORM );
                     NORMARF = CLANHF( NORM, CFORM, UPLO, N, ARF, WORK );
                     NORMA = CLANHE( NORM, UPLO, N, A, LDA, WORK );

                     RESULT[1] = ( NORMA - NORMARF ) / NORMA / EPS;
                     NRUN = NRUN + 1;

                     if ( RESULT(1) >= THRESH ) {
                        if ( NFAIL == 0 && NERRS == 0 ) {
                           WRITE( NOUT, * );
                           WRITE( NOUT, FMT = 9999 );
                        }
                        WRITE( NOUT, FMT = 9997 ) 'CLANHF', N, IIT, UPLO, CFORM, NORM, RESULT(1);
                        NFAIL = NFAIL + 1;
                     }
                  } // 90
               } // 100
            } // 110
         } // 120
      } // 130

      // Print a summary of the results.

      if ( NFAIL == 0 ) {
         WRITE( NOUT, FMT = 9996 )'CLANHF', NRUN;
      } else {
         WRITE( NOUT, FMT = 9995 ) 'CLANHF', NFAIL, NRUN;
      }
      if ( NERRS != 0 ) {
         WRITE( NOUT, FMT = 9994 ) NERRS, 'CLANHF';
      }

 9999 FORMAT('  *** Error(s) or Failure(s) while testing CLANHF ***');
 9998 FORMAT('      Error in ${.a6} with UPLO=''${.a1}'', FORM=''${.a1}'', N=',I5);
 9997 FORMAT('      Failure in ${.a6} N=',I5,' TYPE=',I5,' UPLO=''${.a1}'', FORM =''${.a1}'', NORM=''${.a1}'', test=',G12.5);
 9996 FORMAT(' All tests for ${.a6} auxiliary routine passed the threshold ( ',I5,' tests run)');
 9995 FORMAT(' ${.a6} auxiliary routine: ',I5,' out of ',I5, ' tests failed to pass the threshold');
 9994 FORMAT('${' ' * 26}${.i5} error message recorded (${.a6})');

      }
