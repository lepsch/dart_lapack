      void zdrvrf2(final Nout NOUT, final int NN, final Array<int> NVAL_, final Matrix<double> A_, final int LDA, final int ARF, final int AP, final int ASAV,) {
  final A = A_.having();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                LDA, NN, NOUT;
      int                NVAL( NN );
      Complex         A( LDA, * ), ARF( * ), AP(*), ASAV( LDA, * );
      // ..

// =====================================================================
      bool               LOWER, OK1, OK2;
      String             UPLO, CFORM;
      int                I, IFORM, IIN, INFO, IUPLO, J, N, NRUN;
      String             UPLOS( 2 ), FORMS( 2 );
      final                ISEED=Array<int>( 4 );
      // ..
      // .. External Functions ..
      //- Complex         ZLARND;
      // EXTERNAL ZLARND
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZTFTTR, ZTFTTP, ZTRTTF, ZTRTTP, ZTPTTR, ZTPTTF
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

      // Initialize constants and the random number seed.

      var NRUN = 0;
      var NERRS = Box(0);
      INFO = 0;
      for (I = 1; I <= 4; I++) { // 10
         ISEED[I] = ISEEDY[I - 1];
      } // 10

      for (IIN = 1; IIN <= NN; IIN++) { // 120

         N = NVAL( IIN );

         // Do first for UPLO = 'U', then for UPLO = 'L'

         for (IUPLO = 1; IUPLO <= 2; IUPLO++) { // 110

            final UPLO = UPLOS[IUPLO - 1];
            LOWER = true;
            if (IUPLO == 1) LOWER = false ;

            // Do first for CFORM = 'N', then for CFORM = 'C'

            for (IFORM = 1; IFORM <= 2; IFORM++) { // 100

               CFORM = FORMS( IFORM );

               NRUN++;

               for (J = 1; J <= N; J++) {
                  for (I = 1; I <= N; I++) {
                     A[I][J] = ZLARND( 4, ISEED );
                  }
               }

              srnamc.SRNAMT = 'ZTRTTF';
               ztrttf(CFORM, UPLO, N, A, LDA, ARF, INFO );

              srnamc.SRNAMT = 'ZTFTTP';
               ztfttp(CFORM, UPLO, N, ARF, AP, INFO );

              srnamc.SRNAMT = 'ZTPTTR';
               ztpttr(UPLO, N, AP, ASAV, LDA, INFO );

               OK1 = true;
               if ( LOWER ) {
                  for (J = 1; J <= N; J++) {
                     for (I = J; I <= N; I++) {
                        if ( A(I,J) != ASAV(I,J) ) {
                           OK1 = false;
                        }
                     }
                  }
               } else {
                  for (J = 1; J <= N; J++) {
                     for (I = 1; I <= J; I++) {
                        if ( A(I,J) != ASAV(I,J) ) {
                           OK1 = false;
                        }
                     }
                  }
               }

               NRUN++;

              srnamc.SRNAMT = 'ZTRTTP';
               ztrttp(UPLO, N, A, LDA, AP, INFO );

              srnamc.SRNAMT = 'ZTPTTF';
               ztpttf(CFORM, UPLO, N, AP, ARF, INFO );

              srnamc.SRNAMT = 'ZTFTTR';
               ztfttr(CFORM, UPLO, N, ARF, ASAV, LDA, INFO );

               OK2 = true;
               if ( LOWER ) {
                  for (J = 1; J <= N; J++) {
                     for (I = J; I <= N; I++) {
                        if ( A(I,J) != ASAV(I,J) ) {
                           OK2 = false;
                        }
                     }
                  }
               } else {
                  for (J = 1; J <= N; J++) {
                     for (I = 1; I <= J; I++) {
                        if ( A(I,J) != ASAV(I,J) ) {
                           OK2 = false;
                        }
                     }
                  }
               }

               if (( !OK1 ) || ( !OK2 )) {
                  if ( NERRS == 0 ) {
                     WRITE( NOUT, * );
                     NOUT.println( 9999 );
                  }
                  NOUT.println( 9998 ) N, UPLO, CFORM;
                  NERRS++;
               }

            } // 100
         } // 110
      } // 120

      // Print a summary of the results.

      if ( NERRS == 0 ) {
         NOUT.println( 9997 ) NRUN;
      } else {
         NOUT.println( 9996 ) NERRS, NRUN;
      }

 9999 FORMAT('  *** Error(s) while testing the RFP conversion routines ***');
 9998 FORMAT('      Error in RFP,conversion routines N=',I5, ' UPLO=\'${.a1}\', FORM =\'${.a1}\'');
 9997 FORMAT(' All tests for the RFP conversion routines passed (', I5,' tests run)');
 9996 FORMAT(' RFP conversion routines:',I5,' out of ',I5, ' error message recorded');

      }
