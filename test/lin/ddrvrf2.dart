import 'common.dart';

      void ddrvrf2(final int NOUT, final int NN, final int NVAL, final Matrix<double> A_, final int LDA, final int ARF, final int AP, final int ASAV,) {
  final A = A_.dim();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                LDA, NN, NOUT;
      int                NVAL( NN );
      double             A( LDA, * ), ARF( * ), AP(*), ASAV( LDA, * );
      // ..

// =====================================================================
      bool               LOWER, OK1, OK2;
      String             UPLO, CFORM;
      int                I, IFORM, IIN, INFO, IUPLO, J, N, NERRS, NRUN;
      String             UPLOS( 2 ), FORMS( 2 );
      int                ISEED( 4 ), ISEEDY( 4 );
      // ..
      // .. External Functions ..
      //- double             DLARND;
      // EXTERNAL DLARND
      // ..
      // .. External Subroutines ..
      // EXTERNAL DTFTTR, DTFTTP, DTRTTF, DTRTTP, DTPTTR, DTPTTF
      // ..
      // .. Scalars in Common ..
      String             srnamc.SRNAMT;
      // ..
      // .. Common blocks ..
      // COMMON / SRNAMC / srnamc.SRNAMT
      // ..
      // .. Data statements ..
      const ISEEDY = [ 1988, 1989, 1990, 1991 ];
      const UPLOS = [ 'U', 'L' ];
      const FORMS = [ 'N', 'T' ];

      // Initialize constants and the random number seed.

      NRUN = 0;
      NERRS = 0;
      INFO = 0;
      for (I = 1; I <= 4; I++) { // 10
         ISEED[I] = ISEEDY( I );
      } // 10

      for (IIN = 1; IIN <= NN; IIN++) { // 120

         N = NVAL( IIN );

         // Do first for UPLO = 'U', then for UPLO = 'L'

         for (IUPLO = 1; IUPLO <= 2; IUPLO++) { // 110

            UPLO = UPLOS( IUPLO );
            LOWER = true;
            if (IUPLO == 1) LOWER = false ;

            // Do first for CFORM = 'N', then for CFORM = 'T'

            for (IFORM = 1; IFORM <= 2; IFORM++) { // 100

               CFORM = FORMS( IFORM );

               NRUN = NRUN + 1;

               for (J = 1; J <= N; J++) {
                  for (I = 1; I <= N; I++) {
                     A[I][J] = dlarnd( 2, ISEED );
                  }
               }

               srnamc.SRNAMT = 'DTRTTF';
               dtrttf(CFORM, UPLO, N, A, LDA, ARF, INFO );

               srnamc.SRNAMT = 'DTFTTP';
               dtfttp(CFORM, UPLO, N, ARF, AP, INFO );

               srnamc.SRNAMT = 'DTPTTR';
               dtpttr(UPLO, N, AP, ASAV, LDA, INFO );

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

               NRUN = NRUN + 1;

               srnamc.SRNAMT = 'DTRTTP';
               dtrttp(UPLO, N, A, LDA, AP, INFO );

               srnamc.SRNAMT = 'DTPTTF';
               dtpttf(CFORM, UPLO, N, AP, ARF, INFO );

               srnamc.SRNAMT = 'DTFTTR';
               dtfttr(CFORM, UPLO, N, ARF, ASAV, LDA, INFO );

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
                     WRITE( NOUT, FMT = 9999 );
                  }
                  WRITE( NOUT, FMT = 9998 ) N, UPLO, CFORM;
                  NERRS = NERRS + 1;
               }

            } // 100
         } // 110
      } // 120

      // Print a summary of the results.

      if ( NERRS == 0 ) {
         WRITE( NOUT, FMT = 9997 ) NRUN;
      } else {
         WRITE( NOUT, FMT = 9996 ) NERRS, NRUN;
      }

 9999 FORMAT('  *** Error(s) while testing the RFP conversion routines ***');
 9998 FORMAT('      Error in RFP,conversion routines N=',I5, ' UPLO=\'${.a1}\', FORM =\'${.a1}\'');
 9997 FORMAT(' All tests for the RFP conversion routines passed ( ', I5,' tests run)');
 9996 FORMAT(' RFP conversion routines: ',I5,' out of ',I5, ' error message recorded');

      }
