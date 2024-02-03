      void sdrvrf2(NOUT, NN, NVAL, A, LDA, ARF, AP, ASAV  ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, NN, NOUT;
      // ..
      // .. Array Arguments ..
      int                NVAL( NN );
      REAL               A( LDA, * ), ARF( * ), AP(*), ASAV( LDA, * );
      // ..

// =====================================================================
      // ..
      // .. Local Scalars ..
      bool               LOWER, OK1, OK2;
      String             UPLO, CFORM;
      int                I, IFORM, IIN, INFO, IUPLO, J, N, NERRS, NRUN;
      // ..
      // .. Local Arrays ..
      String             UPLOS( 2 ), FORMS( 2 );
      int                ISEED( 4 ), ISEEDY( 4 );
      // ..
      // .. External Functions ..
      REAL               SLARND;
      // EXTERNAL SLARND
      // ..
      // .. External Subroutines ..
      // EXTERNAL STFTTR, STFTTP, STRTTF, STRTTP, STPTTR, STPTTF
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
      // ..
      // .. Executable Statements ..

      // Initialize constants and the random number seed.

      NRUN = 0;
      NERRS = 0;
      INFO = 0;
      for (I = 1; I <= 4; I++) { // 10
         ISEED( I ) = ISEEDY( I );
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
                     A( I, J) = SLARND( 2, ISEED );
                  }
               }

               SRNAMT = 'DTRTTF';
               strttf(CFORM, UPLO, N, A, LDA, ARF, INFO );

               SRNAMT = 'DTFTTP';
               stfttp(CFORM, UPLO, N, ARF, AP, INFO );

               SRNAMT = 'DTPTTR';
               stpttr(UPLO, N, AP, ASAV, LDA, INFO );

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

               SRNAMT = 'DTRTTP';
               strttp(UPLO, N, A, LDA, AP, INFO );

               SRNAMT = 'DTPTTF';
               stpttf(CFORM, UPLO, N, AP, ARF, INFO );

               SRNAMT = 'DTFTTR';
               stfttr(CFORM, UPLO, N, ARF, ASAV, LDA, INFO );

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

 9999 FORMAT( 1X, ' *** Error(s) while testing the RFP conversion', ' routines ***');
 9998 FORMAT( 1X, '     Error in RFP,conversion routines N=',I5, ' UPLO=''', A1, ''', FORM =''',A1,'''');
 9997 FORMAT( 1X, 'All tests for the RFP conversion routines passed ( ', I5,' tests run)');
 9996 FORMAT( 1X, 'RFP conversion routines: ',I5,' out of ',I5, ' error message recorded');

      return;
      }
