      void sdrvrfp(NOUT, NN, NVAL, NNS, NSVAL, NNT, NTVAL, THRESH, A, ASAV, AFAC, AINV, B, BSAV, XACT, X, ARF, ARFINV, S_WORK_SLATMS, S_WORK_SPOT01, S_TEMP_SPOT02, S_TEMP_SPOT03, S_WORK_SLANSY, S_WORK_SPOT02, S_WORK_SPOT03 ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                NN, NNS, NNT, NOUT;
      REAL               THRESH;
      // ..
      // .. Array Arguments ..
      int                NVAL( NN ), NSVAL( NNS ), NTVAL( NNT );
      REAL               A( * );
      REAL               AINV( * );
      REAL               ASAV( * );
      REAL               B( * );
      REAL               BSAV( * );
      REAL               AFAC( * );
      REAL               ARF( * );
      REAL               ARFINV( * );
      REAL               XACT( * );
      REAL               X( * );
      REAL               S_WORK_SLATMS( * );
      REAL               S_WORK_SPOT01( * );
      REAL               S_TEMP_SPOT02( * );
      REAL               S_TEMP_SPOT03( * );
      REAL               S_WORK_SLANSY( * );
      REAL               S_WORK_SPOT02( * );
      REAL               S_WORK_SPOT03( * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      int                NTESTS;
      const              NTESTS = 4 ;
      // ..
      // .. Local Scalars ..
      bool               ZEROT;
      int                I, INFO, IUPLO, LDA, LDB, IMAT, NERRS, NFAIL, NRHS, NRUN, IZERO, IOFF, K, NT, N, IFORM, IIN, IIT, IIS;
      String             DIST, CTYPE, UPLO, CFORM;
      int                KL, KU, MODE;
      REAL               ANORM, AINVNM, CNDNUM, RCONDC;
      // ..
      // .. Local Arrays ..
      String             UPLOS( 2 ), FORMS( 2 );
      int                ISEED( 4 ), ISEEDY( 4 );
      REAL               RESULT( NTESTS );
      // ..
      // .. External Functions ..
      //- REAL               SLANSY;
      // EXTERNAL SLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALADHD, ALAERH, ALASVM, SGET04, STFTTR, SLACPY, SLARHS, SLATB4, SLATMS, SPFTRI, SPFTRF, SPFTRS, SPOT01, SPOT02, SPOT03, SPOTRI, SPOTRF, STRTTF
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
      NFAIL = 0;
      NERRS = 0;
      for (I = 1; I <= 4; I++) { // 10
         ISEED( I ) = ISEEDY( I );
      } // 10

      for (IIN = 1; IIN <= NN; IIN++) { // 130

         N = NVAL( IIN );
         LDA = max( N, 1 );
         LDB = max( N, 1 );

         for (IIS = 1; IIS <= NNS; IIS++) { // 980

            NRHS = NSVAL( IIS );

            for (IIT = 1; IIT <= NNT; IIT++) { // 120

               IMAT = NTVAL( IIT );

               // If N == 0, only consider the first type

               if (N == 0 && IIT >= 1) GO TO 120;

               // Skip types 3, 4, or 5 if the matrix size is too small.

               if (IMAT == 4 && N <= 1) GO TO 120;
               if (IMAT == 5 && N <= 2) GO TO 120;

               // Do first for UPLO = 'U', then for UPLO = 'L'

               for (IUPLO = 1; IUPLO <= 2; IUPLO++) { // 110
                  UPLO = UPLOS( IUPLO );

                  // Do first for CFORM = 'N', then for CFORM = 'C'

                  for (IFORM = 1; IFORM <= 2; IFORM++) { // 100
                     CFORM = FORMS( IFORM );

                     // Set up parameters with SLATB4 and generate a test
                     // matrix with SLATMS.

                     slatb4('SPO', IMAT, N, N, CTYPE, KL, KU, ANORM, MODE, CNDNUM, DIST );

                     SRNAMT = 'SLATMS';
                     slatms(N, N, DIST, ISEED, CTYPE, S_WORK_SLATMS, MODE, CNDNUM, ANORM, KL, KU, UPLO, A, LDA, S_WORK_SLATMS, INFO );

                     // Check error code from SLATMS.

                     if ( INFO != 0 ) {
                        alaerh('SPF', 'SLATMS', INFO, 0, UPLO, N, N, -1, -1, -1, IIT, NFAIL, NERRS, NOUT );
                        GO TO 100;
                     }

                     // For types 3-5, zero one row and column of the matrix to
                     // test that INFO is returned correctly.

                     ZEROT = IMAT >= 3 && IMAT <= 5;
                     if ( ZEROT ) {
                        if ( IIT == 3 ) {
                           IZERO = 1;
                        } else if ( IIT == 4 ) {
                           IZERO = N;
                        } else {
                           IZERO = N / 2 + 1;
                        }
                        IOFF = ( IZERO-1 )*LDA;

                        // Set row and column IZERO of A to 0.

                        if ( IUPLO == 1 ) {
                           for (I = 1; I <= IZERO - 1; I++) { // 20
                              A( IOFF+I ) = ZERO;
                           } // 20
                           IOFF = IOFF + IZERO;
                           for (I = IZERO; I <= N; I++) { // 30
                              A( IOFF ) = ZERO;
                              IOFF = IOFF + LDA;
                           } // 30
                        } else {
                           IOFF = IZERO;
                           for (I = 1; I <= IZERO - 1; I++) { // 40
                              A( IOFF ) = ZERO;
                              IOFF = IOFF + LDA;
                           } // 40
                           IOFF = IOFF - IZERO;
                           for (I = IZERO; I <= N; I++) { // 50
                              A( IOFF+I ) = ZERO;
                           } // 50
                        }
                     } else {
                        IZERO = 0;
                     }

                     // Save a copy of the matrix A in ASAV.

                     slacpy(UPLO, N, N, A, LDA, ASAV, LDA );

                     // Compute the condition number of A (RCONDC).

                     if ( ZEROT ) {
                        RCONDC = ZERO;
                     } else {

                        // Compute the 1-norm of A.

                        ANORM = SLANSY( '1', UPLO, N, A, LDA, S_WORK_SLANSY );

                        // Factor the matrix A.

                        spotrf(UPLO, N, A, LDA, INFO );

                        // Form the inverse of A.

                        spotri(UPLO, N, A, LDA, INFO );

                        if ( N != 0 ) {

                           // Compute the 1-norm condition number of A.

                           AINVNM = SLANSY( '1', UPLO, N, A, LDA, S_WORK_SLANSY );
                           RCONDC = ( ONE / ANORM ) / AINVNM;

                           // Restore the matrix A.

                           slacpy(UPLO, N, N, ASAV, LDA, A, LDA );
                        }

                     }

                     // Form an exact solution and set the right hand side.

                     SRNAMT = 'SLARHS';
                     slarhs('SPO', 'N', UPLO, ' ', N, N, KL, KU, NRHS, A, LDA, XACT, LDA, B, LDA, ISEED, INFO );
                     slacpy('Full', N, NRHS, B, LDA, BSAV, LDA );

                     // Compute the L*L' or U'*U factorization of the
                     // matrix and solve the system.

                     slacpy(UPLO, N, N, A, LDA, AFAC, LDA );
                     slacpy('Full', N, NRHS, B, LDB, X, LDB );

                     SRNAMT = 'STRTTF';
                     strttf(CFORM, UPLO, N, AFAC, LDA, ARF, INFO );
                     SRNAMT = 'SPFTRF';
                     spftrf(CFORM, UPLO, N, ARF, INFO );

                     // Check error code from SPFTRF.

                     if ( INFO != IZERO ) {

                        // LANGOU: there is a small hick here: IZERO should
                        // always be INFO however if INFO is ZERO, ALAERH does not
                        // complain.

                         alaerh('SPF', 'SPFSV ', INFO, IZERO, UPLO, N, N, -1, -1, NRHS, IIT, NFAIL, NERRS, NOUT );
                         GO TO 100;
                      }

                     // Skip the tests if INFO is not 0.

                     if ( INFO != 0 ) {
                        GO TO 100;
                     }

                     SRNAMT = 'SPFTRS';
                     spftrs(CFORM, UPLO, N, NRHS, ARF, X, LDB, INFO );

                     SRNAMT = 'STFTTR';
                     stfttr(CFORM, UPLO, N, ARF, AFAC, LDA, INFO );

                     // Reconstruct matrix from factors and compute
                     // residual.

                     slacpy(UPLO, N, N, AFAC, LDA, ASAV, LDA );
                     spot01(UPLO, N, A, LDA, AFAC, LDA, S_WORK_SPOT01, RESULT( 1 ) );
                     slacpy(UPLO, N, N, ASAV, LDA, AFAC, LDA );

                     // Form the inverse and compute the residual.

                     if ((N % 2) == 0) {
                        slacpy('A', N+1, N/2, ARF, N+1, ARFINV, N+1 );
                     } else {
                        slacpy('A', N, (N+1)/2, ARF, N, ARFINV, N );
                     }

                     SRNAMT = 'SPFTRI';
                     spftri(CFORM, UPLO, N, ARFINV , INFO );

                     SRNAMT = 'STFTTR';
                     stfttr(CFORM, UPLO, N, ARFINV, AINV, LDA, INFO );

                     // Check error code from SPFTRI.

                     if (INFO != 0) alaerh( 'SPO', 'SPFTRI', INFO, 0, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );

                     spot03(UPLO, N, A, LDA, AINV, LDA, S_TEMP_SPOT03, LDA, S_WORK_SPOT03, RCONDC, RESULT( 2 ) );

                     // Compute residual of the computed solution.

                     slacpy('Full', N, NRHS, B, LDA, S_TEMP_SPOT02, LDA );
                     spot02(UPLO, N, NRHS, A, LDA, X, LDA, S_TEMP_SPOT02, LDA, S_WORK_SPOT02, RESULT( 3 ) );

                     // Check solution from generated exact solution.
                      sget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 4 ) );
                     NT = 4;

                     // Print information about the tests that did not
                     // pass the threshold.

                     for (K = 1; K <= NT; K++) { // 60
                        if ( RESULT( K ) >= THRESH ) {
                           if (NFAIL == 0 && NERRS == 0) aladhd( NOUT, 'SPF' );
                           WRITE( NOUT, FMT = 9999 )'SPFSV ', UPLO, N, IIT, K, RESULT( K );
                           NFAIL = NFAIL + 1;
                        }
                     } // 60
                     NRUN = NRUN + NT;
                  } // 100
               } // 110
            } // 120
         } // 980
      } // 130

      // Print a summary of the results.

      alasvm('SPF', NOUT, NFAIL, NRUN, NERRS );

 9999 FORMAT( 1X, A6, ', UPLO=''', A1, ''', N =', I5, ', type ', I1, ', test(', I1, ')=', G12.5 );

      return;
      }
