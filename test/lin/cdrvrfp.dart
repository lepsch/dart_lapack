      void cdrvrfp(NOUT, NN, NVAL, NNS, NSVAL, NNT, NTVAL, THRESH, A, ASAV, AFAC, AINV, B, BSAV, XACT, X, ARF, ARFINV, C_WORK_CLATMS, C_WORK_CPOT02, C_WORK_CPOT03, S_WORK_CLATMS, S_WORK_CLANHE, S_WORK_CPOT01, S_WORK_CPOT02, S_WORK_CPOT03 ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                NN, NNS, NNT, NOUT;
      double               THRESH;
      int                NVAL( NN ), NSVAL( NNS ), NTVAL( NNT );
      Complex            A( * );
      Complex            AINV( * );
      Complex            ASAV( * );
      Complex            B( * );
      Complex            BSAV( * );
      Complex            AFAC( * );
      Complex            ARF( * );
      Complex            ARFINV( * );
      Complex            XACT( * );
      Complex            X( * );
      Complex            C_WORK_CLATMS( * );
      Complex            C_WORK_CPOT02( * );
      Complex            C_WORK_CPOT03( * );
      double               S_WORK_CLATMS( * );
      double               S_WORK_CLANHE( * );
      double               S_WORK_CPOT01( * );
      double               S_WORK_CPOT02( * );
      double               S_WORK_CPOT03( * );
      // ..

      double               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      int                NTESTS;
      const              NTESTS = 4 ;
      bool               ZEROT;
      int                I, INFO, IUPLO, LDA, LDB, IMAT, NERRS, NFAIL, NRHS, NRUN, IZERO, IOFF, K, NT, N, IFORM, IIN, IIT, IIS;
      String             DIST, CTYPE, UPLO, CFORM;
      int                KL, KU, MODE;
      double               ANORM, AINVNM, CNDNUM, RCONDC;
      String             UPLOS( 2 ), FORMS( 2 );
      int                ISEED( 4 ), ISEEDY( 4 );
      double               RESULT( NTESTS );
      // ..
      // .. External Functions ..
      //- REAL               CLANHE;
      // EXTERNAL CLANHE
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALADHD, ALAERH, ALASVM, CGET04, CTFTTR, CLACPY, CLAIPD, CLARHS, CLATB4, CLATMS, CPFTRI, CPFTRF, CPFTRS, CPOT01, CPOT02, CPOT03, CPOTRI, CPOTRF, CTRTTF
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

      NRUN = 0;
      NFAIL = 0;
      NERRS = 0;
      for (I = 1; I <= 4; I++) { // 10
         ISEED[I] = ISEEDY( I );
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

                     // Set up parameters with CLATB4 and generate a test
                     // matrix with CLATMS.

                     clatb4('CPO', IMAT, N, N, CTYPE, KL, KU, ANORM, MODE, CNDNUM, DIST );

                    srnamc.SRNAMT = 'CLATMS';
                     clatms(N, N, DIST, ISEED, CTYPE, S_WORK_CLATMS, MODE, CNDNUM, ANORM, KL, KU, UPLO, A, LDA, C_WORK_CLATMS, INFO );

                     // Check error code from CLATMS.

                     if ( INFO != 0 ) {
                        alaerh('CPF', 'CLATMS', INFO, 0, UPLO, N, N, -1, -1, -1, IIT, NFAIL, NERRS, NOUT );
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
                              A[IOFF+I] = ZERO;
                           } // 20
                           IOFF = IOFF + IZERO;
                           for (I = IZERO; I <= N; I++) { // 30
                              A[IOFF] = ZERO;
                              IOFF = IOFF + LDA;
                           } // 30
                        } else {
                           IOFF = IZERO;
                           for (I = 1; I <= IZERO - 1; I++) { // 40
                              A[IOFF] = ZERO;
                              IOFF = IOFF + LDA;
                           } // 40
                           IOFF = IOFF - IZERO;
                           for (I = IZERO; I <= N; I++) { // 50
                              A[IOFF+I] = ZERO;
                           } // 50
                        }
                     } else {
                        IZERO = 0;
                     }

                     // Set the imaginary part of the diagonals.

                     claipd(N, A, LDA+1, 0 );

                     // Save a copy of the matrix A in ASAV.

                     clacpy(UPLO, N, N, A, LDA, ASAV, LDA );

                     // Compute the condition number of A (RCONDC).

                     if ( ZEROT ) {
                        RCONDC = ZERO;
                     } else {

                        // Compute the 1-norm of A.

                        ANORM = CLANHE( '1', UPLO, N, A, LDA, S_WORK_CLANHE );

                        // Factor the matrix A.

                        cpotrf(UPLO, N, A, LDA, INFO );

                        // Form the inverse of A.

                        cpotri(UPLO, N, A, LDA, INFO );

                        if ( N != 0 ) {

                           // Compute the 1-norm condition number of A.

                           AINVNM = CLANHE( '1', UPLO, N, A, LDA, S_WORK_CLANHE );
                           RCONDC = ( ONE / ANORM ) / AINVNM;

                           // Restore the matrix A.

                           clacpy(UPLO, N, N, ASAV, LDA, A, LDA );
                        }

                     }

                     // Form an exact solution and set the right hand side.

                    srnamc.SRNAMT = 'CLARHS';
                     clarhs('CPO', 'N', UPLO, ' ', N, N, KL, KU, NRHS, A, LDA, XACT, LDA, B, LDA, ISEED, INFO );
                     clacpy('Full', N, NRHS, B, LDA, BSAV, LDA );

                     // Compute the L*L' or U'*U factorization of the
                     // matrix and solve the system.

                     clacpy(UPLO, N, N, A, LDA, AFAC, LDA );
                     clacpy('Full', N, NRHS, B, LDB, X, LDB );

                    srnamc.SRNAMT = 'CTRTTF';
                     ctrttf(CFORM, UPLO, N, AFAC, LDA, ARF, INFO );
                    srnamc.SRNAMT = 'CPFTRF';
                     cpftrf(CFORM, UPLO, N, ARF, INFO );

                     // Check error code from CPFTRF.

                     if ( INFO != IZERO ) {

                        // LANGOU: there is a small hick here: IZERO should
                        // always be INFO however if INFO is ZERO, ALAERH does not
                        // complain.

                         alaerh('CPF', 'CPFSV ', INFO, IZERO, UPLO, N, N, -1, -1, NRHS, IIT, NFAIL, NERRS, NOUT );
                         GO TO 100;
                      }

                      // Skip the tests if INFO is not 0.

                     if ( INFO != 0 ) {
                        GO TO 100;
                     }

                    srnamc.SRNAMT = 'CPFTRS';
                     cpftrs(CFORM, UPLO, N, NRHS, ARF, X, LDB, INFO );

                    srnamc.SRNAMT = 'CTFTTR';
                     ctfttr(CFORM, UPLO, N, ARF, AFAC, LDA, INFO );

                     // Reconstruct matrix from factors and compute
                     // residual.

                     clacpy(UPLO, N, N, AFAC, LDA, ASAV, LDA );
                     cpot01(UPLO, N, A, LDA, AFAC, LDA, S_WORK_CPOT01, RESULT( 1 ) );
                     clacpy(UPLO, N, N, ASAV, LDA, AFAC, LDA );

                     // Form the inverse and compute the residual.

                    if ((N % 2) == 0) {
                       clacpy('A', N+1, N/2, ARF, N+1, ARFINV, N+1 );
                    } else {
                       clacpy('A', N, (N+1)/2, ARF, N, ARFINV, N );
                    }

                    srnamc.SRNAMT = 'CPFTRI';
                     cpftri(CFORM, UPLO, N, ARFINV , INFO );

                    srnamc.SRNAMT = 'CTFTTR';
                     ctfttr(CFORM, UPLO, N, ARFINV, AINV, LDA, INFO );

                     // Check error code from CPFTRI.

                     if (INFO != 0) alaerh( 'CPO', 'CPFTRI', INFO, 0, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );

                     cpot03(UPLO, N, A, LDA, AINV, LDA, C_WORK_CPOT03, LDA, S_WORK_CPOT03, RCONDC, RESULT( 2 ) );

                     // Compute residual of the computed solution.

                     clacpy('Full', N, NRHS, B, LDA, C_WORK_CPOT02, LDA );
                     cpot02(UPLO, N, NRHS, A, LDA, X, LDA, C_WORK_CPOT02, LDA, S_WORK_CPOT02, RESULT( 3 ) );

                     // Check solution from generated exact solution.

                     cget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 4 ) );
                     NT = 4;

                     // Print information about the tests that did not
                     // pass the threshold.

                     for (K = 1; K <= NT; K++) { // 60
                        if ( RESULT( K ) >= THRESH ) {
                           if (NFAIL == 0 && NERRS == 0) aladhd( NOUT, 'CPF' );
                           WRITE( NOUT, FMT = 9999 )'CPFSV ', UPLO, N, IIT, K, RESULT( K );
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

      alasvm('CPF', NOUT, NFAIL, NRUN, NERRS );

 9999 FORMAT(' ${.a6}, UPLO=''${.a1}'', N =${.i5}, type ${.i1}, test(${.i1})=${.g12_5};

      }
