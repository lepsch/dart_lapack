      void zdrvrfp(final Nout NOUT, NN, final int NVAL, final int NNS, final Array<int> NSVAL_, final int NNT, final int NTVAL, final int THRESH, final int A, final int ASAV, final int AFAC, final int AINV, final Array<Complex> B_, final Array<Complex> BSAV_, final int XACT, final int X, final int ARF, final int ARFINV, final int Z_WORK_ZLATMS, final int Z_WORK_ZPOT02, final int Z_WORK_ZPOT03, final int D_WORK_ZLATMS, final int D_WORK_ZLANHE, final int D_WORK_ZPOT01, final int D_WORK_ZPOT02, final int D_WORK_ZPOT03,) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                NN, NNS, NNT, NOUT;
      double             THRESH;
      int                NVAL( NN ), NSVAL( NNS ), NTVAL( NNT );
      Complex         A( * );
      Complex         AINV( * );
      Complex         ASAV( * );
      Complex         B( * );
      Complex         BSAV( * );
      Complex         AFAC( * );
      Complex         ARF( * );
      Complex         ARFINV( * );
      Complex         XACT( * );
      Complex         X( * );
      Complex         Z_WORK_ZLATMS( * );
      Complex         Z_WORK_ZPOT02( * );
      Complex         Z_WORK_ZPOT03( * );
      double             D_WORK_ZLATMS( * );
      double             D_WORK_ZLANHE( * );
      double             D_WORK_ZPOT01( * );
      double             D_WORK_ZPOT02( * );
      double             D_WORK_ZPOT03( * );
      // ..

      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      int                NTESTS;
      const              NTESTS = 4 ;
      bool               ZEROT;
      int                I, INFO, IUPLO, LDA, LDB, IMAT, NRHS, IZERO, IOFF, K, NT, N, IFORM, IIN, IIT, IIS;
      String             DIST, CTYPE, UPLO, CFORM;
      int                KL, KU, MODE;
      double             ANORM, AINVNM, CNDNUM, RCONDC;
      String             UPLOS( 2 ), FORMS( 2 );
      final                ISEED=Array<int>( 4 );
      final             RESULT=Array<double>( NTESTS );
      // ..
      // .. External Functions ..
      //- double             ZLANHE;
      // EXTERNAL ZLANHE
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALADHD, ALAERH, ALASVM, ZGET04, ZTFTTR, ZLACPY, ZLAIPD, ZLARHS, ZLATB4, ZLATMS, ZPFTRI, ZPFTRF, ZPFTRS, ZPOT01, ZPOT02, ZPOT03, ZPOTRI, ZPOTRF, ZTRTTF
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
      var NFAIL = 0;
      var NERRS = Box(0);
      for (I = 1; I <= 4; I++) { // 10
         ISEED[I] = ISEEDY[I - 1];
      } // 10

      for (IIN = 1; IIN <= NN; IIN++) { // 130

         N = NVAL( IIN );
         final LDA = max( N, 1 );
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
                  final UPLO = UPLOS[IUPLO - 1];

                  // Do first for CFORM = 'N', then for CFORM = 'C'

                  for (IFORM = 1; IFORM <= 2; IFORM++) { // 100
                     CFORM = FORMS( IFORM );

                     // Set up parameters with ZLATB4 and generate a test
                     // matrix with ZLATMS.

                     zlatb4('ZPO', IMAT, N, N, CTYPE, KL, KU, ANORM, MODE, CNDNUM, DIST );

                    srnamc.SRNAMT = 'ZLATMS';
                     zlatms(N, N, DIST, ISEED, CTYPE, D_WORK_ZLATMS, MODE, CNDNUM, ANORM, KL, KU, UPLO, A, LDA, Z_WORK_ZLATMS, INFO );

                     // Check error code from ZLATMS.

                     if ( INFO.value != 0 ) {
                        alaerh('ZPF', 'ZLATMS', INFO, 0, UPLO, N, N, -1, -1, -1, IIT, NFAIL, NERRS, NOUT );
                        GO TO 100;
                     }

                     // For types 3-5, zero one row and column of the matrix to
                     // test that INFO is returned correctly.

                     final ZEROT = IMAT >= 3 && IMAT <= 5;
                     if ( ZEROT ) {
                        if ( IIT == 3 ) {
                           IZERO = 1;
                        } else if ( IIT == 4 ) {
                           IZERO = N;
                        } else {
                           IZERO = N ~/ 2 + 1;
                        }
                        IOFF = ( IZERO-1 )*LDA;

                        // Set row and column IZERO of A to 0.

                        if ( IUPLO == 1 ) {
                           for (I = 1; I <= IZERO - 1; I++) { // 20
                              A[IOFF+I] = ZERO;
                           } // 20
                           IOFF += IZERO;
                           for (I = IZERO; I <= N; I++) { // 30
                              A[IOFF] = ZERO;
                              IOFF += LDA;
                           } // 30
                        } else {
                           IOFF = IZERO;
                           for (I = 1; I <= IZERO - 1; I++) { // 40
                              A[IOFF] = ZERO;
                              IOFF += LDA;
                           } // 40
                           IOFF -= IZERO;
                           for (I = IZERO; I <= N; I++) { // 50
                              A[IOFF+I] = ZERO;
                           } // 50
                        }
                     } else {
                        IZERO = 0;
                     }

                     // Set the imaginary part of the diagonals.

                     zlaipd(N, A, LDA+1, 0 );

                     // Save a copy of the matrix A in ASAV.

                     zlacpy(UPLO, N, N, A, LDA, ASAV, LDA );

                     // Compute the condition number of A (RCONDC).

                     final int IZERO;
                     if ( ZEROT ) {
                        RCONDC = ZERO;
                     } else {

                        // Compute the 1-norm of A.

                        ANORM = ZLANHE( '1', UPLO, N, A, LDA, D_WORK_ZLANHE );

                        // Factor the matrix A.

                        zpotrf(UPLO, N, A, LDA, INFO );

                        // Form the inverse of A.

                        zpotri(UPLO, N, A, LDA, INFO );

                        if ( N != 0 ) {

                           // Compute the 1-norm condition number of A.

                           AINVNM = ZLANHE( '1', UPLO, N, A, LDA, D_WORK_ZLANHE );
                           RCONDC = ( ONE / ANORM ) / AINVNM;

                           // Restore the matrix A.

                           zlacpy(UPLO, N, N, ASAV, LDA, A, LDA );
                        }

                     }

                     // Form an exact solution and set the right hand side.

                    srnamc.SRNAMT = 'ZLARHS';
                     zlarhs('ZPO', 'N', UPLO, ' ', N, N, KL, KU, NRHS, A, LDA, XACT, LDA, B, LDA, ISEED, INFO );
                     zlacpy('Full', N, NRHS, B, LDA, BSAV, LDA );

                     // Compute the L*L' or U'*U factorization of the
                     // matrix and solve the system.

                     zlacpy(UPLO, N, N, A, LDA, AFAC, LDA );
                     zlacpy('Full', N, NRHS, B, LDB, X, LDB );

                    srnamc.SRNAMT = 'ZTRTTF';
                     ztrttf(CFORM, UPLO, N, AFAC, LDA, ARF, INFO );
                    srnamc.SRNAMT = 'ZPFTRF';
                     zpftrf(CFORM, UPLO, N, ARF, INFO );

                     // Check error code from ZPFTRF.

                     if ( INFO.value != IZERO ) {

                        // LANGOU: there is a small hick here: IZERO should
                        // always be INFO however if INFO is ZERO, ALAERH does not
                        // complain.

                         alaerh('ZPF', 'ZPFSV ', INFO, IZERO, UPLO, N, N, -1, -1, NRHS, IIT, NFAIL, NERRS, NOUT );
                         GO TO 100;
                      }

                      // Skip the tests if INFO is not 0.

                     if ( INFO.value != 0 ) {
                        GO TO 100;
                     }

                    srnamc.SRNAMT = 'ZPFTRS';
                     zpftrs(CFORM, UPLO, N, NRHS, ARF, X, LDB, INFO );

                    srnamc.SRNAMT = 'ZTFTTR';
                     ztfttr(CFORM, UPLO, N, ARF, AFAC, LDA, INFO );

                     // Reconstruct matrix from factors and compute
                     // residual.

                     zlacpy(UPLO, N, N, AFAC, LDA, ASAV, LDA );
                     zpot01(UPLO, N, A, LDA, AFAC, LDA, D_WORK_ZPOT01, RESULT( 1 ) );
                     zlacpy(UPLO, N, N, ASAV, LDA, AFAC, LDA );

                     // Form the inverse and compute the residual.

                    if ((N % 2) == 0) {
                       zlacpy('A', N+1, N/2, ARF, N+1, ARFINV, N+1 );
                    } else {
                       zlacpy('A', N, (N+1)/2, ARF, N, ARFINV, N );
                    }

                    srnamc.SRNAMT = 'ZPFTRI';
                     zpftri(CFORM, UPLO, N, ARFINV , INFO );

                    srnamc.SRNAMT = 'ZTFTTR';
                     ztfttr(CFORM, UPLO, N, ARFINV, AINV, LDA, INFO );

                     // Check error code from ZPFTRI.

                     if (INFO != 0) alaerh( 'ZPO', 'ZPFTRI', INFO, 0, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );

                     zpot03(UPLO, N, A, LDA, AINV, LDA, Z_WORK_ZPOT03, LDA, D_WORK_ZPOT03, RCONDC, RESULT( 2 ) );

                     // Compute residual of the computed solution.

                     zlacpy('Full', N, NRHS, B, LDA, Z_WORK_ZPOT02, LDA );
                     zpot02(UPLO, N, NRHS, A, LDA, X, LDA, Z_WORK_ZPOT02, LDA, D_WORK_ZPOT02, RESULT( 3 ) );

                     // Check solution from generated exact solution.

                     zget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 4 ) );
                     NT = 4;

                     // Print information about the tests that did not
                     // pass the threshold.

                     for (K = 1; K <= NT; K++) { // 60
                        if ( RESULT[K] >= THRESH ) {
                           if (NFAIL == 0 && NERRS.value == 0) aladhd( NOUT, 'ZPF' );
                           NOUT.println( 9999 )'ZPFSV ', UPLO, N, IIT, K, RESULT( K );
                           NFAIL++;
                        }
                     } // 60
                     NRUN +=  NT;
                  } // 100
               } // 110
            } // 120
         } // 980
      } // 130

      // Print a summary of the results.

      alasvm('ZPF', NOUT, NFAIL, NRUN, NERRS );

 9999 FORMAT(' ${.a6}, UPLO=\'${.a1}\', N =${N.i5}, type ${IMAT.i1}, test(${.i1})=${RESULT[].g12_5};

      }
