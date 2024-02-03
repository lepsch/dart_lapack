      SUBROUTINE DDRVRFP( NOUT, NN, NVAL, NNS, NSVAL, NNT, NTVAL, THRESH, A, ASAV, AFAC, AINV, B, BSAV, XACT, X, ARF, ARFINV, D_WORK_DLATMS, D_WORK_DPOT01, D_TEMP_DPOT02, D_TEMP_DPOT03, D_WORK_DLANSY, D_WORK_DPOT02, D_WORK_DPOT03 )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                NN, NNS, NNT, NOUT;
      double             THRESH;
      // ..
      // .. Array Arguments ..
      int                NVAL( NN ), NSVAL( NNS ), NTVAL( NNT );
      double             A( * );
      double             AINV( * );
      double             ASAV( * );
      double             B( * );
      double             BSAV( * );
      double             AFAC( * );
      double             ARF( * );
      double             ARFINV( * );
      double             XACT( * );
      double             X( * );
      double             D_WORK_DLATMS( * );
      double             D_WORK_DPOT01( * );
      double             D_TEMP_DPOT02( * );
      double             D_TEMP_DPOT03( * );
      double             D_WORK_DLANSY( * );
      double             D_WORK_DPOT02( * );
      double             D_WORK_DPOT03( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0D+0, ZERO = 0.0D+0 ;
      int                NTESTS;
      const              NTESTS = 4 ;
      // ..
      // .. Local Scalars ..
      bool               ZEROT;
      int                I, INFO, IUPLO, LDA, LDB, IMAT, NERRS, NFAIL, NRHS, NRUN, IZERO, IOFF, K, NT, N, IFORM, IIN, IIT, IIS;
      String             DIST, CTYPE, UPLO, CFORM;
      int                KL, KU, MODE;
      double             ANORM, AINVNM, CNDNUM, RCONDC;
      // ..
      // .. Local Arrays ..
      String             UPLOS( 2 ), FORMS( 2 );
      int                ISEED( 4 ), ISEEDY( 4 );
      double             RESULT( NTESTS );
      // ..
      // .. External Functions ..
      double             DLANSY;
      // EXTERNAL DLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALADHD, ALAERH, ALASVM, DGET04, DTFTTR, DLACPY, DLARHS, DLATB4, DLATMS, DPFTRI, DPFTRF, DPFTRS, DPOT01, DPOT02, DPOT03, DPOTRI, DPOTRF, DTRTTF
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
      DATA               FORMS / 'N', 'T' /
      // ..
      // .. Executable Statements ..

      // Initialize constants and the random number seed.

      NRUN = 0
      NFAIL = 0
      NERRS = 0
      for (I = 1; I <= 4; I++) { // 10
         ISEED( I ) = ISEEDY( I )
      } // 10

      for (IIN = 1; IIN <= NN; IIN++) { // 130

         N = NVAL( IIN )
         LDA = MAX( N, 1 )
         LDB = MAX( N, 1 )

         for (IIS = 1; IIS <= NNS; IIS++) { // 980

            NRHS = NSVAL( IIS )

            for (IIT = 1; IIT <= NNT; IIT++) { // 120

               IMAT = NTVAL( IIT )

               // If N == 0, only consider the first type

               if (N == 0 && IIT.GE.1) GO TO 120;

               // Skip types 3, 4, or 5 if the matrix size is too small.

               if (IMAT == 4 && N.LE.1) GO TO 120;
               if (IMAT == 5 && N.LE.2) GO TO 120;

               // Do first for UPLO = 'U', then for UPLO = 'L'

               for (IUPLO = 1; IUPLO <= 2; IUPLO++) { // 110
                  UPLO = UPLOS( IUPLO )

                  // Do first for CFORM = 'N', then for CFORM = 'C'

                  for (IFORM = 1; IFORM <= 2; IFORM++) { // 100
                     CFORM = FORMS( IFORM )

                     // Set up parameters with DLATB4 and generate a test
                     // matrix with DLATMS.

                     dlatb4('DPO', IMAT, N, N, CTYPE, KL, KU, ANORM, MODE, CNDNUM, DIST );

                     SRNAMT = 'DLATMS'
                     dlatms(N, N, DIST, ISEED, CTYPE, D_WORK_DLATMS, MODE, CNDNUM, ANORM, KL, KU, UPLO, A, LDA, D_WORK_DLATMS, INFO );

                     // Check error code from DLATMS.

                     if ( INFO != 0 ) {
                        alaerh('DPF', 'DLATMS', INFO, 0, UPLO, N, N, -1, -1, -1, IIT, NFAIL, NERRS, NOUT );
                        GO TO 100
                     }

                     // For types 3-5, zero one row and column of the matrix to
                     // test that INFO is returned correctly.

                     ZEROT = IMAT.GE.3 && IMAT.LE.5
                     if ( ZEROT ) {
                        if ( IIT == 3 ) {
                           IZERO = 1
                        } else if ( IIT == 4 ) {
                           IZERO = N
                        } else {
                           IZERO = N / 2 + 1
                        }
                        IOFF = ( IZERO-1 )*LDA

                        // Set row and column IZERO of A to 0.

                        if ( IUPLO == 1 ) {
                           for (I = 1; I <= IZERO - 1; I++) { // 20
                              A( IOFF+I ) = ZERO
                           } // 20
                           IOFF = IOFF + IZERO
                           for (I = IZERO; I <= N; I++) { // 30
                              A( IOFF ) = ZERO
                              IOFF = IOFF + LDA
                           } // 30
                        } else {
                           IOFF = IZERO
                           for (I = 1; I <= IZERO - 1; I++) { // 40
                              A( IOFF ) = ZERO
                              IOFF = IOFF + LDA
                           } // 40
                           IOFF = IOFF - IZERO
                           for (I = IZERO; I <= N; I++) { // 50
                              A( IOFF+I ) = ZERO
                           } // 50
                        }
                     } else {
                        IZERO = 0
                     }

                     // Save a copy of the matrix A in ASAV.

                     dlacpy(UPLO, N, N, A, LDA, ASAV, LDA );

                     // Compute the condition number of A (RCONDC).

                     if ( ZEROT ) {
                        RCONDC = ZERO
                     } else {

                        // Compute the 1-norm of A.

                        ANORM = DLANSY( '1', UPLO, N, A, LDA, D_WORK_DLANSY )

                        // Factor the matrix A.

                        dpotrf(UPLO, N, A, LDA, INFO );

                        // Form the inverse of A.

                        dpotri(UPLO, N, A, LDA, INFO );

                        if ( N != 0 ) {

                           // Compute the 1-norm condition number of A.

                           AINVNM = DLANSY( '1', UPLO, N, A, LDA, D_WORK_DLANSY )
                           RCONDC = ( ONE / ANORM ) / AINVNM

                           // Restore the matrix A.

                           dlacpy(UPLO, N, N, ASAV, LDA, A, LDA );
                        }

                     }

                     // Form an exact solution and set the right hand side.

                     SRNAMT = 'DLARHS'
                     dlarhs('DPO', 'N', UPLO, ' ', N, N, KL, KU, NRHS, A, LDA, XACT, LDA, B, LDA, ISEED, INFO );
                     dlacpy('Full', N, NRHS, B, LDA, BSAV, LDA );

                     // Compute the L*L' or U'*U factorization of the
                     // matrix and solve the system.

                     dlacpy(UPLO, N, N, A, LDA, AFAC, LDA );
                     dlacpy('Full', N, NRHS, B, LDB, X, LDB );

                     SRNAMT = 'DTRTTF'
                     dtrttf(CFORM, UPLO, N, AFAC, LDA, ARF, INFO );
                     SRNAMT = 'DPFTRF'
                     dpftrf(CFORM, UPLO, N, ARF, INFO );

                     // Check error code from DPFTRF.

                     if ( INFO != IZERO ) {

                        // LANGOU: there is a small hick here: IZERO should
                        // always be INFO however if INFO is ZERO, ALAERH does not
                        // complain.

                         alaerh('DPF', 'DPFSV ', INFO, IZERO, UPLO, N, N, -1, -1, NRHS, IIT, NFAIL, NERRS, NOUT );
                         GO TO 100
                      }

                     // Skip the tests if INFO is not 0.

                     if ( INFO != 0 ) {
                        GO TO 100
                     }

                     SRNAMT = 'DPFTRS'
                     dpftrs(CFORM, UPLO, N, NRHS, ARF, X, LDB, INFO );

                     SRNAMT = 'DTFTTR'
                     dtfttr(CFORM, UPLO, N, ARF, AFAC, LDA, INFO );

                     // Reconstruct matrix from factors and compute
                     // residual.

                     dlacpy(UPLO, N, N, AFAC, LDA, ASAV, LDA );
                     dpot01(UPLO, N, A, LDA, AFAC, LDA, D_WORK_DPOT01, RESULT( 1 ) );
                     dlacpy(UPLO, N, N, ASAV, LDA, AFAC, LDA );

                     // Form the inverse and compute the residual.

                     if (MOD(N,2) == 0) {
                        dlacpy('A', N+1, N/2, ARF, N+1, ARFINV, N+1 );
                     } else {
                        dlacpy('A', N, (N+1)/2, ARF, N, ARFINV, N );
                     }

                     SRNAMT = 'DPFTRI'
                     dpftri(CFORM, UPLO, N, ARFINV , INFO );

                     SRNAMT = 'DTFTTR'
                     dtfttr(CFORM, UPLO, N, ARFINV, AINV, LDA, INFO );

                     // Check error code from DPFTRI.

                     if (INFO != 0) CALL ALAERH( 'DPO', 'DPFTRI', INFO, 0, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );

                     dpot03(UPLO, N, A, LDA, AINV, LDA, D_TEMP_DPOT03, LDA, D_WORK_DPOT03, RCONDC, RESULT( 2 ) );

                     // Compute residual of the computed solution.

                     dlacpy('Full', N, NRHS, B, LDA, D_TEMP_DPOT02, LDA );
                     dpot02(UPLO, N, NRHS, A, LDA, X, LDA, D_TEMP_DPOT02, LDA, D_WORK_DPOT02, RESULT( 3 ) );

                     // Check solution from generated exact solution.
                      dget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 4 ) );
                     NT = 4

                     // Print information about the tests that did not
                     // pass the threshold.

                     for (K = 1; K <= NT; K++) { // 60
                        if ( RESULT( K ).GE.THRESH ) {
                           if (NFAIL == 0 && NERRS == 0) CALL ALADHD( NOUT, 'DPF' )                            WRITE( NOUT, FMT = 9999 )'DPFSV ', UPLO, N, IIT, K, RESULT( K );
                           NFAIL = NFAIL + 1
                        }
                     } // 60
                     NRUN = NRUN + NT
                  } // 100
               } // 110
            } // 120
         } // 980
      } // 130

      // Print a summary of the results.

      alasvm('DPF', NOUT, NFAIL, NRUN, NERRS );

 9999 FORMAT( 1X, A6, ', UPLO=''', A1, ''', N =', I5, ', type ', I1, ', test(', I1, ')=', G12.5 )

      RETURN

      // End of DDRVRFP

      }
