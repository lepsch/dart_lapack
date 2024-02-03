      SUBROUTINE ZDRVRFP( NOUT, NN, NVAL, NNS, NSVAL, NNT, NTVAL, THRESH, A, ASAV, AFAC, AINV, B, BSAV, XACT, X, ARF, ARFINV, Z_WORK_ZLATMS, Z_WORK_ZPOT02, Z_WORK_ZPOT03, D_WORK_ZLATMS, D_WORK_ZLANHE, D_WORK_ZPOT01, D_WORK_ZPOT02, D_WORK_ZPOT03 )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                NN, NNS, NNT, NOUT;
      double             THRESH;
      // ..
      // .. Array Arguments ..
      int                NVAL( NN ), NSVAL( NNS ), NTVAL( NNT );
      COMPLEX*16         A( * )
      COMPLEX*16         AINV( * )
      COMPLEX*16         ASAV( * )
      COMPLEX*16         B( * )
      COMPLEX*16         BSAV( * )
      COMPLEX*16         AFAC( * )
      COMPLEX*16         ARF( * )
      COMPLEX*16         ARFINV( * )
      COMPLEX*16         XACT( * )
      COMPLEX*16         X( * )
      COMPLEX*16         Z_WORK_ZLATMS( * )
      COMPLEX*16         Z_WORK_ZPOT02( * )
      COMPLEX*16         Z_WORK_ZPOT03( * )
      double             D_WORK_ZLATMS( * );
      double             D_WORK_ZLANHE( * );
      double             D_WORK_ZPOT01( * );
      double             D_WORK_ZPOT02( * );
      double             D_WORK_ZPOT03( * );
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
      double             ZLANHE;
      // EXTERNAL ZLANHE
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALADHD, ALAERH, ALASVM, ZGET04, ZTFTTR, ZLACPY, ZLAIPD, ZLARHS, ZLATB4, ZLATMS, ZPFTRI, ZPFTRF, ZPFTRS, ZPOT01, ZPOT02, ZPOT03, ZPOTRI, ZPOTRF, ZTRTTF
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

               // If N.EQ.0, only consider the first type

               IF( N.EQ.0 .AND. IIT.GE.1 ) GO TO 120

               // Skip types 3, 4, or 5 if the matrix size is too small.

               IF( IMAT.EQ.4 .AND. N.LE.1 ) GO TO 120
               IF( IMAT.EQ.5 .AND. N.LE.2 ) GO TO 120

               // Do first for UPLO = 'U', then for UPLO = 'L'

               for (IUPLO = 1; IUPLO <= 2; IUPLO++) { // 110
                  UPLO = UPLOS( IUPLO )

                  // Do first for CFORM = 'N', then for CFORM = 'C'

                  for (IFORM = 1; IFORM <= 2; IFORM++) { // 100
                     CFORM = FORMS( IFORM )

                     // Set up parameters with ZLATB4 and generate a test
                     // matrix with ZLATMS.

                     zlatb4('ZPO', IMAT, N, N, CTYPE, KL, KU, ANORM, MODE, CNDNUM, DIST );

                     SRNAMT = 'ZLATMS'
                     zlatms(N, N, DIST, ISEED, CTYPE, D_WORK_ZLATMS, MODE, CNDNUM, ANORM, KL, KU, UPLO, A, LDA, Z_WORK_ZLATMS, INFO );

                     // Check error code from ZLATMS.

                     if ( INFO.NE.0 ) {
                        alaerh('ZPF', 'ZLATMS', INFO, 0, UPLO, N, N, -1, -1, -1, IIT, NFAIL, NERRS, NOUT );
                        GO TO 100
                     }

                     // For types 3-5, zero one row and column of the matrix to
                     // test that INFO is returned correctly.

                     ZEROT = IMAT.GE.3 .AND. IMAT.LE.5
                     if ( ZEROT ) {
                        if ( IIT.EQ.3 ) {
                           IZERO = 1
                        } else if ( IIT.EQ.4 ) {
                           IZERO = N
                        } else {
                           IZERO = N / 2 + 1
                        }
                        IOFF = ( IZERO-1 )*LDA

                        // Set row and column IZERO of A to 0.

                        if ( IUPLO.EQ.1 ) {
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

                     // Set the imaginary part of the diagonals.

                     zlaipd(N, A, LDA+1, 0 );

                     // Save a copy of the matrix A in ASAV.

                     zlacpy(UPLO, N, N, A, LDA, ASAV, LDA );

                     // Compute the condition number of A (RCONDC).

                     if ( ZEROT ) {
                        RCONDC = ZERO
                     } else {

                        // Compute the 1-norm of A.

                        ANORM = ZLANHE( '1', UPLO, N, A, LDA, D_WORK_ZLANHE )

                        // Factor the matrix A.

                        zpotrf(UPLO, N, A, LDA, INFO );

                        // Form the inverse of A.

                        zpotri(UPLO, N, A, LDA, INFO );

                        if ( N .NE. 0 ) {

                           // Compute the 1-norm condition number of A.

                           AINVNM = ZLANHE( '1', UPLO, N, A, LDA, D_WORK_ZLANHE )
                           RCONDC = ( ONE / ANORM ) / AINVNM

                           // Restore the matrix A.

                           zlacpy(UPLO, N, N, ASAV, LDA, A, LDA );
                        }

                     }

                     // Form an exact solution and set the right hand side.

                     SRNAMT = 'ZLARHS'
                     zlarhs('ZPO', 'N', UPLO, ' ', N, N, KL, KU, NRHS, A, LDA, XACT, LDA, B, LDA, ISEED, INFO );
                     zlacpy('Full', N, NRHS, B, LDA, BSAV, LDA );

                     // Compute the L*L' or U'*U factorization of the
                     // matrix and solve the system.

                     zlacpy(UPLO, N, N, A, LDA, AFAC, LDA );
                     zlacpy('Full', N, NRHS, B, LDB, X, LDB );

                     SRNAMT = 'ZTRTTF'
                     ztrttf(CFORM, UPLO, N, AFAC, LDA, ARF, INFO );
                     SRNAMT = 'ZPFTRF'
                     zpftrf(CFORM, UPLO, N, ARF, INFO );

                     // Check error code from ZPFTRF.

                     if ( INFO.NE.IZERO ) {

                        // LANGOU: there is a small hick here: IZERO should
                        // always be INFO however if INFO is ZERO, ALAERH does not
                        // complain.

                         alaerh('ZPF', 'ZPFSV ', INFO, IZERO, UPLO, N, N, -1, -1, NRHS, IIT, NFAIL, NERRS, NOUT );
                         GO TO 100
                      }

                      // Skip the tests if INFO is not 0.

                     if ( INFO.NE.0 ) {
                        GO TO 100
                     }

                     SRNAMT = 'ZPFTRS'
                     zpftrs(CFORM, UPLO, N, NRHS, ARF, X, LDB, INFO );

                     SRNAMT = 'ZTFTTR'
                     ztfttr(CFORM, UPLO, N, ARF, AFAC, LDA, INFO );

                     // Reconstruct matrix from factors and compute
                     // residual.

                     zlacpy(UPLO, N, N, AFAC, LDA, ASAV, LDA );
                     zpot01(UPLO, N, A, LDA, AFAC, LDA, D_WORK_ZPOT01, RESULT( 1 ) );
                     zlacpy(UPLO, N, N, ASAV, LDA, AFAC, LDA );

                     // Form the inverse and compute the residual.

                    if (MOD(N,2).EQ.0) {
                       zlacpy('A', N+1, N/2, ARF, N+1, ARFINV, N+1 );
                    } else {
                       zlacpy('A', N, (N+1)/2, ARF, N, ARFINV, N );
                    }

                     SRNAMT = 'ZPFTRI'
                     zpftri(CFORM, UPLO, N, ARFINV , INFO );

                     SRNAMT = 'ZTFTTR'
                     ztfttr(CFORM, UPLO, N, ARFINV, AINV, LDA, INFO );

                     // Check error code from ZPFTRI.

                     IF( INFO.NE.0 ) CALL ALAERH( 'ZPO', 'ZPFTRI', INFO, 0, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT )

                     zpot03(UPLO, N, A, LDA, AINV, LDA, Z_WORK_ZPOT03, LDA, D_WORK_ZPOT03, RCONDC, RESULT( 2 ) );

                     // Compute residual of the computed solution.

                     zlacpy('Full', N, NRHS, B, LDA, Z_WORK_ZPOT02, LDA )                      CALL ZPOT02( UPLO, N, NRHS, A, LDA, X, LDA, Z_WORK_ZPOT02, LDA, D_WORK_ZPOT02, RESULT( 3 ) );

                     // Check solution from generated exact solution.

                     zget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 4 ) );
                     NT = 4

                     // Print information about the tests that did not
                     // pass the threshold.

                     for (K = 1; K <= NT; K++) { // 60
                        if ( RESULT( K ).GE.THRESH ) {
                           IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALADHD( NOUT, 'ZPF' )                            WRITE( NOUT, FMT = 9999 )'ZPFSV ', UPLO, N, IIT, K, RESULT( K )
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

      alasvm('ZPF', NOUT, NFAIL, NRUN, NERRS );

 9999 FORMAT( 1X, A6, ', UPLO=''', A1, ''', N =', I5, ', type ', I1, ', test(', I1, ')=', G12.5 )

      RETURN

      // End of ZDRVRFP

      }
