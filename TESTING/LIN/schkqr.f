      SUBROUTINE SCHKQR( DOTYPE, NM, MVAL, NN, NVAL, NNB, NBVAL, NXVAL, NRHS, THRESH, TSTERR, NMAX, A, AF, AQ, AR, AC, B, X, XACT, TAU, WORK, RWORK, IWORK, NOUT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               TSTERR;
      int                NM, NMAX, NN, NNB, NOUT, NRHS;
      REAL               THRESH
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      int                IWORK( * ), MVAL( * ), NBVAL( * ), NVAL( * ), NXVAL( * );
      REAL               A( * ), AC( * ), AF( * ), AQ( * ), AR( * ), B( * ), RWORK( * ), TAU( * ), WORK( * ), X( * ), XACT( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      int                NTESTS;
      const              NTESTS = 9 ;
      int                NTYPES;
      const              NTYPES = 8 ;
      REAL               ZERO
      const              ZERO = 0.0E0 ;
      // ..
      // .. Local Scalars ..
      String             DIST, TYPE;
      String             PATH;
      int                I, IK, IM, IMAT, IN, INB, INFO, K, KL, KU, LDA, LWORK, M, MINMN, MODE, N, NB, NERRS, NFAIL, NK, NRUN, NT, NX;
      REAL               ANORM, CNDNUM
      // ..
      // .. Local Arrays ..
      int                ISEED( 4 ), ISEEDY( 4 ), KVAL( 4 );
      REAL               RESULT( NTESTS )
      // ..
      // .. External Functions ..
      bool               SGENND;
      // EXTERNAL SGENND
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAERH, ALAHD, ALASUM, SERRQR, SGELS, SGET02, SLACPY, SLARHS, SLATB4, SLATMS, SQRT01, SQRT01P, SQRT02, SQRT03, XLAENV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Scalars in Common ..
      bool               LERR, OK;
      String             SRNAMT;
      int                INFOT, NUNIT;
      // ..
      // .. Common blocks ..
      COMMON             / INFOC / INFOT, NUNIT, OK, LERR
      COMMON             / SRNAMC / SRNAMT
      // ..
      // .. Data statements ..
      DATA               ISEEDY / 1988, 1989, 1990, 1991 /
      // ..
      // .. Executable Statements ..

      // Initialize constants and the random number seed.

      PATH( 1: 1 ) = 'Single precision'
      PATH( 2: 3 ) = 'QR'
      NRUN = 0
      NFAIL = 0
      NERRS = 0
      DO 10 I = 1, 4
         ISEED( I ) = ISEEDY( I )
   10 CONTINUE

      // Test the error exits

      IF( TSTERR ) CALL SERRQR( PATH, NOUT )
      INFOT = 0
      xlaenv(2, 2 );

      LDA = NMAX
      LWORK = NMAX*MAX( NMAX, NRHS )

      // Do for each value of M in MVAL.

      DO 70 IM = 1, NM
         M = MVAL( IM )

         // Do for each value of N in NVAL.

         DO 60 IN = 1, NN
            N = NVAL( IN )
            MINMN = MIN( M, N )
            DO 50 IMAT = 1, NTYPES

               // Do the tests only if DOTYPE( IMAT ) is true.

               IF( .NOT.DOTYPE( IMAT ) ) GO TO 50

               // Set up parameters with SLATB4 and generate a test matrix
               // with SLATMS.

               slatb4(PATH, IMAT, M, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST );

               SRNAMT = 'SLATMS'
               slatms(M, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU, 'No packing', A, LDA, WORK, INFO );

               // Check error code from SLATMS.

               if ( INFO.NE.0 ) {
                  alaerh(PATH, 'SLATMS', INFO, 0, ' ', M, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );
                  GO TO 50
               }

               // Set some values for K: the first value must be MINMN,
               // corresponding to the call of SQRT01; other values are
               // used in the calls of SQRT02, and must not exceed MINMN.

               KVAL( 1 ) = MINMN
               KVAL( 2 ) = 0
               KVAL( 3 ) = 1
               KVAL( 4 ) = MINMN / 2
               if ( MINMN.EQ.0 ) {
                  NK = 1
               } else if ( MINMN.EQ.1 ) {
                  NK = 2
               } else if ( MINMN.LE.3 ) {
                  NK = 3
               } else {
                  NK = 4
               }

               // Do for each value of K in KVAL

               DO 40 IK = 1, NK
                  K = KVAL( IK )

                  // Do for each pair of values (NB,NX) in NBVAL and NXVAL.

                  DO 30 INB = 1, NNB
                     NB = NBVAL( INB )
                     xlaenv(1, NB );
                     NX = NXVAL( INB )
                     xlaenv(3, NX );
                     DO I = 1, NTESTS
                        RESULT( I ) = ZERO
                     END DO
                     NT = 2
                     if ( IK.EQ.1 ) {

                        // Test SGEQRF

                        sqrt01(M, N, A, AF, AQ, AR, LDA, TAU, WORK, LWORK, RWORK, RESULT( 1 ) );

                        // Test SGEQRFP

                        sqrt01p(M, N, A, AF, AQ, AR, LDA, TAU, WORK, LWORK, RWORK, RESULT( 8 ) );
                          IF( .NOT. SGENND( M, N, AF, LDA ) ) RESULT( 9 ) = 2*THRESH
                         NT = NT + 1
                     } else if ( M.GE.N ) {

                        // Test SORGQR, using factorization
                        // returned by SQRT01

                        sqrt02(M, N, K, A, AF, AQ, AR, LDA, TAU, WORK, LWORK, RWORK, RESULT( 1 ) );
                     }
                     if ( M.GE.K ) {

                        // Test SORMQR, using factorization returned
                        // by SQRT01

                        sqrt03(M, N, K, AF, AC, AR, AQ, LDA, TAU, WORK, LWORK, RWORK, RESULT( 3 ) );
                        NT = NT + 4

                        // If M>=N and K=N, call SGELS to solve a system
                        // with NRHS right hand sides and compute the
                        // residual.

                        if ( K.EQ.N .AND. INB.EQ.1 ) {

                           // Generate a solution and set the right
                           // hand side.

                           SRNAMT = 'SLARHS'
                           slarhs(PATH, 'New', 'Full', 'No transpose', M, N, 0, 0, NRHS, A, LDA, XACT, LDA, B, LDA, ISEED, INFO );

                           slacpy('Full', M, NRHS, B, LDA, X, LDA );

                           // Reset AF to the original matrix. SGELS
                           // factors the matrix before solving the system.

                           slacpy('Full', M, N, A, LDA, AF, LDA );

                           SRNAMT = 'SGELS'
                           sgels('No transpose', M, N, NRHS, AF, LDA, X, LDA, WORK, LWORK, INFO );

                           // Check error code from SGELS.

                           IF( INFO.NE.0 ) CALL ALAERH( PATH, 'SGELS', INFO, 0, 'N', M, N, NRHS, -1, NB, IMAT, NFAIL, NERRS, NOUT )

                           sget02('No transpose', M, N, NRHS, A, LDA, X, LDA, B, LDA, RWORK, RESULT( 7 ) );
                           NT = NT + 1
                        }
                     }

                     // Print information about the tests that did not
                     // pass the threshold.

                     DO 20 I = 1, NTESTS
                        if ( RESULT( I ).GE.THRESH ) {
                           IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALAHD( NOUT, PATH )                            WRITE( NOUT, FMT = 9999 )M, N, K, NB, NX, IMAT, I, RESULT( I )
                           NFAIL = NFAIL + 1
                        }
   20                CONTINUE
                     NRUN = NRUN + NTESTS
   30             CONTINUE
   40          CONTINUE
   50       CONTINUE
   60    CONTINUE
   70 CONTINUE

      // Print a summary of the results.

      alasum(PATH, NOUT, NFAIL, NRUN, NERRS );

 9999 FORMAT( ' M=', I5, ', N=', I5, ', K=', I5, ', NB=', I4, ', NX=', I5, ', type ', I2, ', test(', I2, ')=', G12.5 )
      RETURN

      // End of SCHKQR

      }
