      SUBROUTINE SCHKQL( DOTYPE, NM, MVAL, NN, NVAL, NNB, NBVAL, NXVAL, NRHS, THRESH, TSTERR, NMAX, A, AF, AQ, AL, AC, B, X, XACT, TAU, WORK, RWORK, NOUT )

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
      int                MVAL( * ), NBVAL( * ), NVAL( * ), NXVAL( * );
      REAL               A( * ), AC( * ), AF( * ), AL( * ), AQ( * ), B( * ), RWORK( * ), TAU( * ), WORK( * ), X( * ), XACT( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      int                NTESTS;
      const              NTESTS = 7 ;
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
      // .. External Subroutines ..
      // EXTERNAL ALAERH, ALAHD, ALASUM, SERRQL, SGEQLS, SGET02, SLACPY, SLARHS, SLATB4, SLATMS, SQLT01, SQLT02, SQLT03, XLAENV
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
      PATH( 2: 3 ) = 'QL'
      NRUN = 0
      NFAIL = 0
      NERRS = 0
      for (I = 1; I <= 4; I++) { // 10
         ISEED( I ) = ISEEDY( I )
      } // 10

      // Test the error exits

      IF( TSTERR ) CALL SERRQL( PATH, NOUT )
      INFOT = 0
      xlaenv(2, 2 );

      LDA = NMAX
      LWORK = NMAX*MAX( NMAX, NRHS )

      // Do for each value of M in MVAL.

      for (IM = 1; IM <= NM; IM++) { // 70
         M = MVAL( IM )

         // Do for each value of N in NVAL.

         for (IN = 1; IN <= NN; IN++) { // 60
            N = NVAL( IN )
            MINMN = MIN( M, N )
            for (IMAT = 1; IMAT <= NTYPES; IMAT++) { // 50

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
               // corresponding to the call of SQLT01; other values are
               // used in the calls of SQLT02, and must not exceed MINMN.

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

               for (IK = 1; IK <= NK; IK++) { // 40
                  K = KVAL( IK )

                  // Do for each pair of values (NB,NX) in NBVAL and NXVAL.

                  for (INB = 1; INB <= NNB; INB++) { // 30
                     NB = NBVAL( INB )
                     xlaenv(1, NB );
                     NX = NXVAL( INB )
                     xlaenv(3, NX );
                     for (I = 1; I <= NTESTS; I++) {
                        RESULT( I ) = ZERO
                     END DO
                     NT = 2
                     if ( IK.EQ.1 ) {

                        // Test SGEQLF

                        sqlt01(M, N, A, AF, AQ, AL, LDA, TAU, WORK, LWORK, RWORK, RESULT( 1 ) );
                     } else if ( M.GE.N ) {

                        // Test SORGQL, using factorization
                        // returned by SQLT01

                        sqlt02(M, N, K, A, AF, AQ, AL, LDA, TAU, WORK, LWORK, RWORK, RESULT( 1 ) );
                     }
                     if ( M.GE.K ) {

                        // Test SORMQL, using factorization returned
                        // by SQLT01

                        sqlt03(M, N, K, AF, AC, AL, AQ, LDA, TAU, WORK, LWORK, RWORK, RESULT( 3 ) );
                        NT = NT + 4

                        // If M>=N and K=N, call SGEQLS to solve a system
                        // with NRHS right hand sides and compute the
                        // residual.

                        if ( K.EQ.N .AND. INB.EQ.1 ) {

                           // Generate a solution and set the right
                           // hand side.

                           SRNAMT = 'SLARHS'
                           slarhs(PATH, 'New', 'Full', 'No transpose', M, N, 0, 0, NRHS, A, LDA, XACT, LDA, B, LDA, ISEED, INFO );

                           slacpy('Full', M, NRHS, B, LDA, X, LDA );
                           SRNAMT = 'SGEQLS'
                           sgeqls(M, N, NRHS, AF, LDA, TAU, X, LDA, WORK, LWORK, INFO );

                           // Check error code from SGEQLS.

                           IF( INFO.NE.0 ) CALL ALAERH( PATH, 'SGEQLS', INFO, 0, ' ', M, N, NRHS, -1, NB, IMAT, NFAIL, NERRS, NOUT )

                           sget02('No transpose', M, N, NRHS, A, LDA, X( M-N+1 ), LDA, B, LDA, RWORK, RESULT( 7 ) );
                           NT = NT + 1
                        }
                     }

                     // Print information about the tests that did not
                     // pass the threshold.

                     for (I = 1; I <= NT; I++) { // 20
                        if ( RESULT( I ).GE.THRESH ) {
                           IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALAHD( NOUT, PATH )                            WRITE( NOUT, FMT = 9999 )M, N, K, NB, NX, IMAT, I, RESULT( I )
                           NFAIL = NFAIL + 1
                        }
                     } // 20
                     NRUN = NRUN + NT
                  } // 30
               } // 40
            } // 50
         } // 60
      } // 70

      // Print a summary of the results.

      alasum(PATH, NOUT, NFAIL, NRUN, NERRS );

 9999 FORMAT( ' M=', I5, ', N=', I5, ', K=', I5, ', NB=', I4, ', NX=', I5, ', type ', I2, ', test(', I2, ')=', G12.5 )
      RETURN

      // End of SCHKQL

      }
