      SUBROUTINE ZCHKGE( DOTYPE, NM, MVAL, NN, NVAL, NNB, NBVAL, NNS, NSVAL, THRESH, TSTERR, NMAX, A, AFAC, AINV, B, X, XACT, WORK, RWORK, IWORK, NOUT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               TSTERR;
      int                NM, NMAX, NN, NNB, NNS, NOUT;
      double             THRESH;
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      int                IWORK( * ), MVAL( * ), NBVAL( * ), NSVAL( * ), NVAL( * );
      double             RWORK( * );
      COMPLEX*16         A( * ), AFAC( * ), AINV( * ), B( * ), WORK( * ), X( * ), XACT( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0D+0, ZERO = 0.0D+0 ;
      int                NTYPES;
      const              NTYPES = 11 ;
      int                NTESTS;
      const              NTESTS = 8 ;
      int                NTRAN;
      const              NTRAN = 3 ;
      // ..
      // .. Local Scalars ..
      bool               TRFCON, ZEROT;
      String             DIST, NORM, TRANS, TYPE, XTYPE;
      String             PATH;
      int                I, IM, IMAT, IN, INB, INFO, IOFF, IRHS, ITRAN, IZERO, K, KL, KU, LDA, LWORK, M, MODE, N, NB, NERRS, NFAIL, NIMAT, NRHS, NRUN, NT;
      double             AINVNM, ANORM, ANORMI, ANORMO, CNDNUM, DUMMY, RCOND, RCONDC, RCONDI, RCONDO;
      // ..
      // .. Local Arrays ..
      String             TRANSS( NTRAN );
      int                ISEED( 4 ), ISEEDY( 4 );
      double             RESULT( NTESTS );
      // ..
      // .. External Functions ..
      double             DGET06, ZLANGE;
      // EXTERNAL DGET06, ZLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAERH, ALAHD, ALASUM, XLAENV, ZERRGE, ZGECON, ZGERFS, ZGET01, ZGET02, ZGET03, ZGET04, ZGET07, ZGETRF, ZGETRI, ZGETRS, ZLACPY, ZLARHS, ZLASET, ZLATB4, ZLATMS
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DCMPLX, MAX, MIN
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
      DATA               ISEEDY / 1988, 1989, 1990, 1991 / , TRANSS / 'N', 'T', 'C' /
      // ..
      // .. Executable Statements ..

      // Initialize constants and the random number seed.

      PATH( 1: 1 ) = 'Zomplex precision'
      PATH( 2: 3 ) = 'GE'
      NRUN = 0
      NFAIL = 0
      NERRS = 0
      DO 10 I = 1, 4
         ISEED( I ) = ISEEDY( I )
   10 CONTINUE

      // Test the error exits

      CALL XLAENV( 1, 1 )
      IF( TSTERR ) CALL ZERRGE( PATH, NOUT )
      INFOT = 0
      CALL XLAENV( 2, 2 )

      // Do for each value of M in MVAL

      DO 120 IM = 1, NM
         M = MVAL( IM )
         LDA = MAX( 1, M )

         // Do for each value of N in NVAL

         DO 110 IN = 1, NN
            N = NVAL( IN )
            XTYPE = 'N'
            NIMAT = NTYPES
            IF( M.LE.0 .OR. N.LE.0 ) NIMAT = 1

            DO 100 IMAT = 1, NIMAT

               // Do the tests only if DOTYPE( IMAT ) is true.

               IF( .NOT.DOTYPE( IMAT ) ) GO TO 100

               // Skip types 5, 6, or 7 if the matrix size is too small.

               ZEROT = IMAT.GE.5 .AND. IMAT.LE.7
               IF( ZEROT .AND. N.LT.IMAT-4 ) GO TO 100

               // Set up parameters with ZLATB4 and generate a test matrix
               // with ZLATMS.

               CALL ZLATB4( PATH, IMAT, M, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST )

               SRNAMT = 'ZLATMS'
               CALL ZLATMS( M, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU, 'No packing', A, LDA, WORK, INFO )

               // Check error code from ZLATMS.

               if ( INFO.NE.0 ) {
                  CALL ALAERH( PATH, 'ZLATMS', INFO, 0, ' ', M, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT )
                  GO TO 100
               }

               // For types 5-7, zero one or more columns of the matrix to
              t // est that INFO is returned correctly.

               if ( ZEROT ) {
                  if ( IMAT.EQ.5 ) {
                     IZERO = 1
                  } else if ( IMAT.EQ.6 ) {
                     IZERO = MIN( M, N )
                  } else {
                     IZERO = MIN( M, N ) / 2 + 1
                  }
                  IOFF = ( IZERO-1 )*LDA
                  if ( IMAT.LT.7 ) {
                     DO 20 I = 1, M
                        A( IOFF+I ) = ZERO
   20                CONTINUE
                  } else {
                     CALL ZLASET( 'Full', M, N-IZERO+1, DCMPLX( ZERO ), DCMPLX( ZERO ), A( IOFF+1 ), LDA )
                  }
               } else {
                  IZERO = 0
               }

               // These lines, if used in place of the calls in the DO 60
               // loop, cause the code to bomb on a Sun SPARCstation.

                // ANORMO = ZLANGE( 'O', M, N, A, LDA, RWORK )
                // ANORMI = ZLANGE( 'I', M, N, A, LDA, RWORK )

               // Do for each blocksize in NBVAL

               DO 90 INB = 1, NNB
                  NB = NBVAL( INB )
                  CALL XLAENV( 1, NB )

                  // Compute the LU factorization of the matrix.

                  CALL ZLACPY( 'Full', M, N, A, LDA, AFAC, LDA )
                  SRNAMT = 'ZGETRF'
                  CALL ZGETRF( M, N, AFAC, LDA, IWORK, INFO )

                  // Check error code from ZGETRF.

                  IF( INFO.NE.IZERO ) CALL ALAERH( PATH, 'ZGETRF', INFO, IZERO, ' ', M, N, -1, -1, NB, IMAT, NFAIL, NERRS, NOUT )
                  TRFCON = .FALSE.

*+    TEST 1
                  // Reconstruct matrix from factors and compute residual.

                  CALL ZLACPY( 'Full', M, N, AFAC, LDA, AINV, LDA )
                  CALL ZGET01( M, N, A, LDA, AINV, LDA, IWORK, RWORK, RESULT( 1 ) )
                  NT = 1

*+    TEST 2
                  // Form the inverse if the factorization was successful
                  // and compute the residual.

                  if ( M.EQ.N .AND. INFO.EQ.0 ) {
                     CALL ZLACPY( 'Full', N, N, AFAC, LDA, AINV, LDA )
                     SRNAMT = 'ZGETRI'
                     NRHS = NSVAL( 1 )
                     LWORK = NMAX*MAX( 3, NRHS )
                     CALL ZGETRI( N, AINV, LDA, IWORK, WORK, LWORK, INFO )

                     // Check error code from ZGETRI.

                     IF( INFO.NE.0 ) CALL ALAERH( PATH, 'ZGETRI', INFO, 0, ' ', N, N, -1, -1, NB, IMAT, NFAIL, NERRS, NOUT )

                     // Compute the residual for the matrix times its
                     // inverse.  Also compute the 1-norm condition number
                     // of A.

                     CALL ZGET03( N, A, LDA, AINV, LDA, WORK, LDA, RWORK, RCONDO, RESULT( 2 ) )
                     ANORMO = ZLANGE( 'O', M, N, A, LDA, RWORK )

                     // Compute the infinity-norm condition number of A.

                     ANORMI = ZLANGE( 'I', M, N, A, LDA, RWORK )
                     AINVNM = ZLANGE( 'I', N, N, AINV, LDA, RWORK )
                     if ( ANORMI.LE.ZERO .OR. AINVNM.LE.ZERO ) {
                        RCONDI = ONE
                     } else {
                        RCONDI = ( ONE / ANORMI ) / AINVNM
                     }
                     NT = 2
                  } else {

                     // Do only the condition estimate if INFO > 0.

                     TRFCON = .TRUE.
                     ANORMO = ZLANGE( 'O', M, N, A, LDA, RWORK )
                     ANORMI = ZLANGE( 'I', M, N, A, LDA, RWORK )
                     RCONDO = ZERO
                     RCONDI = ZERO
                  }

                  // Print information about the tests so far that did not
                  // pass the threshold.

                  DO 30 K = 1, NT
                     if ( RESULT( K ).GE.THRESH ) {
                        IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALAHD( NOUT, PATH )                         WRITE( NOUT, FMT = 9999 )M, N, NB, IMAT, K, RESULT( K )
                        NFAIL = NFAIL + 1
                     }
   30             CONTINUE
                  NRUN = NRUN + NT

                  // Skip the remaining tests if this is not the first
                  // block size or if M .ne. N.  Skip the solve tests if
                 t // he matrix is singular.

                  IF( INB.GT.1 .OR. M.NE.N ) GO TO 90                   IF( TRFCON ) GO TO 70

                  DO 60 IRHS = 1, NNS
                     NRHS = NSVAL( IRHS )
                     XTYPE = 'N'

                     DO 50 ITRAN = 1, NTRAN
                        TRANS = TRANSS( ITRAN )
                        if ( ITRAN.EQ.1 ) {
                           RCONDC = RCONDO
                        } else {
                           RCONDC = RCONDI
                        }

*+    TEST 3
                        // Solve and compute residual for A * X = B.

                        SRNAMT = 'ZLARHS'
                        CALL ZLARHS( PATH, XTYPE, ' ', TRANS, N, N, KL, KU, NRHS, A, LDA, XACT, LDA, B, LDA, ISEED, INFO )
                        XTYPE = 'C'

                        CALL ZLACPY( 'Full', N, NRHS, B, LDA, X, LDA )
                        SRNAMT = 'ZGETRS'
                        CALL ZGETRS( TRANS, N, NRHS, AFAC, LDA, IWORK, X, LDA, INFO )

                        // Check error code from ZGETRS.

                        IF( INFO.NE.0 ) CALL ALAERH( PATH, 'ZGETRS', INFO, 0, TRANS, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT )

                        CALL ZLACPY( 'Full', N, NRHS, B, LDA, WORK, LDA )                         CALL ZGET02( TRANS, N, N, NRHS, A, LDA, X, LDA, WORK, LDA, RWORK, RESULT( 3 ) )

*+    TEST 4
                        // Check solution from generated exact solution.

                        CALL ZGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 4 ) )

*+    TESTS 5, 6, and 7
                        // Use iterative refinement to improve the
                        // solution.

                        SRNAMT = 'ZGERFS'
                        CALL ZGERFS( TRANS, N, NRHS, A, LDA, AFAC, LDA, IWORK, B, LDA, X, LDA, RWORK, RWORK( NRHS+1 ), WORK, RWORK( 2*NRHS+1 ), INFO )

                        // Check error code from ZGERFS.

                        IF( INFO.NE.0 ) CALL ALAERH( PATH, 'ZGERFS', INFO, 0, TRANS, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT )

                        CALL ZGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 5 ) )                         CALL ZGET07( TRANS, N, NRHS, A, LDA, B, LDA, X, LDA, XACT, LDA, RWORK, .TRUE., RWORK( NRHS+1 ), RESULT( 6 ) )

                        // Print information about the tests that did not
                        // pass the threshold.

                        DO 40 K = 3, 7
                           if ( RESULT( K ).GE.THRESH ) {
                              IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALAHD( NOUT, PATH )                               WRITE( NOUT, FMT = 9998 )TRANS, N, NRHS, IMAT, K, RESULT( K )
                              NFAIL = NFAIL + 1
                           }
   40                   CONTINUE
                        NRUN = NRUN + 5
   50                CONTINUE
   60             CONTINUE

*+    TEST 8
                     // Get an estimate of RCOND = 1/CNDNUM.

   70             CONTINUE
                  DO 80 ITRAN = 1, 2
                     if ( ITRAN.EQ.1 ) {
                        ANORM = ANORMO
                        RCONDC = RCONDO
                        NORM = 'O'
                     } else {
                        ANORM = ANORMI
                        RCONDC = RCONDI
                        NORM = 'I'
                     }
                     SRNAMT = 'ZGECON'
                     CALL ZGECON( NORM, N, AFAC, LDA, ANORM, RCOND, WORK, RWORK, INFO )

                        // Check error code from ZGECON.

                     IF( INFO.NE.0 ) CALL ALAERH( PATH, 'ZGECON', INFO, 0, NORM, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT )

                        // This line is needed on a Sun SPARCstation.

                     DUMMY = RCOND

                     RESULT( 8 ) = DGET06( RCOND, RCONDC )

                     // Print information about the tests that did not pass
                    t // he threshold.

                     if ( RESULT( 8 ).GE.THRESH ) {
                        IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALAHD( NOUT, PATH )                         WRITE( NOUT, FMT = 9997 )NORM, N, IMAT, 8, RESULT( 8 )
                        NFAIL = NFAIL + 1
                     }
                     NRUN = NRUN + 1
   80             CONTINUE
   90          CONTINUE
  100       CONTINUE

  110    CONTINUE
  120 CONTINUE

      // Print a summary of the results.

      CALL ALASUM( PATH, NOUT, NFAIL, NRUN, NERRS )

 9999 FORMAT( ' M = ', I5, ', N =', I5, ', NB =', I4, ', type ', I2,
     $      ', test(', I2, ') =', G12.5 )
 9998 FORMAT( ' TRANS=''', A1, ''', N =', I5, ', NRHS=', I3, ', type ',
     $      I2, ', test(', I2, ') =', G12.5 )
 9997 FORMAT( ' NORM =''', A1, ''', N =', I5, ',', 10X, ' type ', I2,
     $      ', test(', I2, ') =', G12.5 )
      RETURN

      // End of ZCHKGE

      }
