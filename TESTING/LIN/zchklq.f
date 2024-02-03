      SUBROUTINE ZCHKLQ( DOTYPE, NM, MVAL, NN, NVAL, NNB, NBVAL, NXVAL, NRHS, THRESH, TSTERR, NMAX, A, AF, AQ, AL, AC, B, X, XACT, TAU, WORK, RWORK, NOUT )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      bool               TSTERR;
      int                NM, NMAX, NN, NNB, NOUT, NRHS;
      double             THRESH;
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      int                MVAL( * ), NBVAL( * ), NVAL( * ), NXVAL( * );
      double             RWORK( * );
      COMPLEX*16         A( * ), AC( * ), AF( * ), AL( * ), AQ( * ), B( * ), TAU( * ), WORK( * ), X( * ), XACT( * )
      // ..
*
*  =====================================================================
*
      // .. Parameters ..
      int                NTESTS;
      PARAMETER          ( NTESTS = 7 )
      int                NTYPES;
      PARAMETER          ( NTYPES = 8 )
      double             ZERO;
      PARAMETER          ( ZERO = 0.0D0 )
      // ..
      // .. Local Scalars ..
      String             DIST, TYPE;
      String             PATH;
      int                I, IK, IM, IMAT, IN, INB, INFO, K, KL, KU, LDA, LWORK, M, MINMN, MODE, N, NB, NERRS, NFAIL, NK, NRUN, NT, NX;
      double             ANORM, CNDNUM;
      // ..
      // .. Local Arrays ..
      int                ISEED( 4 ), ISEEDY( 4 ), KVAL( 4 );
      double             RESULT( NTESTS );
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAERH, ALAHD, ALASUM, XLAENV, ZERRLQ, ZGELS, ZGET02, ZLACPY, ZLARHS, ZLATB4, ZLATMS, ZLQT01, ZLQT02, ZLQT03
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
*
      // Initialize constants and the random number seed.
*
      PATH( 1: 1 ) = 'Zomplex precision'
      PATH( 2: 3 ) = 'LQ'
      NRUN = 0
      NFAIL = 0
      NERRS = 0
      DO 10 I = 1, 4
         ISEED( I ) = ISEEDY( I )
   10 CONTINUE
*
      // Test the error exits
*
      IF( TSTERR ) CALL ZERRLQ( PATH, NOUT )
      INFOT = 0
      CALL XLAENV( 2, 2 )
*
      LDA = NMAX
      LWORK = NMAX*MAX( NMAX, NRHS )
*
      // Do for each value of M in MVAL.
*
      DO 70 IM = 1, NM
         M = MVAL( IM )
*
         // Do for each value of N in NVAL.
*
         DO 60 IN = 1, NN
            N = NVAL( IN )
            MINMN = MIN( M, N )
            DO 50 IMAT = 1, NTYPES
*
               // Do the tests only if DOTYPE( IMAT ) is true.
*
               IF( .NOT.DOTYPE( IMAT ) ) GO TO 50
*
               // Set up parameters with ZLATB4 and generate a test matrix
               // with ZLATMS.
*
               CALL ZLATB4( PATH, IMAT, M, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST )
*
               SRNAMT = 'ZLATMS'
               CALL ZLATMS( M, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU, 'No packing', A, LDA, WORK, INFO )
*
               // Check error code from ZLATMS.
*
               IF( INFO.NE.0 ) THEN
                  CALL ALAERH( PATH, 'ZLATMS', INFO, 0, ' ', M, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT )
                  GO TO 50
               END IF
*
               // Set some values for K: the first value must be MINMN,
               // corresponding to the call of ZLQT01; other values are
               // used in the calls of ZLQT02, and must not exceed MINMN.
*
               KVAL( 1 ) = MINMN
               KVAL( 2 ) = 0
               KVAL( 3 ) = 1
               KVAL( 4 ) = MINMN / 2
               IF( MINMN.EQ.0 ) THEN
                  NK = 1
               ELSE IF( MINMN.EQ.1 ) THEN
                  NK = 2
               ELSE IF( MINMN.LE.3 ) THEN
                  NK = 3
               ELSE
                  NK = 4
               END IF
*
               // Do for each value of K in KVAL
*
               DO 40 IK = 1, NK
                  K = KVAL( IK )
*
                  // Do for each pair of values (NB,NX) in NBVAL and NXVAL.
*
                  DO 30 INB = 1, NNB
                     NB = NBVAL( INB )
                     CALL XLAENV( 1, NB )
                     NX = NXVAL( INB )
                     CALL XLAENV( 3, NX )
                     DO I = 1, NTESTS
                        RESULT( I ) = ZERO
                     END DO
                     NT = 2
                     IF( IK.EQ.1 ) THEN
*
                        // Test ZGELQF
*
                        CALL ZLQT01( M, N, A, AF, AQ, AL, LDA, TAU, WORK, LWORK, RWORK, RESULT( 1 ) )
                     ELSE IF( M.LE.N ) THEN
*
                        // Test ZUNGLQ, using factorization
                        // returned by ZLQT01
*
                        CALL ZLQT02( M, N, K, A, AF, AQ, AL, LDA, TAU, WORK, LWORK, RWORK, RESULT( 1 ) )
                     END IF
                     IF( M.GE.K ) THEN
*
                        // Test ZUNMLQ, using factorization returned
                        // by ZLQT01
*
                        CALL ZLQT03( M, N, K, AF, AC, AL, AQ, LDA, TAU, WORK, LWORK, RWORK, RESULT( 3 ) )
                        NT = NT + 4
*
                        // If M<=N and K=M, call ZGELS to solve a system
                        // with NRHS right hand sides and compute the
                        // residual.
*
                        IF( K.EQ.M .AND. INB.EQ.1 ) THEN
*
                           // Generate a solution and set the right
                           // hand side.
*
                           SRNAMT = 'ZLARHS'
                           CALL ZLARHS( PATH, 'New', 'Full', 'No transpose', M, N, 0, 0, NRHS, A, LDA, XACT, LDA, B, LDA, ISEED, INFO )
*
                           CALL ZLACPY( 'Full', M, NRHS, B, LDA, X, LDA )
*
                           // Reset AF to the original matrix. ZGELS
                           // factors the matrix before solving the system.
*
                           CALL ZLACPY( 'Full', M, N, A, LDA, AF, LDA )
*
                           SRNAMT = 'ZGELS'
                           CALL ZGELS( 'No transpose', M, N, NRHS, AF, LDA, X, LDA, WORK, LWORK, INFO )
*
                           // Check error code from ZGELS.
*
                           IF( INFO.NE.0 ) CALL ALAERH( PATH, 'ZGELS', INFO, 0, 'N', M, N, NRHS, -1, NB, IMAT, NFAIL, NERRS, NOUT )
*
                           CALL ZGET02( 'No transpose', M, N, NRHS, A, LDA, X, LDA, B, LDA, RWORK, RESULT( 7 ) )
                           NT = NT + 1
                        END IF
                     END IF
*
                     // Print information about the tests that did not
                     // pass the threshold.
*
                     DO 20 I = 1, NT
                        IF( RESULT( I ).GE.THRESH ) THEN
                           IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALAHD( NOUT, PATH )                            WRITE( NOUT, FMT = 9999 )M, N, K, NB, NX, IMAT, I, RESULT( I )
                           NFAIL = NFAIL + 1
                        END IF
   20                CONTINUE
                     NRUN = NRUN + NT
   30             CONTINUE
   40          CONTINUE
   50       CONTINUE
   60    CONTINUE
   70 CONTINUE
*
      // Print a summary of the results.
*
      CALL ALASUM( PATH, NOUT, NFAIL, NRUN, NERRS )
*
 9999 FORMAT( ' M=', I5, ', N=', I5, ', K=', I5, ', NB=', I4, ', NX=',
     $      I5, ', type ', I2, ', test(', I2, ')=', G12.5 )
      RETURN
*
      // End of ZCHKLQ
*
      END
