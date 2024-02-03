      SUBROUTINE ZDRVPO( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, NMAX, A, AFAC, ASAV, B, BSAV, X, XACT, S, WORK, RWORK, NOUT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               TSTERR;
      int                NMAX, NN, NOUT, NRHS;
      double             THRESH;
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      int                NVAL( * );
      double             RWORK( * ), S( * );
      COMPLEX*16         A( * ), AFAC( * ), ASAV( * ), B( * ), BSAV( * ), WORK( * ), X( * ), XACT( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0D+0, ZERO = 0.0D+0 ;
      int                NTYPES;
      const              NTYPES = 9 ;
      int                NTESTS;
      const              NTESTS = 6 ;
      // ..
      // .. Local Scalars ..
      bool               EQUIL, NOFACT, PREFAC, ZEROT;
      String             DIST, EQUED, FACT, TYPE, UPLO, XTYPE;
      String             PATH;
      int                I, IEQUED, IFACT, IMAT, IN, INFO, IOFF, IUPLO, IZERO, K, K1, KL, KU, LDA, MODE, N, NB, NBMIN, NERRS, NFACT, NFAIL, NIMAT, NRUN, NT;
      double             AINVNM, AMAX, ANORM, CNDNUM, RCOND, RCONDC, ROLDC, SCOND;
      // ..
      // .. Local Arrays ..
      String             EQUEDS( 2 ), FACTS( 3 ), UPLOS( 2 );
      int                ISEED( 4 ), ISEEDY( 4 );
      double             RESULT( NTESTS );
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DGET06, ZLANHE;
      // EXTERNAL LSAME, DGET06, ZLANHE
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALADHD, ALAERH, ALASVM, XLAENV, ZERRVX, ZGET04, ZLACPY, ZLAIPD, ZLAQHE, ZLARHS, ZLASET, ZLATB4, ZLATMS, ZPOEQU, ZPOSV, ZPOSVX, ZPOT01, ZPOT02, ZPOT05, ZPOTRF, ZPOTRI
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
      // .. Intrinsic Functions ..
      // INTRINSIC DCMPLX, MAX
      // ..
      // .. Data statements ..
      DATA               ISEEDY / 1988, 1989, 1990, 1991 /
      DATA               UPLOS / 'U', 'L' /
      DATA               FACTS / 'F', 'N', 'E' /
      DATA               EQUEDS / 'N', 'Y' /
      // ..
      // .. Executable Statements ..

      // Initialize constants and the random number seed.

      PATH( 1: 1 ) = 'Zomplex precision'
      PATH( 2: 3 ) = 'PO'
      NRUN = 0
      NFAIL = 0
      NERRS = 0
      DO 10 I = 1, 4
         ISEED( I ) = ISEEDY( I )
   10 CONTINUE

      // Test the error exits

      IF( TSTERR ) CALL ZERRVX( PATH, NOUT )
      INFOT = 0

      // Set the block size and minimum block size for testing.

      NB = 1
      NBMIN = 2
      CALL XLAENV( 1, NB )
      CALL XLAENV( 2, NBMIN )

      // Do for each value of N in NVAL

      DO 130 IN = 1, NN
         N = NVAL( IN )
         LDA = MAX( N, 1 )
         XTYPE = 'N'
         NIMAT = NTYPES
         IF( N.LE.0 ) NIMAT = 1

         DO 120 IMAT = 1, NIMAT

            // Do the tests only if DOTYPE( IMAT ) is true.

            IF( .NOT.DOTYPE( IMAT ) ) GO TO 120

            // Skip types 3, 4, or 5 if the matrix size is too small.

            ZEROT = IMAT.GE.3 .AND. IMAT.LE.5
            IF( ZEROT .AND. N.LT.IMAT-2 ) GO TO 120

            // Do first for UPLO = 'U', then for UPLO = 'L'

            DO 110 IUPLO = 1, 2
               UPLO = UPLOS( IUPLO )

               // Set up parameters with ZLATB4 and generate a test matrix
               // with ZLATMS.

               CALL ZLATB4( PATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST )

               SRNAMT = 'ZLATMS'
               CALL ZLATMS( N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU, UPLO, A, LDA, WORK, INFO )

               // Check error code from ZLATMS.

               if ( INFO.NE.0 ) {
                  CALL ALAERH( PATH, 'ZLATMS', INFO, 0, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT )
                  GO TO 110
               }

               // For types 3-5, zero one row and column of the matrix to
               // test that INFO is returned correctly.

               if ( ZEROT ) {
                  if ( IMAT.EQ.3 ) {
                     IZERO = 1
                  } else if ( IMAT.EQ.4 ) {
                     IZERO = N
                  } else {
                     IZERO = N / 2 + 1
                  }
                  IOFF = ( IZERO-1 )*LDA

                  // Set row and column IZERO of A to 0.

                  if ( IUPLO.EQ.1 ) {
                     DO 20 I = 1, IZERO - 1
                        A( IOFF+I ) = ZERO
   20                CONTINUE
                     IOFF = IOFF + IZERO
                     DO 30 I = IZERO, N
                        A( IOFF ) = ZERO
                        IOFF = IOFF + LDA
   30                CONTINUE
                  } else {
                     IOFF = IZERO
                     DO 40 I = 1, IZERO - 1
                        A( IOFF ) = ZERO
                        IOFF = IOFF + LDA
   40                CONTINUE
                     IOFF = IOFF - IZERO
                     DO 50 I = IZERO, N
                        A( IOFF+I ) = ZERO
   50                CONTINUE
                  }
               } else {
                  IZERO = 0
               }

               // Set the imaginary part of the diagonals.

               CALL ZLAIPD( N, A, LDA+1, 0 )

               // Save a copy of the matrix A in ASAV.

               CALL ZLACPY( UPLO, N, N, A, LDA, ASAV, LDA )

               DO 100 IEQUED = 1, 2
                  EQUED = EQUEDS( IEQUED )
                  if ( IEQUED.EQ.1 ) {
                     NFACT = 3
                  } else {
                     NFACT = 1
                  }

                  DO 90 IFACT = 1, NFACT
                     FACT = FACTS( IFACT )
                     PREFAC = LSAME( FACT, 'F' )
                     NOFACT = LSAME( FACT, 'N' )
                     EQUIL = LSAME( FACT, 'E' )

                     if ( ZEROT ) {
                        IF( PREFAC ) GO TO 90
                        RCONDC = ZERO

                     } else if ( .NOT.LSAME( FACT, 'N' ) ) {

                        // Compute the condition number for comparison with
                        // the value returned by ZPOSVX (FACT = 'N' reuses
                        // the condition number from the previous iteration
                        // with FACT = 'F').

                        CALL ZLACPY( UPLO, N, N, ASAV, LDA, AFAC, LDA )
                        if ( EQUIL .OR. IEQUED.GT.1 ) {

                           // Compute row and column scale factors to
                           // equilibrate the matrix A.

                           CALL ZPOEQU( N, AFAC, LDA, S, SCOND, AMAX, INFO )
                           if ( INFO.EQ.0 .AND. N.GT.0 ) {
                              IF( IEQUED.GT.1 ) SCOND = ZERO

                              // Equilibrate the matrix.

                              CALL ZLAQHE( UPLO, N, AFAC, LDA, S, SCOND, AMAX, EQUED )
                           }
                        }

                        // Save the condition number of the
                        // non-equilibrated system for use in ZGET04.

                        IF( EQUIL ) ROLDC = RCONDC

                        // Compute the 1-norm of A.

                        ANORM = ZLANHE( '1', UPLO, N, AFAC, LDA, RWORK )

                        // Factor the matrix A.

                        CALL ZPOTRF( UPLO, N, AFAC, LDA, INFO )

                        // Form the inverse of A.

                        CALL ZLACPY( UPLO, N, N, AFAC, LDA, A, LDA )
                        CALL ZPOTRI( UPLO, N, A, LDA, INFO )

                        // Compute the 1-norm condition number of A.

                        AINVNM = ZLANHE( '1', UPLO, N, A, LDA, RWORK )
                        if ( ANORM.LE.ZERO .OR. AINVNM.LE.ZERO ) {
                           RCONDC = ONE
                        } else {
                           RCONDC = ( ONE / ANORM ) / AINVNM
                        }
                     }

                     // Restore the matrix A.

                     CALL ZLACPY( UPLO, N, N, ASAV, LDA, A, LDA )

                     // Form an exact solution and set the right hand side.

                     SRNAMT = 'ZLARHS'
                     CALL ZLARHS( PATH, XTYPE, UPLO, ' ', N, N, KL, KU, NRHS, A, LDA, XACT, LDA, B, LDA, ISEED, INFO )
                     XTYPE = 'C'
                     CALL ZLACPY( 'Full', N, NRHS, B, LDA, BSAV, LDA )

                     if ( NOFACT ) {

                        // --- Test ZPOSV  ---

                        // Compute the L*L' or U'*U factorization of the
                        // matrix and solve the system.

                        CALL ZLACPY( UPLO, N, N, A, LDA, AFAC, LDA )
                        CALL ZLACPY( 'Full', N, NRHS, B, LDA, X, LDA )

                        SRNAMT = 'ZPOSV '
                        CALL ZPOSV( UPLO, N, NRHS, AFAC, LDA, X, LDA, INFO )

                        // Check error code from ZPOSV .

                        if ( INFO.NE.IZERO ) {
                           CALL ALAERH( PATH, 'ZPOSV ', INFO, IZERO, UPLO, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT )
                           GO TO 70
                        } else if ( INFO.NE.0 ) {
                           GO TO 70
                        }

                        // Reconstruct matrix from factors and compute
                        // residual.

                        CALL ZPOT01( UPLO, N, A, LDA, AFAC, LDA, RWORK, RESULT( 1 ) )

                        // Compute residual of the computed solution.

                        CALL ZLACPY( 'Full', N, NRHS, B, LDA, WORK, LDA )                         CALL ZPOT02( UPLO, N, NRHS, A, LDA, X, LDA, WORK, LDA, RWORK, RESULT( 2 ) )

                        // Check solution from generated exact solution.

                        CALL ZGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 3 ) )
                        NT = 3

                        // Print information about the tests that did not
                        // pass the threshold.

                        DO 60 K = 1, NT
                           if ( RESULT( K ).GE.THRESH ) {
                              IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALADHD( NOUT, PATH )                               WRITE( NOUT, FMT = 9999 )'ZPOSV ', UPLO, N, IMAT, K, RESULT( K )
                              NFAIL = NFAIL + 1
                           }
   60                   CONTINUE
                        NRUN = NRUN + NT
   70                   CONTINUE
                     }

                     // --- Test ZPOSVX ---

                     IF( .NOT.PREFAC ) CALL ZLASET( UPLO, N, N, DCMPLX( ZERO ), DCMPLX( ZERO ), AFAC, LDA )
                     CALL ZLASET( 'Full', N, NRHS, DCMPLX( ZERO ), DCMPLX( ZERO ), X, LDA )
                     if ( IEQUED.GT.1 .AND. N.GT.0 ) {

                        // Equilibrate the matrix if FACT='F' and
                        // EQUED='Y'.

                        CALL ZLAQHE( UPLO, N, A, LDA, S, SCOND, AMAX, EQUED )
                     }

                     // Solve the system and compute the condition number
                     // and error bounds using ZPOSVX.

                     SRNAMT = 'ZPOSVX'
                     CALL ZPOSVX( FACT, UPLO, N, NRHS, A, LDA, AFAC, LDA, EQUED, S, B, LDA, X, LDA, RCOND, RWORK, RWORK( NRHS+1 ), WORK, RWORK( 2*NRHS+1 ), INFO )

                     // Check the error code from ZPOSVX.

                     if ( INFO.NE.IZERO ) {
                        CALL ALAERH( PATH, 'ZPOSVX', INFO, IZERO, FACT // UPLO, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT )
                        GO TO 90
                     }

                     if ( INFO.EQ.0 ) {
                        if ( .NOT.PREFAC ) {

                           // Reconstruct matrix from factors and compute
                           // residual.

                           CALL ZPOT01( UPLO, N, A, LDA, AFAC, LDA, RWORK( 2*NRHS+1 ), RESULT( 1 ) )
                           K1 = 1
                        } else {
                           K1 = 2
                        }

                        // Compute residual of the computed solution.

                        CALL ZLACPY( 'Full', N, NRHS, BSAV, LDA, WORK, LDA )                         CALL ZPOT02( UPLO, N, NRHS, ASAV, LDA, X, LDA, WORK, LDA, RWORK( 2*NRHS+1 ), RESULT( 2 ) )

                        // Check solution from generated exact solution.

                        IF( NOFACT .OR. ( PREFAC .AND. LSAME( EQUED, 'N' ) ) ) THEN                            CALL ZGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 3 ) )
                        } else {
                           CALL ZGET04( N, NRHS, X, LDA, XACT, LDA, ROLDC, RESULT( 3 ) )
                        }

                        // Check the error bounds from iterative
                        // refinement.

                        CALL ZPOT05( UPLO, N, NRHS, ASAV, LDA, B, LDA, X, LDA, XACT, LDA, RWORK, RWORK( NRHS+1 ), RESULT( 4 ) )
                     } else {
                        K1 = 6
                     }

                     // Compare RCOND from ZPOSVX with the computed value
                     // in RCONDC.

                     RESULT( 6 ) = DGET06( RCOND, RCONDC )

                     // Print information about the tests that did not pass
                     // the threshold.

                     DO 80 K = K1, 6
                        if ( RESULT( K ).GE.THRESH ) {
                           IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALADHD( NOUT, PATH )
                           if ( PREFAC ) {
                              WRITE( NOUT, FMT = 9997 )'ZPOSVX', FACT, UPLO, N, EQUED, IMAT, K, RESULT( K )
                           } else {
                              WRITE( NOUT, FMT = 9998 )'ZPOSVX', FACT, UPLO, N, IMAT, K, RESULT( K )
                           }
                           NFAIL = NFAIL + 1
                        }
   80                CONTINUE
                     NRUN = NRUN + 7 - K1
   90             CONTINUE
  100          CONTINUE
  110       CONTINUE
  120    CONTINUE
  130 CONTINUE

      // Print a summary of the results.

      CALL ALASVM( PATH, NOUT, NFAIL, NRUN, NERRS )

 9999 FORMAT( 1X, A, ', UPLO=''', A1, ''', N =', I5, ', type ', I1, ', test(', I1, ')=', G12.5 )
 9998 FORMAT( 1X, A, ', FACT=''', A1, ''', UPLO=''', A1, ''', N=', I5, ', type ', I1, ', test(', I1, ')=', G12.5 )
 9997 FORMAT( 1X, A, ', FACT=''', A1, ''', UPLO=''', A1, ''', N=', I5, ', EQUED=''', A1, ''', type ', I1, ', test(', I1, ') =', G12.5 )
      RETURN

      // End of ZDRVPO

      }
