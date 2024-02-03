      SUBROUTINE DDRVPP( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, NMAX, A, AFAC, ASAV, B, BSAV, X, XACT, S, WORK, RWORK, IWORK, NOUT )

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
      int                IWORK( * ), NVAL( * );
      double             A( * ), AFAC( * ), ASAV( * ), B( * ), BSAV( * ), RWORK( * ), S( * ), WORK( * ), X( * ), XACT( * );
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
      String             DIST, EQUED, FACT, PACKIT, TYPE, UPLO, XTYPE;
      String             PATH;
      int                I, IEQUED, IFACT, IMAT, IN, INFO, IOFF, IUPLO, IZERO, K, K1, KL, KU, LDA, MODE, N, NERRS, NFACT, NFAIL, NIMAT, NPP, NRUN, NT;
      double             AINVNM, AMAX, ANORM, CNDNUM, RCOND, RCONDC, ROLDC, SCOND;
      // ..
      // .. Local Arrays ..
      String             EQUEDS( 2 ), FACTS( 3 ), PACKS( 2 ), UPLOS( 2 );
      int                ISEED( 4 ), ISEEDY( 4 );
      double             RESULT( NTESTS );
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DGET06, DLANSP;
      // EXTERNAL LSAME, DGET06, DLANSP
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALADHD, ALAERH, ALASVM, DCOPY, DERRVX, DGET04, DLACPY, DLAQSP, DLARHS, DLASET, DLATB4, DLATMS, DPPEQU, DPPSV, DPPSVX, DPPT01, DPPT02, DPPT05, DPPTRF, DPPTRI
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
      // INTRINSIC MAX
      // ..
      // .. Data statements ..
      DATA               ISEEDY / 1988, 1989, 1990, 1991 /
      DATA               UPLOS / 'U', 'L' / , FACTS / 'F', 'N', 'E' / , PACKS / 'C', 'R' / , EQUEDS / 'N', 'Y' /
      // ..
      // .. Executable Statements ..

      // Initialize constants and the random number seed.

      PATH( 1: 1 ) = 'double          ';
      PATH( 2: 3 ) = 'PP'
      NRUN = 0
      NFAIL = 0
      NERRS = 0
      DO 10 I = 1, 4
         ISEED( I ) = ISEEDY( I )
   10 CONTINUE

      // Test the error exits

      IF( TSTERR ) CALL DERRVX( PATH, NOUT )
      INFOT = 0

      // Do for each value of N in NVAL

      DO 140 IN = 1, NN
         N = NVAL( IN )
         LDA = MAX( N, 1 )
         NPP = N*( N+1 ) / 2
         XTYPE = 'N'
         NIMAT = NTYPES
         IF( N.LE.0 ) NIMAT = 1

         DO 130 IMAT = 1, NIMAT

            // Do the tests only if DOTYPE( IMAT ) is true.

            IF( .NOT.DOTYPE( IMAT ) ) GO TO 130

            // Skip types 3, 4, or 5 if the matrix size is too small.

            ZEROT = IMAT.GE.3 .AND. IMAT.LE.5
            IF( ZEROT .AND. N.LT.IMAT-2 ) GO TO 130

            // Do first for UPLO = 'U', then for UPLO = 'L'

            DO 120 IUPLO = 1, 2
               UPLO = UPLOS( IUPLO )
               PACKIT = PACKS( IUPLO )

               // Set up parameters with DLATB4 and generate a test matrix
               // with DLATMS.

               CALL DLATB4( PATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST )
               RCONDC = ONE / CNDNUM

               SRNAMT = 'DLATMS'
               CALL DLATMS( N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU, PACKIT, A, LDA, WORK, INFO )

               // Check error code from DLATMS.

               if ( INFO.NE.0 ) {
                  CALL ALAERH( PATH, 'DLATMS', INFO, 0, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT )
                  GO TO 120
               }

               // For types 3-5, zero one row and column of the matrix to
              t // est that INFO is returned correctly.

               if ( ZEROT ) {
                  if ( IMAT.EQ.3 ) {
                     IZERO = 1
                  } else if ( IMAT.EQ.4 ) {
                     IZERO = N
                  } else {
                     IZERO = N / 2 + 1
                  }

                  // Set row and column IZERO of A to 0.

                  if ( IUPLO.EQ.1 ) {
                     IOFF = ( IZERO-1 )*IZERO / 2
                     DO 20 I = 1, IZERO - 1
                        A( IOFF+I ) = ZERO
   20                CONTINUE
                     IOFF = IOFF + IZERO
                     DO 30 I = IZERO, N
                        A( IOFF ) = ZERO
                        IOFF = IOFF + I
   30                CONTINUE
                  } else {
                     IOFF = IZERO
                     DO 40 I = 1, IZERO - 1
                        A( IOFF ) = ZERO
                        IOFF = IOFF + N - I
   40                CONTINUE
                     IOFF = IOFF - IZERO
                     DO 50 I = IZERO, N
                        A( IOFF+I ) = ZERO
   50                CONTINUE
                  }
               } else {
                  IZERO = 0
               }

               // Save a copy of the matrix A in ASAV.

               CALL DCOPY( NPP, A, 1, ASAV, 1 )

               DO 110 IEQUED = 1, 2
                  EQUED = EQUEDS( IEQUED )
                  if ( IEQUED.EQ.1 ) {
                     NFACT = 3
                  } else {
                     NFACT = 1
                  }

                  DO 100 IFACT = 1, NFACT
                     FACT = FACTS( IFACT )
                     PREFAC = LSAME( FACT, 'F' )
                     NOFACT = LSAME( FACT, 'N' )
                     EQUIL = LSAME( FACT, 'E' )

                     if ( ZEROT ) {
                        IF( PREFAC ) GO TO 100
                        RCONDC = ZERO

                     } else if ( .NOT.LSAME( FACT, 'N' ) ) {

                        // Compute the condition number for comparison with
                       t // he value returned by DPPSVX (FACT = 'N' reuses
                       t // he condition number from the previous iteration
                        // with FACT = 'F').

                        CALL DCOPY( NPP, ASAV, 1, AFAC, 1 )
                        if ( EQUIL .OR. IEQUED.GT.1 ) {

                           // Compute row and column scale factors to
                           // equilibrate the matrix A.

                           CALL DPPEQU( UPLO, N, AFAC, S, SCOND, AMAX, INFO )
                           if ( INFO.EQ.0 .AND. N.GT.0 ) {
                              IF( IEQUED.GT.1 ) SCOND = ZERO

                              // Equilibrate the matrix.

                              CALL DLAQSP( UPLO, N, AFAC, S, SCOND, AMAX, EQUED )
                           }
                        }

                        // Save the condition number of the
                        // non-equilibrated system for use in DGET04.

                        IF( EQUIL ) ROLDC = RCONDC

                        // Compute the 1-norm of A.

                        ANORM = DLANSP( '1', UPLO, N, AFAC, RWORK )

                        // Factor the matrix A.

                        CALL DPPTRF( UPLO, N, AFAC, INFO )

                        // Form the inverse of A.

                        CALL DCOPY( NPP, AFAC, 1, A, 1 )
                        CALL DPPTRI( UPLO, N, A, INFO )

                        // Compute the 1-norm condition number of A.

                        AINVNM = DLANSP( '1', UPLO, N, A, RWORK )
                        if ( ANORM.LE.ZERO .OR. AINVNM.LE.ZERO ) {
                           RCONDC = ONE
                        } else {
                           RCONDC = ( ONE / ANORM ) / AINVNM
                        }
                     }

                     // Restore the matrix A.

                     CALL DCOPY( NPP, ASAV, 1, A, 1 )

                     // Form an exact solution and set the right hand side.

                     SRNAMT = 'DLARHS'
                     CALL DLARHS( PATH, XTYPE, UPLO, ' ', N, N, KL, KU, NRHS, A, LDA, XACT, LDA, B, LDA, ISEED, INFO )
                     XTYPE = 'C'
                     CALL DLACPY( 'Full', N, NRHS, B, LDA, BSAV, LDA )

                     if ( NOFACT ) {

                        // --- Test DPPSV  ---

                        // Compute the L*L' or U'*U factorization of the
                        // matrix and solve the system.

                        CALL DCOPY( NPP, A, 1, AFAC, 1 )
                        CALL DLACPY( 'Full', N, NRHS, B, LDA, X, LDA )

                        SRNAMT = 'DPPSV '
                        CALL DPPSV( UPLO, N, NRHS, AFAC, X, LDA, INFO )

                        // Check error code from DPPSV .

                        if ( INFO.NE.IZERO ) {
                           CALL ALAERH( PATH, 'DPPSV ', INFO, IZERO, UPLO, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT )
                           GO TO 70
                        } else if ( INFO.NE.0 ) {
                           GO TO 70
                        }

                        // Reconstruct matrix from factors and compute
                        // residual.

                        CALL DPPT01( UPLO, N, A, AFAC, RWORK, RESULT( 1 ) )

                        // Compute residual of the computed solution.

                        CALL DLACPY( 'Full', N, NRHS, B, LDA, WORK, LDA )                         CALL DPPT02( UPLO, N, NRHS, A, X, LDA, WORK, LDA, RWORK, RESULT( 2 ) )

                        // Check solution from generated exact solution.

                        CALL DGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 3 ) )
                        NT = 3

                        // Print information about the tests that did not
                        // pass the threshold.

                        DO 60 K = 1, NT
                           if ( RESULT( K ).GE.THRESH ) {
                              IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALADHD( NOUT, PATH )                               WRITE( NOUT, FMT = 9999 )'DPPSV ', UPLO, N, IMAT, K, RESULT( K )
                              NFAIL = NFAIL + 1
                           }
   60                   CONTINUE
                        NRUN = NRUN + NT
   70                   CONTINUE
                     }

                     // --- Test DPPSVX ---

                     IF( .NOT.PREFAC .AND. NPP.GT.0 ) CALL DLASET( 'Full', NPP, 1, ZERO, ZERO, AFAC, NPP )
                     CALL DLASET( 'Full', N, NRHS, ZERO, ZERO, X, LDA )
                     if ( IEQUED.GT.1 .AND. N.GT.0 ) {

                        // Equilibrate the matrix if FACT='F' and
                        // EQUED='Y'.

                        CALL DLAQSP( UPLO, N, A, S, SCOND, AMAX, EQUED )
                     }

                     // Solve the system and compute the condition number
                     // and error bounds using DPPSVX.

                     SRNAMT = 'DPPSVX'
                     CALL DPPSVX( FACT, UPLO, N, NRHS, A, AFAC, EQUED, S, B, LDA, X, LDA, RCOND, RWORK, RWORK( NRHS+1 ), WORK, IWORK, INFO )

                     // Check the error code from DPPSVX.

                     if ( INFO.NE.IZERO ) {
                        CALL ALAERH( PATH, 'DPPSVX', INFO, IZERO, FACT // UPLO, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT )
                        GO TO 90
                     }

                     if ( INFO.EQ.0 ) {
                        if ( .NOT.PREFAC ) {

                           // Reconstruct matrix from factors and compute
                           // residual.

                           CALL DPPT01( UPLO, N, A, AFAC, RWORK( 2*NRHS+1 ), RESULT( 1 ) )
                           K1 = 1
                        } else {
                           K1 = 2
                        }

                        // Compute residual of the computed solution.

                        CALL DLACPY( 'Full', N, NRHS, BSAV, LDA, WORK, LDA )                         CALL DPPT02( UPLO, N, NRHS, ASAV, X, LDA, WORK, LDA, RWORK( 2*NRHS+1 ), RESULT( 2 ) )

                        // Check solution from generated exact solution.

                        IF( NOFACT .OR. ( PREFAC .AND. LSAME( EQUED, 'N' ) ) ) THEN                            CALL DGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 3 ) )
                        } else {
                           CALL DGET04( N, NRHS, X, LDA, XACT, LDA, ROLDC, RESULT( 3 ) )
                        }

                        // Check the error bounds from iterative
                        // refinement.

                        CALL DPPT05( UPLO, N, NRHS, ASAV, B, LDA, X, LDA, XACT, LDA, RWORK, RWORK( NRHS+1 ), RESULT( 4 ) )
                     } else {
                        K1 = 6
                     }

                     // Compare RCOND from DPPSVX with the computed value
                     // in RCONDC.

                     RESULT( 6 ) = DGET06( RCOND, RCONDC )

                     // Print information about the tests that did not pass
                    t // he threshold.

                     DO 80 K = K1, 6
                        if ( RESULT( K ).GE.THRESH ) {
                           IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALADHD( NOUT, PATH )
                           if ( PREFAC ) {
                              WRITE( NOUT, FMT = 9997 )'DPPSVX', FACT, UPLO, N, EQUED, IMAT, K, RESULT( K )
                           } else {
                              WRITE( NOUT, FMT = 9998 )'DPPSVX', FACT, UPLO, N, IMAT, K, RESULT( K )
                           }
                           NFAIL = NFAIL + 1
                        }
   80                CONTINUE
                     NRUN = NRUN + 7 - K1
   90                CONTINUE
  100             CONTINUE
  110          CONTINUE
  120       CONTINUE
  130    CONTINUE
  140 CONTINUE

      // Print a summary of the results.

      CALL ALASVM( PATH, NOUT, NFAIL, NRUN, NERRS )

 9999 FORMAT( 1X, A, ', UPLO=''', A1, ''', N =', I5, ', type ', I1,
     $      ', test(', I1, ')=', G12.5 )
 9998 FORMAT( 1X, A, ', FACT=''', A1, ''', UPLO=''', A1, ''', N=', I5,
     $      ', type ', I1, ', test(', I1, ')=', G12.5 )
 9997 FORMAT( 1X, A, ', FACT=''', A1, ''', UPLO=''', A1, ''', N=', I5,
     $      ', EQUED=''', A1, ''', type ', I1, ', test(', I1, ')=',
     $      G12.5 )
      RETURN

      // End of DDRVPP

      }
