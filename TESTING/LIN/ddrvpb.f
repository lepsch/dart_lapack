      SUBROUTINE DDRVPB( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, NMAX, A, AFAC, ASAV, B, BSAV, X, XACT, S, WORK, RWORK, IWORK, NOUT )

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
      int                NTYPES, NTESTS;
      const              NTYPES = 8, NTESTS = 6 ;
      int                NBW;
      const              NBW = 4 ;
      // ..
      // .. Local Scalars ..
      bool               EQUIL, NOFACT, PREFAC, ZEROT;
      String             DIST, EQUED, FACT, PACKIT, TYPE, UPLO, XTYPE;
      String             PATH;
      int                I, I1, I2, IEQUED, IFACT, IKD, IMAT, IN, INFO, IOFF, IUPLO, IW, IZERO, K, K1, KD, KL, KOFF, KU, LDA, LDAB, MODE, N, NB, NBMIN, NERRS, NFACT, NFAIL, NIMAT, NKD, NRUN, NT;
      double             AINVNM, AMAX, ANORM, CNDNUM, RCOND, RCONDC, ROLDC, SCOND;
      // ..
      // .. Local Arrays ..
      String             EQUEDS( 2 ), FACTS( 3 );
      int                ISEED( 4 ), ISEEDY( 4 ), KDVAL( NBW );
      double             RESULT( NTESTS );
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DGET06, DLANGE, DLANSB;
      // EXTERNAL LSAME, DGET06, DLANGE, DLANSB
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALADHD, ALAERH, ALASVM, DCOPY, DERRVX, DGET04, DLACPY, DLAQSB, DLARHS, DLASET, DLATB4, DLATMS, DPBEQU, DPBSV, DPBSVX, DPBT01, DPBT02, DPBT05, DPBTRF, DPBTRS, DSWAP, XLAENV
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
      DATA               FACTS / 'F', 'N', 'E' /
      DATA               EQUEDS / 'N', 'Y' /
      // ..
      // .. Executable Statements ..

      // Initialize constants and the random number seed.

      PATH( 1: 1 ) = 'double          ';
      PATH( 2: 3 ) = 'PB'
      NRUN = 0
      NFAIL = 0
      NERRS = 0
      DO 10 I = 1, 4
         ISEED( I ) = ISEEDY( I )
   10 CONTINUE

      // Test the error exits

      IF( TSTERR ) CALL DERRVX( PATH, NOUT )
      INFOT = 0
      KDVAL( 1 ) = 0

      // Set the block size and minimum block size for testing.

      NB = 1
      NBMIN = 2
      CALL XLAENV( 1, NB )
      CALL XLAENV( 2, NBMIN )

      // Do for each value of N in NVAL

      DO 110 IN = 1, NN
         N = NVAL( IN )
         LDA = MAX( N, 1 )
         XTYPE = 'N'

         // Set limits on the number of loop iterations.

         NKD = MAX( 1, MIN( N, 4 ) )
         NIMAT = NTYPES
         IF( N.EQ.0 ) NIMAT = 1

         KDVAL( 2 ) = N + ( N+1 ) / 4
         KDVAL( 3 ) = ( 3*N-1 ) / 4
         KDVAL( 4 ) = ( N+1 ) / 4

         DO 100 IKD = 1, NKD

            // Do for KD = 0, (5*N+1)/4, (3N-1)/4, and (N+1)/4. This order
            // makes it easier to skip redundant values for small values
            // of N.

            KD = KDVAL( IKD )
            LDAB = KD + 1

            // Do first for UPLO = 'U', then for UPLO = 'L'

            DO 90 IUPLO = 1, 2
               KOFF = 1
               if ( IUPLO.EQ.1 ) {
                  UPLO = 'U'
                  PACKIT = 'Q'
                  KOFF = MAX( 1, KD+2-N )
               } else {
                  UPLO = 'L'
                  PACKIT = 'B'
               }

               DO 80 IMAT = 1, NIMAT

                  // Do the tests only if DOTYPE( IMAT ) is true.

                  IF( .NOT.DOTYPE( IMAT ) ) GO TO 80

                  // Skip types 2, 3, or 4 if the matrix size is too small.

                  ZEROT = IMAT.GE.2 .AND. IMAT.LE.4
                  IF( ZEROT .AND. N.LT.IMAT-1 ) GO TO 80

                  if ( .NOT.ZEROT .OR. .NOT.DOTYPE( 1 ) ) {

                     // Set up parameters with DLATB4 and generate a test
                     // matrix with DLATMS.

                     CALL DLATB4( PATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST )

                     SRNAMT = 'DLATMS'
                     CALL DLATMS( N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KD, KD, PACKIT, A( KOFF ), LDAB, WORK, INFO )

                     // Check error code from DLATMS.

                     if ( INFO.NE.0 ) {
                        CALL ALAERH( PATH, 'DLATMS', INFO, 0, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT )
                        GO TO 80
                     }
                  } else if ( IZERO.GT.0 ) {

                     // Use the same matrix for types 3 and 4 as for type
                     // 2 by copying back the zeroed out column,

                     IW = 2*LDA + 1
                     if ( IUPLO.EQ.1 ) {
                        IOFF = ( IZERO-1 )*LDAB + KD + 1
                        CALL DCOPY( IZERO-I1, WORK( IW ), 1, A( IOFF-IZERO+I1 ), 1 )
                        IW = IW + IZERO - I1
                        CALL DCOPY( I2-IZERO+1, WORK( IW ), 1, A( IOFF ), MAX( LDAB-1, 1 ) )
                     } else {
                        IOFF = ( I1-1 )*LDAB + 1
                        CALL DCOPY( IZERO-I1, WORK( IW ), 1, A( IOFF+IZERO-I1 ), MAX( LDAB-1, 1 ) )
                        IOFF = ( IZERO-1 )*LDAB + 1
                        IW = IW + IZERO - I1
                        CALL DCOPY( I2-IZERO+1, WORK( IW ), 1, A( IOFF ), 1 )
                     }
                  }

                  // For types 2-4, zero one row and column of the matrix
                 t // o test that INFO is returned correctly.

                  IZERO = 0
                  if ( ZEROT ) {
                     if ( IMAT.EQ.2 ) {
                        IZERO = 1
                     } else if ( IMAT.EQ.3 ) {
                        IZERO = N
                     } else {
                        IZERO = N / 2 + 1
                     }

                     // Save the zeroed out row and column in WORK(*,3)

                     IW = 2*LDA
                     DO 20 I = 1, MIN( 2*KD+1, N )
                        WORK( IW+I ) = ZERO
   20                CONTINUE
                     IW = IW + 1
                     I1 = MAX( IZERO-KD, 1 )
                     I2 = MIN( IZERO+KD, N )

                     if ( IUPLO.EQ.1 ) {
                        IOFF = ( IZERO-1 )*LDAB + KD + 1
                        CALL DSWAP( IZERO-I1, A( IOFF-IZERO+I1 ), 1, WORK( IW ), 1 )
                        IW = IW + IZERO - I1
                        CALL DSWAP( I2-IZERO+1, A( IOFF ), MAX( LDAB-1, 1 ), WORK( IW ), 1 )
                     } else {
                        IOFF = ( I1-1 )*LDAB + 1
                        CALL DSWAP( IZERO-I1, A( IOFF+IZERO-I1 ), MAX( LDAB-1, 1 ), WORK( IW ), 1 )
                        IOFF = ( IZERO-1 )*LDAB + 1
                        IW = IW + IZERO - I1
                        CALL DSWAP( I2-IZERO+1, A( IOFF ), 1, WORK( IW ), 1 )
                     }
                  }

                  // Save a copy of the matrix A in ASAV.

                  CALL DLACPY( 'Full', KD+1, N, A, LDAB, ASAV, LDAB )

                  DO 70 IEQUED = 1, 2
                     EQUED = EQUEDS( IEQUED )
                     if ( IEQUED.EQ.1 ) {
                        NFACT = 3
                     } else {
                        NFACT = 1
                     }

                     DO 60 IFACT = 1, NFACT
                        FACT = FACTS( IFACT )
                        PREFAC = LSAME( FACT, 'F' )
                        NOFACT = LSAME( FACT, 'N' )
                        EQUIL = LSAME( FACT, 'E' )

                        if ( ZEROT ) {
                           IF( PREFAC ) GO TO 60
                           RCONDC = ZERO

                        } else if ( .NOT.LSAME( FACT, 'N' ) ) {

                           // Compute the condition number for comparison
                           // with the value returned by DPBSVX (FACT =
                           // 'N' reuses the condition number from the
                           // previous iteration with FACT = 'F').

                           CALL DLACPY( 'Full', KD+1, N, ASAV, LDAB, AFAC, LDAB )
                           if ( EQUIL .OR. IEQUED.GT.1 ) {

                              // Compute row and column scale factors to
                              // equilibrate the matrix A.

                              CALL DPBEQU( UPLO, N, KD, AFAC, LDAB, S, SCOND, AMAX, INFO )
                              if ( INFO.EQ.0 .AND. N.GT.0 ) {
                                 IF( IEQUED.GT.1 ) SCOND = ZERO

                                 // Equilibrate the matrix.

                                 CALL DLAQSB( UPLO, N, KD, AFAC, LDAB, S, SCOND, AMAX, EQUED )
                              }
                           }

                           // Save the condition number of the
                           // non-equilibrated system for use in DGET04.

                           IF( EQUIL ) ROLDC = RCONDC

                           // Compute the 1-norm of A.

                           ANORM = DLANSB( '1', UPLO, N, KD, AFAC, LDAB, RWORK )

                           // Factor the matrix A.

                           CALL DPBTRF( UPLO, N, KD, AFAC, LDAB, INFO )

                           // Form the inverse of A.

                           CALL DLASET( 'Full', N, N, ZERO, ONE, A, LDA )
                           SRNAMT = 'DPBTRS'
                           CALL DPBTRS( UPLO, N, KD, N, AFAC, LDAB, A, LDA, INFO )

                           // Compute the 1-norm condition number of A.

                           AINVNM = DLANGE( '1', N, N, A, LDA, RWORK )
                           if ( ANORM.LE.ZERO .OR. AINVNM.LE.ZERO ) {
                              RCONDC = ONE
                           } else {
                              RCONDC = ( ONE / ANORM ) / AINVNM
                           }
                        }

                        // Restore the matrix A.

                        CALL DLACPY( 'Full', KD+1, N, ASAV, LDAB, A, LDAB )

                        // Form an exact solution and set the right hand
                        // side.

                        SRNAMT = 'DLARHS'
                        CALL DLARHS( PATH, XTYPE, UPLO, ' ', N, N, KD, KD, NRHS, A, LDAB, XACT, LDA, B, LDA, ISEED, INFO )
                        XTYPE = 'C'
                        CALL DLACPY( 'Full', N, NRHS, B, LDA, BSAV, LDA )

                        if ( NOFACT ) {

                           // --- Test DPBSV  ---

                           // Compute the L*L' or U'*U factorization of the
                           // matrix and solve the system.

                           CALL DLACPY( 'Full', KD+1, N, A, LDAB, AFAC, LDAB )                            CALL DLACPY( 'Full', N, NRHS, B, LDA, X, LDA )

                           SRNAMT = 'DPBSV '
                           CALL DPBSV( UPLO, N, KD, NRHS, AFAC, LDAB, X, LDA, INFO )

                           // Check error code from DPBSV .

                           if ( INFO.NE.IZERO ) {
                              CALL ALAERH( PATH, 'DPBSV ', INFO, IZERO, UPLO, N, N, KD, KD, NRHS, IMAT, NFAIL, NERRS, NOUT )
                              GO TO 40
                           } else if ( INFO.NE.0 ) {
                              GO TO 40
                           }

                           // Reconstruct matrix from factors and compute
                           // residual.

                           CALL DPBT01( UPLO, N, KD, A, LDAB, AFAC, LDAB, RWORK, RESULT( 1 ) )

                           // Compute residual of the computed solution.

                           CALL DLACPY( 'Full', N, NRHS, B, LDA, WORK, LDA )                            CALL DPBT02( UPLO, N, KD, NRHS, A, LDAB, X, LDA, WORK, LDA, RWORK, RESULT( 2 ) )

                           // Check solution from generated exact solution.

                           CALL DGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 3 ) )
                           NT = 3

                           // Print information about the tests that did
                           // not pass the threshold.

                           DO 30 K = 1, NT
                              if ( RESULT( K ).GE.THRESH ) {
                                 IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALADHD( NOUT, PATH )                                  WRITE( NOUT, FMT = 9999 )'DPBSV ', UPLO, N, KD, IMAT, K, RESULT( K )
                                 NFAIL = NFAIL + 1
                              }
   30                      CONTINUE
                           NRUN = NRUN + NT
   40                      CONTINUE
                        }

                        // --- Test DPBSVX ---

                        IF( .NOT.PREFAC ) CALL DLASET( 'Full', KD+1, N, ZERO, ZERO, AFAC, LDAB )
                        CALL DLASET( 'Full', N, NRHS, ZERO, ZERO, X, LDA )
                        if ( IEQUED.GT.1 .AND. N.GT.0 ) {

                           // Equilibrate the matrix if FACT='F' and
                           // EQUED='Y'

                           CALL DLAQSB( UPLO, N, KD, A, LDAB, S, SCOND, AMAX, EQUED )
                        }

                        // Solve the system and compute the condition
                        // number and error bounds using DPBSVX.

                        SRNAMT = 'DPBSVX'
                        CALL DPBSVX( FACT, UPLO, N, KD, NRHS, A, LDAB, AFAC, LDAB, EQUED, S, B, LDA, X, LDA, RCOND, RWORK, RWORK( NRHS+1 ), WORK, IWORK, INFO )

                        // Check the error code from DPBSVX.

                        if ( INFO.NE.IZERO ) {
                           CALL ALAERH( PATH, 'DPBSVX', INFO, IZERO, FACT // UPLO, N, N, KD, KD, NRHS, IMAT, NFAIL, NERRS, NOUT )
                           GO TO 60
                        }

                        if ( INFO.EQ.0 ) {
                           if ( .NOT.PREFAC ) {

                              // Reconstruct matrix from factors and
                              // compute residual.

                              CALL DPBT01( UPLO, N, KD, A, LDAB, AFAC, LDAB, RWORK( 2*NRHS+1 ), RESULT( 1 ) )
                              K1 = 1
                           } else {
                              K1 = 2
                           }

                           // Compute residual of the computed solution.

                           CALL DLACPY( 'Full', N, NRHS, BSAV, LDA, WORK, LDA )                            CALL DPBT02( UPLO, N, KD, NRHS, ASAV, LDAB, X, LDA, WORK, LDA, RWORK( 2*NRHS+1 ), RESULT( 2 ) )

                           // Check solution from generated exact solution.

                           IF( NOFACT .OR. ( PREFAC .AND. LSAME( EQUED, 'N' ) ) ) THEN                               CALL DGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 3 ) )
                           } else {
                              CALL DGET04( N, NRHS, X, LDA, XACT, LDA, ROLDC, RESULT( 3 ) )
                           }

                           // Check the error bounds from iterative
                           // refinement.

                           CALL DPBT05( UPLO, N, KD, NRHS, ASAV, LDAB, B, LDA, X, LDA, XACT, LDA, RWORK, RWORK( NRHS+1 ), RESULT( 4 ) )
                        } else {
                           K1 = 6
                        }

                        // Compare RCOND from DPBSVX with the computed
                        // value in RCONDC.

                        RESULT( 6 ) = DGET06( RCOND, RCONDC )

                        // Print information about the tests that did not
                        // pass the threshold.

                        DO 50 K = K1, 6
                           if ( RESULT( K ).GE.THRESH ) {
                              IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALADHD( NOUT, PATH )
                              if ( PREFAC ) {
                                 WRITE( NOUT, FMT = 9997 )'DPBSVX', FACT, UPLO, N, KD, EQUED, IMAT, K, RESULT( K )
                              } else {
                                 WRITE( NOUT, FMT = 9998 )'DPBSVX', FACT, UPLO, N, KD, IMAT, K, RESULT( K )
                              }
                              NFAIL = NFAIL + 1
                           }
   50                   CONTINUE
                        NRUN = NRUN + 7 - K1
   60                CONTINUE
   70             CONTINUE
   80          CONTINUE
   90       CONTINUE
  100    CONTINUE
  110 CONTINUE

      // Print a summary of the results.

      CALL ALASVM( PATH, NOUT, NFAIL, NRUN, NERRS )

 9999 FORMAT( 1X, A, ', UPLO=''', A1, ''', N =', I5, ', KD =', I5,
     $      ', type ', I1, ', test(', I1, ')=', G12.5 )
 9998 FORMAT( 1X, A, '( ''', A1, ''', ''', A1, ''', ', I5, ', ', I5,
     $      ', ... ), type ', I1, ', test(', I1, ')=', G12.5 )
 9997 FORMAT( 1X, A, '( ''', A1, ''', ''', A1, ''', ', I5, ', ', I5,
     $      ', ... ), EQUED=''', A1, ''', type ', I1, ', test(', I1,
     $      ')=', G12.5 )
      RETURN

      // End of DDRVPB

      }
