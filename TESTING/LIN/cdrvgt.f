      SUBROUTINE CDRVGT( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, A, AF, B, X, XACT, WORK, RWORK, IWORK, NOUT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               TSTERR;
      int                NN, NOUT, NRHS;
      REAL               THRESH
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      int                IWORK( * ), NVAL( * );
      REAL               RWORK( * )
      COMPLEX            A( * ), AF( * ), B( * ), WORK( * ), X( * ), XACT( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO
      const              ONE = 1.0E+0, ZERO = 0.0E+0 ;
      int                NTYPES;
      const              NTYPES = 12 ;
      int                NTESTS;
      const              NTESTS = 6 ;
      // ..
      // .. Local Scalars ..
      bool               TRFCON, ZEROT;
      String             DIST, FACT, TRANS, TYPE;
      String             PATH;
      int                I, IFACT, IMAT, IN, INFO, ITRAN, IX, IZERO, J, K, K1, KL, KOFF, KU, LDA, M, MODE, N, NERRS, NFAIL, NIMAT, NRUN, NT;
      REAL               AINVNM, ANORM, ANORMI, ANORMO, COND, RCOND, RCONDC, RCONDI, RCONDO
      // ..
      // .. Local Arrays ..
      String             TRANSS( 3 );
      int                ISEED( 4 ), ISEEDY( 4 );
      REAL               RESULT( NTESTS ), Z( 3 )
      // ..
      // .. External Functions ..
      REAL               CLANGT, SCASUM, SGET06
      // EXTERNAL CLANGT, SCASUM, SGET06
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALADHD, ALAERH, ALASVM, CCOPY, CERRVX, CGET04, CGTSV, CGTSVX, CGTT01, CGTT02, CGTT05, CGTTRF, CGTTRS, CLACPY, CLAGTM, CLARNV, CLASET, CLATB4, CLATMS, CSSCAL
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CMPLX, MAX
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
      DATA               ISEEDY / 0, 0, 0, 1 / , TRANSS / 'N', 'T', 'C' /
      // ..
      // .. Executable Statements ..

      PATH( 1: 1 ) = 'Complex precision'
      PATH( 2: 3 ) = 'GT'
      NRUN = 0
      NFAIL = 0
      NERRS = 0
      DO 10 I = 1, 4
         ISEED( I ) = ISEEDY( I )
   10 CONTINUE

      // Test the error exits

      IF( TSTERR ) CALL CERRVX( PATH, NOUT )
      INFOT = 0

      DO 140 IN = 1, NN

         // Do for each value of N in NVAL.

         N = NVAL( IN )
         M = MAX( N-1, 0 )
         LDA = MAX( 1, N )
         NIMAT = NTYPES
         IF( N.LE.0 ) NIMAT = 1

         DO 130 IMAT = 1, NIMAT

            // Do the tests only if DOTYPE( IMAT ) is true.

            IF( .NOT.DOTYPE( IMAT ) ) GO TO 130

            // Set up parameters with CLATB4.

            CALL CLATB4( PATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, COND, DIST )

            ZEROT = IMAT.GE.8 .AND. IMAT.LE.10
            if ( IMAT.LE.6 ) {

               // Types 1-6:  generate matrices of known condition number.

               KOFF = MAX( 2-KU, 3-MAX( 1, N ) )
               SRNAMT = 'CLATMS'
               CALL CLATMS( N, N, DIST, ISEED, TYPE, RWORK, MODE, COND, ANORM, KL, KU, 'Z', AF( KOFF ), 3, WORK, INFO )

               // Check the error code from CLATMS.

               if ( INFO.NE.0 ) {
                  CALL ALAERH( PATH, 'CLATMS', INFO, 0, ' ', N, N, KL, KU, -1, IMAT, NFAIL, NERRS, NOUT )
                  GO TO 130
               }
               IZERO = 0

               if ( N.GT.1 ) {
                  CALL CCOPY( N-1, AF( 4 ), 3, A, 1 )
                  CALL CCOPY( N-1, AF( 3 ), 3, A( N+M+1 ), 1 )
               }
               CALL CCOPY( N, AF( 2 ), 3, A( M+1 ), 1 )
            } else {

               // Types 7-12:  generate tridiagonal matrices with
               // unknown condition numbers.

               if ( .NOT.ZEROT .OR. .NOT.DOTYPE( 7 ) ) {

                  // Generate a matrix with elements from [-1,1].

                  CALL CLARNV( 2, ISEED, N+2*M, A )
                  IF( ANORM.NE.ONE ) CALL CSSCAL( N+2*M, ANORM, A, 1 )
               } else if ( IZERO.GT.0 ) {

                  // Reuse the last matrix by copying back the zeroed out
                  // elements.

                  if ( IZERO.EQ.1 ) {
                     A( N ) = Z( 2 )
                     IF( N.GT.1 ) A( 1 ) = Z( 3 )
                  } else if ( IZERO.EQ.N ) {
                     A( 3*N-2 ) = Z( 1 )
                     A( 2*N-1 ) = Z( 2 )
                  } else {
                     A( 2*N-2+IZERO ) = Z( 1 )
                     A( N-1+IZERO ) = Z( 2 )
                     A( IZERO ) = Z( 3 )
                  }
               }

               // If IMAT > 7, set one column of the matrix to 0.

               if ( .NOT.ZEROT ) {
                  IZERO = 0
               } else if ( IMAT.EQ.8 ) {
                  IZERO = 1
                  Z( 2 ) = REAL( A( N ) )
                  A( N ) = ZERO
                  if ( N.GT.1 ) {
                     Z( 3 ) = REAL( A( 1 ) )
                     A( 1 ) = ZERO
                  }
               } else if ( IMAT.EQ.9 ) {
                  IZERO = N
                  Z( 1 ) = REAL( A( 3*N-2 ) )
                  Z( 2 ) = REAL( A( 2*N-1 ) )
                  A( 3*N-2 ) = ZERO
                  A( 2*N-1 ) = ZERO
               } else {
                  IZERO = ( N+1 ) / 2
                  DO 20 I = IZERO, N - 1
                     A( 2*N-2+I ) = ZERO
                     A( N-1+I ) = ZERO
                     A( I ) = ZERO
   20             CONTINUE
                  A( 3*N-2 ) = ZERO
                  A( 2*N-1 ) = ZERO
               }
            }

            DO 120 IFACT = 1, 2
               if ( IFACT.EQ.1 ) {
                  FACT = 'F'
               } else {
                  FACT = 'N'
               }

               // Compute the condition number for comparison with
              t // he value returned by CGTSVX.

               if ( ZEROT ) {
                  IF( IFACT.EQ.1 ) GO TO 120
                  RCONDO = ZERO
                  RCONDI = ZERO

               } else if ( IFACT.EQ.1 ) {
                  CALL CCOPY( N+2*M, A, 1, AF, 1 )

                  // Compute the 1-norm and infinity-norm of A.

                  ANORMO = CLANGT( '1', N, A, A( M+1 ), A( N+M+1 ) )
                  ANORMI = CLANGT( 'I', N, A, A( M+1 ), A( N+M+1 ) )

                  // Factor the matrix A.

                  CALL CGTTRF( N, AF, AF( M+1 ), AF( N+M+1 ), AF( N+2*M+1 ), IWORK, INFO )

                  // Use CGTTRS to solve for one column at a time of
                  // inv(A), computing the maximum column sum as we go.

                  AINVNM = ZERO
                  DO 40 I = 1, N
                     DO 30 J = 1, N
                        X( J ) = ZERO
   30                CONTINUE
                     X( I ) = ONE
                     CALL CGTTRS( 'No transpose', N, 1, AF, AF( M+1 ), AF( N+M+1 ), AF( N+2*M+1 ), IWORK, X, LDA, INFO )
                     AINVNM = MAX( AINVNM, SCASUM( N, X, 1 ) )
   40             CONTINUE

                  // Compute the 1-norm condition number of A.

                  if ( ANORMO.LE.ZERO .OR. AINVNM.LE.ZERO ) {
                     RCONDO = ONE
                  } else {
                     RCONDO = ( ONE / ANORMO ) / AINVNM
                  }

                  // Use CGTTRS to solve for one column at a time of
                  // inv(A'), computing the maximum column sum as we go.

                  AINVNM = ZERO
                  DO 60 I = 1, N
                     DO 50 J = 1, N
                        X( J ) = ZERO
   50                CONTINUE
                     X( I ) = ONE
                     CALL CGTTRS( 'Conjugate transpose', N, 1, AF, AF( M+1 ), AF( N+M+1 ), AF( N+2*M+1 ), IWORK, X, LDA, INFO )
                     AINVNM = MAX( AINVNM, SCASUM( N, X, 1 ) )
   60             CONTINUE

                  // Compute the infinity-norm condition number of A.

                  if ( ANORMI.LE.ZERO .OR. AINVNM.LE.ZERO ) {
                     RCONDI = ONE
                  } else {
                     RCONDI = ( ONE / ANORMI ) / AINVNM
                  }
               }

               DO 110 ITRAN = 1, 3
                  TRANS = TRANSS( ITRAN )
                  if ( ITRAN.EQ.1 ) {
                     RCONDC = RCONDO
                  } else {
                     RCONDC = RCONDI
                  }

                  // Generate NRHS random solution vectors.

                  IX = 1
                  DO 70 J = 1, NRHS
                     CALL CLARNV( 2, ISEED, N, XACT( IX ) )
                     IX = IX + LDA
   70             CONTINUE

                  // Set the right hand side.

                  CALL CLAGTM( TRANS, N, NRHS, ONE, A, A( M+1 ), A( N+M+1 ), XACT, LDA, ZERO, B, LDA )

                  if ( IFACT.EQ.2 .AND. ITRAN.EQ.1 ) {

                     // --- Test CGTSV  ---

                     // Solve the system using Gaussian elimination with
                     // partial pivoting.

                     CALL CCOPY( N+2*M, A, 1, AF, 1 )
                     CALL CLACPY( 'Full', N, NRHS, B, LDA, X, LDA )

                     SRNAMT = 'CGTSV '
                     CALL CGTSV( N, NRHS, AF, AF( M+1 ), AF( N+M+1 ), X, LDA, INFO )

                     // Check error code from CGTSV .

                     IF( INFO.NE.IZERO ) CALL ALAERH( PATH, 'CGTSV ', INFO, IZERO, ' ', N, N, 1, 1, NRHS, IMAT, NFAIL, NERRS, NOUT )
                     NT = 1
                     if ( IZERO.EQ.0 ) {

                        // Check residual of computed solution.

                        CALL CLACPY( 'Full', N, NRHS, B, LDA, WORK, LDA )                         CALL CGTT02( TRANS, N, NRHS, A, A( M+1 ), A( N+M+1 ), X, LDA, WORK, LDA, RESULT( 2 ) )

                        // Check solution from generated exact solution.

                        CALL CGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 3 ) )
                        NT = 3
                     }

                     // Print information about the tests that did not pass
                    t // he threshold.

                     DO 80 K = 2, NT
                        if ( RESULT( K ).GE.THRESH ) {
                           IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALADHD( NOUT, PATH )                            WRITE( NOUT, FMT = 9999 )'CGTSV ', N, IMAT, K, RESULT( K )
                           NFAIL = NFAIL + 1
                        }
   80                CONTINUE
                     NRUN = NRUN + NT - 1
                  }

                  // --- Test CGTSVX ---

                  if ( IFACT.GT.1 ) {

                     // Initialize AF to zero.

                     DO 90 I = 1, 3*N - 2
                        AF( I ) = ZERO
   90                CONTINUE
                  }
                  CALL CLASET( 'Full', N, NRHS, CMPLX( ZERO ), CMPLX( ZERO ), X, LDA )

                  // Solve the system and compute the condition number and
                  // error bounds using CGTSVX.

                  SRNAMT = 'CGTSVX'
                  CALL CGTSVX( FACT, TRANS, N, NRHS, A, A( M+1 ), A( N+M+1 ), AF, AF( M+1 ), AF( N+M+1 ), AF( N+2*M+1 ), IWORK, B, LDA, X, LDA, RCOND, RWORK, RWORK( NRHS+1 ), WORK, RWORK( 2*NRHS+1 ), INFO )

                  // Check the error code from CGTSVX.

                  IF( INFO.NE.IZERO ) CALL ALAERH( PATH, 'CGTSVX', INFO, IZERO, FACT // TRANS, N, N, 1, 1, NRHS, IMAT, NFAIL, NERRS, NOUT )

                  if ( IFACT.GE.2 ) {

                     // Reconstruct matrix from factors and compute
                     // residual.

                     CALL CGTT01( N, A, A( M+1 ), A( N+M+1 ), AF, AF( M+1 ), AF( N+M+1 ), AF( N+2*M+1 ), IWORK, WORK, LDA, RWORK, RESULT( 1 ) )
                     K1 = 1
                  } else {
                     K1 = 2
                  }

                  if ( INFO.EQ.0 ) {
                     TRFCON = .FALSE.

                     // Check residual of computed solution.

                     CALL CLACPY( 'Full', N, NRHS, B, LDA, WORK, LDA )
                     CALL CGTT02( TRANS, N, NRHS, A, A( M+1 ), A( N+M+1 ), X, LDA, WORK, LDA, RESULT( 2 ) )

                     // Check solution from generated exact solution.

                     CALL CGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 3 ) )

                     // Check the error bounds from iterative refinement.

                     CALL CGTT05( TRANS, N, NRHS, A, A( M+1 ), A( N+M+1 ), B, LDA, X, LDA, XACT, LDA, RWORK, RWORK( NRHS+1 ), RESULT( 4 ) )
                     NT = 5
                  }

                  // Print information about the tests that did not pass
                 t // he threshold.

                  DO 100 K = K1, NT
                     if ( RESULT( K ).GE.THRESH ) {
                        IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALADHD( NOUT, PATH )                         WRITE( NOUT, FMT = 9998 )'CGTSVX', FACT, TRANS, N, IMAT, K, RESULT( K )
                        NFAIL = NFAIL + 1
                     }
  100             CONTINUE

                  // Check the reciprocal of the condition number.

                  RESULT( 6 ) = SGET06( RCOND, RCONDC )
                  if ( RESULT( 6 ).GE.THRESH ) {
                     IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALADHD( NOUT, PATH )                      WRITE( NOUT, FMT = 9998 )'CGTSVX', FACT, TRANS, N, IMAT, K, RESULT( K )
                     NFAIL = NFAIL + 1
                  }
                  NRUN = NRUN + NT - K1 + 2

  110          CONTINUE
  120       CONTINUE
  130    CONTINUE
  140 CONTINUE

      // Print a summary of the results.

      CALL ALASVM( PATH, NOUT, NFAIL, NRUN, NERRS )

 9999 FORMAT( 1X, A, ', N =', I5, ', type ', I2, ', test ', I2,
     $      ', ratio = ', G12.5 )
 9998 FORMAT( 1X, A, ', FACT=''', A1, ''', TRANS=''', A1, ''', N =',
     $      I5, ', type ', I2, ', test ', I2, ', ratio = ', G12.5 )
      RETURN

      // End of CDRVGT

      }
