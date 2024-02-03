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
      COMMON             / SRNAMC / SRNAMT
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
      DO 10 I = 1, 4
         ISEED( I ) = ISEEDY( I )
   10 CONTINUE

      DO 130 IIN = 1, NN

         N = NVAL( IIN )
         LDA = MAX( N, 1 )
         LDB = MAX( N, 1 )

         DO 980 IIS = 1, NNS

            NRHS = NSVAL( IIS )

            DO 120 IIT = 1, NNT

               IMAT = NTVAL( IIT )

               // If N.EQ.0, only consider the first type

               IF( N.EQ.0 .AND. IIT.GE.1 ) GO TO 120

               // Skip types 3, 4, or 5 if the matrix size is too small.

               IF( IMAT.EQ.4 .AND. N.LE.1 ) GO TO 120
               IF( IMAT.EQ.5 .AND. N.LE.2 ) GO TO 120

               // Do first for UPLO = 'U', then for UPLO = 'L'

               DO 110 IUPLO = 1, 2
                  UPLO = UPLOS( IUPLO )

                  // Do first for CFORM = 'N', then for CFORM = 'C'

                  DO 100 IFORM = 1, 2
                     CFORM = FORMS( IFORM )

                     // Set up parameters with DLATB4 and generate a test
                     // matrix with DLATMS.

                     CALL DLATB4( 'DPO', IMAT, N, N, CTYPE, KL, KU, ANORM, MODE, CNDNUM, DIST )

                     SRNAMT = 'DLATMS'
                     CALL DLATMS( N, N, DIST, ISEED, CTYPE, D_WORK_DLATMS, MODE, CNDNUM, ANORM, KL, KU, UPLO, A, LDA, D_WORK_DLATMS, INFO )

                     // Check error code from DLATMS.

                     if ( INFO.NE.0 ) {
                        CALL ALAERH( 'DPF', 'DLATMS', INFO, 0, UPLO, N, N, -1, -1, -1, IIT, NFAIL, NERRS, NOUT )
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
                           DO 20 I = 1, IZERO - 1
                              A( IOFF+I ) = ZERO
   20                      CONTINUE
                           IOFF = IOFF + IZERO
                           DO 30 I = IZERO, N
                              A( IOFF ) = ZERO
                              IOFF = IOFF + LDA
   30                      CONTINUE
                        } else {
                           IOFF = IZERO
                           DO 40 I = 1, IZERO - 1
                              A( IOFF ) = ZERO
                              IOFF = IOFF + LDA
   40                      CONTINUE
                           IOFF = IOFF - IZERO
                           DO 50 I = IZERO, N
                              A( IOFF+I ) = ZERO
   50                      CONTINUE
                        }
                     } else {
                        IZERO = 0
                     }

                     // Save a copy of the matrix A in ASAV.

                     CALL DLACPY( UPLO, N, N, A, LDA, ASAV, LDA )

                     // Compute the condition number of A (RCONDC).

                     if ( ZEROT ) {
                        RCONDC = ZERO
                     } else {

                        // Compute the 1-norm of A.

                        ANORM = DLANSY( '1', UPLO, N, A, LDA, D_WORK_DLANSY )

                        // Factor the matrix A.

                        CALL DPOTRF( UPLO, N, A, LDA, INFO )

                        // Form the inverse of A.

                        CALL DPOTRI( UPLO, N, A, LDA, INFO )

                        if ( N .NE. 0 ) {

                           // Compute the 1-norm condition number of A.

                           AINVNM = DLANSY( '1', UPLO, N, A, LDA, D_WORK_DLANSY )
                           RCONDC = ( ONE / ANORM ) / AINVNM

                           // Restore the matrix A.

                           CALL DLACPY( UPLO, N, N, ASAV, LDA, A, LDA )
                        }

                     }

                     // Form an exact solution and set the right hand side.

                     SRNAMT = 'DLARHS'
                     CALL DLARHS( 'DPO', 'N', UPLO, ' ', N, N, KL, KU, NRHS, A, LDA, XACT, LDA, B, LDA, ISEED, INFO )
                     CALL DLACPY( 'Full', N, NRHS, B, LDA, BSAV, LDA )

                     // Compute the L*L' or U'*U factorization of the
                     // matrix and solve the system.

                     CALL DLACPY( UPLO, N, N, A, LDA, AFAC, LDA )
                     CALL DLACPY( 'Full', N, NRHS, B, LDB, X, LDB )

                     SRNAMT = 'DTRTTF'
                     CALL DTRTTF( CFORM, UPLO, N, AFAC, LDA, ARF, INFO )
                     SRNAMT = 'DPFTRF'
                     CALL DPFTRF( CFORM, UPLO, N, ARF, INFO )

                     // Check error code from DPFTRF.

                     if ( INFO.NE.IZERO ) {

                        // LANGOU: there is a small hick here: IZERO should
                        // always be INFO however if INFO is ZERO, ALAERH does not
                        // complain.

                         CALL ALAERH( 'DPF', 'DPFSV ', INFO, IZERO, UPLO, N, N, -1, -1, NRHS, IIT, NFAIL, NERRS, NOUT )
                         GO TO 100
                      }

                     // Skip the tests if INFO is not 0.

                     if ( INFO.NE.0 ) {
                        GO TO 100
                     }

                     SRNAMT = 'DPFTRS'
                     CALL DPFTRS( CFORM, UPLO, N, NRHS, ARF, X, LDB, INFO )

                     SRNAMT = 'DTFTTR'
                     CALL DTFTTR( CFORM, UPLO, N, ARF, AFAC, LDA, INFO )

                     // Reconstruct matrix from factors and compute
                     // residual.

                     CALL DLACPY( UPLO, N, N, AFAC, LDA, ASAV, LDA )
                     CALL DPOT01( UPLO, N, A, LDA, AFAC, LDA, D_WORK_DPOT01, RESULT( 1 ) )
                     CALL DLACPY( UPLO, N, N, ASAV, LDA, AFAC, LDA )

                     // Form the inverse and compute the residual.

                     if (MOD(N,2).EQ.0) {
                        CALL DLACPY( 'A', N+1, N/2, ARF, N+1, ARFINV, N+1 )
                     } else {
                        CALL DLACPY( 'A', N, (N+1)/2, ARF, N, ARFINV, N )
                     }

                     SRNAMT = 'DPFTRI'
                     CALL DPFTRI( CFORM, UPLO, N, ARFINV , INFO )

                     SRNAMT = 'DTFTTR'
                     CALL DTFTTR( CFORM, UPLO, N, ARFINV, AINV, LDA, INFO )

                     // Check error code from DPFTRI.

                     IF( INFO.NE.0 ) CALL ALAERH( 'DPO', 'DPFTRI', INFO, 0, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT )

                     CALL DPOT03( UPLO, N, A, LDA, AINV, LDA, D_TEMP_DPOT03, LDA, D_WORK_DPOT03, RCONDC, RESULT( 2 ) )

                     // Compute residual of the computed solution.

                     CALL DLACPY( 'Full', N, NRHS, B, LDA, D_TEMP_DPOT02, LDA )                      CALL DPOT02( UPLO, N, NRHS, A, LDA, X, LDA, D_TEMP_DPOT02, LDA, D_WORK_DPOT02, RESULT( 3 ) )

                     // Check solution from generated exact solution.
                      CALL DGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 4 ) )
                     NT = 4

                     // Print information about the tests that did not
                     // pass the threshold.

                     DO 60 K = 1, NT
                        if ( RESULT( K ).GE.THRESH ) {
                           IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALADHD( NOUT, 'DPF' )                            WRITE( NOUT, FMT = 9999 )'DPFSV ', UPLO, N, IIT, K, RESULT( K )
                           NFAIL = NFAIL + 1
                        }
   60                CONTINUE
                     NRUN = NRUN + NT
  100             CONTINUE
  110          CONTINUE
  120       CONTINUE
  980    CONTINUE
  130 CONTINUE

      // Print a summary of the results.

      CALL ALASVM( 'DPF', NOUT, NFAIL, NRUN, NERRS )

 9999 FORMAT( 1X, A6, ', UPLO=''', A1, ''', N =', I5, ', type ', I1,
     +      ', test(', I1, ')=', G12.5 )

      RETURN

      // End of DDRVRFP

      }
