      SUBROUTINE DCHKSY_RK( DOTYPE, NN, NVAL, NNB, NBVAL, NNS, NSVAL, THRESH, TSTERR, NMAX, A, AFAC, E, AINV, B, X, XACT, WORK, RWORK, IWORK, NOUT )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      LOGICAL            TSTERR
      int                NMAX, NN, NNB, NNS, NOUT
      DOUBLE PRECISION   THRESH
*     ..
*     .. Array Arguments ..
      LOGICAL            DOTYPE( * )
      int                IWORK( * ), NBVAL( * ), NSVAL( * ), NVAL( * )
      DOUBLE PRECISION   A( * ), AFAC( * ), AINV( * ), B( * ), E( * ), RWORK( * ), WORK( * ), X( * ), XACT( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      DOUBLE PRECISION   EIGHT, SEVTEN
      PARAMETER          ( EIGHT = 8.0D+0, SEVTEN = 17.0D+0 )
      int                NTYPES
      PARAMETER          ( NTYPES = 10 )
      int                NTESTS
      PARAMETER          ( NTESTS = 7 )
*     ..
*     .. Local Scalars ..
      LOGICAL            TRFCON, ZEROT
      String             DIST, TYPE, UPLO, XTYPE;
      String             PATH, MATPATH;
      int                I, I1, I2, IMAT, IN, INB, INFO, IOFF, IRHS, ITEMP, IUPLO, IZERO, J, K, KL, KU, LDA, LWORK, MODE, N, NB, NERRS, NFAIL, NIMAT, NRHS, NRUN, NT
      DOUBLE PRECISION   ALPHA, ANORM, CNDNUM, CONST, DTEMP, SING_MAX, SING_MIN, RCOND, RCONDC
*     ..
*     .. Local Arrays ..
      String             UPLOS( 2 );
      int                IDUMMY( 1 ), ISEED( 4 ), ISEEDY( 4 )
      DOUBLE PRECISION   BLOCK( 2, 2 ), DDUMMY( 1 ), RESULT( NTESTS )
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DGET06, DLANGE, DLANSY
      EXTERNAL           DGET06, DLANGE, DLANSY
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALAERH, ALAHD, ALASUM, DERRSY, DGESVD, DGET04, DLACPY, DLARHS, DLATB4, DLATMS, DPOT02, DPOT03, DSYCON_3, DSYT01_3, DSYTRF_RK, DSYTRI_3, DSYTRS_3, XLAENV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, SQRT
*     ..
*     .. Scalars in Common ..
      LOGICAL            LERR, OK
      String             SRNAMT;
      int                INFOT, NUNIT
*     ..
*     .. Common blocks ..
      COMMON             / INFOC / INFOT, NUNIT, OK, LERR
      COMMON             / SRNAMC / SRNAMT
*     ..
*     .. Data statements ..
      DATA               ISEEDY / 1988, 1989, 1990, 1991 /
      DATA               UPLOS / 'U', 'L' /
*     ..
*     .. Executable Statements ..
*
*     Initialize constants and the random number seed.
*
      ALPHA = ( ONE+SQRT( SEVTEN ) ) / EIGHT
*
*     Test path
*
      PATH( 1: 1 ) = 'Double precision'
      PATH( 2: 3 ) = 'SK'
*
*     Path to generate matrices
*
      MATPATH( 1: 1 ) = 'Double precision'
      MATPATH( 2: 3 ) = 'SY'
*
      NRUN = 0
      NFAIL = 0
      NERRS = 0
      DO 10 I = 1, 4
         ISEED( I ) = ISEEDY( I )
   10 CONTINUE
*
*     Test the error exits
*
      IF( TSTERR ) CALL DERRSY( PATH, NOUT )
      INFOT = 0
*
*     Set the minimum block size for which the block routine should
*     be used, which will be later returned by ILAENV
*
      CALL XLAENV( 2, 2 )
*
*     Do for each value of N in NVAL
*
      DO 270 IN = 1, NN
         N = NVAL( IN )
         LDA = MAX( N, 1 )
         XTYPE = 'N'
         NIMAT = NTYPES
         IF( N.LE.0 ) NIMAT = 1
*
         IZERO = 0
*
*        Do for each value of matrix type IMAT
*
         DO 260 IMAT = 1, NIMAT
*
*           Do the tests only if DOTYPE( IMAT ) is true.
*
            IF( .NOT.DOTYPE( IMAT ) ) GO TO 260
*
*           Skip types 3, 4, 5, or 6 if the matrix size is too small.
*
            ZEROT = IMAT.GE.3 .AND. IMAT.LE.6
            IF( ZEROT .AND. N.LT.IMAT-2 ) GO TO 260
*
*           Do first for UPLO = 'U', then for UPLO = 'L'
*
            DO 250 IUPLO = 1, 2
               UPLO = UPLOS( IUPLO )
*
*              Begin generate the test matrix A.
*
*              Set up parameters with DLATB4 for the matrix generator
*              based on the type of matrix to be generated.
*
               CALL DLATB4( MATPATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST )
*
*              Generate a matrix with DLATMS.
*
               SRNAMT = 'DLATMS'
               CALL DLATMS( N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU, UPLO, A, LDA, WORK, INFO )
*
*              Check error code from DLATMS and handle error.
*
               IF( INFO.NE.0 ) THEN
                  CALL ALAERH( PATH, 'DLATMS', INFO, 0, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT )
*
*                 Skip all tests for this generated matrix
*
                  GO TO 250
               END IF
*
*              For matrix types 3-6, zero one or more rows and
*              columns of the matrix to test that INFO is returned
*              correctly.
*
               IF( ZEROT ) THEN
                  IF( IMAT.EQ.3 ) THEN
                     IZERO = 1
                  ELSE IF( IMAT.EQ.4 ) THEN
                     IZERO = N
                  ELSE
                     IZERO = N / 2 + 1
                  END IF
*
                  IF( IMAT.LT.6 ) THEN
*
*                    Set row and column IZERO to zero.
*
                     IF( IUPLO.EQ.1 ) THEN
                        IOFF = ( IZERO-1 )*LDA
                        DO 20 I = 1, IZERO - 1
                           A( IOFF+I ) = ZERO
   20                   CONTINUE
                        IOFF = IOFF + IZERO
                        DO 30 I = IZERO, N
                           A( IOFF ) = ZERO
                           IOFF = IOFF + LDA
   30                   CONTINUE
                     ELSE
                        IOFF = IZERO
                        DO 40 I = 1, IZERO - 1
                           A( IOFF ) = ZERO
                           IOFF = IOFF + LDA
   40                   CONTINUE
                        IOFF = IOFF - IZERO
                        DO 50 I = IZERO, N
                           A( IOFF+I ) = ZERO
   50                   CONTINUE
                     END IF
                  ELSE
                     IF( IUPLO.EQ.1 ) THEN
*
*                       Set the first IZERO rows and columns to zero.
*
                        IOFF = 0
                        DO 70 J = 1, N
                           I2 = MIN( J, IZERO )
                           DO 60 I = 1, I2
                              A( IOFF+I ) = ZERO
   60                      CONTINUE
                           IOFF = IOFF + LDA
   70                   CONTINUE
                     ELSE
*
*                       Set the last IZERO rows and columns to zero.
*
                        IOFF = 0
                        DO 90 J = 1, N
                           I1 = MAX( J, IZERO )
                           DO 80 I = I1, N
                              A( IOFF+I ) = ZERO
   80                      CONTINUE
                           IOFF = IOFF + LDA
   90                   CONTINUE
                     END IF
                  END IF
               ELSE
                  IZERO = 0
               END IF
*
*              End generate the test matrix A.
*
*
*              Do for each value of NB in NBVAL
*
               DO 240 INB = 1, NNB
*
*                 Set the optimal blocksize, which will be later
*                 returned by ILAENV.
*
                  NB = NBVAL( INB )
                  CALL XLAENV( 1, NB )
*
*                 Copy the test matrix A into matrix AFAC which
*                 will be factorized in place. This is needed to
*                 preserve the test matrix A for subsequent tests.
*
                  CALL DLACPY( UPLO, N, N, A, LDA, AFAC, LDA )
*
*                 Compute the L*D*L**T or U*D*U**T factorization of the
*                 matrix. IWORK stores details of the interchanges and
*                 the block structure of D. AINV is a work array for
*                 block factorization, LWORK is the length of AINV.
*
                  LWORK = MAX( 2, NB )*LDA
                  SRNAMT = 'DSYTRF_RK'
                  CALL DSYTRF_RK( UPLO, N, AFAC, LDA, E, IWORK, AINV, LWORK, INFO )
*
*                 Adjust the expected value of INFO to account for
*                 pivoting.
*
                  K = IZERO
                  IF( K.GT.0 ) THEN
  100                CONTINUE
                     IF( IWORK( K ).LT.0 ) THEN
                        IF( IWORK( K ).NE.-K ) THEN
                           K = -IWORK( K )
                           GO TO 100
                        END IF
                     ELSE IF( IWORK( K ).NE.K ) THEN
                        K = IWORK( K )
                        GO TO 100
                     END IF
                  END IF
*
*                 Check error code from DSYTRF_RK and handle error.
*
                  IF( INFO.NE.K) CALL ALAERH( PATH, 'DSYTRF_RK', INFO, K, UPLO, N, N, -1, -1, NB, IMAT, NFAIL, NERRS, NOUT )
*
*                 Set the condition estimate flag if the INFO is not 0.
*
                  IF( INFO.NE.0 ) THEN
                     TRFCON = .TRUE.
                  ELSE
                     TRFCON = .FALSE.
                  END IF
*
*+    TEST 1
*                 Reconstruct matrix from factors and compute residual.
*
                  CALL DSYT01_3( UPLO, N, A, LDA, AFAC, LDA, E, IWORK, AINV, LDA, RWORK, RESULT( 1 ) )
                  NT = 1
*
*+    TEST 2
*                 Form the inverse and compute the residual,
*                 if the factorization was competed without INFO > 0
*                 (i.e. there is no zero rows and columns).
*                 Do it only for the first block size.
*
                  IF( INB.EQ.1 .AND. .NOT.TRFCON ) THEN
                     CALL DLACPY( UPLO, N, N, AFAC, LDA, AINV, LDA )
                     SRNAMT = 'DSYTRI_3'
*
*                    Another reason that we need to compute the inverse
*                    is that DPOT03 produces RCONDC which is used later
*                    in TEST6 and TEST7.
*
                     LWORK = (N+NB+1)*(NB+3)
                     CALL DSYTRI_3( UPLO, N, AINV, LDA, E, IWORK, WORK, LWORK, INFO )
*
*                    Check error code from DSYTRI_3 and handle error.
*
                     IF( INFO.NE.0 ) CALL ALAERH( PATH, 'DSYTRI_3', INFO, -1, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT )
*
*                    Compute the residual for a symmetric matrix times
*                    its inverse.
*
                     CALL DPOT03( UPLO, N, A, LDA, AINV, LDA, WORK, LDA, RWORK, RCONDC, RESULT( 2 ) )
                     NT = 2
                  END IF
*
*                 Print information about the tests that did not pass
*                 the threshold.
*
                  DO 110 K = 1, NT
                     IF( RESULT( K ).GE.THRESH ) THEN
                        IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALAHD( NOUT, PATH )                         WRITE( NOUT, FMT = 9999 )UPLO, N, NB, IMAT, K, RESULT( K )
                        NFAIL = NFAIL + 1
                     END IF
  110             CONTINUE
                  NRUN = NRUN + NT
*
*+    TEST 3
*                 Compute largest element in U or L
*
                  RESULT( 3 ) = ZERO
                  DTEMP = ZERO
*
                  CONST = ONE / ( ONE-ALPHA )
*
                  IF( IUPLO.EQ.1 ) THEN
*
*                 Compute largest element in U
*
                     K = N
  120                CONTINUE
                     IF( K.LE.1 ) GO TO 130
*
                     IF( IWORK( K ).GT.ZERO ) THEN
*
*                       Get max absolute value from elements
*                       in column k in in U
*
                        DTEMP = DLANGE( 'M', K-1, 1, AFAC( ( K-1 )*LDA+1 ), LDA, RWORK )
                     ELSE
*
*                       Get max absolute value from elements
*                       in columns k and k-1 in U
*
                        DTEMP = DLANGE( 'M', K-2, 2, AFAC( ( K-2 )*LDA+1 ), LDA, RWORK )
                        K = K - 1
*
                     END IF
*
*                    DTEMP should be bounded by CONST
*
                     DTEMP = DTEMP - CONST + THRESH
                     IF( DTEMP.GT.RESULT( 3 ) ) RESULT( 3 ) = DTEMP
*
                     K = K - 1
*
                     GO TO 120
  130                CONTINUE
*
                  ELSE
*
*                 Compute largest element in L
*
                     K = 1
  140                CONTINUE
                     IF( K.GE.N ) GO TO 150
*
                     IF( IWORK( K ).GT.ZERO ) THEN
*
*                       Get max absolute value from elements
*                       in column k in in L
*
                        DTEMP = DLANGE( 'M', N-K, 1, AFAC( ( K-1 )*LDA+K+1 ), LDA, RWORK )
                     ELSE
*
*                       Get max absolute value from elements
*                       in columns k and k+1 in L
*
                        DTEMP = DLANGE( 'M', N-K-1, 2, AFAC( ( K-1 )*LDA+K+2 ), LDA, RWORK )
                        K = K + 1
*
                     END IF
*
*                    DTEMP should be bounded by CONST
*
                     DTEMP = DTEMP - CONST + THRESH
                     IF( DTEMP.GT.RESULT( 3 ) ) RESULT( 3 ) = DTEMP
*
                     K = K + 1
*
                     GO TO 140
  150                CONTINUE
                  END IF
*
*+    TEST 4
*                 Compute largest 2-Norm (condition number)
*                 of 2-by-2 diag blocks
*
                  RESULT( 4 ) = ZERO
                  DTEMP = ZERO
*
                  CONST = ( ONE+ALPHA ) / ( ONE-ALPHA )
                  CALL DLACPY( UPLO, N, N, AFAC, LDA, AINV, LDA )
*
                  IF( IUPLO.EQ.1 ) THEN
*
*                    Loop backward for UPLO = 'U'
*
                     K = N
  160                CONTINUE
                     IF( K.LE.1 ) GO TO 170
*
                     IF( IWORK( K ).LT.ZERO ) THEN
*
*                       Get the two singular values
*                       (real and non-negative) of a 2-by-2 block,
*                       store them in RWORK array
*
                        BLOCK( 1, 1 ) = AFAC( ( K-2 )*LDA+K-1 )
                        BLOCK( 1, 2 ) = E( K )
                        BLOCK( 2, 1 ) = BLOCK( 1, 2 )
                        BLOCK( 2, 2 ) = AFAC( (K-1)*LDA+K )
*
                        CALL DGESVD( 'N', 'N', 2, 2, BLOCK, 2, RWORK, DDUMMY, 1, DDUMMY, 1, WORK, 10, INFO )
*
                        SING_MAX = RWORK( 1 )
                        SING_MIN = RWORK( 2 )
*
                        DTEMP = SING_MAX / SING_MIN
*
*                       DTEMP should be bounded by CONST
*
                        DTEMP = DTEMP - CONST + THRESH
                        IF( DTEMP.GT.RESULT( 4 ) ) RESULT( 4 ) = DTEMP
                        K = K - 1
*
                     END IF
*
                     K = K - 1
*
                     GO TO 160
  170                CONTINUE
*
                  ELSE
*
*                    Loop forward for UPLO = 'L'
*
                     K = 1
  180                CONTINUE
                     IF( K.GE.N ) GO TO 190
*
                     IF( IWORK( K ).LT.ZERO ) THEN
*
*                       Get the two singular values
*                       (real and non-negative) of a 2-by-2 block,
*                       store them in RWORK array
*
                        BLOCK( 1, 1 ) = AFAC( ( K-1 )*LDA+K )
                        BLOCK( 2, 1 ) = E( K )
                        BLOCK( 1, 2 ) = BLOCK( 2, 1 )
                        BLOCK( 2, 2 ) = AFAC( K*LDA+K+1 )
*
                        CALL DGESVD( 'N', 'N', 2, 2, BLOCK, 2, RWORK, DDUMMY, 1, DDUMMY, 1, WORK, 10, INFO )
*
*
                        SING_MAX = RWORK( 1 )
                        SING_MIN = RWORK( 2 )
*
                        DTEMP = SING_MAX / SING_MIN
*
*                       DTEMP should be bounded by CONST
*
                        DTEMP = DTEMP - CONST + THRESH
                        IF( DTEMP.GT.RESULT( 4 ) ) RESULT( 4 ) = DTEMP
                        K = K + 1
*
                     END IF
*
                     K = K + 1
*
                     GO TO 180
  190                CONTINUE
                  END IF
*
*                 Print information about the tests that did not pass
*                 the threshold.
*
                  DO 200 K = 3, 4
                     IF( RESULT( K ).GE.THRESH ) THEN
                        IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALAHD( NOUT, PATH )                         WRITE( NOUT, FMT = 9999 )UPLO, N, NB, IMAT, K, RESULT( K )
                        NFAIL = NFAIL + 1
                     END IF
  200             CONTINUE
                  NRUN = NRUN + 2
*
*                 Skip the other tests if this is not the first block
*                 size.
*
                  IF( INB.GT.1 ) GO TO 240
*
*                 Do only the condition estimate if INFO is not 0.
*
                  IF( TRFCON ) THEN
                     RCONDC = ZERO
                     GO TO 230
                  END IF
*
*                 Do for each value of NRHS in NSVAL.
*
                  DO 220 IRHS = 1, NNS
                     NRHS = NSVAL( IRHS )
*
*+    TEST 5 ( Using TRS_3)
*                 Solve and compute residual for  A * X = B.
*
*                    Choose a set of NRHS random solution vectors
*                    stored in XACT and set up the right hand side B
*
                     SRNAMT = 'DLARHS'
                     CALL DLARHS( MATPATH, XTYPE, UPLO, ' ', N, N, KL, KU, NRHS, A, LDA, XACT, LDA, B, LDA, ISEED, INFO )
                     CALL DLACPY( 'Full', N, NRHS, B, LDA, X, LDA )
*
                     SRNAMT = 'DSYTRS_3'
                     CALL DSYTRS_3( UPLO, N, NRHS, AFAC, LDA, E, IWORK, X, LDA, INFO )
*
*                    Check error code from DSYTRS_3 and handle error.
*
                     IF( INFO.NE.0 ) CALL ALAERH( PATH, 'DSYTRS_3', INFO, 0, UPLO, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT )
*
                     CALL DLACPY( 'Full', N, NRHS, B, LDA, WORK, LDA )
*
*                    Compute the residual for the solution
*
                     CALL DPOT02( UPLO, N, NRHS, A, LDA, X, LDA, WORK, LDA, RWORK, RESULT( 5 ) )
*
*+    TEST 6
*                    Check solution from generated exact solution.
*
                     CALL DGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 6 ) )
*
*                    Print information about the tests that did not pass
*                    the threshold.
*
                     DO 210 K = 5, 6
                        IF( RESULT( K ).GE.THRESH ) THEN
                           IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALAHD( NOUT, PATH )                            WRITE( NOUT, FMT = 9998 )UPLO, N, NRHS, IMAT, K, RESULT( K )
                           NFAIL = NFAIL + 1
                        END IF
  210                CONTINUE
                     NRUN = NRUN + 2
*
*                 End do for each value of NRHS in NSVAL.
*
  220             CONTINUE
*
*+    TEST 7
*                 Get an estimate of RCOND = 1/CNDNUM.
*
  230             CONTINUE
                  ANORM = DLANSY( '1', UPLO, N, A, LDA, RWORK )
                  SRNAMT = 'DSYCON_3'
                  CALL DSYCON_3( UPLO, N, AFAC, LDA, E, IWORK, ANORM, RCOND, WORK, IWORK( N+1 ), INFO )
*
*                 Check error code from DSYCON_3 and handle error.
*
                  IF( INFO.NE.0 ) CALL ALAERH( PATH, 'DSYCON_3', INFO, 0, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT )
*
*                 Compute the test ratio to compare to values of RCOND
*
                  RESULT( 7 ) = DGET06( RCOND, RCONDC )
*
*                 Print information about the tests that did not pass
*                 the threshold.
*
                  IF( RESULT( 7 ).GE.THRESH ) THEN
                     IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALAHD( NOUT, PATH )                      WRITE( NOUT, FMT = 9997 ) UPLO, N, IMAT, 7, RESULT( 7 )
                     NFAIL = NFAIL + 1
                  END IF
                  NRUN = NRUN + 1
  240          CONTINUE
*
  250       CONTINUE
  260    CONTINUE
  270 CONTINUE
*
*     Print a summary of the results.
*
      CALL ALASUM( PATH, NOUT, NFAIL, NRUN, NERRS )
*
 9999 FORMAT( ' UPLO = ''', A1, ''', N =', I5, ', NB =', I4, ', type ',
     $      I2, ', test ', I2, ', ratio =', G12.5 )
 9998 FORMAT( ' UPLO = ''', A1, ''', N =', I5, ', NRHS=', I3, ', type ',
     $      I2, ', test(', I2, ') =', G12.5 )
 9997 FORMAT( ' UPLO = ''', A1, ''', N =', I5, ',', 10X, ' type ', I2,
     $      ', test(', I2, ') =', G12.5 )
      RETURN
*
*     End of DCHKSY_RK
*
      END
