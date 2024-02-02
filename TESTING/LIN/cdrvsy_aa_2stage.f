      SUBROUTINE CDRVSY_AA_2STAGE(
     $                         DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR,
     $                         NMAX, A, AFAC, AINV, B, X, XACT, WORK,
     $                         RWORK, IWORK, NOUT )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      LOGICAL            TSTERR
      INTEGER            NMAX, NN, NOUT, NRHS
      REAL               THRESH
*     ..
*     .. Array Arguments ..
      LOGICAL            DOTYPE( * )
      INTEGER            IWORK( * ), NVAL( * )
      REAL               RWORK( * )
      COMPLEX            A( * ), AFAC( * ), AINV( * ), B( * ),
     $                   WORK( * ), X( * ), XACT( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX            CZERO
      PARAMETER          ( CZERO = ( 0.0E+0, 0.0E+0 ) )
      INTEGER            NTYPES, NTESTS
      PARAMETER          ( NTYPES = 10, NTESTS = 3 )
      INTEGER            NFACT
      PARAMETER          ( NFACT = 2 )
*     ..
*     .. Local Scalars ..
      LOGICAL            ZEROT
      CHARACTER          DIST, FACT, TYPE, UPLO, XTYPE
      CHARACTER*3        MATPATH, PATH
      INTEGER            I, I1, I2, IFACT, IMAT, IN, INFO, IOFF, IUPLO,
     $                   IZERO, J, K, KL, KU, LDA, LWORK, MODE, N,
     $                   NB, NBMIN, NERRS, NFAIL, NIMAT, NRUN, NT
      REAL               ANORM, CNDNUM
*     ..
*     .. Local Arrays ..
      CHARACTER          FACTS( NFACT ), UPLOS( 2 )
      INTEGER            ISEED( 4 ), ISEEDY( 4 )
      REAL               RESULT( NTESTS )
*     ..
*     .. External Functions ..
      COMPLEX            CLANSY, SGET06
      EXTERNAL           CLANSY, SGET06
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALADHD, ALAERH, ALASVM, XLAENV, CERRVX,
     $                   CGET04, CLACPY, CLARHS, CLATB4, CLATMS,
     $                   CSYSV_AA_2STAGE, CSYT01_AA, CSYT02,
     $                   CSYTRF_AA_2STAGE
*     ..
*     .. Scalars in Common ..
      LOGICAL            LERR, OK
      CHARACTER*32       SRNAMT
      INTEGER            INFOT, NUNIT
*     ..
*     .. Common blocks ..
      COMMON             / INFOC / INFOT, NUNIT, OK, LERR
      COMMON             / SRNAMC / SRNAMT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CMPLX, MAX, MIN
*     ..
*     .. Data statements ..
      DATA               ISEEDY / 1988, 1989, 1990, 1991 /
      DATA               UPLOS / 'U', 'L' / , FACTS / 'F', 'N' /
*     ..
*     .. Executable Statements ..
*
*     Initialize constants and the random number seed.
*
*     Test path
*
      PATH( 1: 1 ) = 'Complex precision'
      PATH( 2: 3 ) = 'S2'
*
*     Path to generate matrices
*
      MATPATH( 1: 1 ) = 'Complex precision'
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
      IF( TSTERR )
     $   CALL CERRVX( PATH, NOUT )
      INFOT = 0
*
*     Set the block size and minimum block size for testing.
*
      NB = 1
      NBMIN = 2
      CALL XLAENV( 1, NB )
      CALL XLAENV( 2, NBMIN )
*
*     Do for each value of N in NVAL
*
      DO 180 IN = 1, NN
         N = NVAL( IN )
         LDA = MAX( N, 1 )
         XTYPE = 'N'
         NIMAT = NTYPES
         IF( N.LE.0 )
     $      NIMAT = 1
*
         DO 170 IMAT = 1, NIMAT
*
*           Do the tests only if DOTYPE( IMAT ) is true.
*
            IF( .NOT.DOTYPE( IMAT ) )
     $         GO TO 170
*
*           Skip types 3, 4, 5, or 6 if the matrix size is too small.
*
            ZEROT = IMAT.GE.3 .AND. IMAT.LE.6
            IF( ZEROT .AND. N.LT.IMAT-2 )
     $         GO TO 170
*
*           Do first for UPLO = 'U', then for UPLO = 'L'
*
            DO 160 IUPLO = 1, 2
               UPLO = UPLOS( IUPLO )
*
*              Begin generate the test matrix A.
*
*              Set up parameters with CLATB4 for the matrix generator
*              based on the type of matrix to be generated.
*
              CALL CLATB4( MATPATH, IMAT, N, N, TYPE, KL, KU, ANORM,
     $                         MODE, CNDNUM, DIST )
*
*              Generate a matrix with CLATMS.
*
                  SRNAMT = 'CLATMS'
                  CALL CLATMS( N, N, DIST, ISEED, TYPE, RWORK, MODE,
     $                         CNDNUM, ANORM, KL, KU, UPLO, A, LDA,
     $                         WORK, INFO )
*
*                 Check error code from CLATMS and handle error.
*
                  IF( INFO.NE.0 ) THEN
                     CALL ALAERH( PATH, 'CLATMS', INFO, 0, UPLO, N, N,
     $                            -1, -1, -1, IMAT, NFAIL, NERRS, NOUT )
                     GO TO 160
                  END IF
*
*                 For types 3-6, zero one or more rows and columns of
*                 the matrix to test that INFO is returned correctly.
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
*                       Set row and column IZERO to zero.
*
                        IF( IUPLO.EQ.1 ) THEN
                           IOFF = ( IZERO-1 )*LDA
                           DO 20 I = 1, IZERO - 1
                              A( IOFF+I ) = CZERO
   20                      CONTINUE
                           IOFF = IOFF + IZERO
                           DO 30 I = IZERO, N
                              A( IOFF ) = CZERO
                              IOFF = IOFF + LDA
   30                      CONTINUE
                        ELSE
                           IOFF = IZERO
                           DO 40 I = 1, IZERO - 1
                              A( IOFF ) = CZERO
                              IOFF = IOFF + LDA
   40                      CONTINUE
                           IOFF = IOFF - IZERO
                           DO 50 I = IZERO, N
                              A( IOFF+I ) = CZERO
   50                      CONTINUE
                        END IF
                     ELSE
                        IOFF = 0
                        IF( IUPLO.EQ.1 ) THEN
*
*                       Set the first IZERO rows and columns to zero.
*
                           DO 70 J = 1, N
                              I2 = MIN( J, IZERO )
                              DO 60 I = 1, I2
                                 A( IOFF+I ) = CZERO
   60                         CONTINUE
                              IOFF = IOFF + LDA
   70                      CONTINUE
                           IZERO = 1
                        ELSE
*
*                       Set the first IZERO rows and columns to zero.
*
                           IOFF = 0
                           DO 90 J = 1, N
                              I1 = MAX( J, IZERO )
                              DO 80 I = I1, N
                                 A( IOFF+I ) = CZERO
   80                         CONTINUE
                              IOFF = IOFF + LDA
   90                      CONTINUE
                        END IF
                     END IF
                  ELSE
                     IZERO = 0
                  END IF
*
*                 End generate the test matrix A.
*
*
               DO 150 IFACT = 1, NFACT
*
*                 Do first for FACT = 'F', then for other values.
*
                  FACT = FACTS( IFACT )
*
*                 Form an exact solution and set the right hand side.
*
                  SRNAMT = 'CLARHS'
                  CALL CLARHS( MATPATH, XTYPE, UPLO, ' ', N, N, KL, KU,
     $                         NRHS, A, LDA, XACT, LDA, B, LDA, ISEED,
     $                         INFO )
                  XTYPE = 'C'
*
*                 --- Test CSYSV_AA_2STAGE  ---
*
                  IF( IFACT.EQ.2 ) THEN
                     CALL CLACPY( UPLO, N, N, A, LDA, AFAC, LDA )
                     CALL CLACPY( 'Full', N, NRHS, B, LDA, X, LDA )
*
*                    Factor the matrix and solve the system using CSYSV_AA.
*
                     SRNAMT = 'CSYSV_AA_2STAGE '
                     LWORK = MIN(N*NB, 3*NMAX*NMAX)
                     CALL CSYSV_AA_2STAGE( UPLO, N, NRHS, AFAC, LDA,
     $                                 AINV, (3*NB+1)*N, 
     $                                 IWORK, IWORK( 1+N ),
     $                                 X, LDA, WORK, LWORK, INFO )
*
*                    Adjust the expected value of INFO to account for
*                    pivoting.
*
                     IF( IZERO.GT.0 ) THEN
                        J = 1
                        K = IZERO
  100                   CONTINUE
                        IF( J.EQ.K ) THEN
                           K = IWORK( J )
                        ELSE IF( IWORK( J ).EQ.K ) THEN
                           K = J
                        END IF
                        IF( J.LT.K ) THEN
                           J = J + 1
                           GO TO 100
                        END IF
                     ELSE
                        K = 0
                     END IF
*
*                    Check error code from CSYSV_AA .
*
                     IF( INFO.NE.K ) THEN
                        CALL ALAERH( PATH, 'CSYSV_AA', INFO, K,
     $                               UPLO, N, N, -1, -1, NRHS,
     $                               IMAT, NFAIL, NERRS, NOUT )
                        GO TO 120
                     ELSE IF( INFO.NE.0 ) THEN
                        GO TO 120
                     END IF
*
*                    Compute residual of the computed solution.
*
                     CALL CLACPY( 'Full', N, NRHS, B, LDA, WORK, LDA )
                     CALL CSYT02( UPLO, N, NRHS, A, LDA, X, LDA, WORK,
     $                            LDA, RWORK, RESULT( 1 ) )
*
*                    Reconstruct matrix from factors and compute
*                    residual.
*
c                     CALL CSY01_AA( UPLO, N, A, LDA, AFAC, LDA,
c     $                                  IWORK, AINV, LDA, RWORK,
c     $                                  RESULT( 2 ) )
c                     NT = 2
                     NT = 1
*
*                    Print information about the tests that did not pass
*                    the threshold.
*
                     DO 110 K = 1, NT
                        IF( RESULT( K ).GE.THRESH ) THEN
                           IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 )
     $                        CALL ALADHD( NOUT, PATH )
                           WRITE( NOUT, FMT = 9999 )'CSYSV_AA_2STAGE ',
     $                         UPLO, N, IMAT, K, RESULT( K )
                           NFAIL = NFAIL + 1
                        END IF
  110                CONTINUE
                     NRUN = NRUN + NT
  120                CONTINUE
                  END IF
*
  150          CONTINUE
*
  160       CONTINUE
  170    CONTINUE
  180 CONTINUE
*
*     Print a summary of the results.
*
      CALL ALASVM( PATH, NOUT, NFAIL, NRUN, NERRS )
*
 9999 FORMAT( 1X, A, ', UPLO=''', A1, ''', N =', I5, ', type ', I2,
     $      ', test ', I2, ', ratio =', G12.5 )
      RETURN
*
*     End of CSDRVSY_AA_2STAGE
*
      END