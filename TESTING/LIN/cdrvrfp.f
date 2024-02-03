      SUBROUTINE CDRVRFP( NOUT, NN, NVAL, NNS, NSVAL, NNT, NTVAL, THRESH, A, ASAV, AFAC, AINV, B, BSAV, XACT, X, ARF, ARFINV, C_WORK_CLATMS, C_WORK_CPOT02, C_WORK_CPOT03, S_WORK_CLATMS, S_WORK_CLANHE, S_WORK_CPOT01, S_WORK_CPOT02, S_WORK_CPOT03 )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                NN, NNS, NNT, NOUT
      REAL               THRESH
*     ..
*     .. Array Arguments ..
      int                NVAL( NN ), NSVAL( NNS ), NTVAL( NNT )
      COMPLEX            A( * )
      COMPLEX            AINV( * )
      COMPLEX            ASAV( * )
      COMPLEX            B( * )
      COMPLEX            BSAV( * )
      COMPLEX            AFAC( * )
      COMPLEX            ARF( * )
      COMPLEX            ARFINV( * )
      COMPLEX            XACT( * )
      COMPLEX            X( * )
      COMPLEX            C_WORK_CLATMS( * )
      COMPLEX            C_WORK_CPOT02( * )
      COMPLEX            C_WORK_CPOT03( * )
      REAL               S_WORK_CLATMS( * )
      REAL               S_WORK_CLANHE( * )
      REAL               S_WORK_CPOT01( * )
      REAL               S_WORK_CPOT02( * )
      REAL               S_WORK_CPOT03( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
      int                NTESTS
      PARAMETER          ( NTESTS = 4 )
*     ..
*     .. Local Scalars ..
      LOGICAL            ZEROT
      int                I, INFO, IUPLO, LDA, LDB, IMAT, NERRS, NFAIL, NRHS, NRUN, IZERO, IOFF, K, NT, N, IFORM, IIN, IIT, IIS
      String             DIST, CTYPE, UPLO, CFORM;
      int                KL, KU, MODE
      REAL               ANORM, AINVNM, CNDNUM, RCONDC
*     ..
*     .. Local Arrays ..
      String             UPLOS( 2 ), FORMS( 2 );
      int                ISEED( 4 ), ISEEDY( 4 )
      REAL               RESULT( NTESTS )
*     ..
*     .. External Functions ..
      REAL               CLANHE
      EXTERNAL           CLANHE
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALADHD, ALAERH, ALASVM, CGET04, CTFTTR, CLACPY, CLAIPD, CLARHS, CLATB4, CLATMS, CPFTRI, CPFTRF, CPFTRS, CPOT01, CPOT02, CPOT03, CPOTRI, CPOTRF, CTRTTF
*     ..
*     .. Scalars in Common ..
      String             SRNAMT;
*     ..
*     .. Common blocks ..
      COMMON             / SRNAMC / SRNAMT
*     ..
*     .. Data statements ..
      DATA               ISEEDY / 1988, 1989, 1990, 1991 /
      DATA               UPLOS / 'U', 'L' /
      DATA               FORMS / 'N', 'C' /
*     ..
*     .. Executable Statements ..
*
*     Initialize constants and the random number seed.
*
      NRUN = 0
      NFAIL = 0
      NERRS = 0
      DO 10 I = 1, 4
         ISEED( I ) = ISEEDY( I )
   10 CONTINUE
*
      DO 130 IIN = 1, NN
*
         N = NVAL( IIN )
         LDA = MAX( N, 1 )
         LDB = MAX( N, 1 )
*
         DO 980 IIS = 1, NNS
*
            NRHS = NSVAL( IIS )
*
            DO 120 IIT = 1, NNT
*
               IMAT = NTVAL( IIT )
*
*              If N.EQ.0, only consider the first type
*
               IF( N.EQ.0 .AND. IIT.GE.1 ) GO TO 120
*
*              Skip types 3, 4, or 5 if the matrix size is too small.
*
               IF( IMAT.EQ.4 .AND. N.LE.1 ) GO TO 120
               IF( IMAT.EQ.5 .AND. N.LE.2 ) GO TO 120
*
*              Do first for UPLO = 'U', then for UPLO = 'L'
*
               DO 110 IUPLO = 1, 2
                  UPLO = UPLOS( IUPLO )
*
*                 Do first for CFORM = 'N', then for CFORM = 'C'
*
                  DO 100 IFORM = 1, 2
                     CFORM = FORMS( IFORM )
*
*                    Set up parameters with CLATB4 and generate a test
*                    matrix with CLATMS.
*
                     CALL CLATB4( 'CPO', IMAT, N, N, CTYPE, KL, KU, ANORM, MODE, CNDNUM, DIST )
*
                     SRNAMT = 'CLATMS'
                     CALL CLATMS( N, N, DIST, ISEED, CTYPE, S_WORK_CLATMS, MODE, CNDNUM, ANORM, KL, KU, UPLO, A, LDA, C_WORK_CLATMS, INFO )
*
*                    Check error code from CLATMS.
*
                     IF( INFO.NE.0 ) THEN
                        CALL ALAERH( 'CPF', 'CLATMS', INFO, 0, UPLO, N, N, -1, -1, -1, IIT, NFAIL, NERRS, NOUT )
                        GO TO 100
                     END IF
*
*                    For types 3-5, zero one row and column of the matrix to
*                    test that INFO is returned correctly.
*
                     ZEROT = IMAT.GE.3 .AND. IMAT.LE.5
                     IF( ZEROT ) THEN
                        IF( IIT.EQ.3 ) THEN
                           IZERO = 1
                        ELSE IF( IIT.EQ.4 ) THEN
                           IZERO = N
                        ELSE
                           IZERO = N / 2 + 1
                        END IF
                        IOFF = ( IZERO-1 )*LDA
*
*                       Set row and column IZERO of A to 0.
*
                        IF( IUPLO.EQ.1 ) THEN
                           DO 20 I = 1, IZERO - 1
                              A( IOFF+I ) = ZERO
   20                      CONTINUE
                           IOFF = IOFF + IZERO
                           DO 30 I = IZERO, N
                              A( IOFF ) = ZERO
                              IOFF = IOFF + LDA
   30                      CONTINUE
                        ELSE
                           IOFF = IZERO
                           DO 40 I = 1, IZERO - 1
                              A( IOFF ) = ZERO
                              IOFF = IOFF + LDA
   40                      CONTINUE
                           IOFF = IOFF - IZERO
                           DO 50 I = IZERO, N
                              A( IOFF+I ) = ZERO
   50                      CONTINUE
                        END IF
                     ELSE
                        IZERO = 0
                     END IF
*
*                    Set the imaginary part of the diagonals.
*
                     CALL CLAIPD( N, A, LDA+1, 0 )
*
*                    Save a copy of the matrix A in ASAV.
*
                     CALL CLACPY( UPLO, N, N, A, LDA, ASAV, LDA )
*
*                    Compute the condition number of A (RCONDC).
*
                     IF( ZEROT ) THEN
                        RCONDC = ZERO
                     ELSE
*
*                       Compute the 1-norm of A.
*
                        ANORM = CLANHE( '1', UPLO, N, A, LDA, S_WORK_CLANHE )
*
*                       Factor the matrix A.
*
                        CALL CPOTRF( UPLO, N, A, LDA, INFO )
*
*                       Form the inverse of A.
*
                        CALL CPOTRI( UPLO, N, A, LDA, INFO )

                        IF ( N .NE. 0 ) THEN
*
*                          Compute the 1-norm condition number of A.
*
                           AINVNM = CLANHE( '1', UPLO, N, A, LDA, S_WORK_CLANHE )
                           RCONDC = ( ONE / ANORM ) / AINVNM
*
*                          Restore the matrix A.
*
                           CALL CLACPY( UPLO, N, N, ASAV, LDA, A, LDA )
                        END IF
*
                     END IF
*
*                    Form an exact solution and set the right hand side.
*
                     SRNAMT = 'CLARHS'
                     CALL CLARHS( 'CPO', 'N', UPLO, ' ', N, N, KL, KU, NRHS, A, LDA, XACT, LDA, B, LDA, ISEED, INFO )
                     CALL CLACPY( 'Full', N, NRHS, B, LDA, BSAV, LDA )
*
*                    Compute the L*L' or U'*U factorization of the
*                    matrix and solve the system.
*
                     CALL CLACPY( UPLO, N, N, A, LDA, AFAC, LDA )
                     CALL CLACPY( 'Full', N, NRHS, B, LDB, X, LDB )
*
                     SRNAMT = 'CTRTTF'
                     CALL CTRTTF( CFORM, UPLO, N, AFAC, LDA, ARF, INFO )
                     SRNAMT = 'CPFTRF'
                     CALL CPFTRF( CFORM, UPLO, N, ARF, INFO )
*
*                    Check error code from CPFTRF.
*
                     IF( INFO.NE.IZERO ) THEN
*
*                       LANGOU: there is a small hick here: IZERO should
*                       always be INFO however if INFO is ZERO, ALAERH does not
*                       complain.
*
                         CALL ALAERH( 'CPF', 'CPFSV ', INFO, IZERO, UPLO, N, N, -1, -1, NRHS, IIT, NFAIL, NERRS, NOUT )
                         GO TO 100
                      END IF
*
*                     Skip the tests if INFO is not 0.
*
                     IF( INFO.NE.0 ) THEN
                        GO TO 100
                     END IF
*
                     SRNAMT = 'CPFTRS'
                     CALL CPFTRS( CFORM, UPLO, N, NRHS, ARF, X, LDB, INFO )
*
                     SRNAMT = 'CTFTTR'
                     CALL CTFTTR( CFORM, UPLO, N, ARF, AFAC, LDA, INFO )
*
*                    Reconstruct matrix from factors and compute
*                    residual.
*
                     CALL CLACPY( UPLO, N, N, AFAC, LDA, ASAV, LDA )
                     CALL CPOT01( UPLO, N, A, LDA, AFAC, LDA, S_WORK_CPOT01, RESULT( 1 ) )
                     CALL CLACPY( UPLO, N, N, ASAV, LDA, AFAC, LDA )
*
*                    Form the inverse and compute the residual.
*
                    IF(MOD(N,2).EQ.0)THEN
                       CALL CLACPY( 'A', N+1, N/2, ARF, N+1, ARFINV, N+1 )
                    ELSE
                       CALL CLACPY( 'A', N, (N+1)/2, ARF, N, ARFINV, N )
                    END IF
*
                     SRNAMT = 'CPFTRI'
                     CALL CPFTRI( CFORM, UPLO, N, ARFINV , INFO )
*
                     SRNAMT = 'CTFTTR'
                     CALL CTFTTR( CFORM, UPLO, N, ARFINV, AINV, LDA, INFO )
*
*                    Check error code from CPFTRI.
*
                     IF( INFO.NE.0 ) CALL ALAERH( 'CPO', 'CPFTRI', INFO, 0, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT )
*
                     CALL CPOT03( UPLO, N, A, LDA, AINV, LDA, C_WORK_CPOT03, LDA, S_WORK_CPOT03, RCONDC, RESULT( 2 ) )
*
*                    Compute residual of the computed solution.
*
                     CALL CLACPY( 'Full', N, NRHS, B, LDA, C_WORK_CPOT02, LDA )                      CALL CPOT02( UPLO, N, NRHS, A, LDA, X, LDA, C_WORK_CPOT02, LDA, S_WORK_CPOT02, RESULT( 3 ) )
*
*                    Check solution from generated exact solution.
*
                     CALL CGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 4 ) )
                     NT = 4
*
*                    Print information about the tests that did not
*                    pass the threshold.
*
                     DO 60 K = 1, NT
                        IF( RESULT( K ).GE.THRESH ) THEN
                           IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALADHD( NOUT, 'CPF' )                            WRITE( NOUT, FMT = 9999 )'CPFSV ', UPLO, N, IIT, K, RESULT( K )
                           NFAIL = NFAIL + 1
                        END IF
   60                CONTINUE
                     NRUN = NRUN + NT
  100             CONTINUE
  110          CONTINUE
  120       CONTINUE
  980    CONTINUE
  130 CONTINUE
*
*     Print a summary of the results.
*
      CALL ALASVM( 'CPF', NOUT, NFAIL, NRUN, NERRS )
*
 9999 FORMAT( 1X, A6, ', UPLO=''', A1, ''', N =', I5, ', type ', I1,
     +      ', test(', I1, ')=', G12.5 )
*
      RETURN
*
*     End of CDRVRFP
*
      END
