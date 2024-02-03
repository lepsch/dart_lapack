      SUBROUTINE DDRVLS( DOTYPE, NM, MVAL, NN, NVAL, NNS, NSVAL, NNB, NBVAL, NXVAL, THRESH, TSTERR, A, COPYA, B, COPYB, C, S, COPYS, NOUT )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      LOGICAL            TSTERR
      int                NM, NN, NNB, NNS, NOUT
      DOUBLE PRECISION   THRESH
*     ..
*     .. Array Arguments ..
      LOGICAL            DOTYPE( * )
      int                MVAL( * ), NBVAL( * ), NSVAL( * ), NVAL( * ), NXVAL( * )       DOUBLE PRECISION   A( * ), B( * ), C( * ), COPYA( * ), COPYB( * ), COPYS( * ), S( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      int                NTESTS
      PARAMETER          ( NTESTS = 18 )
      int                SMLSIZ
      PARAMETER          ( SMLSIZ = 25 )
      DOUBLE PRECISION   ONE, TWO, ZERO
      PARAMETER          ( ONE = 1.0D0, TWO = 2.0D0, ZERO = 0.0D0 )
*     ..
*     .. Local Scalars ..
      CHARACTER          TRANS
      String             PATH;
      int                CRANK, I, IM, IMB, IN, INB, INFO, INS, IRANK, ISCALE, ITRAN, ITYPE, J, K, LDA, LDB, LDWORK, LWLSY, LWORK, M, MNMIN, N, NB, NCOLS, NERRS, NFAIL, NRHS, NROWS, NRUN, RANK, MB, MMAX, NMAX, NSMAX, LIWORK, LWORK_DGELS, LWORK_DGELST, LWORK_DGETSLS, LWORK_DGELSS, LWORK_DGELSY, LWORK_DGELSD
      DOUBLE PRECISION   EPS, NORMA, NORMB, RCOND
*     ..
*     .. Local Arrays ..
      int                ISEED( 4 ), ISEEDY( 4 ), IWQ( 1 )
      DOUBLE PRECISION   RESULT( NTESTS ), WQ( 1 )
*     ..
*     .. Allocatable Arrays ..
      DOUBLE PRECISION, ALLOCATABLE :: WORK (:)
      int    , ALLOCATABLE :: IWORK (:)
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DASUM, DLAMCH, DQRT12, DQRT14, DQRT17
      EXTERNAL           DASUM, DLAMCH, DQRT12, DQRT14, DQRT17
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALAERH, ALAHD, ALASVM, DAXPY, DERRLS, DGELS, DGELSD, DGELSS, DGELST, DGELSY, DGEMM, DGETSLS, DLACPY, DLARNV, DQRT13, DQRT15, DQRT16, DSCAL, XLAENV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, INT, MAX, MIN, SQRT
*     ..
*     .. Scalars in Common ..
      LOGICAL            LERR, OK
      String             SRNAMT;
      int                INFOT, IOUNIT
*     ..
*     .. Common blocks ..
      COMMON             / INFOC / INFOT, IOUNIT, OK, LERR
      COMMON             / SRNAMC / SRNAMT
*     ..
*     .. Data statements ..
      DATA               ISEEDY / 1988, 1989, 1990, 1991 /
*     ..
*     .. Executable Statements ..
*
*     Initialize constants and the random number seed.
*
      PATH( 1: 1 ) = 'Double precision'
      PATH( 2: 3 ) = 'LS'
      NRUN = 0
      NFAIL = 0
      NERRS = 0
      DO 10 I = 1, 4
         ISEED( I ) = ISEEDY( I )
   10 CONTINUE
      EPS = DLAMCH( 'Epsilon' )
*
*     Threshold for rank estimation
*
      RCOND = SQRT( EPS ) - ( SQRT( EPS )-EPS ) / 2
*
*     Test the error exits
*
      CALL XLAENV( 2, 2 )
      CALL XLAENV( 9, SMLSIZ )
      IF( TSTERR ) CALL DERRLS( PATH, NOUT )
*
*     Print the header if NM = 0 or NN = 0 and THRESH = 0.
*
      IF( ( NM.EQ.0 .OR. NN.EQ.0 ) .AND. THRESH.EQ.ZERO ) CALL ALAHD( NOUT, PATH )
      INFOT = 0
      CALL XLAENV( 2, 2 )
      CALL XLAENV( 9, SMLSIZ )
*
*     Compute maximal workspace needed for all routines
*
      NMAX = 0
      MMAX = 0
      NSMAX = 0
      DO I = 1, NM
         IF ( MVAL( I ).GT.MMAX ) THEN
            MMAX = MVAL( I )
         END IF
      ENDDO
      DO I = 1, NN
         IF ( NVAL( I ).GT.NMAX ) THEN
            NMAX = NVAL( I )
         END IF
      ENDDO
      DO I = 1, NNS
         IF ( NSVAL( I ).GT.NSMAX ) THEN
            NSMAX = NSVAL( I )
         END IF
      ENDDO
      M = MMAX
      N = NMAX
      NRHS = NSMAX
      MNMIN = MAX( MIN( M, N ), 1 )
*
*     Compute workspace needed for routines
*     DQRT14, DQRT17 (two side cases), DQRT15 and DQRT12
*
      LWORK = MAX( 1, ( M+N )*NRHS, ( N+NRHS )*( M+2 ), ( M+NRHS )*( N+2 ), MAX( M+MNMIN, NRHS*MNMIN,2*N+M ), MAX( M*N+4*MNMIN+MAX(M,N), M*N+2*MNMIN+4*N ) )
      LIWORK = 1
*
*     Iterate through all test cases and compute necessary workspace
*     sizes for ?GELS, ?GELST, ?GETSLS, ?GELSY, ?GELSS and ?GELSD
*     routines.
*
      DO IM = 1, NM
         M = MVAL( IM )
         LDA = MAX( 1, M )
         DO IN = 1, NN
            N = NVAL( IN )
            MNMIN = MAX(MIN( M, N ),1)
            LDB = MAX( 1, M, N )
            DO INS = 1, NNS
               NRHS = NSVAL( INS )
               DO IRANK = 1, 2
                  DO ISCALE = 1, 3
                     ITYPE = ( IRANK-1 )*3 + ISCALE
                     IF( DOTYPE( ITYPE ) ) THEN
                        IF( IRANK.EQ.1 ) THEN
                           DO ITRAN = 1, 2
                              IF( ITRAN.EQ.1 ) THEN
                                 TRANS = 'N'
                              ELSE
                                 TRANS = 'T'
                              END IF
*
*                             Compute workspace needed for DGELS
                              CALL DGELS( TRANS, M, N, NRHS, A, LDA, B, LDB, WQ, -1, INFO )
                              LWORK_DGELS = INT ( WQ ( 1 ) )
*                             Compute workspace needed for DGELST
                              CALL DGELST( TRANS, M, N, NRHS, A, LDA, B, LDB, WQ, -1, INFO )
                              LWORK_DGELST = INT ( WQ ( 1 ) )
*                             Compute workspace needed for DGETSLS
                              CALL DGETSLS( TRANS, M, N, NRHS, A, LDA, B, LDB, WQ, -1, INFO )
                              LWORK_DGETSLS = INT( WQ ( 1 ) )
                           ENDDO
                        END IF
*                       Compute workspace needed for DGELSY
                        CALL DGELSY( M, N, NRHS, A, LDA, B, LDB, IWQ, RCOND, CRANK, WQ, -1, INFO )
                        LWORK_DGELSY = INT( WQ ( 1 ) )
*                       Compute workspace needed for DGELSS
                        CALL DGELSS( M, N, NRHS, A, LDA, B, LDB, S, RCOND, CRANK, WQ, -1 , INFO )
                        LWORK_DGELSS = INT( WQ ( 1 ) )
*                       Compute workspace needed for DGELSD
                        CALL DGELSD( M, N, NRHS, A, LDA, B, LDB, S, RCOND, CRANK, WQ, -1, IWQ, INFO )
                        LWORK_DGELSD = INT( WQ ( 1 ) )
*                       Compute LIWORK workspace needed for DGELSY and DGELSD
                        LIWORK = MAX( LIWORK, N, IWQ( 1 ) )
*                       Compute LWORK workspace needed for all functions
                        LWORK = MAX( LWORK, LWORK_DGELS, LWORK_DGELST, LWORK_DGETSLS, LWORK_DGELSY, LWORK_DGELSS, LWORK_DGELSD )
                     END IF
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
*
      LWLSY = LWORK
*
      ALLOCATE( WORK( LWORK ) )
      ALLOCATE( IWORK( LIWORK ) )
*
      DO 150 IM = 1, NM
         M = MVAL( IM )
         LDA = MAX( 1, M )
*
         DO 140 IN = 1, NN
            N = NVAL( IN )
            MNMIN = MAX(MIN( M, N ),1)
            LDB = MAX( 1, M, N )
            MB = (MNMIN+1)
*
            DO 130 INS = 1, NNS
               NRHS = NSVAL( INS )
*
               DO 120 IRANK = 1, 2
                  DO 110 ISCALE = 1, 3
                     ITYPE = ( IRANK-1 )*3 + ISCALE
                     IF( .NOT.DOTYPE( ITYPE ) ) GO TO 110
*                 =====================================================
*                       Begin test DGELS
*                 =====================================================
                     IF( IRANK.EQ.1 ) THEN
*
*                       Generate a matrix of scaling type ISCALE
*
                        CALL DQRT13( ISCALE, M, N, COPYA, LDA, NORMA, ISEED )
*
*                       Loop for testing different block sizes.
*
                        DO INB = 1, NNB
                           NB = NBVAL( INB )
                           CALL XLAENV( 1, NB )
                           CALL XLAENV( 3, NXVAL( INB ) )
*
*                          Loop for testing non-transposed and transposed.
*
                           DO ITRAN = 1, 2
                              IF( ITRAN.EQ.1 ) THEN
                                 TRANS = 'N'
                                 NROWS = M
                                 NCOLS = N
                              ELSE
                                 TRANS = 'T'
                                 NROWS = N
                                 NCOLS = M
                              END IF
                              LDWORK = MAX( 1, NCOLS )
*
*                             Set up a consistent rhs
*
                              IF( NCOLS.GT.0 ) THEN
                                 CALL DLARNV( 2, ISEED, NCOLS*NRHS, WORK )                                  CALL DSCAL( NCOLS*NRHS, ONE / DBLE( NCOLS ), WORK, 1 )
                              END IF
                              CALL DGEMM( TRANS, 'No transpose', NROWS, NRHS, NCOLS, ONE, COPYA, LDA, WORK, LDWORK, ZERO, B, LDB )
                              CALL DLACPY( 'Full', NROWS, NRHS, B, LDB, COPYB, LDB )
*
*                             Solve LS or overdetermined system
*
                              IF( M.GT.0 .AND. N.GT.0 ) THEN
                                 CALL DLACPY( 'Full', M, N, COPYA, LDA, A, LDA )                                  CALL DLACPY( 'Full', NROWS, NRHS, COPYB, LDB, B, LDB )
                              END IF
                              SRNAMT = 'DGELS '
                              CALL DGELS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO )                               IF( INFO.NE.0 ) CALL ALAERH( PATH, 'DGELS ', INFO, 0, TRANS, M, N, NRHS, -1, NB, ITYPE, NFAIL, NERRS, NOUT )
*
*                             Test 1: Check correctness of results
*                             for DGELS, compute the residual:
*                             RESID = norm(B - A*X) /
*                             / ( max(m,n) * norm(A) * norm(X) * EPS )
*
                              IF( NROWS.GT.0 .AND. NRHS.GT.0 ) CALL DLACPY( 'Full', NROWS, NRHS, COPYB, LDB, C, LDB )                               CALL DQRT16( TRANS, M, N, NRHS, COPYA, LDA, B, LDB, C, LDB, WORK, RESULT( 1 ) )
*
*                             Test 2: Check correctness of results
*                             for DGELS.
*
                              IF( ( ITRAN.EQ.1 .AND. M.GE.N ) .OR. ( ITRAN.EQ.2 .AND. M.LT.N ) ) THEN
*
*                                Solving LS system, compute:
*                                r = norm((B- A*X)**T * A) /
*                                 / (norm(A)*norm(B)*max(M,N,NRHS)*EPS)
*
                                 RESULT( 2 ) = DQRT17( TRANS, 1, M, N, NRHS, COPYA, LDA, B, LDB, COPYB, LDB, C, WORK, LWORK )
                              ELSE
*
*                                Solving overdetermined system
*
                                 RESULT( 2 ) = DQRT14( TRANS, M, N, NRHS, COPYA, LDA, B, LDB, WORK, LWORK )
                              END IF
*
*                             Print information about the tests that
*                             did not pass the threshold.
*
                              DO K = 1, 2
                                 IF( RESULT( K ).GE.THRESH ) THEN
                                    IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALAHD( NOUT, PATH )                                     WRITE( NOUT, FMT = 9999 ) TRANS, M, N, NRHS, NB, ITYPE, K, RESULT( K )
                                    NFAIL = NFAIL + 1
                                 END IF
                              END DO
                              NRUN = NRUN + 2
                           END DO
                        END DO
                     END IF
*                 =====================================================
*                       End test DGELS
*                 =====================================================
*                 =====================================================
*                       Begin test DGELST
*                 =====================================================
                     IF( IRANK.EQ.1 ) THEN
*
*                       Generate a matrix of scaling type ISCALE
*
                        CALL DQRT13( ISCALE, M, N, COPYA, LDA, NORMA, ISEED )
*
*                       Loop for testing different block sizes.
*
                        DO INB = 1, NNB
                           NB = NBVAL( INB )
                           CALL XLAENV( 1, NB )
*
*                          Loop for testing non-transposed and transposed.
*
                           DO ITRAN = 1, 2
                              IF( ITRAN.EQ.1 ) THEN
                                 TRANS = 'N'
                                 NROWS = M
                                 NCOLS = N
                              ELSE
                                 TRANS = 'T'
                                 NROWS = N
                                 NCOLS = M
                              END IF
                              LDWORK = MAX( 1, NCOLS )
*
*                             Set up a consistent rhs
*
                              IF( NCOLS.GT.0 ) THEN
                                 CALL DLARNV( 2, ISEED, NCOLS*NRHS, WORK )                                  CALL DSCAL( NCOLS*NRHS, ONE / DBLE( NCOLS ), WORK, 1 )
                              END IF
                              CALL DGEMM( TRANS, 'No transpose', NROWS, NRHS, NCOLS, ONE, COPYA, LDA, WORK, LDWORK, ZERO, B, LDB )
                              CALL DLACPY( 'Full', NROWS, NRHS, B, LDB, COPYB, LDB )
*
*                             Solve LS or overdetermined system
*
                              IF( M.GT.0 .AND. N.GT.0 ) THEN
                                 CALL DLACPY( 'Full', M, N, COPYA, LDA, A, LDA )                                  CALL DLACPY( 'Full', NROWS, NRHS, COPYB, LDB, B, LDB )
                              END IF
                              SRNAMT = 'DGELST'
                              CALL DGELST( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO )                               IF( INFO.NE.0 ) CALL ALAERH( PATH, 'DGELST', INFO, 0, TRANS, M, N, NRHS, -1, NB, ITYPE, NFAIL, NERRS, NOUT )
*
*                             Test 3: Check correctness of results
*                             for DGELST, compute the residual:
*                             RESID = norm(B - A*X) /
*                             / ( max(m,n) * norm(A) * norm(X) * EPS )
*
                              IF( NROWS.GT.0 .AND. NRHS.GT.0 ) CALL DLACPY( 'Full', NROWS, NRHS, COPYB, LDB, C, LDB )                               CALL DQRT16( TRANS, M, N, NRHS, COPYA, LDA, B, LDB, C, LDB, WORK, RESULT( 3 ) )
*
*                             Test 4: Check correctness of results
*                             for DGELST.
*
                              IF( ( ITRAN.EQ.1 .AND. M.GE.N ) .OR. ( ITRAN.EQ.2 .AND. M.LT.N ) ) THEN
*
*                                Solving LS system, compute:
*                                r = norm((B- A*X)**T * A) /
*                                 / (norm(A)*norm(B)*max(M,N,NRHS)*EPS)
*
                                 RESULT( 4 ) = DQRT17( TRANS, 1, M, N, NRHS, COPYA, LDA, B, LDB, COPYB, LDB, C, WORK, LWORK )
                              ELSE
*
*                                Solving overdetermined system
*
                                 RESULT( 4 ) = DQRT14( TRANS, M, N, NRHS, COPYA, LDA, B, LDB, WORK, LWORK )
                              END IF
*
*                             Print information about the tests that
*                             did not pass the threshold.
*
                              DO K = 3, 4
                                 IF( RESULT( K ).GE.THRESH ) THEN
                                    IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALAHD( NOUT, PATH )                                     WRITE( NOUT, FMT = 9999 ) TRANS, M, N, NRHS, NB, ITYPE, K, RESULT( K )
                                    NFAIL = NFAIL + 1
                                 END IF
                              END DO
                              NRUN = NRUN + 2
                           END DO
                        END DO
                     END IF
*                 =====================================================
*                       End test DGELST
*                 =====================================================
*                 =====================================================
*                       Begin test DGETSLS
*                 =====================================================
                     IF( IRANK.EQ.1 ) THEN
*
*                       Generate a matrix of scaling type ISCALE
*
                        CALL DQRT13( ISCALE, M, N, COPYA, LDA, NORMA, ISEED )
*
*                       Loop for testing different block sizes MB.
*
                        DO IMB = 1, NNB
                           MB = NBVAL( IMB )
                           CALL XLAENV( 1, MB )
*
*                          Loop for testing different block sizes NB.
*
                           DO INB = 1, NNB
                              NB = NBVAL( INB )
                              CALL XLAENV( 2, NB )
*
*                             Loop for testing non-transposed
*                             and transposed.
*
                              DO ITRAN = 1, 2
                                 IF( ITRAN.EQ.1 ) THEN
                                    TRANS = 'N'
                                    NROWS = M
                                    NCOLS = N
                                 ELSE
                                    TRANS = 'T'
                                    NROWS = N
                                    NCOLS = M
                                 END IF
                                 LDWORK = MAX( 1, NCOLS )
*
*                                Set up a consistent rhs
*
                                 IF( NCOLS.GT.0 ) THEN
                                    CALL DLARNV( 2, ISEED, NCOLS*NRHS, WORK )                                     CALL DSCAL( NCOLS*NRHS, ONE / DBLE( NCOLS ), WORK, 1 )
                                 END IF
                                 CALL DGEMM( TRANS, 'No transpose', NROWS, NRHS, NCOLS, ONE, COPYA, LDA, WORK, LDWORK, ZERO, B, LDB )
                                 CALL DLACPY( 'Full', NROWS, NRHS, B, LDB, COPYB, LDB )
*
*                                Solve LS or overdetermined system
*
                                 IF( M.GT.0 .AND. N.GT.0 ) THEN
                                    CALL DLACPY( 'Full', M, N, COPYA, LDA, A, LDA )                                     CALL DLACPY( 'Full', NROWS, NRHS, COPYB, LDB, B, LDB )
                                 END IF
                                 SRNAMT = 'DGETSLS'
                                 CALL DGETSLS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO )                                  IF( INFO.NE.0 ) CALL ALAERH( PATH, 'DGETSLS', INFO, 0, TRANS, M, N, NRHS, -1, NB, ITYPE, NFAIL, NERRS, NOUT )
*
*                             Test 5: Check correctness of results
*                             for DGETSLS, compute the residual:
*                             RESID = norm(B - A*X) /
*                             / ( max(m,n) * norm(A) * norm(X) * EPS )
*
                                 IF( NROWS.GT.0 .AND. NRHS.GT.0 ) CALL DLACPY( 'Full', NROWS, NRHS, COPYB, LDB, C, LDB )                                  CALL DQRT16( TRANS, M, N, NRHS, COPYA, LDA, B, LDB, C, LDB, WORK, RESULT( 5 ) )
*
*                             Test 6: Check correctness of results
*                             for DGETSLS.
*
                                 IF( ( ITRAN.EQ.1 .AND. M.GE.N ) .OR. ( ITRAN.EQ.2 .AND. M.LT.N ) ) THEN
*
*                                   Solving LS system, compute:
*                                   r = norm((B- A*X)**T * A) /
*                                 / (norm(A)*norm(B)*max(M,N,NRHS)*EPS)
*
                                    RESULT( 6 ) = DQRT17( TRANS, 1, M, N, NRHS, COPYA, LDA, B, LDB, COPYB, LDB, C, WORK, LWORK )
                                 ELSE
*
*                                   Solving overdetermined system
*
                                    RESULT( 6 ) = DQRT14( TRANS, M, N, NRHS, COPYA, LDA, B, LDB, WORK, LWORK )
                                 END IF
*
*                                Print information about the tests that
*                                did not pass the threshold.
*
                                 DO K = 5, 6
                                    IF( RESULT( K ).GE.THRESH ) THEN
                                       IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALAHD( NOUT, PATH )                                        WRITE( NOUT, FMT = 9997 ) TRANS, M, N, NRHS, MB, NB, ITYPE, K, RESULT( K )
                                       NFAIL = NFAIL + 1
                                    END IF
                                 END DO
                                 NRUN = NRUN + 2
                              END DO
                           END DO
                        END DO
                     END IF
*                 =====================================================
*                       End test DGETSLS
*                 =====================================================
*
*                    Generate a matrix of scaling type ISCALE and rank
*                    type IRANK.
*
                     CALL DQRT15( ISCALE, IRANK, M, N, NRHS, COPYA, LDA, COPYB, LDB, COPYS, RANK, NORMA, NORMB, ISEED, WORK, LWORK )
*
*                    workspace used: MAX(M+MIN(M,N),NRHS*MIN(M,N),2*N+M)
*
                     LDWORK = MAX( 1, M )
*
*                    Loop for testing different block sizes.
*
                     DO 100 INB = 1, NNB
                        NB = NBVAL( INB )
                        CALL XLAENV( 1, NB )
                        CALL XLAENV( 3, NXVAL( INB ) )
*
*                       Test DGELSY
*
*                       DGELSY:  Compute the minimum-norm solution X
*                       to min( norm( A * X - B ) )
*                       using the rank-revealing orthogonal
*                       factorization.
*
*                       Initialize vector IWORK.
*
                        DO 70 J = 1, N
                           IWORK( J ) = 0
   70                   CONTINUE
*
                        CALL DLACPY( 'Full', M, N, COPYA, LDA, A, LDA )
                        CALL DLACPY( 'Full', M, NRHS, COPYB, LDB, B, LDB )
*
                        SRNAMT = 'DGELSY'
                        CALL DGELSY( M, N, NRHS, A, LDA, B, LDB, IWORK, RCOND, CRANK, WORK, LWLSY, INFO )                         IF( INFO.NE.0 ) CALL ALAERH( PATH, 'DGELSY', INFO, 0, ' ', M, N, NRHS, -1, NB, ITYPE, NFAIL, NERRS, NOUT )
*
*                       Test 7:  Compute relative error in svd
*                                workspace: M*N + 4*MIN(M,N) + MAX(M,N)
*
                        RESULT( 7 ) = DQRT12( CRANK, CRANK, A, LDA, COPYS, WORK, LWORK )
*
*                       Test 8:  Compute error in solution
*                                workspace:  M*NRHS + M
*
                        CALL DLACPY( 'Full', M, NRHS, COPYB, LDB, WORK, LDWORK )                         CALL DQRT16( 'No transpose', M, N, NRHS, COPYA, LDA, B, LDB, WORK, LDWORK, WORK( M*NRHS+1 ), RESULT( 8 ) )
*
*                       Test 9:  Check norm of r'*A
*                                workspace: NRHS*(M+N)
*
                        RESULT( 9 ) = ZERO
                        IF( M.GT.CRANK ) RESULT( 9 ) = DQRT17( 'No transpose', 1, M, N, NRHS, COPYA, LDA, B, LDB, COPYB, LDB, C, WORK, LWORK )
*
*                       Test 10:  Check if x is in the rowspace of A
*                                workspace: (M+NRHS)*(N+2)
*
                        RESULT( 10 ) = ZERO
*
                        IF( N.GT.CRANK ) RESULT( 10 ) = DQRT14( 'No transpose', M, N, NRHS, COPYA, LDA, B, LDB, WORK, LWORK )
*
*                       Test DGELSS
*
*                       DGELSS:  Compute the minimum-norm solution X
*                       to min( norm( A * X - B ) )
*                       using the SVD.
*
                        CALL DLACPY( 'Full', M, N, COPYA, LDA, A, LDA )
                        CALL DLACPY( 'Full', M, NRHS, COPYB, LDB, B, LDB )
                        SRNAMT = 'DGELSS'
                        CALL DGELSS( M, N, NRHS, A, LDA, B, LDB, S, RCOND, CRANK, WORK, LWORK, INFO )                         IF( INFO.NE.0 ) CALL ALAERH( PATH, 'DGELSS', INFO, 0, ' ', M, N, NRHS, -1, NB, ITYPE, NFAIL, NERRS, NOUT )
*
*                       workspace used: 3*min(m,n) +
*                                       max(2*min(m,n),nrhs,max(m,n))
*
*                       Test 11:  Compute relative error in svd
*
                        IF( RANK.GT.0 ) THEN
                           CALL DAXPY( MNMIN, -ONE, COPYS, 1, S, 1 )
                           RESULT( 11 ) = DASUM( MNMIN, S, 1 ) / DASUM( MNMIN, COPYS, 1 ) / ( EPS*DBLE( MNMIN ) )
                        ELSE
                           RESULT( 11 ) = ZERO
                        END IF
*
*                       Test 12:  Compute error in solution
*
                        CALL DLACPY( 'Full', M, NRHS, COPYB, LDB, WORK, LDWORK )                         CALL DQRT16( 'No transpose', M, N, NRHS, COPYA, LDA, B, LDB, WORK, LDWORK, WORK( M*NRHS+1 ), RESULT( 12 ) )
*
*                       Test 13:  Check norm of r'*A
*
                        RESULT( 13 ) = ZERO
                        IF( M.GT.CRANK ) RESULT( 13 ) = DQRT17( 'No transpose', 1, M, N, NRHS, COPYA, LDA, B, LDB, COPYB, LDB, C, WORK, LWORK )
*
*                       Test 14:  Check if x is in the rowspace of A
*
                        RESULT( 14 ) = ZERO
                        IF( N.GT.CRANK ) RESULT( 14 ) = DQRT14( 'No transpose', M, N, NRHS, COPYA, LDA, B, LDB, WORK, LWORK )
*
*                       Test DGELSD
*
*                       DGELSD:  Compute the minimum-norm solution X
*                       to min( norm( A * X - B ) ) using a
*                       divide and conquer SVD.
*
*                       Initialize vector IWORK.
*
                        DO 80 J = 1, N
                           IWORK( J ) = 0
   80                   CONTINUE
*
                        CALL DLACPY( 'Full', M, N, COPYA, LDA, A, LDA )
                        CALL DLACPY( 'Full', M, NRHS, COPYB, LDB, B, LDB )
*
                        SRNAMT = 'DGELSD'
                        CALL DGELSD( M, N, NRHS, A, LDA, B, LDB, S, RCOND, CRANK, WORK, LWORK, IWORK, INFO )                         IF( INFO.NE.0 ) CALL ALAERH( PATH, 'DGELSD', INFO, 0, ' ', M, N, NRHS, -1, NB, ITYPE, NFAIL, NERRS, NOUT )
*
*                       Test 15:  Compute relative error in svd
*
                        IF( RANK.GT.0 ) THEN
                           CALL DAXPY( MNMIN, -ONE, COPYS, 1, S, 1 )
                           RESULT( 15 ) = DASUM( MNMIN, S, 1 ) / DASUM( MNMIN, COPYS, 1 ) / ( EPS*DBLE( MNMIN ) )
                        ELSE
                           RESULT( 15 ) = ZERO
                        END IF
*
*                       Test 16:  Compute error in solution
*
                        CALL DLACPY( 'Full', M, NRHS, COPYB, LDB, WORK, LDWORK )                         CALL DQRT16( 'No transpose', M, N, NRHS, COPYA, LDA, B, LDB, WORK, LDWORK, WORK( M*NRHS+1 ), RESULT( 16 ) )
*
*                       Test 17:  Check norm of r'*A
*
                        RESULT( 17 ) = ZERO
                        IF( M.GT.CRANK ) RESULT( 17 ) = DQRT17( 'No transpose', 1, M, N, NRHS, COPYA, LDA, B, LDB, COPYB, LDB, C, WORK, LWORK )
*
*                       Test 18:  Check if x is in the rowspace of A
*
                        RESULT( 18 ) = ZERO
                        IF( N.GT.CRANK ) RESULT( 18 ) = DQRT14( 'No transpose', M, N, NRHS, COPYA, LDA, B, LDB, WORK, LWORK )
*
*                       Print information about the tests that did not
*                       pass the threshold.
*
                        DO 90 K = 7, 18
                           IF( RESULT( K ).GE.THRESH ) THEN
                              IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALAHD( NOUT, PATH )                               WRITE( NOUT, FMT = 9998 )M, N, NRHS, NB, ITYPE, K, RESULT( K )
                              NFAIL = NFAIL + 1
                           END IF
   90                   CONTINUE
                        NRUN = NRUN + 12
*
  100                CONTINUE






  110             CONTINUE
  120          CONTINUE
  130       CONTINUE
  140    CONTINUE
  150 CONTINUE
*
*     Print a summary of the results.
*
      CALL ALASVM( PATH, NOUT, NFAIL, NRUN, NERRS )
*
 9999 FORMAT( ' TRANS=''', A1, ''', M=', I5, ', N=', I5, ', NRHS=', I4,
     $      ', NB=', I4, ', type', I2, ', test(', I2, ')=', G12.5 )
 9998 FORMAT( ' M=', I5, ', N=', I5, ', NRHS=', I4, ', NB=', I4,
     $      ', type', I2, ', test(', I2, ')=', G12.5 )
 9997 FORMAT( ' TRANS=''', A1,' M=', I5, ', N=', I5, ', NRHS=', I4,
     $      ', MB=', I4,', NB=', I4,', type', I2,
     $      ', test(', I2, ')=', G12.5 )
*
      DEALLOCATE( WORK )
      DEALLOCATE( IWORK )
      RETURN
*
*     End of DDRVLS
*
      END
