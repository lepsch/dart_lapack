      SUBROUTINE SSYGVD( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, W, WORK, LWORK, IWORK, LIWORK, INFO )
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      String             JOBZ, UPLO;
      int                INFO, ITYPE, LDA, LDB, LIWORK, LWORK, N;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      REAL               A( LDA, * ), B( LDB, * ), W( * ), WORK( * )
      // ..
*
*  =====================================================================
*
      // .. Parameters ..
      REAL               ONE
      PARAMETER          ( ONE = 1.0E+0 )
      // ..
      // .. Local Scalars ..
      bool               LQUERY, UPPER, WANTZ;
      String             TRANS;
      int                LIOPT, LIWMIN, LOPT, LWMIN;
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               SROUNDUP_LWORK
      // EXTERNAL LSAME, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL SPOTRF, SSYEVD, SSYGST, STRMM, STRSM, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, REAL
      // ..
      // .. Executable Statements ..
*
      // Test the input parameters.
*
      WANTZ = LSAME( JOBZ, 'V' )
      UPPER = LSAME( UPLO, 'U' )
      LQUERY = ( LWORK.EQ.-1 .OR. LIWORK.EQ.-1 )
*
      INFO = 0
      IF( N.LE.1 ) THEN
         LIWMIN = 1
         LWMIN = 1
      ELSE IF( WANTZ ) THEN
         LIWMIN = 3 + 5*N
         LWMIN = 1 + 6*N + 2*N**2
      ELSE
         LIWMIN = 1
         LWMIN = 2*N + 1
      END IF
      LOPT = LWMIN
      LIOPT = LIWMIN
      IF( ITYPE.LT.1 .OR. ITYPE.GT.3 ) THEN
         INFO = -1
      ELSE IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
         INFO = -2
      ELSE IF( .NOT.( UPPER .OR. LSAME( UPLO, 'L' ) ) ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
*
      IF( INFO.EQ.0 ) THEN
         WORK( 1 ) = SROUNDUP_LWORK(LOPT)
         IWORK( 1 ) = LIOPT
*
         IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
            INFO = -11
         ELSE IF( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) THEN
            INFO = -13
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SSYGVD', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
      // Quick return if possible
*
      IF( N.EQ.0 ) RETURN
*
      // Form a Cholesky factorization of B.
*
      CALL SPOTRF( UPLO, N, B, LDB, INFO )
      IF( INFO.NE.0 ) THEN
         INFO = N + INFO
         RETURN
      END IF
*
      // Transform problem to standard eigenvalue problem and solve.
*
      CALL SSYGST( ITYPE, UPLO, N, A, LDA, B, LDB, INFO )
      CALL SSYEVD( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, IWORK, LIWORK, INFO )
      LOPT = INT( MAX( REAL( LOPT ), REAL( WORK( 1 ) ) ) )
      LIOPT = INT( MAX( REAL( LIOPT ), REAL( IWORK( 1 ) ) ) )
*
      IF( WANTZ .AND. INFO.EQ.0 ) THEN
*
         // Backtransform eigenvectors to the original problem.
*
         IF( ITYPE.EQ.1 .OR. ITYPE.EQ.2 ) THEN
*
            // For A*x=(lambda)*B*x and A*B*x=(lambda)*x;
            // backtransform eigenvectors: x = inv(L)**T*y or inv(U)*y
*
            IF( UPPER ) THEN
               TRANS = 'N'
            ELSE
               TRANS = 'T'
            END IF
*
            CALL STRSM( 'Left', UPLO, TRANS, 'Non-unit', N, N, ONE, B, LDB, A, LDA )
*
         ELSE IF( ITYPE.EQ.3 ) THEN
*
            // For B*A*x=(lambda)*x;
            // backtransform eigenvectors: x = L*y or U**T*y
*
            IF( UPPER ) THEN
               TRANS = 'T'
            ELSE
               TRANS = 'N'
            END IF
*
            CALL STRMM( 'Left', UPLO, TRANS, 'Non-unit', N, N, ONE, B, LDB, A, LDA )
         END IF
      END IF
*
      WORK( 1 ) = SROUNDUP_LWORK(LOPT)
      IWORK( 1 ) = LIOPT
*
      RETURN
*
      // End of SSYGVD
*
      END
