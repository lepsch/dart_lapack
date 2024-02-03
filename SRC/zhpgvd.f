      SUBROUTINE ZHPGVD( ITYPE, JOBZ, UPLO, N, AP, BP, W, Z, LDZ, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO )
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             JOBZ, UPLO;
      int                INFO, ITYPE, LDZ, LIWORK, LRWORK, LWORK, N;
*     ..
*     .. Array Arguments ..
      int                IWORK( * );
      double             RWORK( * ), W( * );
      COMPLEX*16         AP( * ), BP( * ), WORK( * ), Z( LDZ, * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      bool               LQUERY, UPPER, WANTZ;
      String             TRANS;
      int                J, LIWMIN, LRWMIN, LWMIN, NEIG;
*     ..
*     .. External Functions ..
      bool               LSAME;
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZHPEVD, ZHPGST, ZPPTRF, ZTPMV, ZTPSV
*     ..
*     .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      WANTZ = LSAME( JOBZ, 'V' )
      UPPER = LSAME( UPLO, 'U' )
      LQUERY = ( LWORK.EQ.-1 .OR. LRWORK.EQ.-1 .OR. LIWORK.EQ.-1 )
*
      INFO = 0
      IF( ITYPE.LT.1 .OR. ITYPE.GT.3 ) THEN
         INFO = -1
      ELSE IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
         INFO = -2
      ELSE IF( .NOT.( UPPER .OR. LSAME( UPLO, 'L' ) ) ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.N ) ) THEN
         INFO = -9
      END IF
*
      IF( INFO.EQ.0 ) THEN
         IF( N.LE.1 ) THEN
            LWMIN = 1
            LIWMIN = 1
            LRWMIN = 1
         ELSE
            IF( WANTZ ) THEN
               LWMIN = 2*N
               LRWMIN = 1 + 5*N + 2*N**2
               LIWMIN = 3 + 5*N
            ELSE
               LWMIN = N
               LRWMIN = N
               LIWMIN = 1
            END IF
         END IF
*
         WORK( 1 ) = LWMIN
         RWORK( 1 ) = LRWMIN
         IWORK( 1 ) = LIWMIN
         IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
            INFO = -11
         ELSE IF( LRWORK.LT.LRWMIN .AND. .NOT.LQUERY ) THEN
            INFO = -13
         ELSE IF( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) THEN
            INFO = -15
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZHPGVD', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 ) RETURN
*
*     Form a Cholesky factorization of B.
*
      CALL ZPPTRF( UPLO, N, BP, INFO )
      IF( INFO.NE.0 ) THEN
         INFO = N + INFO
         RETURN
      END IF
*
*     Transform problem to standard eigenvalue problem and solve.
*
      CALL ZHPGST( ITYPE, UPLO, N, AP, BP, INFO )
      CALL ZHPEVD( JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO )
      LWMIN = INT( MAX( DBLE( LWMIN ), DBLE( WORK( 1 ) ) ) )
      LRWMIN = INT( MAX( DBLE( LRWMIN ), DBLE( RWORK( 1 ) ) ) )
      LIWMIN = INT( MAX( DBLE( LIWMIN ), DBLE( IWORK( 1 ) ) ) )
*
      IF( WANTZ ) THEN
*
*        Backtransform eigenvectors to the original problem.
*
         NEIG = N
         IF( INFO.GT.0 ) NEIG = INFO - 1
         IF( ITYPE.EQ.1 .OR. ITYPE.EQ.2 ) THEN
*
*           For A*x=(lambda)*B*x and A*B*x=(lambda)*x;
*           backtransform eigenvectors: x = inv(L)**H *y or inv(U)*y
*
            IF( UPPER ) THEN
               TRANS = 'N'
            ELSE
               TRANS = 'C'
            END IF
*
            DO 10 J = 1, NEIG
               CALL ZTPSV( UPLO, TRANS, 'Non-unit', N, BP, Z( 1, J ), 1 )
   10       CONTINUE
*
         ELSE IF( ITYPE.EQ.3 ) THEN
*
*           For B*A*x=(lambda)*x;
*           backtransform eigenvectors: x = L*y or U**H *y
*
            IF( UPPER ) THEN
               TRANS = 'C'
            ELSE
               TRANS = 'N'
            END IF
*
            DO 20 J = 1, NEIG
               CALL ZTPMV( UPLO, TRANS, 'Non-unit', N, BP, Z( 1, J ), 1 )
   20       CONTINUE
         END IF
      END IF
*
      WORK( 1 ) = LWMIN
      RWORK( 1 ) = LRWMIN
      IWORK( 1 ) = LIWMIN
      RETURN
*
*     End of ZHPGVD
*
      END
