      SUBROUTINE ZHPGV( ITYPE, JOBZ, UPLO, N, AP, BP, W, Z, LDZ, WORK, RWORK, INFO )
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          JOBZ, UPLO
      INTEGER            INFO, ITYPE, LDZ, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   RWORK( * ), W( * )
      COMPLEX*16         AP( * ), BP( * ), WORK( * ), Z( LDZ, * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            UPPER, WANTZ
      CHARACTER          TRANS
      INTEGER            J, NEIG
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZHPEV, ZHPGST, ZPPTRF, ZTPMV, ZTPSV
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      WANTZ = LSAME( JOBZ, 'V' )
      UPPER = LSAME( UPLO, 'U' )
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
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZHPGV ', -INFO )
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
      CALL ZHPEV( JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, RWORK, INFO )
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
      RETURN
*
*     End of ZHPGV
*
      END
