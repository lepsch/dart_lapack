      SUBROUTINE STPLQT( M, N, L, MB, A, LDA, B, LDB, T, LDT, WORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int     INFO, LDA, LDB, LDT, N, M, L, MB;
*     ..
*     .. Array Arguments ..
      REAL A( LDA, * ), B( LDB, * ), T( LDT, * ), WORK( * )
*     ..
*
* =====================================================================
*
*     ..
*     .. Local Scalars ..
      int        I, IB, LB, NB, IINFO;
*     ..
*     .. External Subroutines ..
      // EXTERNAL STPLQT2, STPRFB, XERBLA
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( L.LT.0 .OR. (L.GT.MIN(M,N) .AND. MIN(M,N).GE.0)) THEN
         INFO = -3
      ELSE IF( MB.LT.1 .OR. (MB.GT.M .AND. M.GT.0)) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -6
      ELSE IF( LDB.LT.MAX( 1, M ) ) THEN
         INFO = -8
      ELSE IF( LDT.LT.MB ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'STPLQT', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 ) RETURN
*
      DO I = 1, M, MB
*
*     Compute the QR factorization of the current block
*
         IB = MIN( M-I+1, MB )
         NB = MIN( N-L+I+IB-1, N )
         IF( I.GE.L ) THEN
            LB = 0
         ELSE
            LB = NB-N+L-I+1
         END IF
*
         CALL STPLQT2( IB, NB, LB, A(I,I), LDA, B( I, 1 ), LDB, T(1, I ), LDT, IINFO )
*
*     Update by applying H**T to B(I+IB:M,:) from the right
*
         IF( I+IB.LE.M ) THEN
            CALL STPRFB( 'R', 'N', 'F', 'R', M-I-IB+1, NB, IB, LB, B( I, 1 ), LDB, T( 1, I ), LDT, A( I+IB, I ), LDA, B( I+IB, 1 ), LDB, WORK, M-I-IB+1)
         END IF
      END DO
      RETURN
*
*     End of STPLQT
*
      END
