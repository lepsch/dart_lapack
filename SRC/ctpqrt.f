      SUBROUTINE CTPQRT( M, N, L, NB, A, LDA, B, LDB, T, LDT, WORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      int     INFO, LDA, LDB, LDT, N, M, L, NB;
      // ..
      // .. Array Arguments ..
      COMPLEX A( LDA, * ), B( LDB, * ), T( LDT, * ), WORK( * )
      // ..
*
* =====================================================================
*
      // ..
      // .. Local Scalars ..
      int        I, IB, LB, MB, IINFO;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CTPQRT2, CTPRFB, XERBLA
      // ..
      // .. Executable Statements ..
*
      // Test the input arguments
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( L.LT.0 .OR. (L.GT.MIN(M,N) .AND. MIN(M,N).GE.0)) THEN
         INFO = -3
      ELSE IF( NB.LT.1 .OR. (NB.GT.N .AND. N.GT.0)) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDB.LT.MAX( 1, M ) ) THEN
         INFO = -8
      ELSE IF( LDT.LT.NB ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CTPQRT', -INFO )
         RETURN
      END IF
*
      // Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 ) RETURN
*
      DO I = 1, N, NB
*
      // Compute the QR factorization of the current block
*
         IB = MIN( N-I+1, NB )
         MB = MIN( M-L+I+IB-1, M )
         IF( I.GE.L ) THEN
            LB = 0
         ELSE
            LB = MB-M+L-I+1
         END IF
*
         CALL CTPQRT2( MB, IB, LB, A(I,I), LDA, B( 1, I ), LDB, T(1, I ), LDT, IINFO )
*
      // Update by applying H**H to B(:,I+IB:N) from the left
*
         IF( I+IB.LE.N ) THEN
            CALL CTPRFB( 'L', 'C', 'F', 'C', MB, N-I-IB+1, IB, LB, B( 1, I ), LDB, T( 1, I ), LDT, A( I, I+IB ), LDA, B( 1, I+IB ), LDB, WORK, IB )
         END IF
      END DO
      RETURN
*
      // End of CTPQRT
*
      END
