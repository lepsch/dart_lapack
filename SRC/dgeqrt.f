      SUBROUTINE DGEQRT( M, N, NB, A, LDA, T, LDT, WORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int     INFO, LDA, LDT, M, N, NB;
      // ..
      // .. Array Arguments ..
      double           A( LDA, * ), T( LDT, * ), WORK( * );
      // ..

* =====================================================================

      // ..
      // .. Local Scalars ..
      int        I, IB, IINFO, K;
      bool       USE_RECURSIVE_QR;
      const    USE_RECURSIVE_QR=.TRUE. ;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEQRT2, DGEQRT3, DLARFB, XERBLA
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NB.LT.1 .OR. ( NB.GT.MIN(M,N) .AND. MIN(M,N).GT.0 ) )THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LDT.LT.NB ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGEQRT', -INFO )
         RETURN
      END IF

      // Quick return if possible

      K = MIN( M, N )
      IF( K.EQ.0 ) RETURN

      // Blocked loop of length K

      DO I = 1, K,  NB
         IB = MIN( K-I+1, NB )

      // Compute the QR factorization of the current block A(I:M,I:I+IB-1)

         IF( USE_RECURSIVE_QR ) THEN
            CALL DGEQRT3( M-I+1, IB, A(I,I), LDA, T(1,I), LDT, IINFO )
         ELSE
            CALL DGEQRT2( M-I+1, IB, A(I,I), LDA, T(1,I), LDT, IINFO )
         END IF
         IF( I+IB.LE.N ) THEN

      // Update by applying H**T to A(I:M,I+IB:N) from the left

            CALL DLARFB( 'L', 'T', 'F', 'C', M-I+1, N-I-IB+1, IB, A( I, I ), LDA, T( 1, I ), LDT, A( I, I+IB ), LDA, WORK , N-I-IB+1 )
         END IF
      END DO
      RETURN

      // End of DGEQRT

      }
