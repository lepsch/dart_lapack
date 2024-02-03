      SUBROUTINE ZGEBD2( M, N, A, LDA, D, E, TAUQ, TAUP, WORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                INFO, LDA, M, N;
*     ..
*     .. Array Arguments ..
      double             D( * ), E( * );
      COMPLEX*16         A( LDA, * ), TAUP( * ), TAUQ( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ZERO, ONE
      PARAMETER          ( ZERO = ( 0.0D+0, 0.0D+0 ), ONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      int                I;
      COMPLEX*16         ALPHA
*     ..
*     .. External Subroutines ..
      // EXTERNAL XERBLA, ZLACGV, ZLARF, ZLARFG
*     ..
*     .. Intrinsic Functions ..
      // INTRINSIC DCONJG, MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.LT.0 ) THEN
         CALL XERBLA( 'ZGEBD2', -INFO )
         RETURN
      END IF
*
      IF( M.GE.N ) THEN
*
*        Reduce to upper bidiagonal form
*
         DO 10 I = 1, N
*
*           Generate elementary reflector H(i) to annihilate A(i+1:m,i)
*
            ALPHA = A( I, I )
            CALL ZLARFG( M-I+1, ALPHA, A( MIN( I+1, M ), I ), 1, TAUQ( I ) )
            D( I ) = DBLE( ALPHA )
            A( I, I ) = ONE
*
*           Apply H(i)**H to A(i:m,i+1:n) from the left
*
            IF( I.LT.N ) CALL ZLARF( 'Left', M-I+1, N-I, A( I, I ), 1, DCONJG( TAUQ( I ) ), A( I, I+1 ), LDA, WORK )
            A( I, I ) = D( I )
*
            IF( I.LT.N ) THEN
*
*              Generate elementary reflector G(i) to annihilate
*              A(i,i+2:n)
*
               CALL ZLACGV( N-I, A( I, I+1 ), LDA )
               ALPHA = A( I, I+1 )
               CALL ZLARFG( N-I, ALPHA, A( I, MIN( I+2, N ) ), LDA, TAUP( I ) )
               E( I ) = DBLE( ALPHA )
               A( I, I+1 ) = ONE
*
*              Apply G(i) to A(i+1:m,i+1:n) from the right
*
               CALL ZLARF( 'Right', M-I, N-I, A( I, I+1 ), LDA, TAUP( I ), A( I+1, I+1 ), LDA, WORK )
               CALL ZLACGV( N-I, A( I, I+1 ), LDA )
               A( I, I+1 ) = E( I )
            ELSE
               TAUP( I ) = ZERO
            END IF
   10    CONTINUE
      ELSE
*
*        Reduce to lower bidiagonal form
*
         DO 20 I = 1, M
*
*           Generate elementary reflector G(i) to annihilate A(i,i+1:n)
*
            CALL ZLACGV( N-I+1, A( I, I ), LDA )
            ALPHA = A( I, I )
            CALL ZLARFG( N-I+1, ALPHA, A( I, MIN( I+1, N ) ), LDA, TAUP( I ) )
            D( I ) = DBLE( ALPHA )
            A( I, I ) = ONE
*
*           Apply G(i) to A(i+1:m,i:n) from the right
*
            IF( I.LT.M ) CALL ZLARF( 'Right', M-I, N-I+1, A( I, I ), LDA, TAUP( I ), A( I+1, I ), LDA, WORK )
            CALL ZLACGV( N-I+1, A( I, I ), LDA )
            A( I, I ) = D( I )
*
            IF( I.LT.M ) THEN
*
*              Generate elementary reflector H(i) to annihilate
*              A(i+2:m,i)
*
               ALPHA = A( I+1, I )
               CALL ZLARFG( M-I, ALPHA, A( MIN( I+2, M ), I ), 1, TAUQ( I ) )
               E( I ) = DBLE( ALPHA )
               A( I+1, I ) = ONE
*
*              Apply H(i)**H to A(i+1:m,i+1:n) from the left
*
               CALL ZLARF( 'Left', M-I, N-I, A( I+1, I ), 1, DCONJG( TAUQ( I ) ), A( I+1, I+1 ), LDA, WORK )
               A( I+1, I ) = E( I )
            ELSE
               TAUQ( I ) = ZERO
            END IF
   20    CONTINUE
      END IF
      RETURN
*
*     End of ZGEBD2
*
      END
