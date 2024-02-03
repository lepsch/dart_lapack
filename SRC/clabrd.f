      SUBROUTINE CLABRD( M, N, NB, A, LDA, D, E, TAUQ, TAUP, X, LDX, Y, LDY )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                LDA, LDX, LDY, M, N, NB;
*     ..
*     .. Array Arguments ..
      REAL               D( * ), E( * )
      COMPLEX            A( LDA, * ), TAUP( * ), TAUQ( * ), X( LDX, * ), Y( LDY, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX            ZERO, ONE
      PARAMETER          ( ZERO = ( 0.0E+0, 0.0E+0 ), ONE = ( 1.0E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      int                I;
      COMPLEX            ALPHA
*     ..
*     .. External Subroutines ..
      // EXTERNAL CGEMV, CLACGV, CLARFG, CSCAL
*     ..
*     .. Intrinsic Functions ..
      // INTRINSIC MIN
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( M.LE.0 .OR. N.LE.0 ) RETURN
*
      IF( M.GE.N ) THEN
*
*        Reduce to upper bidiagonal form
*
         DO 10 I = 1, NB
*
*           Update A(i:m,i)
*
            CALL CLACGV( I-1, Y( I, 1 ), LDY )
            CALL CGEMV( 'No transpose', M-I+1, I-1, -ONE, A( I, 1 ), LDA, Y( I, 1 ), LDY, ONE, A( I, I ), 1 )
            CALL CLACGV( I-1, Y( I, 1 ), LDY )
            CALL CGEMV( 'No transpose', M-I+1, I-1, -ONE, X( I, 1 ), LDX, A( 1, I ), 1, ONE, A( I, I ), 1 )
*
*           Generate reflection Q(i) to annihilate A(i+1:m,i)
*
            ALPHA = A( I, I )
            CALL CLARFG( M-I+1, ALPHA, A( MIN( I+1, M ), I ), 1, TAUQ( I ) )
            D( I ) = REAL( ALPHA )
            IF( I.LT.N ) THEN
               A( I, I ) = ONE
*
*              Compute Y(i+1:n,i)
*
               CALL CGEMV( 'Conjugate transpose', M-I+1, N-I, ONE, A( I, I+1 ), LDA, A( I, I ), 1, ZERO, Y( I+1, I ), 1 )                CALL CGEMV( 'Conjugate transpose', M-I+1, I-1, ONE, A( I, 1 ), LDA, A( I, I ), 1, ZERO, Y( 1, I ), 1 )                CALL CGEMV( 'No transpose', N-I, I-1, -ONE, Y( I+1, 1 ), LDY, Y( 1, I ), 1, ONE, Y( I+1, I ), 1 )                CALL CGEMV( 'Conjugate transpose', M-I+1, I-1, ONE, X( I, 1 ), LDX, A( I, I ), 1, ZERO, Y( 1, I ), 1 )                CALL CGEMV( 'Conjugate transpose', I-1, N-I, -ONE, A( 1, I+1 ), LDA, Y( 1, I ), 1, ONE, Y( I+1, I ), 1 )
               CALL CSCAL( N-I, TAUQ( I ), Y( I+1, I ), 1 )
*
*              Update A(i,i+1:n)
*
               CALL CLACGV( N-I, A( I, I+1 ), LDA )
               CALL CLACGV( I, A( I, 1 ), LDA )
               CALL CGEMV( 'No transpose', N-I, I, -ONE, Y( I+1, 1 ), LDY, A( I, 1 ), LDA, ONE, A( I, I+1 ), LDA )
               CALL CLACGV( I, A( I, 1 ), LDA )
               CALL CLACGV( I-1, X( I, 1 ), LDX )
               CALL CGEMV( 'Conjugate transpose', I-1, N-I, -ONE, A( 1, I+1 ), LDA, X( I, 1 ), LDX, ONE, A( I, I+1 ), LDA )
               CALL CLACGV( I-1, X( I, 1 ), LDX )
*
*              Generate reflection P(i) to annihilate A(i,i+2:n)
*
               ALPHA = A( I, I+1 )
               CALL CLARFG( N-I, ALPHA, A( I, MIN( I+2, N ) ), LDA, TAUP( I ) )
               E( I ) = REAL( ALPHA )
               A( I, I+1 ) = ONE
*
*              Compute X(i+1:m,i)
*
               CALL CGEMV( 'No transpose', M-I, N-I, ONE, A( I+1, I+1 ), LDA, A( I, I+1 ), LDA, ZERO, X( I+1, I ), 1 )                CALL CGEMV( 'Conjugate transpose', N-I, I, ONE, Y( I+1, 1 ), LDY, A( I, I+1 ), LDA, ZERO, X( 1, I ), 1 )
               CALL CGEMV( 'No transpose', M-I, I, -ONE, A( I+1, 1 ), LDA, X( 1, I ), 1, ONE, X( I+1, I ), 1 )                CALL CGEMV( 'No transpose', I-1, N-I, ONE, A( 1, I+1 ), LDA, A( I, I+1 ), LDA, ZERO, X( 1, I ), 1 )                CALL CGEMV( 'No transpose', M-I, I-1, -ONE, X( I+1, 1 ), LDX, X( 1, I ), 1, ONE, X( I+1, I ), 1 )
               CALL CSCAL( M-I, TAUP( I ), X( I+1, I ), 1 )
               CALL CLACGV( N-I, A( I, I+1 ), LDA )
            END IF
   10    CONTINUE
      ELSE
*
*        Reduce to lower bidiagonal form
*
         DO 20 I = 1, NB
*
*           Update A(i,i:n)
*
            CALL CLACGV( N-I+1, A( I, I ), LDA )
            CALL CLACGV( I-1, A( I, 1 ), LDA )
            CALL CGEMV( 'No transpose', N-I+1, I-1, -ONE, Y( I, 1 ), LDY, A( I, 1 ), LDA, ONE, A( I, I ), LDA )
            CALL CLACGV( I-1, A( I, 1 ), LDA )
            CALL CLACGV( I-1, X( I, 1 ), LDX )
            CALL CGEMV( 'Conjugate transpose', I-1, N-I+1, -ONE, A( 1, I ), LDA, X( I, 1 ), LDX, ONE, A( I, I ), LDA )
            CALL CLACGV( I-1, X( I, 1 ), LDX )
*
*           Generate reflection P(i) to annihilate A(i,i+1:n)
*
            ALPHA = A( I, I )
            CALL CLARFG( N-I+1, ALPHA, A( I, MIN( I+1, N ) ), LDA, TAUP( I ) )
            D( I ) = REAL( ALPHA )
            IF( I.LT.M ) THEN
               A( I, I ) = ONE
*
*              Compute X(i+1:m,i)
*
               CALL CGEMV( 'No transpose', M-I, N-I+1, ONE, A( I+1, I ), LDA, A( I, I ), LDA, ZERO, X( I+1, I ), 1 )                CALL CGEMV( 'Conjugate transpose', N-I+1, I-1, ONE, Y( I, 1 ), LDY, A( I, I ), LDA, ZERO, X( 1, I ), 1 )
               CALL CGEMV( 'No transpose', M-I, I-1, -ONE, A( I+1, 1 ), LDA, X( 1, I ), 1, ONE, X( I+1, I ), 1 )                CALL CGEMV( 'No transpose', I-1, N-I+1, ONE, A( 1, I ), LDA, A( I, I ), LDA, ZERO, X( 1, I ), 1 )                CALL CGEMV( 'No transpose', M-I, I-1, -ONE, X( I+1, 1 ), LDX, X( 1, I ), 1, ONE, X( I+1, I ), 1 )
               CALL CSCAL( M-I, TAUP( I ), X( I+1, I ), 1 )
               CALL CLACGV( N-I+1, A( I, I ), LDA )
*
*              Update A(i+1:m,i)
*
               CALL CLACGV( I-1, Y( I, 1 ), LDY )
               CALL CGEMV( 'No transpose', M-I, I-1, -ONE, A( I+1, 1 ), LDA, Y( I, 1 ), LDY, ONE, A( I+1, I ), 1 )
               CALL CLACGV( I-1, Y( I, 1 ), LDY )
               CALL CGEMV( 'No transpose', M-I, I, -ONE, X( I+1, 1 ), LDX, A( 1, I ), 1, ONE, A( I+1, I ), 1 )
*
*              Generate reflection Q(i) to annihilate A(i+2:m,i)
*
               ALPHA = A( I+1, I )
               CALL CLARFG( M-I, ALPHA, A( MIN( I+2, M ), I ), 1, TAUQ( I ) )
               E( I ) = REAL( ALPHA )
               A( I+1, I ) = ONE
*
*              Compute Y(i+1:n,i)
*
               CALL CGEMV( 'Conjugate transpose', M-I, N-I, ONE, A( I+1, I+1 ), LDA, A( I+1, I ), 1, ZERO, Y( I+1, I ), 1 )                CALL CGEMV( 'Conjugate transpose', M-I, I-1, ONE, A( I+1, 1 ), LDA, A( I+1, I ), 1, ZERO, Y( 1, I ), 1 )                CALL CGEMV( 'No transpose', N-I, I-1, -ONE, Y( I+1, 1 ), LDY, Y( 1, I ), 1, ONE, Y( I+1, I ), 1 )                CALL CGEMV( 'Conjugate transpose', M-I, I, ONE, X( I+1, 1 ), LDX, A( I+1, I ), 1, ZERO, Y( 1, I ), 1 )                CALL CGEMV( 'Conjugate transpose', I, N-I, -ONE, A( 1, I+1 ), LDA, Y( 1, I ), 1, ONE, Y( I+1, I ), 1 )
               CALL CSCAL( N-I, TAUQ( I ), Y( I+1, I ), 1 )
            ELSE
               CALL CLACGV( N-I+1, A( I, I ), LDA )
            END IF
   20    CONTINUE
      END IF
      RETURN
*
*     End of CLABRD
*
      END
