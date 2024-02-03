      SUBROUTINE ZLATRD( UPLO, N, NB, A, LDA, E, TAU, W, LDW )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                LDA, LDW, N, NB;
      // ..
      // .. Array Arguments ..
      double             E( * );
      COMPLEX*16         A( LDA, * ), TAU( * ), W( LDW, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX*16         ZERO, ONE, HALF
      PARAMETER          ( ZERO = ( 0.0D+0, 0.0D+0 ), ONE = ( 1.0D+0, 0.0D+0 ), HALF = ( 0.5D+0, 0.0D+0 ) )
      // ..
      // .. Local Scalars ..
      int                I, IW;
      COMPLEX*16         ALPHA
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZAXPY, ZGEMV, ZHEMV, ZLACGV, ZLARFG, ZSCAL
      // ..
      // .. External Functions ..
      bool               LSAME;
      COMPLEX*16         ZDOTC
      // EXTERNAL LSAME, ZDOTC
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MIN
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      IF( N.LE.0 ) RETURN

      IF( LSAME( UPLO, 'U' ) ) THEN

         // Reduce last NB columns of upper triangle

         DO 10 I = N, N - NB + 1, -1
            IW = I - N + NB
            IF( I.LT.N ) THEN

               // Update A(1:i,i)

               A( I, I ) = DBLE( A( I, I ) )
               CALL ZLACGV( N-I, W( I, IW+1 ), LDW )
               CALL ZGEMV( 'No transpose', I, N-I, -ONE, A( 1, I+1 ), LDA, W( I, IW+1 ), LDW, ONE, A( 1, I ), 1 )
               CALL ZLACGV( N-I, W( I, IW+1 ), LDW )
               CALL ZLACGV( N-I, A( I, I+1 ), LDA )
               CALL ZGEMV( 'No transpose', I, N-I, -ONE, W( 1, IW+1 ), LDW, A( I, I+1 ), LDA, ONE, A( 1, I ), 1 )
               CALL ZLACGV( N-I, A( I, I+1 ), LDA )
               A( I, I ) = DBLE( A( I, I ) )
            END IF
            IF( I.GT.1 ) THEN

               // Generate elementary reflector H(i) to annihilate
               // A(1:i-2,i)

               ALPHA = A( I-1, I )
               CALL ZLARFG( I-1, ALPHA, A( 1, I ), 1, TAU( I-1 ) )
               E( I-1 ) = DBLE( ALPHA )
               A( I-1, I ) = ONE

               // Compute W(1:i-1,i)

               CALL ZHEMV( 'Upper', I-1, ONE, A, LDA, A( 1, I ), 1, ZERO, W( 1, IW ), 1 )
               IF( I.LT.N ) THEN
                  CALL ZGEMV( 'Conjugate transpose', I-1, N-I, ONE, W( 1, IW+1 ), LDW, A( 1, I ), 1, ZERO, W( I+1, IW ), 1 )                   CALL ZGEMV( 'No transpose', I-1, N-I, -ONE, A( 1, I+1 ), LDA, W( I+1, IW ), 1, ONE, W( 1, IW ), 1 )                   CALL ZGEMV( 'Conjugate transpose', I-1, N-I, ONE, A( 1, I+1 ), LDA, A( 1, I ), 1, ZERO, W( I+1, IW ), 1 )                   CALL ZGEMV( 'No transpose', I-1, N-I, -ONE, W( 1, IW+1 ), LDW, W( I+1, IW ), 1, ONE, W( 1, IW ), 1 )
               END IF
               CALL ZSCAL( I-1, TAU( I-1 ), W( 1, IW ), 1 )
               ALPHA = -HALF*TAU( I-1 )*ZDOTC( I-1, W( 1, IW ), 1, A( 1, I ), 1 )
               CALL ZAXPY( I-1, ALPHA, A( 1, I ), 1, W( 1, IW ), 1 )
            END IF

   10    CONTINUE
      ELSE

         // Reduce first NB columns of lower triangle

         DO 20 I = 1, NB

            // Update A(i:n,i)

            A( I, I ) = DBLE( A( I, I ) )
            CALL ZLACGV( I-1, W( I, 1 ), LDW )
            CALL ZGEMV( 'No transpose', N-I+1, I-1, -ONE, A( I, 1 ), LDA, W( I, 1 ), LDW, ONE, A( I, I ), 1 )
            CALL ZLACGV( I-1, W( I, 1 ), LDW )
            CALL ZLACGV( I-1, A( I, 1 ), LDA )
            CALL ZGEMV( 'No transpose', N-I+1, I-1, -ONE, W( I, 1 ), LDW, A( I, 1 ), LDA, ONE, A( I, I ), 1 )
            CALL ZLACGV( I-1, A( I, 1 ), LDA )
            A( I, I ) = DBLE( A( I, I ) )
            IF( I.LT.N ) THEN

               // Generate elementary reflector H(i) to annihilate
               // A(i+2:n,i)

               ALPHA = A( I+1, I )
               CALL ZLARFG( N-I, ALPHA, A( MIN( I+2, N ), I ), 1, TAU( I ) )
               E( I ) = DBLE( ALPHA )
               A( I+1, I ) = ONE

               // Compute W(i+1:n,i)

               CALL ZHEMV( 'Lower', N-I, ONE, A( I+1, I+1 ), LDA, A( I+1, I ), 1, ZERO, W( I+1, I ), 1 )                CALL ZGEMV( 'Conjugate transpose', N-I, I-1, ONE, W( I+1, 1 ), LDW, A( I+1, I ), 1, ZERO, W( 1, I ), 1 )                CALL ZGEMV( 'No transpose', N-I, I-1, -ONE, A( I+1, 1 ), LDA, W( 1, I ), 1, ONE, W( I+1, I ), 1 )                CALL ZGEMV( 'Conjugate transpose', N-I, I-1, ONE, A( I+1, 1 ), LDA, A( I+1, I ), 1, ZERO, W( 1, I ), 1 )
               CALL ZGEMV( 'No transpose', N-I, I-1, -ONE, W( I+1, 1 ), LDW, W( 1, I ), 1, ONE, W( I+1, I ), 1 )
               CALL ZSCAL( N-I, TAU( I ), W( I+1, I ), 1 )
               ALPHA = -HALF*TAU( I )*ZDOTC( N-I, W( I+1, I ), 1, A( I+1, I ), 1 )
               CALL ZAXPY( N-I, ALPHA, A( I+1, I ), 1, W( I+1, I ), 1 )
            END IF

   20    CONTINUE
      END IF

      RETURN

      // End of ZLATRD

      END
