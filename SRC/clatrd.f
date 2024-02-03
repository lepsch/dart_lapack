      SUBROUTINE CLATRD( UPLO, N, NB, A, LDA, E, TAU, W, LDW )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                LDA, LDW, N, NB;
      // ..
      // .. Array Arguments ..
      REAL               E( * )
      COMPLEX            A( LDA, * ), TAU( * ), W( LDW, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX            ZERO, ONE, HALF
      PARAMETER          ( ZERO = ( 0.0E+0, 0.0E+0 ), ONE = ( 1.0E+0, 0.0E+0 ), HALF = ( 0.5E+0, 0.0E+0 ) )
      // ..
      // .. Local Scalars ..
      int                I, IW;
      COMPLEX            ALPHA
      // ..
      // .. External Subroutines ..
      // EXTERNAL CAXPY, CGEMV, CHEMV, CLACGV, CLARFG, CSCAL
      // ..
      // .. External Functions ..
      bool               LSAME;
      COMPLEX            CDOTC
      // EXTERNAL LSAME, CDOTC
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MIN, REAL
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

               A( I, I ) = REAL( A( I, I ) )
               CALL CLACGV( N-I, W( I, IW+1 ), LDW )
               CALL CGEMV( 'No transpose', I, N-I, -ONE, A( 1, I+1 ), LDA, W( I, IW+1 ), LDW, ONE, A( 1, I ), 1 )
               CALL CLACGV( N-I, W( I, IW+1 ), LDW )
               CALL CLACGV( N-I, A( I, I+1 ), LDA )
               CALL CGEMV( 'No transpose', I, N-I, -ONE, W( 1, IW+1 ), LDW, A( I, I+1 ), LDA, ONE, A( 1, I ), 1 )
               CALL CLACGV( N-I, A( I, I+1 ), LDA )
               A( I, I ) = REAL( A( I, I ) )
            END IF
            IF( I.GT.1 ) THEN

               // Generate elementary reflector H(i) to annihilate
               // A(1:i-2,i)

               ALPHA = A( I-1, I )
               CALL CLARFG( I-1, ALPHA, A( 1, I ), 1, TAU( I-1 ) )
               E( I-1 ) = REAL( ALPHA )
               A( I-1, I ) = ONE

               // Compute W(1:i-1,i)

               CALL CHEMV( 'Upper', I-1, ONE, A, LDA, A( 1, I ), 1, ZERO, W( 1, IW ), 1 )
               IF( I.LT.N ) THEN
                  CALL CGEMV( 'Conjugate transpose', I-1, N-I, ONE, W( 1, IW+1 ), LDW, A( 1, I ), 1, ZERO, W( I+1, IW ), 1 )                   CALL CGEMV( 'No transpose', I-1, N-I, -ONE, A( 1, I+1 ), LDA, W( I+1, IW ), 1, ONE, W( 1, IW ), 1 )                   CALL CGEMV( 'Conjugate transpose', I-1, N-I, ONE, A( 1, I+1 ), LDA, A( 1, I ), 1, ZERO, W( I+1, IW ), 1 )                   CALL CGEMV( 'No transpose', I-1, N-I, -ONE, W( 1, IW+1 ), LDW, W( I+1, IW ), 1, ONE, W( 1, IW ), 1 )
               END IF
               CALL CSCAL( I-1, TAU( I-1 ), W( 1, IW ), 1 )
               ALPHA = -HALF*TAU( I-1 )*CDOTC( I-1, W( 1, IW ), 1, A( 1, I ), 1 )
               CALL CAXPY( I-1, ALPHA, A( 1, I ), 1, W( 1, IW ), 1 )
            END IF

   10    CONTINUE
      ELSE

         // Reduce first NB columns of lower triangle

         DO 20 I = 1, NB

            // Update A(i:n,i)

            A( I, I ) = REAL( A( I, I ) )
            CALL CLACGV( I-1, W( I, 1 ), LDW )
            CALL CGEMV( 'No transpose', N-I+1, I-1, -ONE, A( I, 1 ), LDA, W( I, 1 ), LDW, ONE, A( I, I ), 1 )
            CALL CLACGV( I-1, W( I, 1 ), LDW )
            CALL CLACGV( I-1, A( I, 1 ), LDA )
            CALL CGEMV( 'No transpose', N-I+1, I-1, -ONE, W( I, 1 ), LDW, A( I, 1 ), LDA, ONE, A( I, I ), 1 )
            CALL CLACGV( I-1, A( I, 1 ), LDA )
            A( I, I ) = REAL( A( I, I ) )
            IF( I.LT.N ) THEN

               // Generate elementary reflector H(i) to annihilate
               // A(i+2:n,i)

               ALPHA = A( I+1, I )
               CALL CLARFG( N-I, ALPHA, A( MIN( I+2, N ), I ), 1, TAU( I ) )
               E( I ) = REAL( ALPHA )
               A( I+1, I ) = ONE

               // Compute W(i+1:n,i)

               CALL CHEMV( 'Lower', N-I, ONE, A( I+1, I+1 ), LDA, A( I+1, I ), 1, ZERO, W( I+1, I ), 1 )                CALL CGEMV( 'Conjugate transpose', N-I, I-1, ONE, W( I+1, 1 ), LDW, A( I+1, I ), 1, ZERO, W( 1, I ), 1 )                CALL CGEMV( 'No transpose', N-I, I-1, -ONE, A( I+1, 1 ), LDA, W( 1, I ), 1, ONE, W( I+1, I ), 1 )                CALL CGEMV( 'Conjugate transpose', N-I, I-1, ONE, A( I+1, 1 ), LDA, A( I+1, I ), 1, ZERO, W( 1, I ), 1 )
               CALL CGEMV( 'No transpose', N-I, I-1, -ONE, W( I+1, 1 ), LDW, W( 1, I ), 1, ONE, W( I+1, I ), 1 )
               CALL CSCAL( N-I, TAU( I ), W( I+1, I ), 1 )
               ALPHA = -HALF*TAU( I )*CDOTC( N-I, W( I+1, I ), 1, A( I+1, I ), 1 )
               CALL CAXPY( N-I, ALPHA, A( I+1, I ), 1, W( I+1, I ), 1 )
            END IF

   20    CONTINUE
      END IF

      RETURN

      // End of CLATRD

      END
