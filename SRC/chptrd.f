      SUBROUTINE CHPTRD( UPLO, N, AP, D, E, TAU, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, N;
      // ..
      // .. Array Arguments ..
      REAL               D( * ), E( * )
      COMPLEX            AP( * ), TAU( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX            ONE, ZERO, HALF
      PARAMETER          ( ONE = ( 1.0E+0, 0.0E+0 ), ZERO = ( 0.0E+0, 0.0E+0 ), HALF = ( 0.5E+0, 0.0E+0 ) )
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                I, I1, I1I1, II;
      COMPLEX            ALPHA, TAUI
      // ..
      // .. External Subroutines ..
      // EXTERNAL CAXPY, CHPMV, CHPR2, CLARFG, XERBLA
      // ..
      // .. External Functions ..
      bool               LSAME;
      COMPLEX            CDOTC
      // EXTERNAL LSAME, CDOTC
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC REAL
      // ..
      // .. Executable Statements ..

      // Test the input parameters

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CHPTRD', -INFO )
         RETURN
      END IF

      // Quick return if possible

      IF( N.LE.0 ) RETURN

      IF( UPPER ) THEN

         // Reduce the upper triangle of A.
         // I1 is the index in AP of A(1,I+1).

         I1 = N*( N-1 ) / 2 + 1
         AP( I1+N-1 ) = REAL( AP( I1+N-1 ) )
         DO 10 I = N - 1, 1, -1

            // Generate elementary reflector H(i) = I - tau * v * v**H
           t // o annihilate A(1:i-1,i+1)

            ALPHA = AP( I1+I-1 )
            CALL CLARFG( I, ALPHA, AP( I1 ), 1, TAUI )
            E( I ) = REAL( ALPHA )

            IF( TAUI.NE.ZERO ) THEN

               // Apply H(i) from both sides to A(1:i,1:i)

               AP( I1+I-1 ) = ONE

               // Compute  y := tau * A * v  storing y in TAU(1:i)

               CALL CHPMV( UPLO, I, TAUI, AP, AP( I1 ), 1, ZERO, TAU, 1 )

               // Compute  w := y - 1/2 * tau * (y**H *v) * v

               ALPHA = -HALF*TAUI*CDOTC( I, TAU, 1, AP( I1 ), 1 )
               CALL CAXPY( I, ALPHA, AP( I1 ), 1, TAU, 1 )

               // Apply the transformation as a rank-2 update:
                  // A := A - v * w**H - w * v**H

               CALL CHPR2( UPLO, I, -ONE, AP( I1 ), 1, TAU, 1, AP )

            END IF
            AP( I1+I-1 ) = E( I )
            D( I+1 ) = REAL( AP( I1+I ) )
            TAU( I ) = TAUI
            I1 = I1 - I
   10    CONTINUE
         D( 1 ) = REAL( AP( 1 ) )
      ELSE

         // Reduce the lower triangle of A. II is the index in AP of
         // A(i,i) and I1I1 is the index of A(i+1,i+1).

         II = 1
         AP( 1 ) = REAL( AP( 1 ) )
         DO 20 I = 1, N - 1
            I1I1 = II + N - I + 1

            // Generate elementary reflector H(i) = I - tau * v * v**H
           t // o annihilate A(i+2:n,i)

            ALPHA = AP( II+1 )
            CALL CLARFG( N-I, ALPHA, AP( II+2 ), 1, TAUI )
            E( I ) = REAL( ALPHA )

            IF( TAUI.NE.ZERO ) THEN

               // Apply H(i) from both sides to A(i+1:n,i+1:n)

               AP( II+1 ) = ONE

               // Compute  y := tau * A * v  storing y in TAU(i:n-1)

               CALL CHPMV( UPLO, N-I, TAUI, AP( I1I1 ), AP( II+1 ), 1, ZERO, TAU( I ), 1 )

               // Compute  w := y - 1/2 * tau * (y**H *v) * v

               ALPHA = -HALF*TAUI*CDOTC( N-I, TAU( I ), 1, AP( II+1 ), 1 )
               CALL CAXPY( N-I, ALPHA, AP( II+1 ), 1, TAU( I ), 1 )

               // Apply the transformation as a rank-2 update:
                  // A := A - v * w**H - w * v**H

               CALL CHPR2( UPLO, N-I, -ONE, AP( II+1 ), 1, TAU( I ), 1, AP( I1I1 ) )

            END IF
            AP( II+1 ) = E( I )
            D( I ) = REAL( AP( II ) )
            TAU( I ) = TAUI
            II = I1I1
   20    CONTINUE
         D( N ) = REAL( AP( II ) )
      END IF

      RETURN

      // End of CHPTRD

      }
