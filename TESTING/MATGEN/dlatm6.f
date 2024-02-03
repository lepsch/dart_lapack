      SUBROUTINE DLATM6( TYPE, N, A, LDA, B, X, LDX, Y, LDY, ALPHA, BETA, WX, WY, S, DIF )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, LDX, LDY, N, TYPE;
      double             ALPHA, BETA, WX, WY;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), B( LDA, * ), DIF( * ), S( * ), X( LDX, * ), Y( LDY, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE, TWO, THREE;
      const              ZERO = 0.0, ONE = 1.0, TWO = 2.0, THREE = 3.0 ;
      // ..
      // .. Local Scalars ..
      int                I, INFO, J;
      // ..
      // .. Local Arrays ..
      double             WORK( 100 ), Z( 12, 12 );
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, SQRT
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGESVD, DLACPY, DLAKF2
      // ..
      // .. Executable Statements ..

      // Generate test problem ...
      // (Da, Db) ...

      for (I = 1; I <= N; I++) { // 20
         for (J = 1; J <= N; J++) { // 10

            if ( I == J ) {
               A( I, I ) = DBLE( I ) + ALPHA
               B( I, I ) = ONE
            } else {
               A( I, J ) = ZERO
               B( I, J ) = ZERO
            }

         } // 10
      } // 20

      // Form X and Y

      dlacpy('F', N, N, B, LDA, Y, LDY );
      Y( 3, 1 ) = -WY
      Y( 4, 1 ) = WY
      Y( 5, 1 ) = -WY
      Y( 3, 2 ) = -WY
      Y( 4, 2 ) = WY
      Y( 5, 2 ) = -WY

      dlacpy('F', N, N, B, LDA, X, LDX );
      X( 1, 3 ) = -WX
      X( 1, 4 ) = -WX
      X( 1, 5 ) = WX
      X( 2, 3 ) = WX
      X( 2, 4 ) = -WX
      X( 2, 5 ) = -WX

      // Form (A, B)

      B( 1, 3 ) = WX + WY
      B( 2, 3 ) = -WX + WY
      B( 1, 4 ) = WX - WY
      B( 2, 4 ) = WX - WY
      B( 1, 5 ) = -WX + WY
      B( 2, 5 ) = WX + WY
      if ( TYPE == 1 ) {
         A( 1, 3 ) = WX*A( 1, 1 ) + WY*A( 3, 3 )
         A( 2, 3 ) = -WX*A( 2, 2 ) + WY*A( 3, 3 )
         A( 1, 4 ) = WX*A( 1, 1 ) - WY*A( 4, 4 )
         A( 2, 4 ) = WX*A( 2, 2 ) - WY*A( 4, 4 )
         A( 1, 5 ) = -WX*A( 1, 1 ) + WY*A( 5, 5 )
         A( 2, 5 ) = WX*A( 2, 2 ) + WY*A( 5, 5 )
      } else if ( TYPE == 2 ) {
         A( 1, 3 ) = TWO*WX + WY
         A( 2, 3 ) = WY
         A( 1, 4 ) = -WY*( TWO+ALPHA+BETA )
         A( 2, 4 ) = TWO*WX - WY*( TWO+ALPHA+BETA )
         A( 1, 5 ) = -TWO*WX + WY*( ALPHA-BETA )
         A( 2, 5 ) = WY*( ALPHA-BETA )
         A( 1, 1 ) = ONE
         A( 1, 2 ) = -ONE
         A( 2, 1 ) = ONE
         A( 2, 2 ) = A( 1, 1 )
         A( 3, 3 ) = ONE
         A( 4, 4 ) = ONE + ALPHA
         A( 4, 5 ) = ONE + BETA
         A( 5, 4 ) = -A( 4, 5 )
         A( 5, 5 ) = A( 4, 4 )
      }

      // Compute condition numbers

      if ( TYPE == 1 ) {

         S( 1 ) = ONE / SQRT( ( ONE+THREE*WY*WY ) / ( ONE+A( 1, 1 )*A( 1, 1 ) ) )          S( 2 ) = ONE / SQRT( ( ONE+THREE*WY*WY ) / ( ONE+A( 2, 2 )*A( 2, 2 ) ) )          S( 3 ) = ONE / SQRT( ( ONE+TWO*WX*WX ) / ( ONE+A( 3, 3 )*A( 3, 3 ) ) )          S( 4 ) = ONE / SQRT( ( ONE+TWO*WX*WX ) / ( ONE+A( 4, 4 )*A( 4, 4 ) ) )          S( 5 ) = ONE / SQRT( ( ONE+TWO*WX*WX ) / ( ONE+A( 5, 5 )*A( 5, 5 ) ) )

         dlakf2(1, 4, A, LDA, A( 2, 2 ), B, B( 2, 2 ), Z, 12 );
         dgesvd('N', 'N', 8, 8, Z, 12, WORK, WORK( 9 ), 1, WORK( 10 ), 1, WORK( 11 ), 40, INFO );
         DIF( 1 ) = WORK( 8 )

         dlakf2(4, 1, A, LDA, A( 5, 5 ), B, B( 5, 5 ), Z, 12 );
         dgesvd('N', 'N', 8, 8, Z, 12, WORK, WORK( 9 ), 1, WORK( 10 ), 1, WORK( 11 ), 40, INFO );
         DIF( 5 ) = WORK( 8 )

      } else if ( TYPE == 2 ) {

         S( 1 ) = ONE / SQRT( ONE / THREE+WY*WY )
         S( 2 ) = S( 1 )
         S( 3 ) = ONE / SQRT( ONE / TWO+WX*WX )
         S( 4 ) = ONE / SQRT( ( ONE+TWO*WX*WX ) / ( ONE+( ONE+ALPHA )*( ONE+ALPHA )+( ONE+BETA )*( ONE+ BETA ) ) )
         S( 5 ) = S( 4 )

         dlakf2(2, 3, A, LDA, A( 3, 3 ), B, B( 3, 3 ), Z, 12 );
         dgesvd('N', 'N', 12, 12, Z, 12, WORK, WORK( 13 ), 1, WORK( 14 ), 1, WORK( 15 ), 60, INFO );
         DIF( 1 ) = WORK( 12 )

         dlakf2(3, 2, A, LDA, A( 4, 4 ), B, B( 4, 4 ), Z, 12 );
         dgesvd('N', 'N', 12, 12, Z, 12, WORK, WORK( 13 ), 1, WORK( 14 ), 1, WORK( 15 ), 60, INFO );
         DIF( 5 ) = WORK( 12 )

      }

      RETURN

      // End of DLATM6

      }
