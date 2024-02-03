      SUBROUTINE SLATRD( UPLO, N, NB, A, LDA, E, TAU, W, LDW );

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                LDA, LDW, N, NB;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * ), E( * ), TAU( * ), W( LDW, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE, HALF;
      const              ZERO = 0.0, ONE = 1.0, HALF = 0.5 ;
      // ..
      // .. Local Scalars ..
      int                I, IW;
      REAL               ALPHA;
      // ..
      // .. External Subroutines ..
      // EXTERNAL SAXPY, SGEMV, SLARFG, SSCAL, SSYMV
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               SDOT;
      // EXTERNAL LSAME, SDOT
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MIN
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      if (N <= 0) RETURN;

      if ( LSAME( UPLO, 'U' ) ) {

         // Reduce last NB columns of upper triangle

         DO 10 I = N, N - NB + 1, -1;
            IW = I - N + NB;
            if ( I < N ) {

               // Update A(1:i,i)

               sgemv('No transpose', I, N-I, -ONE, A( 1, I+1 ), LDA, W( I, IW+1 ), LDW, ONE, A( 1, I ), 1 );
               sgemv('No transpose', I, N-I, -ONE, W( 1, IW+1 ), LDW, A( I, I+1 ), LDA, ONE, A( 1, I ), 1 );
            }
            if ( I > 1 ) {

               // Generate elementary reflector H(i) to annihilate
               // A(1:i-2,i)

               slarfg(I-1, A( I-1, I ), A( 1, I ), 1, TAU( I-1 ) );
               E( I-1 ) = A( I-1, I );
               A( I-1, I ) = ONE;

               // Compute W(1:i-1,i)

               ssymv('Upper', I-1, ONE, A, LDA, A( 1, I ), 1, ZERO, W( 1, IW ), 1 );
               if ( I < N ) {
                  sgemv('Transpose', I-1, N-I, ONE, W( 1, IW+1 ), LDW, A( 1, I ), 1, ZERO, W( I+1, IW ), 1 );
                  sgemv('No transpose', I-1, N-I, -ONE, A( 1, I+1 ), LDA, W( I+1, IW ), 1, ONE, W( 1, IW ), 1 );
                  sgemv('Transpose', I-1, N-I, ONE, A( 1, I+1 ), LDA, A( 1, I ), 1, ZERO, W( I+1, IW ), 1 );
                  sgemv('No transpose', I-1, N-I, -ONE, W( 1, IW+1 ), LDW, W( I+1, IW ), 1, ONE, W( 1, IW ), 1 );
               }
               sscal(I-1, TAU( I-1 ), W( 1, IW ), 1 );
               ALPHA = -HALF*TAU( I-1 )*SDOT( I-1, W( 1, IW ), 1, A( 1, I ), 1 );
               saxpy(I-1, ALPHA, A( 1, I ), 1, W( 1, IW ), 1 );
            }

         } // 10
      } else {

         // Reduce first NB columns of lower triangle

         for (I = 1; I <= NB; I++) { // 20

            // Update A(i:n,i)

            sgemv('No transpose', N-I+1, I-1, -ONE, A( I, 1 ), LDA, W( I, 1 ), LDW, ONE, A( I, I ), 1 );
            sgemv('No transpose', N-I+1, I-1, -ONE, W( I, 1 ), LDW, A( I, 1 ), LDA, ONE, A( I, I ), 1 );
            if ( I < N ) {

               // Generate elementary reflector H(i) to annihilate
               // A(i+2:n,i)

               slarfg(N-I, A( I+1, I ), A( MIN( I+2, N ), I ), 1, TAU( I ) );
               E( I ) = A( I+1, I );
               A( I+1, I ) = ONE;

               // Compute W(i+1:n,i)

               ssymv('Lower', N-I, ONE, A( I+1, I+1 ), LDA, A( I+1, I ), 1, ZERO, W( I+1, I ), 1 );
               sgemv('Transpose', N-I, I-1, ONE, W( I+1, 1 ), LDW, A( I+1, I ), 1, ZERO, W( 1, I ), 1 );
               sgemv('No transpose', N-I, I-1, -ONE, A( I+1, 1 ), LDA, W( 1, I ), 1, ONE, W( I+1, I ), 1 );
               sgemv('Transpose', N-I, I-1, ONE, A( I+1, 1 ), LDA, A( I+1, I ), 1, ZERO, W( 1, I ), 1 );
               sgemv('No transpose', N-I, I-1, -ONE, W( I+1, 1 ), LDW, W( 1, I ), 1, ONE, W( I+1, I ), 1 );
               sscal(N-I, TAU( I ), W( I+1, I ), 1 );
               ALPHA = -HALF*TAU( I )*SDOT( N-I, W( I+1, I ), 1, A( I+1, I ), 1 );
               saxpy(N-I, ALPHA, A( I+1, I ), 1, W( I+1, I ), 1 );
            }

         } // 20
      }

      RETURN;

      // End of SLATRD

      }
