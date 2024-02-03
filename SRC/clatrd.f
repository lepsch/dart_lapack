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
      const              ZERO = ( 0.0E+0, 0.0E+0 ), ONE = ( 1.0E+0, 0.0E+0 ), HALF = ( 0.5E+0, 0.0E+0 ) ;
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

      if (N.LE.0) RETURN;

      if ( LSAME( UPLO, 'U' ) ) {

         // Reduce last NB columns of upper triangle

         DO 10 I = N, N - NB + 1, -1
            IW = I - N + NB
            if ( I.LT.N ) {

               // Update A(1:i,i)

               A( I, I ) = REAL( A( I, I ) )
               clacgv(N-I, W( I, IW+1 ), LDW );
               cgemv('No transpose', I, N-I, -ONE, A( 1, I+1 ), LDA, W( I, IW+1 ), LDW, ONE, A( 1, I ), 1 );
               clacgv(N-I, W( I, IW+1 ), LDW );
               clacgv(N-I, A( I, I+1 ), LDA );
               cgemv('No transpose', I, N-I, -ONE, W( 1, IW+1 ), LDW, A( I, I+1 ), LDA, ONE, A( 1, I ), 1 );
               clacgv(N-I, A( I, I+1 ), LDA );
               A( I, I ) = REAL( A( I, I ) )
            }
            if ( I.GT.1 ) {

               // Generate elementary reflector H(i) to annihilate
               // A(1:i-2,i)

               ALPHA = A( I-1, I )
               clarfg(I-1, ALPHA, A( 1, I ), 1, TAU( I-1 ) );
               E( I-1 ) = REAL( ALPHA )
               A( I-1, I ) = ONE

               // Compute W(1:i-1,i)

               chemv('Upper', I-1, ONE, A, LDA, A( 1, I ), 1, ZERO, W( 1, IW ), 1 );
               if ( I.LT.N ) {
                  cgemv('Conjugate transpose', I-1, N-I, ONE, W( 1, IW+1 ), LDW, A( 1, I ), 1, ZERO, W( I+1, IW ), 1 );
                  cgemv('No transpose', I-1, N-I, -ONE, A( 1, I+1 ), LDA, W( I+1, IW ), 1, ONE, W( 1, IW ), 1 );
                  cgemv('Conjugate transpose', I-1, N-I, ONE, A( 1, I+1 ), LDA, A( 1, I ), 1, ZERO, W( I+1, IW ), 1 );
                  cgemv('No transpose', I-1, N-I, -ONE, W( 1, IW+1 ), LDW, W( I+1, IW ), 1, ONE, W( 1, IW ), 1 );
               }
               cscal(I-1, TAU( I-1 ), W( 1, IW ), 1 );
               ALPHA = -HALF*TAU( I-1 )*CDOTC( I-1, W( 1, IW ), 1, A( 1, I ), 1 )
               caxpy(I-1, ALPHA, A( 1, I ), 1, W( 1, IW ), 1 );
            }

         } // 10
      } else {

         // Reduce first NB columns of lower triangle

         for (I = 1; I <= NB; I++) { // 20

            // Update A(i:n,i)

            A( I, I ) = REAL( A( I, I ) )
            clacgv(I-1, W( I, 1 ), LDW );
            cgemv('No transpose', N-I+1, I-1, -ONE, A( I, 1 ), LDA, W( I, 1 ), LDW, ONE, A( I, I ), 1 );
            clacgv(I-1, W( I, 1 ), LDW );
            clacgv(I-1, A( I, 1 ), LDA );
            cgemv('No transpose', N-I+1, I-1, -ONE, W( I, 1 ), LDW, A( I, 1 ), LDA, ONE, A( I, I ), 1 );
            clacgv(I-1, A( I, 1 ), LDA );
            A( I, I ) = REAL( A( I, I ) )
            if ( I.LT.N ) {

               // Generate elementary reflector H(i) to annihilate
               // A(i+2:n,i)

               ALPHA = A( I+1, I )
               clarfg(N-I, ALPHA, A( MIN( I+2, N ), I ), 1, TAU( I ) );
               E( I ) = REAL( ALPHA )
               A( I+1, I ) = ONE

               // Compute W(i+1:n,i)

               chemv('Lower', N-I, ONE, A( I+1, I+1 ), LDA, A( I+1, I ), 1, ZERO, W( I+1, I ), 1 );
               cgemv('Conjugate transpose', N-I, I-1, ONE, W( I+1, 1 ), LDW, A( I+1, I ), 1, ZERO, W( 1, I ), 1 );
               cgemv('No transpose', N-I, I-1, -ONE, A( I+1, 1 ), LDA, W( 1, I ), 1, ONE, W( I+1, I ), 1 );
               cgemv('Conjugate transpose', N-I, I-1, ONE, A( I+1, 1 ), LDA, A( I+1, I ), 1, ZERO, W( 1, I ), 1 );
               cgemv('No transpose', N-I, I-1, -ONE, W( I+1, 1 ), LDW, W( 1, I ), 1, ONE, W( I+1, I ), 1 );
               cscal(N-I, TAU( I ), W( I+1, I ), 1 );
               ALPHA = -HALF*TAU( I )*CDOTC( N-I, W( I+1, I ), 1, A( I+1, I ), 1 )
               caxpy(N-I, ALPHA, A( I+1, I ), 1, W( I+1, I ), 1 );
            }

         } // 20
      }

      RETURN

      // End of CLATRD

      }
