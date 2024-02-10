      void clatrd(UPLO, N, NB, A, LDA, E, TAU, W, LDW ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                LDA, LDW, N, NB;
      double               E( * );
      Complex            A( LDA, * ), TAU( * ), W( LDW, * );
      // ..

      Complex            ZERO, ONE, HALF;
      const              ZERO = ( 0.0, 0.0 ), ONE = ( 1.0, 0.0 ), HALF = ( 0.5, 0.0 ) ;
      int                I, IW;
      Complex            ALPHA;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CAXPY, CGEMV, CHEMV, CLACGV, CLARFG, CSCAL
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- COMPLEX            CDOTC;
      // EXTERNAL lsame, CDOTC
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MIN, REAL

      // Quick return if possible

      if (N <= 0) return;

      if ( lsame( UPLO, 'U' ) ) {

         // Reduce last NB columns of upper triangle

         for (I = N; I >= N - NB + 1; I--) { // 10
            IW = I - N + NB;
            if ( I < N ) {

               // Update A(1:i,i)

               A[I][I] = double( A( I, I ) );
               clacgv(N-I, W( I, IW+1 ), LDW );
               cgemv('No transpose', I, N-I, -ONE, A( 1, I+1 ), LDA, W( I, IW+1 ), LDW, ONE, A( 1, I ), 1 );
               clacgv(N-I, W( I, IW+1 ), LDW );
               clacgv(N-I, A( I, I+1 ), LDA );
               cgemv('No transpose', I, N-I, -ONE, W( 1, IW+1 ), LDW, A( I, I+1 ), LDA, ONE, A( 1, I ), 1 );
               clacgv(N-I, A( I, I+1 ), LDA );
               A[I][I] = double( A( I, I ) );
            }
            if ( I > 1 ) {

               // Generate elementary reflector H(i) to annihilate
               // A(1:i-2,i)

               ALPHA = A( I-1, I );
               clarfg(I-1, ALPHA, A( 1, I ), 1, TAU( I-1 ) );
               E[I-1] = double( ALPHA );
               A[I-1][I] = ONE;

               // Compute W(1:i-1,i)

               chemv('Upper', I-1, ONE, A, LDA, A( 1, I ), 1, ZERO, W( 1, IW ), 1 );
               if ( I < N ) {
                  cgemv('Conjugate transpose', I-1, N-I, ONE, W( 1, IW+1 ), LDW, A( 1, I ), 1, ZERO, W( I+1, IW ), 1 );
                  cgemv('No transpose', I-1, N-I, -ONE, A( 1, I+1 ), LDA, W( I+1, IW ), 1, ONE, W( 1, IW ), 1 );
                  cgemv('Conjugate transpose', I-1, N-I, ONE, A( 1, I+1 ), LDA, A( 1, I ), 1, ZERO, W( I+1, IW ), 1 );
                  cgemv('No transpose', I-1, N-I, -ONE, W( 1, IW+1 ), LDW, W( I+1, IW ), 1, ONE, W( 1, IW ), 1 );
               }
               cscal(I-1, TAU( I-1 ), W( 1, IW ), 1 );
               ALPHA = -HALF*TAU( I-1 )*CDOTC( I-1, W( 1, IW ), 1, A( 1, I ), 1 );
               caxpy(I-1, ALPHA, A( 1, I ), 1, W( 1, IW ), 1 );
            }

         } // 10
      } else {

         // Reduce first NB columns of lower triangle

         for (I = 1; I <= NB; I++) { // 20

            // Update A(i:n,i)

            A[I][I] = double( A( I, I ) );
            clacgv(I-1, W( I, 1 ), LDW );
            cgemv('No transpose', N-I+1, I-1, -ONE, A( I, 1 ), LDA, W( I, 1 ), LDW, ONE, A( I, I ), 1 );
            clacgv(I-1, W( I, 1 ), LDW );
            clacgv(I-1, A( I, 1 ), LDA );
            cgemv('No transpose', N-I+1, I-1, -ONE, W( I, 1 ), LDW, A( I, 1 ), LDA, ONE, A( I, I ), 1 );
            clacgv(I-1, A( I, 1 ), LDA );
            A[I][I] = double( A( I, I ) );
            if ( I < N ) {

               // Generate elementary reflector H(i) to annihilate
               // A(i+2:n,i)

               ALPHA = A( I+1, I );
               clarfg(N-I, ALPHA, A( min( I+2, N ), I ), 1, TAU( I ) );
               E[I] = double( ALPHA );
               A[I+1][I] = ONE;

               // Compute W(i+1:n,i)

               chemv('Lower', N-I, ONE, A( I+1, I+1 ), LDA, A( I+1, I ), 1, ZERO, W( I+1, I ), 1 );
               cgemv('Conjugate transpose', N-I, I-1, ONE, W( I+1, 1 ), LDW, A( I+1, I ), 1, ZERO, W( 1, I ), 1 );
               cgemv('No transpose', N-I, I-1, -ONE, A( I+1, 1 ), LDA, W( 1, I ), 1, ONE, W( I+1, I ), 1 );
               cgemv('Conjugate transpose', N-I, I-1, ONE, A( I+1, 1 ), LDA, A( I+1, I ), 1, ZERO, W( 1, I ), 1 );
               cgemv('No transpose', N-I, I-1, -ONE, W( I+1, 1 ), LDW, W( 1, I ), 1, ONE, W( I+1, I ), 1 );
               cscal(N-I, TAU( I ), W( I+1, I ), 1 );
               ALPHA = -HALF*TAU( I )*CDOTC( N-I, W( I+1, I ), 1, A( I+1, I ), 1 );
               caxpy(N-I, ALPHA, A( I+1, I ), 1, W( I+1, I ), 1 );
            }

         } // 20
      }

      }
