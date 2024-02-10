      void zgebal(final int JOB, final int N, final Matrix<double> A, final int LDA, final int ILO, final int IHI, final int SCALE, final Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             JOB;
      int                IHI, ILO, INFO, LDA, N;
      double             SCALE( * );
      Complex         A( LDA, * );
      // ..

      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      double             SCLFAC;
      const              SCLFAC = 2.0 ;
      double             FACTOR;
      const              FACTOR = 0.95 ;
      bool               NOCONV, CANSWAP;
      int                I, ICA, IRA, J, K, L;
      double             C, CA, F, G, R, RA, S, SFMAX1, SFMAX2, SFMIN1, SFMIN2;
      // ..
      // .. External Functions ..
      //- bool               DISNAN, lsame;
      //- int                IZAMAX;
      //- double             DLAMCH, DZNRM2;
      // EXTERNAL DISNAN, lsame, IZAMAX, DLAMCH, DZNRM2
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZDSCAL, ZSWAP
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DIMAG, MAX, MIN

      // Test the input parameters

      INFO = 0;
      if ( !lsame( JOB, 'N' ) && !lsame( JOB, 'P' ) && !lsame( JOB, 'S' ) && !lsame( JOB, 'B' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -4;
      }
      if ( INFO != 0 ) {
         xerbla('ZGEBAL', -INFO );
         return;
      }

      // Quick returns.

      if ( N == 0 ) {
         ILO = 1;
         IHI = 0;
         return;
      }

      if ( lsame( JOB, 'N' ) ) {
         for (I = 1; I <= N; I++) {
            SCALE[I] = ONE;
         }
         ILO = 1;
         IHI = N;
         return;
      }

      // Permutation to isolate eigenvalues if possible.

      K = 1;
      L = N;

      if ( !lsame( JOB, 'S' ) ) {

         // Row and column exchange.

         NOCONV = true;
         while (NOCONV) {

            // Search for rows isolating an eigenvalue and push them down.

            NOCONV = false;
            for (I = L; I >= 1; I--) {
               CANSWAP = true;
               for (J = 1; J <= L; J++) {
                  if ( I != J && ( (A( I, J )).toDouble() != ZERO || DIMAG( A( I, J ) ) != ZERO ) ) {
                     CANSWAP = false;
                     break;
                  }
               }

               if ( CANSWAP ) {
                  SCALE[L] = I;
                  if ( I != L ) {
                     zswap(L, A( 1, I ), 1, A( 1, L ), 1 );
                     zswap(N-K+1, A( I, K ), LDA, A( L, K ), LDA );
                  }
                  NOCONV = true;

                  if ( L == 1 ) {
                     ILO = 1;
                     IHI = 1;
                     return;
                  }

                  L = L - 1;
               }
            }

         }

         NOCONV = true;
         while (NOCONV) {

            // Search for columns isolating an eigenvalue and push them left.

            NOCONV = false;
            for (J = K; J <= L; J++) {
               CANSWAP = true;
               for (I = K; I <= L; I++) {
                  if ( I != J && ( (A( I, J )).toDouble() != ZERO || DIMAG( A( I, J ) ) != ZERO ) ) {
                     CANSWAP = false;
                     break;
                  }
               }

               if ( CANSWAP ) {
                  SCALE[K] = J;
                  if ( J != K ) {
                     zswap(L, A( 1, J ), 1, A( 1, K ), 1 );
                     zswap(N-K+1, A( J, K ), LDA, A( K, K ), LDA );
                  }
                  NOCONV = true;

                  K = K + 1;
               }
            }

         }

      }

      // Initialize SCALE for non-permuted submatrix.

      for (I = K; I <= L; I++) {
         SCALE[I] = ONE;
      }

      // If we only had to permute, we are done.

      if ( lsame( JOB, 'P' ) ) {
         ILO = K;
         IHI = L;
         return;
      }

      // Balance the submatrix in rows K to L.

      // Iterative loop for norm reduction.

      SFMIN1 = dlamch( 'S' ) / dlamch( 'P' );
      SFMAX1 = ONE / SFMIN1;
      SFMIN2 = SFMIN1*SCLFAC;
      SFMAX2 = ONE / SFMIN2;

      NOCONV = true;
      while (NOCONV) {
         NOCONV = false;

         for (I = K; I <= L; I++) {

            C = DZNRM2( L-K+1, A( K, I ), 1 );
            R = DZNRM2( L-K+1, A( I, K ), LDA );
            ICA = IZAMAX( L, A( 1, I ), 1 );
            CA = ( A( ICA, I ) ).abs();
            IRA = IZAMAX( N-K+1, A( I, K ), LDA );
            RA = ( A( I, IRA+K-1 ) ).abs();

            // Guard against zero C or R due to underflow.

            if (C == ZERO || R == ZERO) continue;

            // Exit if NaN to avoid infinite loop

            if ( disnan( C+CA+R+RA ) ) {
               INFO = -3;
               xerbla('ZGEBAL', -INFO );
               return;
            }

            G = R / SCLFAC;
            F = ONE;
            S = C + R;

            while (C < G && max( F, C, CA ) < SFMAX2 && min( R, G, RA ) > SFMIN2) {
               F = F*SCLFAC;
               C = C*SCLFAC;
               CA = CA*SCLFAC;
               R = R / SCLFAC;
               G = G / SCLFAC;
               RA = RA / SCLFAC;
            }

            G = C / SCLFAC;

            while (G >= R && max( R, RA ) < SFMAX2 && min( F, C, G, CA ) > SFMIN2) {
               F = F / SCLFAC;
               C = C / SCLFAC;
               G = G / SCLFAC;
               CA = CA / SCLFAC;
               R = R*SCLFAC;
               RA = RA*SCLFAC;
            }

            // Now balance.

            if( ( C+R ) >= FACTOR*S ) continue;
            if ( F < ONE && SCALE( I ) < ONE ) {
               if( F*SCALE( I ) <= SFMIN1 ) continue;
            }
            if ( F > ONE && SCALE( I ) > ONE ) {
               if( SCALE( I ) >= SFMAX1 / F ) continue;
            }
            G = ONE / F;
            SCALE[I] = SCALE( I )*F;
            NOCONV = true;

            zdscal(N-K+1, G, A( I, K ), LDA );
            zdscal(L, F, A( 1, I ), 1 );

         }

      }

      ILO = K;
      IHI = L;

      }
