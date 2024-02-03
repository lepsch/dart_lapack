      SUBROUTINE DGEBAL( JOB, N, A, LDA, ILO, IHI, SCALE, INFO );

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOB;
      int                IHI, ILO, INFO, LDA, N;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), SCALE( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      double             SCLFAC;
      const              SCLFAC = 2.0 ;
      double             FACTOR;
      const              FACTOR = 0.95 ;
      // ..
      // .. Local Scalars ..
      bool               NOCONV, CANSWAP;
      int                I, ICA, IRA, J, K, L;
      double             C, CA, F, G, R, RA, S, SFMAX1, SFMAX2, SFMIN1, SFMIN2;
      // ..
      // .. External Functions ..
      bool               DISNAN, LSAME;
      int                IDAMAX;
      double             DLAMCH, DNRM2;
      // EXTERNAL DISNAN, LSAME, IDAMAX, DLAMCH, DNRM2
      // ..
      // .. External Subroutines ..
      // EXTERNAL DSCAL, DSWAP, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN
      // ..
      // Test the input parameters

      INFO = 0;
      if ( !LSAME( JOB, 'N' ) && !LSAME( JOB, 'P' ) && !LSAME( JOB, 'S' ) && !LSAME( JOB, 'B' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( LDA < MAX( 1, N ) ) {
         INFO = -4;
      }
      if ( INFO != 0 ) {
         xerbla('DGEBAL', -INFO );
         return;
      }

      // Quick returns.

      if ( N == 0 ) {
         ILO = 1;
         IHI = 0;
         return;
      }

      if ( LSAME( JOB, 'N' ) ) {
         for (I = 1; I <= N; I++) {
            SCALE( I ) = ONE;
         }
         ILO = 1;
         IHI = N;
         return;
      }

      // Permutation to isolate eigenvalues if possible.

      K = 1;
      L = N;

      if ( !LSAME( JOB, 'S' ) ) {

         // Row and column exchange.

         NOCONV = true;
         DO WHILE( NOCONV );

            // Search for rows isolating an eigenvalue and push them down.

            NOCONV = false;
            DO I = L, 1, -1;
               CANSWAP = true;
               for (J = 1; J <= L; J++) {
                  if ( I != J && A( I, J ) != ZERO ) {
                     CANSWAP = false;
                     EXIT;
                  }
               }

               if ( CANSWAP ) {
                  SCALE( L ) = I;
                  if ( I != L ) {
                     dswap(L, A( 1, I ), 1, A( 1, L ), 1 );
                     dswap(N-K+1, A( I, K ), LDA, A( L, K ), LDA );
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
         DO WHILE( NOCONV );

            // Search for columns isolating an eigenvalue and push them left.

            NOCONV = false;
            for (J = K; J <= L; J++) {
               CANSWAP = true;
               for (I = K; I <= L; I++) {
                  if ( I != J && A( I, J ) != ZERO ) {
                     CANSWAP = false;
                     EXIT;
                  }
               }

               if ( CANSWAP ) {
                  SCALE( K ) = J;
                  if ( J != K ) {
                     dswap(L, A( 1, J ), 1, A( 1, K ), 1 );
                     dswap(N-K+1, A( J, K ), LDA, A( K, K ), LDA );
                  }
                  NOCONV = true;

                  K = K + 1;
               }
            }

         }

      }

      // Initialize SCALE for non-permuted submatrix.

      for (I = K; I <= L; I++) {
         SCALE( I ) = ONE;
      }

      // If we only had to permute, we are done.

      if ( LSAME( JOB, 'P' ) ) {
         ILO = K;
         IHI = L;
         return;
      }

      // Balance the submatrix in rows K to L.

      // Iterative loop for norm reduction.

      SFMIN1 = DLAMCH( 'S' ) / DLAMCH( 'P' );
      SFMAX1 = ONE / SFMIN1;
      SFMIN2 = SFMIN1*SCLFAC;
      SFMAX2 = ONE / SFMIN2;

      NOCONV = true;
      DO WHILE( NOCONV );
         NOCONV = false;

         for (I = K; I <= L; I++) {

            C = DNRM2( L-K+1, A( K, I ), 1 );
            R = DNRM2( L-K+1, A( I, K ), LDA );
            ICA = IDAMAX( L, A( 1, I ), 1 );
            CA = ABS( A( ICA, I ) );
            IRA = IDAMAX( N-K+1, A( I, K ), LDA );
            RA = ABS( A( I, IRA+K-1 ) );

            // Guard against zero C or R due to underflow.

            if (C == ZERO || R == ZERO) CYCLE;

            // Exit if NaN to avoid infinite loop

            if ( DISNAN( C+CA+R+RA ) ) {
               INFO = -3;
               xerbla('DGEBAL', -INFO );
               return;
            }

            G = R / SCLFAC;
            F = ONE;
            S = C + R;

            DO WHILE( C < G && MAX( F, C, CA ) < SFMAX2 && MIN( R, G, RA ) > SFMIN2 );
               F = F*SCLFAC;
               C = C*SCLFAC;
               CA = CA*SCLFAC;
               R = R / SCLFAC;
               G = G / SCLFAC;
               RA = RA / SCLFAC;
            }

            G = C / SCLFAC;

            DO WHILE( G >= R && MAX( R, RA ) < SFMAX2 && MIN( F, C, G, CA ) > SFMIN2 );
               F = F / SCLFAC;
               C = C / SCLFAC;
               G = G / SCLFAC;
               CA = CA / SCLFAC;
               R = R*SCLFAC;
               RA = RA*SCLFAC;
            }

            // Now balance.

            IF( ( C+R ) >= FACTOR*S ) CYCLE;
            if ( F < ONE && SCALE( I ) < ONE ) {
               IF( F*SCALE( I ) <= SFMIN1 ) CYCLE;
            }
            if ( F > ONE && SCALE( I ) > ONE ) {
               IF( SCALE( I ) >= SFMAX1 / F ) CYCLE;
            }
            G = ONE / F;
            SCALE( I ) = SCALE( I )*F;
            NOCONV = true;

            dscal(N-K+1, G, A( I, K ), LDA );
            dscal(L, F, A( 1, I ), 1 );

         }

      }

      ILO = K;
      IHI = L;

      return;

      // End of DGEBAL

      }
