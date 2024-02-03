      SUBROUTINE SGEBAL( JOB, N, A, LDA, ILO, IHI, SCALE, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOB;
      int                IHI, ILO, INFO, LDA, N;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * ), SCALE( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      REAL               SCLFAC
      const              SCLFAC = 2.0E+0 ;
      REAL               FACTOR
      const              FACTOR = 0.95E+0 ;
      // ..
      // .. Local Scalars ..
      bool               NOCONV, CANSWAP;
      int                I, ICA, IRA, J, K, L;
      REAL               C, CA, F, G, R, RA, S, SFMAX1, SFMAX2, SFMIN1, SFMIN2
      // ..
      // .. External Functions ..
      bool               SISNAN, LSAME;
      int                ISAMAX;
      REAL               SLAMCH, SNRM2
      // EXTERNAL SISNAN, LSAME, ISAMAX, SLAMCH, SNRM2
      // ..
      // .. External Subroutines ..
      // EXTERNAL SSCAL, SSWAP, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN
      // ..
      // Test the input parameters

      INFO = 0
      if ( .NOT.LSAME( JOB, 'N' ) .AND. .NOT.LSAME( JOB, 'P' ) .AND. .NOT.LSAME( JOB, 'S' ) .AND. .NOT.LSAME( JOB, 'B' ) ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -4
      }
      if ( INFO.NE.0 ) {
         xerbla('SGEBAL', -INFO );
         RETURN
      }

      // Quick returns.

      if ( N == 0 ) {
         ILO = 1
         IHI = 0
         RETURN
      }

      if ( LSAME( JOB, 'N' ) ) {
         for (I = 1; I <= N; I++) {
            SCALE( I ) = ONE
         }
         ILO = 1
         IHI = N
         RETURN
      }

      // Permutation to isolate eigenvalues if possible.

      K = 1
      L = N

      if ( .NOT.LSAME( JOB, 'S' ) ) {

         // Row and column exchange.

         NOCONV = true;
         DO WHILE( NOCONV )

            // Search for rows isolating an eigenvalue and push them down.

            NOCONV = false;
            DO I = L, 1, -1
               CANSWAP = true;
               for (J = 1; J <= L; J++) {
                  if ( I.NE.J .AND. A( I, J ).NE.ZERO ) {
                     CANSWAP = false;
                     EXIT
                  }
               }

               if ( CANSWAP ) {
                  SCALE( L ) = I
                  if ( I.NE.L ) {
                     sswap(L, A( 1, I ), 1, A( 1, L ), 1 );
                     sswap(N-K+1, A( I, K ), LDA, A( L, K ), LDA );
                  }
                  NOCONV = true;

                  if ( L == 1 ) {
                     ILO = 1
                     IHI = 1
                     RETURN
                  }

                  L = L - 1
               }
            }

         }

         NOCONV = true;
         DO WHILE( NOCONV )

            // Search for columns isolating an eigenvalue and push them left.

            NOCONV = false;
            for (J = K; J <= L; J++) {
               CANSWAP = true;
               for (I = K; I <= L; I++) {
                  if ( I.NE.J .AND. A( I, J ).NE.ZERO ) {
                     CANSWAP = false;
                     EXIT
                  }
               }

               if ( CANSWAP ) {
                  SCALE( K ) = J
                  if ( J.NE.K ) {
                     sswap(L, A( 1, J ), 1, A( 1, K ), 1 );
                     sswap(N-K+1, A( J, K ), LDA, A( K, K ), LDA );
                  }
                  NOCONV = true;

                  K = K + 1
               }
            }

         }

      }

      // Initialize SCALE for non-permuted submatrix.

      for (I = K; I <= L; I++) {
         SCALE( I ) = ONE
      }

      // If we only had to permute, we are done.

      if ( LSAME( JOB, 'P' ) ) {
         ILO = K
         IHI = L
         RETURN
      }

      // Balance the submatrix in rows K to L.

      // Iterative loop for norm reduction.

      SFMIN1 = SLAMCH( 'S' ) / SLAMCH( 'P' )
      SFMAX1 = ONE / SFMIN1
      SFMIN2 = SFMIN1*SCLFAC
      SFMAX2 = ONE / SFMIN2

      NOCONV = true;
      DO WHILE( NOCONV )
         NOCONV = false;

         for (I = K; I <= L; I++) {

            C = SNRM2( L-K+1, A( K, I ), 1 )
            R = SNRM2( L-K+1, A( I, K ), LDA )
            ICA = ISAMAX( L, A( 1, I ), 1 )
            CA = ABS( A( ICA, I ) )
            IRA = ISAMAX( N-K+1, A( I, K ), LDA )
            RA = ABS( A( I, IRA+K-1 ) )

            // Guard against zero C or R due to underflow.

            if (C == ZERO .OR. R == ZERO) CYCLE;

            // Exit if NaN to avoid infinite loop

            if ( SISNAN( C+CA+R+RA ) ) {
               INFO = -3
               xerbla('SGEBAL', -INFO );
               RETURN
            }

            G = R / SCLFAC
            F = ONE
            S = C + R

            DO WHILE( C.LT.G .AND. MAX( F, C, CA ).LT.SFMAX2 .AND. MIN( R, G, RA ).GT.SFMIN2 )
               F = F*SCLFAC
               C = C*SCLFAC
               CA = CA*SCLFAC
               R = R / SCLFAC
               G = G / SCLFAC
               RA = RA / SCLFAC
            }

            G = C / SCLFAC

            DO WHILE( G.GE.R .AND. MAX( R, RA ).LT.SFMAX2 .AND. MIN( F, C, G, CA ).GT.SFMIN2 )
               F = F / SCLFAC
               C = C / SCLFAC
               G = G / SCLFAC
               CA = CA / SCLFAC
               R = R*SCLFAC
               RA = RA*SCLFAC
            }

            // Now balance.

            IF( ( C+R ).GE.FACTOR*S ) CYCLE
            if ( F.LT.ONE .AND. SCALE( I ).LT.ONE ) {
               IF( F*SCALE( I ).LE.SFMIN1 ) CYCLE
            }
            if ( F.GT.ONE .AND. SCALE( I ).GT.ONE ) {
               IF( SCALE( I ).GE.SFMAX1 / F ) CYCLE
            }
            G = ONE / F
            SCALE( I ) = SCALE( I )*F
            NOCONV = true;

            sscal(N-K+1, G, A( I, K ), LDA );
            sscal(L, F, A( 1, I ), 1 );

         }

      }

      ILO = K
      IHI = L

      RETURN

      // End of SGEBAL

      }
