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
         CALL XERBLA( 'SGEBAL', -INFO )
         RETURN
      }

      // Quick returns.

      if ( N.EQ.0 ) {
         ILO = 1
         IHI = 0
         RETURN
      }

      if ( LSAME( JOB, 'N' ) ) {
         DO I = 1, N
            SCALE( I ) = ONE
         END DO
         ILO = 1
         IHI = N
         RETURN
      }

      // Permutation to isolate eigenvalues if possible.

      K = 1
      L = N

      if ( .NOT.LSAME( JOB, 'S' ) ) {

         // Row and column exchange.

         NOCONV = .TRUE.
         DO WHILE( NOCONV )

            // Search for rows isolating an eigenvalue and push them down.

            NOCONV = .FALSE.
            DO I = L, 1, -1
               CANSWAP = .TRUE.
               DO J = 1, L
                  if ( I.NE.J .AND. A( I, J ).NE.ZERO ) {
                     CANSWAP = .FALSE.
                     EXIT
                  }
               END DO

               if ( CANSWAP ) {
                  SCALE( L ) = I
                  if ( I.NE.L ) {
                     CALL SSWAP( L, A( 1, I ), 1, A( 1, L ), 1 )
                     CALL SSWAP( N-K+1, A( I, K ), LDA, A( L, K ), LDA )
                  }
                  NOCONV = .TRUE.

                  if ( L.EQ.1 ) {
                     ILO = 1
                     IHI = 1
                     RETURN
                  }

                  L = L - 1
               }
            END DO

         END DO

         NOCONV = .TRUE.
         DO WHILE( NOCONV )

            // Search for columns isolating an eigenvalue and push them left.

            NOCONV = .FALSE.
            DO J = K, L
               CANSWAP = .TRUE.
               DO I = K, L
                  if ( I.NE.J .AND. A( I, J ).NE.ZERO ) {
                     CANSWAP = .FALSE.
                     EXIT
                  }
               END DO

               if ( CANSWAP ) {
                  SCALE( K ) = J
                  if ( J.NE.K ) {
                     CALL SSWAP( L, A( 1, J ), 1, A( 1, K ), 1 )
                     CALL SSWAP( N-K+1, A( J, K ), LDA, A( K, K ), LDA )
                  }
                  NOCONV = .TRUE.

                  K = K + 1
               }
            END DO

         END DO

      }

      // Initialize SCALE for non-permuted submatrix.

      DO I = K, L
         SCALE( I ) = ONE
      END DO

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

      NOCONV = .TRUE.
      DO WHILE( NOCONV )
         NOCONV = .FALSE.

         DO I = K, L

            C = SNRM2( L-K+1, A( K, I ), 1 )
            R = SNRM2( L-K+1, A( I, K ), LDA )
            ICA = ISAMAX( L, A( 1, I ), 1 )
            CA = ABS( A( ICA, I ) )
            IRA = ISAMAX( N-K+1, A( I, K ), LDA )
            RA = ABS( A( I, IRA+K-1 ) )

            // Guard against zero C or R due to underflow.

            IF( C.EQ.ZERO .OR. R.EQ.ZERO ) CYCLE

            // Exit if NaN to avoid infinite loop

            if ( SISNAN( C+CA+R+RA ) ) {
               INFO = -3
               CALL XERBLA( 'SGEBAL', -INFO )
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
            END DO

            G = C / SCLFAC

            DO WHILE( G.GE.R .AND. MAX( R, RA ).LT.SFMAX2 .AND. MIN( F, C, G, CA ).GT.SFMIN2 )
               F = F / SCLFAC
               C = C / SCLFAC
               G = G / SCLFAC
               CA = CA / SCLFAC
               R = R*SCLFAC
               RA = RA*SCLFAC
            END DO

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
            NOCONV = .TRUE.

            CALL SSCAL( N-K+1, G, A( I, K ), LDA )
            CALL SSCAL( L, F, A( 1, I ), 1 )

         END DO

      END DO

      ILO = K
      IHI = L

      RETURN

      // End of SGEBAL

      }
