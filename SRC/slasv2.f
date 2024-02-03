      SUBROUTINE SLASV2( F, G, H, SSMIN, SSMAX, SNR, CSR, SNL, CSL );

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      REAL               CSL, CSR, F, G, H, SNL, SNR, SSMAX, SSMIN;
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ZERO;
      const              ZERO = 0.0 ;
      REAL               HALF;
      const              HALF = 0.5 ;
      REAL               ONE;
      const              ONE = 1.0 ;
      REAL               TWO;
      const              TWO = 2.0 ;
      REAL               FOUR;
      const              FOUR = 4.0 ;
      // ..
      // .. Local Scalars ..
      bool               GASMAL, SWAP;
      int                PMAX;
      REAL               A, CLT, CRT, D, FA, FT, GA, GT, HA, HT, L, M, MM, R, S, SLT, SRT, T, TEMP, TSIGN, TT;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, SIGN, SQRT
      // ..
      // .. External Functions ..
      REAL               SLAMCH;
      // EXTERNAL SLAMCH
      // ..
      // .. Executable Statements ..

      FT = F;
      FA = ABS( FT );
      HT = H;
      HA = ABS( H );

      // PMAX points to the maximum absolute element of matrix
        // PMAX = 1 if F largest in absolute values
        // PMAX = 2 if G largest in absolute values
        // PMAX = 3 if H largest in absolute values

      PMAX = 1;
      SWAP = ( HA > FA );
      if ( SWAP ) {
         PMAX = 3;
         TEMP = FT;
         FT = HT;
         HT = TEMP;
         TEMP = FA;
         FA = HA;
         HA = TEMP;

         // Now FA >= HA

      }
      GT = G;
      GA = ABS( GT );
      if ( GA == ZERO ) {

         // Diagonal matrix

         SSMIN = HA;
         SSMAX = FA;
         CLT = ONE;
         CRT = ONE;
         SLT = ZERO;
         SRT = ZERO;
      } else {
         GASMAL = true;
         if ( GA > FA ) {
            PMAX = 2;
            if ( ( FA / GA ) < SLAMCH( 'EPS' ) ) {

               // Case of very large GA

               GASMAL = false;
               SSMAX = GA;
               if ( HA > ONE ) {
                  SSMIN = FA / ( GA / HA );
               } else {
                  SSMIN = ( FA / GA )*HA;
               }
               CLT = ONE;
               SLT = HT / GT;
               SRT = ONE;
               CRT = FT / GT;
            }
         }
         if ( GASMAL ) {

            // Normal case

            D = FA - HA;
            if ( D == FA ) {

               // Copes with infinite F or H

               L = ONE;
            } else {
               L = D / FA;
            }

            // Note that 0 <= L <= 1

            M = GT / FT;

            // Note that abs(M) <= 1/macheps

            T = TWO - L;

            // Note that T >= 1

            MM = M*M;
            TT = T*T;
            S = SQRT( TT+MM );

            // Note that 1 <= S <= 1 + 1/macheps

            if ( L == ZERO ) {
               R = ABS( M );
            } else {
               R = SQRT( L*L+MM );
            }

            // Note that 0 <= R <= 1 + 1/macheps

            A = HALF*( S+R );

            // Note that 1 <= A <= 1 + abs(M)

            SSMIN = HA / A;
            SSMAX = FA*A;
            if ( MM == ZERO ) {

               // Note that M is very tiny

               if ( L == ZERO ) {
                  T = SIGN( TWO, FT )*SIGN( ONE, GT );
               } else {
                  T = GT / SIGN( D, FT ) + M / T;
               }
            } else {
               T = ( M / ( S+T )+M / ( R+L ) )*( ONE+A );
            }
            L = SQRT( T*T+FOUR );
            CRT = TWO / L;
            SRT = T / L;
            CLT = ( CRT+SRT*M ) / A;
            SLT = ( HT / FT )*SRT / A;
         }
      }
      if ( SWAP ) {
         CSL = SRT;
         SNL = CRT;
         CSR = SLT;
         SNR = CLT;
      } else {
         CSL = CLT;
         SNL = SLT;
         CSR = CRT;
         SNR = SRT;
      }

      // Correct signs of SSMAX and SSMIN

      if (PMAX == 1) TSIGN = SIGN( ONE, CSR )*SIGN( ONE, CSL )*SIGN( ONE, F )       IF( PMAX == 2 ) TSIGN = SIGN( ONE, SNR )*SIGN( ONE, CSL )*SIGN( ONE, G )       IF( PMAX == 3 ) TSIGN = SIGN( ONE, SNR )*SIGN( ONE, SNL )*SIGN( ONE, H );
      SSMAX = SIGN( SSMAX, TSIGN );
      SSMIN = SIGN( SSMIN, TSIGN*SIGN( ONE, F )*SIGN( ONE, H ) );
      return;

      // End of SLASV2

      }
