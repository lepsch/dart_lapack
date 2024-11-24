// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      void slasv2(final int F, final int G, final int H, final int SSMIN, final int SSMAX, final int SNR, final int CSR, final int SNL, final int CSL,) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      double               CSL, CSR, F, G, H, SNL, SNR, SSMAX, SSMIN;
      // ..

      double               ZERO;
      const              ZERO = 0.0 ;
      double               HALF;
      const              HALF = 0.5 ;
      double               ONE;
      const              ONE = 1.0 ;
      double               TWO;
      const              TWO = 2.0 ;
      double               FOUR;
      const              FOUR = 4.0 ;
      bool               GASMAL, SWAP;
      int                PMAX;
      double               A, CLT, CRT, D, FA, FT, GA, GT, HA, HT, L, M, MM, R, S, SLT, SRT, T, TEMP, TSIGN, TT;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, SIGN, SQRT
      // ..
      // .. External Functions ..
      //- REAL               SLAMCH;
      // EXTERNAL SLAMCH

      FT = F;
      FA = ( FT ).abs();
      HT = H;
      HA = ( H ).abs();

      // PMAX points to the maximum absolute element of matrix
      //   PMAX = 1 if F largest in absolute values
      //   PMAX = 2 if G largest in absolute values
      //   PMAX = 3 if H largest in absolute values

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
      GA = ( GT ).abs();
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
            S = sqrt( TT+MM );

            // Note that 1 <= S <= 1 + 1/macheps

            if ( L == ZERO ) {
               R = ( M ).abs();
            } else {
               R = sqrt( L*L+MM );
            }

            // Note that 0 <= R <= 1 + 1/macheps

            A = HALF*( S+R );

            // Note that 1 <= A <= 1 + abs(M)

            SSMIN = HA / A;
            SSMAX = FA*A;
            if ( MM == ZERO ) {

               // Note that M is very tiny

               if ( L == ZERO ) {
                  T = sign( TWO, FT )*sign( ONE, GT );
               } else {
                  T = GT / sign( D, FT ) + M / T;
               }
            } else {
               T = ( M / ( S+T )+M / ( R+L ) )*( ONE+A );
            }
            L = sqrt( T*T+FOUR );
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

      if (PMAX == 1) TSIGN = sign( ONE, CSR )*sign( ONE, CSL )*sign( ONE, F );
      if( PMAX == 2 ) TSIGN = sign( ONE, SNR )*sign( ONE, CSL )*sign( ONE, G );
      IF( PMAX == 3 ) TSIGN = sign( ONE, SNR )*sign( ONE, SNL )*sign( ONE, H );
      SSMAX = sign( SSMAX, TSIGN );
      SSMIN = sign( SSMIN, TSIGN*sign( ONE, F )*sign( ONE, H ) );
      }
