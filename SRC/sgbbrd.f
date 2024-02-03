      void sgbbrd(VECT, M, N, NCC, KL, KU, AB, LDAB, D, E, Q, LDQ, PT, LDPT, C, LDC, WORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             VECT;
      int                INFO, KL, KU, LDAB, LDC, LDPT, LDQ, M, N, NCC;
      // ..
      // .. Array Arguments ..
      REAL               AB( LDAB, * ), C( LDC, * ), D( * ), E( * ), PT( LDPT, * ), Q( LDQ, * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      bool               WANTB, WANTC, WANTPT, WANTQ;
      int                I, INCA, J, J1, J2, KB, KB1, KK, KLM, KLU1, KUN, L, MINMN, ML, ML0, MN, MU, MU0, NR, NRT;
      REAL               RA, RB, RC, RS;
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLARGV, SLARTG, SLARTV, SLASET, SROT, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. Executable Statements ..

      // Test the input parameters

      WANTB = LSAME( VECT, 'B' );
      WANTQ = LSAME( VECT, 'Q' ) || WANTB;
      WANTPT = LSAME( VECT, 'P' ) || WANTB;
      WANTC = NCC > 0;
      KLU1 = KL + KU + 1;
      INFO = 0;
      if ( !WANTQ && !WANTPT && !LSAME( VECT, 'N' ) ) {
         INFO = -1;
      } else if ( M < 0 ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -3;
      } else if ( NCC < 0 ) {
         INFO = -4;
      } else if ( KL < 0 ) {
         INFO = -5;
      } else if ( KU < 0 ) {
         INFO = -6;
      } else if ( LDAB < KLU1 ) {
         INFO = -8;
      } else if ( LDQ < 1 || WANTQ && LDQ < max( 1, M ) ) {
         INFO = -12;
      } else if ( LDPT < 1 || WANTPT && LDPT < max( 1, N ) ) {
         INFO = -14;
      } else if ( LDC < 1 || WANTC && LDC < max( 1, M ) ) {
         INFO = -16;
      }
      if ( INFO != 0 ) {
         xerbla('SGBBRD', -INFO );
         return;
      }

      // Initialize Q and P**T to the unit matrix, if needed

      if (WANTQ) CALL SLASET( 'Full', M, M, ZERO, ONE, Q, LDQ );
      IF( WANTPT ) CALL SLASET( 'Full', N, N, ZERO, ONE, PT, LDPT );

      // Quick return if possible.

      if (M == 0 || N == 0) return;

      MINMN = min( M, N );

      if ( KL+KU > 1 ) {

         // Reduce to upper bidiagonal form if KU > 0; if KU = 0, reduce
         // first to lower bidiagonal form and then transform to upper
         // bidiagonal

         if ( KU > 0 ) {
            ML0 = 1;
            MU0 = 2;
         } else {
            ML0 = 2;
            MU0 = 1;
         }

         // Wherever possible, plane rotations are generated and applied in
         // vector operations of length NR over the index set J1:J2:KLU1.

         // The sines of the plane rotations are stored in WORK(1:max(m,n))
         // and the cosines in WORK(max(m,n)+1:2*max(m,n)).

         MN = max( M, N );
         KLM = min( M-1, KL );
         KUN = min( N-1, KU );
         KB = KLM + KUN;
         KB1 = KB + 1;
         INCA = KB1*LDAB;
         NR = 0;
         J1 = KLM + 2;
         J2 = 1 - KUN;

         for (I = 1; I <= MINMN; I++) { // 90

            // Reduce i-th column and i-th row of matrix to bidiagonal form

            ML = KLM + 1;
            MU = KUN + 1;
            for (KK = 1; KK <= KB; KK++) { // 80
               J1 = J1 + KB;
               J2 = J2 + KB;

               // generate plane rotations to annihilate nonzero elements
               // which have been created below the band

               if (NR > 0) CALL SLARGV( NR, AB( KLU1, J1-KLM-1 ), INCA, WORK( J1 ), KB1, WORK( MN+J1 ), KB1 );

               // apply plane rotations from the left

               for (L = 1; L <= KB; L++) { // 10
                  if ( J2-KLM+L-1 > N ) {
                     NRT = NR - 1;
                  } else {
                     NRT = NR;
                  }
                  if (NRT > 0) CALL SLARTV( NRT, AB( KLU1-L, J1-KLM+L-1 ), INCA, AB( KLU1-L+1, J1-KLM+L-1 ), INCA, WORK( MN+J1 ), WORK( J1 ), KB1 );
               } // 10

               if ( ML > ML0 ) {
                  if ( ML <= M-I+1 ) {

                     // generate plane rotation to annihilate a(i+ml-1,i)
                     // within the band, and apply rotation from the left

                     slartg(AB( KU+ML-1, I ), AB( KU+ML, I ), WORK( MN+I+ML-1 ), WORK( I+ML-1 ), RA );
                     AB( KU+ML-1, I ) = RA;
                     if (I < N) CALL SROT( min( KU+ML-2, N-I ), AB( KU+ML-2, I+1 ), LDAB-1, AB( KU+ML-1, I+1 ), LDAB-1, WORK( MN+I+ML-1 ), WORK( I+ML-1 ) );
                  }
                  NR = NR + 1;
                  J1 = J1 - KB1;
               }

               if ( WANTQ ) {

                  // accumulate product of plane rotations in Q

                  DO 20 J = J1, J2, KB1;
                     srot(M, Q( 1, J-1 ), 1, Q( 1, J ), 1, WORK( MN+J ), WORK( J ) );
                  } // 20
               }

               if ( WANTC ) {

                  // apply plane rotations to C

                  DO 30 J = J1, J2, KB1;
                     srot(NCC, C( J-1, 1 ), LDC, C( J, 1 ), LDC, WORK( MN+J ), WORK( J ) );
                  } // 30
               }

               if ( J2+KUN > N ) {

                  // adjust J2 to keep within the bounds of the matrix

                  NR = NR - 1;
                  J2 = J2 - KB1;
               }

               DO 40 J = J1, J2, KB1;

                  // create nonzero element a(j-1,j+ku) above the band
                  // and store it in WORK(n+1:2*n)

                  WORK( J+KUN ) = WORK( J )*AB( 1, J+KUN );
                  AB( 1, J+KUN ) = WORK( MN+J )*AB( 1, J+KUN );
               } // 40

               // generate plane rotations to annihilate nonzero elements
               // which have been generated above the band

               if (NR > 0) CALL SLARGV( NR, AB( 1, J1+KUN-1 ), INCA, WORK( J1+KUN ), KB1, WORK( MN+J1+KUN ), KB1 );

               // apply plane rotations from the right

               for (L = 1; L <= KB; L++) { // 50
                  if ( J2+L-1 > M ) {
                     NRT = NR - 1;
                  } else {
                     NRT = NR;
                  }
                  if (NRT > 0) CALL SLARTV( NRT, AB( L+1, J1+KUN-1 ), INCA, AB( L, J1+KUN ), INCA, WORK( MN+J1+KUN ), WORK( J1+KUN ), KB1 );
               } // 50

               if ( ML == ML0 && MU > MU0 ) {
                  if ( MU <= N-I+1 ) {

                     // generate plane rotation to annihilate a(i,i+mu-1)
                     // within the band, and apply rotation from the right

                     slartg(AB( KU-MU+3, I+MU-2 ), AB( KU-MU+2, I+MU-1 ), WORK( MN+I+MU-1 ), WORK( I+MU-1 ), RA );
                     AB( KU-MU+3, I+MU-2 ) = RA;
                     srot(min( KL+MU-2, M-I ), AB( KU-MU+4, I+MU-2 ), 1, AB( KU-MU+3, I+MU-1 ), 1, WORK( MN+I+MU-1 ), WORK( I+MU-1 ) );
                  }
                  NR = NR + 1;
                  J1 = J1 - KB1;
               }

               if ( WANTPT ) {

                  // accumulate product of plane rotations in P**T

                  DO 60 J = J1, J2, KB1;
                     srot(N, PT( J+KUN-1, 1 ), LDPT, PT( J+KUN, 1 ), LDPT, WORK( MN+J+KUN ), WORK( J+KUN ) );
                  } // 60
               }

               if ( J2+KB > M ) {

                  // adjust J2 to keep within the bounds of the matrix

                  NR = NR - 1;
                  J2 = J2 - KB1;
               }

               DO 70 J = J1, J2, KB1;

                  // create nonzero element a(j+kl+ku,j+ku-1) below the
                  // band and store it in WORK(1:n)

                  WORK( J+KB ) = WORK( J+KUN )*AB( KLU1, J+KUN );
                  AB( KLU1, J+KUN ) = WORK( MN+J+KUN )*AB( KLU1, J+KUN );
               } // 70

               if ( ML > ML0 ) {
                  ML = ML - 1;
               } else {
                  MU = MU - 1;
               }
            } // 80
         } // 90
      }

      if ( KU == 0 && KL > 0 ) {

         // A has been reduced to lower bidiagonal form

         // Transform lower bidiagonal form to upper bidiagonal by applying
         // plane rotations from the left, storing diagonal elements in D
         // and off-diagonal elements in E

         DO 100 I = 1, min( M-1, N );
            slartg(AB( 1, I ), AB( 2, I ), RC, RS, RA );
            D( I ) = RA;
            if ( I < N ) {
               E( I ) = RS*AB( 1, I+1 );
               AB( 1, I+1 ) = RC*AB( 1, I+1 );
            }
            if (WANTQ) CALL SROT( M, Q( 1, I ), 1, Q( 1, I+1 ), 1, RC, RS );
            IF( WANTC ) CALL SROT( NCC, C( I, 1 ), LDC, C( I+1, 1 ), LDC, RC, RS );
         } // 100
         if (M <= N) D( M ) = AB( 1, M );
      } else if ( KU > 0 ) {

         // A has been reduced to upper bidiagonal form

         if ( M < N ) {

            // Annihilate a(m,m+1) by applying plane rotations from the
            // right, storing diagonal elements in D and off-diagonal
            // elements in E

            RB = AB( KU, M+1 );
            DO 110 I = M, 1, -1;
               slartg(AB( KU+1, I ), RB, RC, RS, RA );
               D( I ) = RA;
               if ( I > 1 ) {
                  RB = -RS*AB( KU, I );
                  E( I-1 ) = RC*AB( KU, I );
               }
               if (WANTPT) CALL SROT( N, PT( I, 1 ), LDPT, PT( M+1, 1 ), LDPT, RC, RS );
            } // 110
         } else {

            // Copy off-diagonal elements to E and diagonal elements to D

            for (I = 1; I <= MINMN - 1; I++) { // 120
               E( I ) = AB( KU, I+1 );
            } // 120
            for (I = 1; I <= MINMN; I++) { // 130
               D( I ) = AB( KU+1, I );
            } // 130
         }
      } else {

         // A is diagonal. Set elements of E to zero and copy diagonal
         // elements to D.

         for (I = 1; I <= MINMN - 1; I++) { // 140
            E( I ) = ZERO;
         } // 140
         for (I = 1; I <= MINMN; I++) { // 150
            D( I ) = AB( 1, I );
         } // 150
      }
      return;

      // End of SGBBRD

      }
