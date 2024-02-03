      void cgbbrd(VECT, M, N, NCC, KL, KU, AB, LDAB, D, E, Q, LDQ, PT, LDPT, C, LDC, WORK, RWORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             VECT;
      int                INFO, KL, KU, LDAB, LDC, LDPT, LDQ, M, N, NCC;
      // ..
      // .. Array Arguments ..
      REAL               D( * ), E( * ), RWORK( * );
      COMPLEX            AB( LDAB, * ), C( LDC, * ), PT( LDPT, * ), Q( LDQ, * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ZERO;
      const              ZERO = 0.0 ;
      COMPLEX            CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      bool               WANTB, WANTC, WANTPT, WANTQ;
      int                I, INCA, J, J1, J2, KB, KB1, KK, KLM, KLU1, KUN, L, MINMN, ML, ML0, MU, MU0, NR, NRT;
      REAL               ABST, RC;
      COMPLEX            RA, RB, RS, T;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLARGV, CLARTG, CLARTV, CLASET, CROT, CSCAL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, CONJG, MAX, MIN
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
         xerbla('CGBBRD', -INFO );
         return;
      }

      // Initialize Q and P**H to the unit matrix, if needed

      if (WANTQ) claset( 'Full', M, M, CZERO, CONE, Q, LDQ );
      IF( WANTPT ) claset( 'Full', N, N, CZERO, CONE, PT, LDPT );

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

         // The complex sines of the plane rotations are stored in WORK,
         // and the real cosines in RWORK.

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

               if (NR > 0) clargv( NR, AB( KLU1, J1-KLM-1 ), INCA, WORK( J1 ), KB1, RWORK( J1 ), KB1 );

               // apply plane rotations from the left

               for (L = 1; L <= KB; L++) { // 10
                  if ( J2-KLM+L-1 > N ) {
                     NRT = NR - 1;
                  } else {
                     NRT = NR;
                  }
                  if (NRT > 0) clartv( NRT, AB( KLU1-L, J1-KLM+L-1 ), INCA, AB( KLU1-L+1, J1-KLM+L-1 ), INCA, RWORK( J1 ), WORK( J1 ), KB1 );
               } // 10

               if ( ML > ML0 ) {
                  if ( ML <= M-I+1 ) {

                     // generate plane rotation to annihilate a(i+ml-1,i)
                     // within the band, and apply rotation from the left

                     clartg(AB( KU+ML-1, I ), AB( KU+ML, I ), RWORK( I+ML-1 ), WORK( I+ML-1 ), RA );
                     AB( KU+ML-1, I ) = RA;
                     if (I < N) crot( min( KU+ML-2, N-I ), AB( KU+ML-2, I+1 ), LDAB-1, AB( KU+ML-1, I+1 ), LDAB-1, RWORK( I+ML-1 ), WORK( I+ML-1 ) );
                  }
                  NR = NR + 1;
                  J1 = J1 - KB1;
               }

               if ( WANTQ ) {

                  // accumulate product of plane rotations in Q

                  DO 20 J = J1, J2, KB1;
                     crot(M, Q( 1, J-1 ), 1, Q( 1, J ), 1, RWORK( J ), CONJG( WORK( J ) ) );
                  } // 20
               }

               if ( WANTC ) {

                  // apply plane rotations to C

                  DO 30 J = J1, J2, KB1;
                     crot(NCC, C( J-1, 1 ), LDC, C( J, 1 ), LDC, RWORK( J ), WORK( J ) );
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
                  AB( 1, J+KUN ) = RWORK( J )*AB( 1, J+KUN );
               } // 40

               // generate plane rotations to annihilate nonzero elements
               // which have been generated above the band

               if (NR > 0) clargv( NR, AB( 1, J1+KUN-1 ), INCA, WORK( J1+KUN ), KB1, RWORK( J1+KUN ), KB1 );

               // apply plane rotations from the right

               for (L = 1; L <= KB; L++) { // 50
                  if ( J2+L-1 > M ) {
                     NRT = NR - 1;
                  } else {
                     NRT = NR;
                  }
                  if (NRT > 0) clartv( NRT, AB( L+1, J1+KUN-1 ), INCA, AB( L, J1+KUN ), INCA, RWORK( J1+KUN ), WORK( J1+KUN ), KB1 );
               } // 50

               if ( ML == ML0 && MU > MU0 ) {
                  if ( MU <= N-I+1 ) {

                     // generate plane rotation to annihilate a(i,i+mu-1)
                     // within the band, and apply rotation from the right

                     clartg(AB( KU-MU+3, I+MU-2 ), AB( KU-MU+2, I+MU-1 ), RWORK( I+MU-1 ), WORK( I+MU-1 ), RA );
                     AB( KU-MU+3, I+MU-2 ) = RA;
                     crot(min( KL+MU-2, M-I ), AB( KU-MU+4, I+MU-2 ), 1, AB( KU-MU+3, I+MU-1 ), 1, RWORK( I+MU-1 ), WORK( I+MU-1 ) );
                  }
                  NR = NR + 1;
                  J1 = J1 - KB1;
               }

               if ( WANTPT ) {

                  // accumulate product of plane rotations in P**H

                  DO 60 J = J1, J2, KB1;
                     crot(N, PT( J+KUN-1, 1 ), LDPT, PT( J+KUN, 1 ), LDPT, RWORK( J+KUN ), CONJG( WORK( J+KUN ) ) );
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
                  AB( KLU1, J+KUN ) = RWORK( J+KUN )*AB( KLU1, J+KUN );
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

         // A has been reduced to complex lower bidiagonal form

         // Transform lower bidiagonal form to upper bidiagonal by applying
         // plane rotations from the left, overwriting superdiagonal
         // elements on subdiagonal elements

         DO 100 I = 1, min( M-1, N );
            clartg(AB( 1, I ), AB( 2, I ), RC, RS, RA );
            AB( 1, I ) = RA;
            if ( I < N ) {
               AB( 2, I ) = RS*AB( 1, I+1 );
               AB( 1, I+1 ) = RC*AB( 1, I+1 );
            }
            if (WANTQ) crot( M, Q( 1, I ), 1, Q( 1, I+1 ), 1, RC, CONJG( RS ) );
            IF( WANTC ) crot( NCC, C( I, 1 ), LDC, C( I+1, 1 ), LDC, RC, RS );
         } // 100
      } else {

         // A has been reduced to complex upper bidiagonal form or is
         // diagonal

         if ( KU > 0 && M < N ) {

            // Annihilate a(m,m+1) by applying plane rotations from the
            // right

            RB = AB( KU, M+1 );
            DO 110 I = M, 1, -1;
               clartg(AB( KU+1, I ), RB, RC, RS, RA );
               AB( KU+1, I ) = RA;
               if ( I > 1 ) {
                  RB = -CONJG( RS )*AB( KU, I );
                  AB( KU, I ) = RC*AB( KU, I );
               }
               if (WANTPT) crot( N, PT( I, 1 ), LDPT, PT( M+1, 1 ), LDPT, RC, CONJG( RS ) );
            } // 110
         }
      }

      // Make diagonal and superdiagonal elements real, storing them in D
      // and E

      T = AB( KU+1, 1 );
      for (I = 1; I <= MINMN; I++) { // 120
         ABST = ABS( T );
         D( I ) = ABST;
         if ( ABST != ZERO ) {
            T = T / ABST;
         } else {
            T = CONE;
         }
         if (WANTQ) cscal( M, T, Q( 1, I ), 1 );
         IF( WANTC ) cscal( NCC, CONJG( T ), C( I, 1 ), LDC );
         if ( I < MINMN ) {
            if ( KU == 0 && KL == 0 ) {
               E( I ) = ZERO;
               T = AB( 1, I+1 );
            } else {
               if ( KU == 0 ) {
                  T = AB( 2, I )*CONJG( T );
               } else {
                  T = AB( KU, I+1 )*CONJG( T );
               }
               ABST = ABS( T );
               E( I ) = ABST;
               if ( ABST != ZERO ) {
                  T = T / ABST;
               } else {
                  T = CONE;
               }
               if (WANTPT) cscal( N, T, PT( I+1, 1 ), LDPT );
               T = AB( KU+1, I+1 )*CONJG( T );
            }
         }
      } // 120
      return;
      }
