// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      void ssbtrd(final int VECT, final int UPLO, final int N, final int KD, final Matrix<double> AB_, final int LDAB, final int D, final int E, final Matrix<double> Q_, final int LDQ, final Array<double> _WORK_, final Box<int> INFO,) {
  final AB = AB_.dim();
  final Q = Q_.dim();
  final _WORK = _WORK_.dim();

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO, VECT;
      int                INFO, KD, LDAB, LDQ, N;
      double               AB( LDAB, * ), D( * ), E( * ), Q( LDQ, * ), WORK( * );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      bool               INITQ, UPPER, WANTQ;
      int                I, I2, IBL, INCA, INCX, IQAEND, IQB, IQEND, J, J1, J1END, J1INC, J2, JEND, JIN, JINC, K, KD1, KDM1, KDN, L, LAST, LEND, NQ, NR, NRT;
      double               TEMP;
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLAR2V, SLARGV, SLARTG, SLARTV, SLASET, SROT, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame

      // Test the input parameters

      INITQ = lsame( VECT, 'V' );
      WANTQ = INITQ || lsame( VECT, 'U' );
      UPPER = lsame( UPLO, 'U' );
      KD1 = KD + 1;
      KDM1 = KD - 1;
      INCX = LDAB - 1;
      IQEND = 1;

      INFO = 0;
      if ( !WANTQ && !lsame( VECT, 'N' ) ) {
         INFO = -1;
      } else if ( !UPPER && !lsame( UPLO, 'L' ) ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -3;
      } else if ( KD < 0 ) {
         INFO = -4;
      } else if ( LDAB < KD1 ) {
         INFO = -6;
      } else if ( LDQ < max( 1, N ) && WANTQ ) {
         INFO = -10;
      }
      if ( INFO != 0 ) {
         xerbla('SSBTRD', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      // Initialize Q to the unit matrix, if needed

      if (INITQ) slaset( 'Full', N, N, ZERO, ONE, Q, LDQ );

      // Wherever possible, plane rotations are generated and applied in
      // vector operations of length NR over the index set J1:J2:KD1.

      // The cosines and sines of the plane rotations are stored in the
      // arrays D and WORK.

      INCA = KD1*LDAB;
      KDN = min( N-1, KD );
      if ( UPPER ) {

         if ( KD > 1 ) {

            // Reduce to tridiagonal form, working with upper triangle

            NR = 0;
            J1 = KDN + 2;
            J2 = 1;

            for (I = 1; I <= N - 2; I++) { // 90

               // Reduce i-th row of matrix to tridiagonal form

               for (K = KDN + 1; K >= 2; K--) { // 80
                  J1 = J1 + KDN;
                  J2 = J2 + KDN;

                  if ( NR > 0 ) {

                     // generate plane rotations to annihilate nonzero
                     // elements which have been created outside the band

                     slargv(NR, AB( 1, J1-1 ), INCA, WORK( J1 ), KD1, D( J1 ), KD1 );

                     // apply rotations from the right


                     // Dependent on the the number of diagonals either
                     // SLARTV or SROT is used

                     if ( NR >= 2*KD-1 ) {
                        for (L = 1; L <= KD - 1; L++) { // 10
                           slartv(NR, AB( L+1, J1-1 ), INCA, AB( L, J1 ), INCA, D( J1 ), WORK( J1 ), KD1 );
                        } // 10

                     } else {
                        JEND = J1 + ( NR-1 )*KD1;
                        for (JINC = J1; KD1 < 0 ? JINC >= JEND : JINC <= JEND; JINC += KD1) { // 20
                           srot(KDM1, AB( 2, JINC-1 ), 1, AB( 1, JINC ), 1, D( JINC ), WORK( JINC ) );
                        } // 20
                     }
                  }


                  if ( K > 2 ) {
                     if ( K <= N-I+1 ) {

                        // generate plane rotation to annihilate a(i,i+k-1)
                        // within the band

                        slartg(AB( KD-K+3, I+K-2 ), AB( KD-K+2, I+K-1 ), D( I+K-1 ), WORK( I+K-1 ), TEMP );
                        AB[KD-K+3][I+K-2] = TEMP;

                        // apply rotation from the right

                        srot(K-3, AB( KD-K+4, I+K-2 ), 1, AB( KD-K+3, I+K-1 ), 1, D( I+K-1 ), WORK( I+K-1 ) );
                     }
                     NR = NR + 1;
                     J1 = J1 - KDN - 1;
                  }

                  // apply plane rotations from both sides to diagonal
                  // blocks

                  if (NR > 0) slar2v( NR, AB( KD1, J1-1 ), AB( KD1, J1 ), AB( KD, J1 ), INCA, D( J1 ), WORK( J1 ), KD1 );

                  // apply plane rotations from the left

                  if ( NR > 0 ) {
                     if ( 2*KD-1 < NR ) {

                     // Dependent on the the number of diagonals either
                     // SLARTV or SROT is used

                        for (L = 1; L <= KD - 1; L++) { // 30
                           if ( J2+L > N ) {
                              NRT = NR - 1;
                           } else {
                              NRT = NR;
                           }
                           if (NRT > 0) slartv( NRT, AB( KD-L, J1+L ), INCA, AB( KD-L+1, J1+L ), INCA, D( J1 ), WORK( J1 ), KD1 );
                        } // 30
                     } else {
                        J1END = J1 + KD1*( NR-2 );
                        if ( J1END >= J1 ) {
                           for (JIN = J1; KD1 < 0 ? JIN >= J1END : JIN <= J1END; JIN += KD1) { // 40
                              srot(KD-1, AB( KD-1, JIN+1 ), INCX, AB( KD, JIN+1 ), INCX, D( JIN ), WORK( JIN ) );
                           } // 40
                        }
                        LEND = min( KDM1, N-J2 );
                        LAST = J1END + KD1;
                        if (LEND > 0) srot( LEND, AB( KD-1, LAST+1 ), INCX, AB( KD, LAST+1 ), INCX, D( LAST ), WORK( LAST ) );
                     }
                  }

                  if ( WANTQ ) {

                     // accumulate product of plane rotations in Q

                     if ( INITQ ) {

                  // take advantage of the fact that Q was
                  // initially the Identity matrix

                        IQEND = max( IQEND, J2 );
                        I2 = max( 0, K-3 );
                        IQAEND = 1 + I*KD;
                        if (K == 2) IQAEND = IQAEND + KD;
                        IQAEND = min( IQAEND, IQEND );
                        for (J = J1; KD1 < 0 ? J >= J2 : J <= J2; J += KD1) { // 50
                           IBL = I - I2 / KDM1;
                           I2 = I2 + 1;
                           IQB = max( 1, J-IBL );
                           NQ = 1 + IQAEND - IQB;
                           IQAEND = min( IQAEND+KD, IQEND );
                           srot(NQ, Q( IQB, J-1 ), 1, Q( IQB, J ), 1, D( J ), WORK( J ) );
                        } // 50
                     } else {

                        for (J = J1; KD1 < 0 ? J >= J2 : J <= J2; J += KD1) { // 60
                           srot(N, Q( 1, J-1 ), 1, Q( 1, J ), 1, D( J ), WORK( J ) );
                        } // 60
                     }

                  }

                  if ( J2+KDN > N ) {

                     // adjust J2 to keep within the bounds of the matrix

                     NR = NR - 1;
                     J2 = J2 - KDN - 1;
                  }

                  for (J = J1; KD1 < 0 ? J >= J2 : J <= J2; J += KD1) { // 70

                     // create nonzero element a(j-1,j+kd) outside the band
                     // and store it in WORK

                     WORK[J+KD] = WORK( J )*AB( 1, J+KD );
                     AB[1][J+KD] = D( J )*AB( 1, J+KD );
                  } // 70
               } // 80
            } // 90
         }

         if ( KD > 0 ) {

            // copy off-diagonal elements to E

            for (I = 1; I <= N - 1; I++) { // 100
               E[I] = AB( KD, I+1 );
            } // 100
         } else {

            // set E to zero if original matrix was diagonal

            for (I = 1; I <= N - 1; I++) { // 110
               E[I] = ZERO;
            } // 110
         }

         // copy diagonal elements to D

         for (I = 1; I <= N; I++) { // 120
            D[I] = AB( KD1, I );
         } // 120

      } else {

         if ( KD > 1 ) {

            // Reduce to tridiagonal form, working with lower triangle

            NR = 0;
            J1 = KDN + 2;
            J2 = 1;

            for (I = 1; I <= N - 2; I++) { // 210

               // Reduce i-th column of matrix to tridiagonal form

               for (K = KDN + 1; K >= 2; K--) { // 200
                  J1 = J1 + KDN;
                  J2 = J2 + KDN;

                  if ( NR > 0 ) {

                     // generate plane rotations to annihilate nonzero
                     // elements which have been created outside the band

                     slargv(NR, AB( KD1, J1-KD1 ), INCA, WORK( J1 ), KD1, D( J1 ), KD1 );

                     // apply plane rotations from one side


                     // Dependent on the the number of diagonals either
                     // SLARTV or SROT is used

                     if ( NR > 2*KD-1 ) {
                        for (L = 1; L <= KD - 1; L++) { // 130
                           slartv(NR, AB( KD1-L, J1-KD1+L ), INCA, AB( KD1-L+1, J1-KD1+L ), INCA, D( J1 ), WORK( J1 ), KD1 );
                        } // 130
                     } else {
                        JEND = J1 + KD1*( NR-1 );
                        for (JINC = J1; KD1 < 0 ? JINC >= JEND : JINC <= JEND; JINC += KD1) { // 140
                           srot(KDM1, AB( KD, JINC-KD ), INCX, AB( KD1, JINC-KD ), INCX, D( JINC ), WORK( JINC ) );
                        } // 140
                     }

                  }

                  if ( K > 2 ) {
                     if ( K <= N-I+1 ) {

                        // generate plane rotation to annihilate a(i+k-1,i)
                        // within the band

                        slartg(AB( K-1, I ), AB( K, I ), D( I+K-1 ), WORK( I+K-1 ), TEMP );
                        AB[K-1][I] = TEMP;

                        // apply rotation from the left

                        srot(K-3, AB( K-2, I+1 ), LDAB-1, AB( K-1, I+1 ), LDAB-1, D( I+K-1 ), WORK( I+K-1 ) );
                     }
                     NR = NR + 1;
                     J1 = J1 - KDN - 1;
                  }

                  // apply plane rotations from both sides to diagonal
                  // blocks

                  if (NR > 0) slar2v( NR, AB( 1, J1-1 ), AB( 1, J1 ), AB( 2, J1-1 ), INCA, D( J1 ), WORK( J1 ), KD1 );

                  // apply plane rotations from the right


                     // Dependent on the the number of diagonals either
                     // SLARTV or SROT is used

                  if ( NR > 0 ) {
                     if ( NR > 2*KD-1 ) {
                        for (L = 1; L <= KD - 1; L++) { // 150
                           if ( J2+L > N ) {
                              NRT = NR - 1;
                           } else {
                              NRT = NR;
                           }
                           if (NRT > 0) slartv( NRT, AB( L+2, J1-1 ), INCA, AB( L+1, J1 ), INCA, D( J1 ), WORK( J1 ), KD1 );
                        } // 150
                     } else {
                        J1END = J1 + KD1*( NR-2 );
                        if ( J1END >= J1 ) {
                           for (J1INC = J1; KD1 < 0 ? J1INC >= J1END : J1INC <= J1END; J1INC += KD1) { // 160
                              srot(KDM1, AB( 3, J1INC-1 ), 1, AB( 2, J1INC ), 1, D( J1INC ), WORK( J1INC ) );
                           } // 160
                        }
                        LEND = min( KDM1, N-J2 );
                        LAST = J1END + KD1;
                        if (LEND > 0) srot( LEND, AB( 3, LAST-1 ), 1, AB( 2, LAST ), 1, D( LAST ), WORK( LAST ) );
                     }
                  }



                  if ( WANTQ ) {

                     // accumulate product of plane rotations in Q

                     if ( INITQ ) {

                  // take advantage of the fact that Q was
                  // initially the Identity matrix

                        IQEND = max( IQEND, J2 );
                        I2 = max( 0, K-3 );
                        IQAEND = 1 + I*KD;
                        if (K == 2) IQAEND = IQAEND + KD;
                        IQAEND = min( IQAEND, IQEND );
                        for (J = J1; KD1 < 0 ? J >= J2 : J <= J2; J += KD1) { // 170
                           IBL = I - I2 / KDM1;
                           I2 = I2 + 1;
                           IQB = max( 1, J-IBL );
                           NQ = 1 + IQAEND - IQB;
                           IQAEND = min( IQAEND+KD, IQEND );
                           srot(NQ, Q( IQB, J-1 ), 1, Q( IQB, J ), 1, D( J ), WORK( J ) );
                        } // 170
                     } else {

                        for (J = J1; KD1 < 0 ? J >= J2 : J <= J2; J += KD1) { // 180
                           srot(N, Q( 1, J-1 ), 1, Q( 1, J ), 1, D( J ), WORK( J ) );
                        } // 180
                     }
                  }

                  if ( J2+KDN > N ) {

                     // adjust J2 to keep within the bounds of the matrix

                     NR = NR - 1;
                     J2 = J2 - KDN - 1;
                  }

                  for (J = J1; KD1 < 0 ? J >= J2 : J <= J2; J += KD1) { // 190

                     // create nonzero element a(j+kd,j-1) outside the
                     // band and store it in WORK

                     WORK[J+KD] = WORK( J )*AB( KD1, J );
                     AB[KD1][J] = D( J )*AB( KD1, J );
                  } // 190
               } // 200
            } // 210
         }

         if ( KD > 0 ) {

            // copy off-diagonal elements to E

            for (I = 1; I <= N - 1; I++) { // 220
               E[I] = AB( 2, I );
            } // 220
         } else {

            // set E to zero if original matrix was diagonal

            for (I = 1; I <= N - 1; I++) { // 230
               E[I] = ZERO;
            } // 230
         }

         // copy diagonal elements to D

         for (I = 1; I <= N; I++) { // 240
            D[I] = AB( 1, I );
         } // 240
      }

      }
