      void zhbtrd(VECT, UPLO, N, KD, AB, LDAB, D, E, Q, LDQ, WORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO, VECT;
      int                INFO, KD, LDAB, LDQ, N;
      // ..
      // .. Array Arguments ..
      double             D( * ), E( * );
      Complex         AB( LDAB, * ), Q( LDQ, * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO;
      const              ZERO = 0.0 ;
      Complex         CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      bool               INITQ, UPPER, WANTQ;
      int                I, I2, IBL, INCA, INCX, IQAEND, IQB, IQEND, J, J1, J1END, J1INC, J2, JEND, JIN, JINC, K, KD1, KDM1, KDN, L, LAST, LEND, NQ, NR, NRT;
      double             ABST;
      Complex         T, TEMP;
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZLACGV, ZLAR2V, ZLARGV, ZLARTG, ZLARTV, ZLASET, ZROT, ZSCAL
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DCONJG, MAX, MIN
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. Executable Statements ..

      // Test the input parameters

      INITQ = LSAME( VECT, 'V' );
      WANTQ = INITQ || LSAME( VECT, 'U' );
      UPPER = LSAME( UPLO, 'U' );
      KD1 = KD + 1;
      KDM1 = KD - 1;
      INCX = LDAB - 1;
      IQEND = 1;

      INFO = 0;
      if ( !WANTQ && !LSAME( VECT, 'N' ) ) {
         INFO = -1;
      } else if ( !UPPER && !LSAME( UPLO, 'L' ) ) {
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
         xerbla('ZHBTRD', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      // Initialize Q to the unit matrix, if needed

      if (INITQ) CALL ZLASET( 'Full', N, N, CZERO, CONE, Q, LDQ );

      // Wherever possible, plane rotations are generated and applied in
      // vector operations of length NR over the index set J1:J2:KD1.

      // The real cosines and complex sines of the plane rotations are
      // stored in the arrays D and WORK.

      INCA = KD1*LDAB;
      KDN = min( N-1, KD );
      if ( UPPER ) {

         if ( KD > 1 ) {

            // Reduce to complex Hermitian tridiagonal form, working with
            // the upper triangle

            NR = 0;
            J1 = KDN + 2;
            J2 = 1;

            AB( KD1, 1 ) = DBLE( AB( KD1, 1 ) );
            for (I = 1; I <= N - 2; I++) { // 90

               // Reduce i-th row of matrix to tridiagonal form

               DO 80 K = KDN + 1, 2, -1;
                  J1 = J1 + KDN;
                  J2 = J2 + KDN;

                  if ( NR > 0 ) {

                     // generate plane rotations to annihilate nonzero
                     // elements which have been created outside the band

                     zlargv(NR, AB( 1, J1-1 ), INCA, WORK( J1 ), KD1, D( J1 ), KD1 );

                     // apply rotations from the right


                     // Dependent on the the number of diagonals either
                     // ZLARTV or ZROT is used

                     if ( NR >= 2*KD-1 ) {
                        for (L = 1; L <= KD - 1; L++) { // 10
                           zlartv(NR, AB( L+1, J1-1 ), INCA, AB( L, J1 ), INCA, D( J1 ), WORK( J1 ), KD1 );
                        } // 10

                     } else {
                        JEND = J1 + ( NR-1 )*KD1;
                        DO 20 JINC = J1, JEND, KD1;
                           zrot(KDM1, AB( 2, JINC-1 ), 1, AB( 1, JINC ), 1, D( JINC ), WORK( JINC ) );
                        } // 20
                     }
                  }


                  if ( K > 2 ) {
                     if ( K <= N-I+1 ) {

                        // generate plane rotation to annihilate a(i,i+k-1)
                        // within the band

                        zlartg(AB( KD-K+3, I+K-2 ), AB( KD-K+2, I+K-1 ), D( I+K-1 ), WORK( I+K-1 ), TEMP );
                        AB( KD-K+3, I+K-2 ) = TEMP;

                        // apply rotation from the right

                        zrot(K-3, AB( KD-K+4, I+K-2 ), 1, AB( KD-K+3, I+K-1 ), 1, D( I+K-1 ), WORK( I+K-1 ) );
                     }
                     NR = NR + 1;
                     J1 = J1 - KDN - 1;
                  }

                  // apply plane rotations from both sides to diagonal
                  // blocks

                  if (NR > 0) CALL ZLAR2V( NR, AB( KD1, J1-1 ), AB( KD1, J1 ), AB( KD, J1 ), INCA, D( J1 ), WORK( J1 ), KD1 );

                  // apply plane rotations from the left

                  if ( NR > 0 ) {
                     zlacgv(NR, WORK( J1 ), KD1 );
                     if ( 2*KD-1 < NR ) {

                     // Dependent on the the number of diagonals either
                     // ZLARTV or ZROT is used

                        for (L = 1; L <= KD - 1; L++) { // 30
                           if ( J2+L > N ) {
                              NRT = NR - 1;
                           } else {
                              NRT = NR;
                           }
                           if (NRT > 0) CALL ZLARTV( NRT, AB( KD-L, J1+L ), INCA, AB( KD-L+1, J1+L ), INCA, D( J1 ), WORK( J1 ), KD1 );
                        } // 30
                     } else {
                        J1END = J1 + KD1*( NR-2 );
                        if ( J1END >= J1 ) {
                           DO 40 JIN = J1, J1END, KD1;
                              zrot(KD-1, AB( KD-1, JIN+1 ), INCX, AB( KD, JIN+1 ), INCX, D( JIN ), WORK( JIN ) );
                           } // 40
                        }
                        LEND = min( KDM1, N-J2 );
                        LAST = J1END + KD1;
                        if (LEND > 0) CALL ZROT( LEND, AB( KD-1, LAST+1 ), INCX, AB( KD, LAST+1 ), INCX, D( LAST ), WORK( LAST ) );
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
                        DO 50 J = J1, J2, KD1;
                           IBL = I - I2 / KDM1;
                           I2 = I2 + 1;
                           IQB = max( 1, J-IBL );
                           NQ = 1 + IQAEND - IQB;
                           IQAEND = min( IQAEND+KD, IQEND );
                           zrot(NQ, Q( IQB, J-1 ), 1, Q( IQB, J ), 1, D( J ), DCONJG( WORK( J ) ) );
                        } // 50
                     } else {

                        DO 60 J = J1, J2, KD1;
                           zrot(N, Q( 1, J-1 ), 1, Q( 1, J ), 1, D( J ), DCONJG( WORK( J ) ) );
                        } // 60
                     }

                  }

                  if ( J2+KDN > N ) {

                     // adjust J2 to keep within the bounds of the matrix

                     NR = NR - 1;
                     J2 = J2 - KDN - 1;
                  }

                  DO 70 J = J1, J2, KD1;

                     // create nonzero element a(j-1,j+kd) outside the band
                     // and store it in WORK

                     WORK( J+KD ) = WORK( J )*AB( 1, J+KD );
                     AB( 1, J+KD ) = D( J )*AB( 1, J+KD );
                  } // 70
               } // 80
            } // 90
         }

         if ( KD > 0 ) {

            // make off-diagonal elements real and copy them to E

            for (I = 1; I <= N - 1; I++) { // 100
               T = AB( KD, I+1 );
               ABST = ABS( T );
               AB( KD, I+1 ) = ABST;
               E( I ) = ABST;
               if ( ABST != ZERO ) {
                  T = T / ABST;
               } else {
                  T = CONE;
               }
               if (I < N-1) AB( KD, I+2 ) = AB( KD, I+2 )*T;
               if ( WANTQ ) {
                  zscal(N, DCONJG( T ), Q( 1, I+1 ), 1 );
               }
            } // 100
         } else {

            // set E to zero if original matrix was diagonal

            for (I = 1; I <= N - 1; I++) { // 110
               E( I ) = ZERO;
            } // 110
         }

         // copy diagonal elements to D

         for (I = 1; I <= N; I++) { // 120
            D( I ) = DBLE( AB( KD1, I ) );
         } // 120

      } else {

         if ( KD > 1 ) {

            // Reduce to complex Hermitian tridiagonal form, working with
            // the lower triangle

            NR = 0;
            J1 = KDN + 2;
            J2 = 1;

            AB( 1, 1 ) = DBLE( AB( 1, 1 ) );
            for (I = 1; I <= N - 2; I++) { // 210

               // Reduce i-th column of matrix to tridiagonal form

               DO 200 K = KDN + 1, 2, -1;
                  J1 = J1 + KDN;
                  J2 = J2 + KDN;

                  if ( NR > 0 ) {

                     // generate plane rotations to annihilate nonzero
                     // elements which have been created outside the band

                     zlargv(NR, AB( KD1, J1-KD1 ), INCA, WORK( J1 ), KD1, D( J1 ), KD1 );

                     // apply plane rotations from one side


                     // Dependent on the the number of diagonals either
                     // ZLARTV or ZROT is used

                     if ( NR > 2*KD-1 ) {
                        for (L = 1; L <= KD - 1; L++) { // 130
                           zlartv(NR, AB( KD1-L, J1-KD1+L ), INCA, AB( KD1-L+1, J1-KD1+L ), INCA, D( J1 ), WORK( J1 ), KD1 );
                        } // 130
                     } else {
                        JEND = J1 + KD1*( NR-1 );
                        DO 140 JINC = J1, JEND, KD1;
                           zrot(KDM1, AB( KD, JINC-KD ), INCX, AB( KD1, JINC-KD ), INCX, D( JINC ), WORK( JINC ) );
                        } // 140
                     }

                  }

                  if ( K > 2 ) {
                     if ( K <= N-I+1 ) {

                        // generate plane rotation to annihilate a(i+k-1,i)
                        // within the band

                        zlartg(AB( K-1, I ), AB( K, I ), D( I+K-1 ), WORK( I+K-1 ), TEMP );
                        AB( K-1, I ) = TEMP;

                        // apply rotation from the left

                        zrot(K-3, AB( K-2, I+1 ), LDAB-1, AB( K-1, I+1 ), LDAB-1, D( I+K-1 ), WORK( I+K-1 ) );
                     }
                     NR = NR + 1;
                     J1 = J1 - KDN - 1;
                  }

                  // apply plane rotations from both sides to diagonal
                  // blocks

                  if (NR > 0) CALL ZLAR2V( NR, AB( 1, J1-1 ), AB( 1, J1 ), AB( 2, J1-1 ), INCA, D( J1 ), WORK( J1 ), KD1 );

                  // apply plane rotations from the right


                     // Dependent on the the number of diagonals either
                     // ZLARTV or ZROT is used

                  if ( NR > 0 ) {
                     zlacgv(NR, WORK( J1 ), KD1 );
                     if ( NR > 2*KD-1 ) {
                        for (L = 1; L <= KD - 1; L++) { // 150
                           if ( J2+L > N ) {
                              NRT = NR - 1;
                           } else {
                              NRT = NR;
                           }
                           if (NRT > 0) CALL ZLARTV( NRT, AB( L+2, J1-1 ), INCA, AB( L+1, J1 ), INCA, D( J1 ), WORK( J1 ), KD1 );
                        } // 150
                     } else {
                        J1END = J1 + KD1*( NR-2 );
                        if ( J1END >= J1 ) {
                           DO 160 J1INC = J1, J1END, KD1;
                              zrot(KDM1, AB( 3, J1INC-1 ), 1, AB( 2, J1INC ), 1, D( J1INC ), WORK( J1INC ) );
                           } // 160
                        }
                        LEND = min( KDM1, N-J2 );
                        LAST = J1END + KD1;
                        if (LEND > 0) CALL ZROT( LEND, AB( 3, LAST-1 ), 1, AB( 2, LAST ), 1, D( LAST ), WORK( LAST ) );
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
                        DO 170 J = J1, J2, KD1;
                           IBL = I - I2 / KDM1;
                           I2 = I2 + 1;
                           IQB = max( 1, J-IBL );
                           NQ = 1 + IQAEND - IQB;
                           IQAEND = min( IQAEND+KD, IQEND );
                           zrot(NQ, Q( IQB, J-1 ), 1, Q( IQB, J ), 1, D( J ), WORK( J ) );
                        } // 170
                     } else {

                        DO 180 J = J1, J2, KD1;
                           zrot(N, Q( 1, J-1 ), 1, Q( 1, J ), 1, D( J ), WORK( J ) );
                        } // 180
                     }
                  }

                  if ( J2+KDN > N ) {

                     // adjust J2 to keep within the bounds of the matrix

                     NR = NR - 1;
                     J2 = J2 - KDN - 1;
                  }

                  DO 190 J = J1, J2, KD1;

                     // create nonzero element a(j+kd,j-1) outside the
                     // band and store it in WORK

                     WORK( J+KD ) = WORK( J )*AB( KD1, J );
                     AB( KD1, J ) = D( J )*AB( KD1, J );
                  } // 190
               } // 200
            } // 210
         }

         if ( KD > 0 ) {

            // make off-diagonal elements real and copy them to E

            for (I = 1; I <= N - 1; I++) { // 220
               T = AB( 2, I );
               ABST = ABS( T );
               AB( 2, I ) = ABST;
               E( I ) = ABST;
               if ( ABST != ZERO ) {
                  T = T / ABST;
               } else {
                  T = CONE;
               }
               if (I < N-1) AB( 2, I+1 ) = AB( 2, I+1 )*T;
               if ( WANTQ ) {
                  zscal(N, T, Q( 1, I+1 ), 1 );
               }
            } // 220
         } else {

            // set E to zero if original matrix was diagonal

            for (I = 1; I <= N - 1; I++) { // 230
               E( I ) = ZERO;
            } // 230
         }

         // copy diagonal elements to D

         for (I = 1; I <= N; I++) { // 240
            D( I ) = DBLE( AB( 1, I ) );
         } // 240
      }

      return;

      // End of ZHBTRD

      }
