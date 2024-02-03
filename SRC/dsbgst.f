      void dsbgst(VECT, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, X, LDX, WORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO, VECT;
      int                INFO, KA, KB, LDAB, LDBB, LDX, N;
      // ..
      // .. Array Arguments ..
      double             AB( LDAB, * ), BB( LDBB, * ), WORK( * ), X( LDX, * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      bool               UPDATE, UPPER, WANTX;
      int                I, I0, I1, I2, INCA, J, J1, J1T, J2, J2T, K, KA1, KB1, KBT, L, M, NR, NRT, NX;
      double             BII, RA, RA1, T;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGER, DLAR2V, DLARGV, DLARTG, DLARTV, DLASET, DROT, DSCAL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..

      // Test the input parameters

      WANTX = LSAME( VECT, 'V' );
      UPPER = LSAME( UPLO, 'U' );
      KA1 = KA + 1;
      KB1 = KB + 1;
      INFO = 0;
      if ( !WANTX && !LSAME( VECT, 'N' ) ) {
         INFO = -1;
      } else if ( !UPPER && !LSAME( UPLO, 'L' ) ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -3;
      } else if ( KA < 0 ) {
         INFO = -4;
      } else if ( KB < 0 || KB > KA ) {
         INFO = -5;
      } else if ( LDAB < KA+1 ) {
         INFO = -7;
      } else if ( LDBB < KB+1 ) {
         INFO = -9;
      } else if ( LDX < 1 || WANTX && LDX < max( 1, N ) ) {
         INFO = -11;
      }
      if ( INFO != 0 ) {
         xerbla('DSBGST', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      INCA = LDAB*KA1;

      // Initialize X to the unit matrix, if needed

      if (WANTX) dlaset( 'Full', N, N, ZERO, ONE, X, LDX );

      // Set M to the splitting point m. It must be the same value as is
      // used in DPBSTF. The chosen value allows the arrays WORK and RWORK
      // to be of dimension (N).

      M = ( N+KB ) / 2;

      // The routine works in two phases, corresponding to the two halves
      // of the split Cholesky factorization of B as S**T*S where

      // S = ( U    )
          // ( M  L )

      // with U upper triangular of order m, and L lower triangular of
      // order n-m. S has the same bandwidth as B.

      // S is treated as a product of elementary matrices:

      // S = S(m)*S(m-1)*...*S(2)*S(1)*S(m+1)*S(m+2)*...*S(n-1)*S(n)

      // where S(i) is determined by the i-th row of S.

      // In phase 1, the index i takes the values n, n-1, ... , m+1;
      // in phase 2, it takes the values 1, 2, ... , m.

      // For each value of i, the current matrix A is updated by forming
      // inv(S(i))**T*A*inv(S(i)). This creates a triangular bulge outside
      // the band of A. The bulge is then pushed down toward the bottom of
      // A in phase 1, and up toward the top of A in phase 2, by applying
      // plane rotations.

      // There are kb*(kb+1)/2 elements in the bulge, but at most 2*kb-1
      // of them are linearly independent, so annihilating a bulge requires
      // only 2*kb-1 plane rotations. The rotations are divided into a 1st
      // set of kb-1 rotations, and a 2nd set of kb rotations.

      // Wherever possible, rotations are generated and applied in vector
      // operations of length NR between the indices J1 and J2 (sometimes
      // replaced by modified values NRT, J1T or J2T).

      // The cosines and sines of the rotations are stored in the array
      // WORK. The cosines of the 1st set of rotations are stored in
      // elements n+2:n+m-kb-1 and the sines of the 1st set in elements
      // 2:m-kb-1; the cosines of the 2nd set are stored in elements
      // n+m-kb+1:2*n and the sines of the second set in elements m-kb+1:n.

      // The bulges are not formed explicitly; nonzero elements outside the
      // band are created only when they are required for generating new
      // rotations; they are stored in the array WORK, in positions where
      // they are later overwritten by the sines of the rotations which
      // annihilate them.

      // **************************** Phase 1 *****************************

      // The logical structure of this phase is:

      // UPDATE = true;
      // DO I = N, M + 1, -1
         // use S(i) to update A and create a new bulge
         // apply rotations to push all bulges KA positions downward
      // END DO
      // UPDATE = false;
      // DO I = M + KA + 1, N - 1
         // apply rotations to push all bulges KA positions downward
      // END DO

      // To avoid duplicating code, the two loops are merged.

      UPDATE = true;
      I = N + 1;
      } // 10
      if ( UPDATE ) {
         I = I - 1;
         KBT = min( KB, I-1 );
         I0 = I - 1;
         I1 = min( N, I+KA );
         I2 = I - KBT + KA1;
         if ( I < M+1 ) {
            UPDATE = false;
            I = I + 1;
            I0 = M;
            if (KA == 0) GO TO 480;
            GO TO 10;
         }
      } else {
         I = I + KA;
         if (I > N-1) GO TO 480;
      }

      if ( UPPER ) {

         // Transform A, working with the upper triangle

         if ( UPDATE ) {

            // Form  inv(S(i))**T * A * inv(S(i))

            BII = BB( KB1, I );
            for (J = I; J <= I1; J++) { // 20
               AB( I-J+KA1, J ) = AB( I-J+KA1, J ) / BII;
            } // 20
            DO 30 J = max( 1, I-KA ), I;
               AB( J-I+KA1, I ) = AB( J-I+KA1, I ) / BII;
            } // 30
            for (K = I - KBT; K <= I - 1; K++) { // 60
               for (J = I - KBT; J <= K; J++) { // 40
                  AB( J-K+KA1, K ) = AB( J-K+KA1, K ) - BB( J-I+KB1, I )*AB( K-I+KA1, I ) - BB( K-I+KB1, I )*AB( J-I+KA1, I ) + AB( KA1, I )*BB( J-I+KB1, I )* BB( K-I+KB1, I );
               } // 40
               DO 50 J = max( 1, I-KA ), I - KBT - 1;
                  AB( J-K+KA1, K ) = AB( J-K+KA1, K ) - BB( K-I+KB1, I )*AB( J-I+KA1, I );
               } // 50
            } // 60
            for (J = I; J <= I1; J++) { // 80
               DO 70 K = max( J-KA, I-KBT ), I - 1;
                  AB( K-J+KA1, J ) = AB( K-J+KA1, J ) - BB( K-I+KB1, I )*AB( I-J+KA1, J );
               } // 70
            } // 80

            if ( WANTX ) {

               // post-multiply X by inv(S(i))

               dscal(N-M, ONE / BII, X( M+1, I ), 1 );
               if (KBT > 0) dger( N-M, KBT, -ONE, X( M+1, I ), 1, BB( KB1-KBT, I ), 1, X( M+1, I-KBT ), LDX );
            }

            // store a(i,i1) in RA1 for use in next loop over K

            RA1 = AB( I-I1+KA1, I1 );
         }

         // Generate and apply vectors of rotations to chase all the
         // existing bulges KA positions down toward the bottom of the
         // band

         for (K = 1; K <= KB - 1; K++) { // 130
            if ( UPDATE ) {

               // Determine the rotations which would annihilate the bulge
               // which has in theory just been created

               if ( I-K+KA < N && I-K > 1 ) {

                  // generate rotation to annihilate a(i,i-k+ka+1)

                  dlartg(AB( K+1, I-K+KA ), RA1, WORK( N+I-K+KA-M ), WORK( I-K+KA-M ), RA );

                  // create nonzero element a(i-k,i-k+ka+1) outside the
                  // band and store it in WORK(i-k)

                  T = -BB( KB1-K, I )*RA1;
                  WORK( I-K ) = WORK( N+I-K+KA-M )*T - WORK( I-K+KA-M )*AB( 1, I-K+KA )                   AB( 1, I-K+KA ) = WORK( I-K+KA-M )*T + WORK( N+I-K+KA-M )*AB( 1, I-K+KA );
                  RA1 = RA;
               }
            }
            J2 = I - K - 1 + max( 1, K-I0+2 )*KA1;
            NR = ( N-J2+KA ) / KA1;
            J1 = J2 + ( NR-1 )*KA1;
            if ( UPDATE ) {
               J2T = max( J2, I+2*KA-K+1 );
            } else {
               J2T = J2;
            }
            NRT = ( N-J2T+KA ) / KA1;
            DO 90 J = J2T, J1, KA1;

               // create nonzero element a(j-ka,j+1) outside the band
               // and store it in WORK(j-m)

               WORK( J-M ) = WORK( J-M )*AB( 1, J+1 );
               AB( 1, J+1 ) = WORK( N+J-M )*AB( 1, J+1 );
            } // 90

            // generate rotations in 1st set to annihilate elements which
            // have been created outside the band

            if (NRT > 0) dlargv( NRT, AB( 1, J2T ), INCA, WORK( J2T-M ), KA1, WORK( N+J2T-M ), KA1 );
            if ( NR > 0 ) {

               // apply rotations in 1st set from the right

               for (L = 1; L <= KA - 1; L++) { // 100
                  dlartv(NR, AB( KA1-L, J2 ), INCA, AB( KA-L, J2+1 ), INCA, WORK( N+J2-M ), WORK( J2-M ), KA1 );
               } // 100

               // apply rotations in 1st set from both sides to diagonal
               // blocks

               dlar2v(NR, AB( KA1, J2 ), AB( KA1, J2+1 ), AB( KA, J2+1 ), INCA, WORK( N+J2-M ), WORK( J2-M ), KA1 );

            }

            // start applying rotations in 1st set from the left

            DO 110 L = KA - 1, KB - K + 1, -1;
               NRT = ( N-J2+L ) / KA1;
               if (NRT > 0) dlartv( NRT, AB( L, J2+KA1-L ), INCA, AB( L+1, J2+KA1-L ), INCA, WORK( N+J2-M ), WORK( J2-M ), KA1 );
            } // 110

            if ( WANTX ) {

               // post-multiply X by product of rotations in 1st set

               DO 120 J = J2, J1, KA1;
                  drot(N-M, X( M+1, J ), 1, X( M+1, J+1 ), 1, WORK( N+J-M ), WORK( J-M ) );
               } // 120
            }
         } // 130

         if ( UPDATE ) {
            if ( I2 <= N && KBT > 0 ) {

               // create nonzero element a(i-kbt,i-kbt+ka+1) outside the
               // band and store it in WORK(i-kbt)

               WORK( I-KBT ) = -BB( KB1-KBT, I )*RA1;
            }
         }

         DO 170 K = KB, 1, -1;
            if ( UPDATE ) {
               J2 = I - K - 1 + max( 2, K-I0+1 )*KA1;
            } else {
               J2 = I - K - 1 + max( 1, K-I0+1 )*KA1;
            }

            // finish applying rotations in 2nd set from the left

            DO 140 L = KB - K, 1, -1;
               NRT = ( N-J2+KA+L ) / KA1;
               if (NRT > 0) dlartv( NRT, AB( L, J2-L+1 ), INCA, AB( L+1, J2-L+1 ), INCA, WORK( N+J2-KA ), WORK( J2-KA ), KA1 );
            } // 140
            NR = ( N-J2+KA ) / KA1;
            J1 = J2 + ( NR-1 )*KA1;
            DO 150 J = J1, J2, -KA1;
               WORK( J ) = WORK( J-KA );
               WORK( N+J ) = WORK( N+J-KA );
            } // 150
            DO 160 J = J2, J1, KA1;

               // create nonzero element a(j-ka,j+1) outside the band
               // and store it in WORK(j)

               WORK( J ) = WORK( J )*AB( 1, J+1 );
               AB( 1, J+1 ) = WORK( N+J )*AB( 1, J+1 );
            } // 160
            if ( UPDATE ) {
               if (I-K < N-KA && K <= KBT) WORK( I-K+KA ) = WORK( I-K );
            }
         } // 170

         DO 210 K = KB, 1, -1;
            J2 = I - K - 1 + max( 1, K-I0+1 )*KA1;
            NR = ( N-J2+KA ) / KA1;
            J1 = J2 + ( NR-1 )*KA1;
            if ( NR > 0 ) {

               // generate rotations in 2nd set to annihilate elements
               // which have been created outside the band

               dlargv(NR, AB( 1, J2 ), INCA, WORK( J2 ), KA1, WORK( N+J2 ), KA1 );

               // apply rotations in 2nd set from the right

               for (L = 1; L <= KA - 1; L++) { // 180
                  dlartv(NR, AB( KA1-L, J2 ), INCA, AB( KA-L, J2+1 ), INCA, WORK( N+J2 ), WORK( J2 ), KA1 );
               } // 180

               // apply rotations in 2nd set from both sides to diagonal
               // blocks

               dlar2v(NR, AB( KA1, J2 ), AB( KA1, J2+1 ), AB( KA, J2+1 ), INCA, WORK( N+J2 ), WORK( J2 ), KA1 );

            }

            // start applying rotations in 2nd set from the left

            DO 190 L = KA - 1, KB - K + 1, -1;
               NRT = ( N-J2+L ) / KA1;
               if (NRT > 0) dlartv( NRT, AB( L, J2+KA1-L ), INCA, AB( L+1, J2+KA1-L ), INCA, WORK( N+J2 ), WORK( J2 ), KA1 );
            } // 190

            if ( WANTX ) {

               // post-multiply X by product of rotations in 2nd set

               DO 200 J = J2, J1, KA1;
                  drot(N-M, X( M+1, J ), 1, X( M+1, J+1 ), 1, WORK( N+J ), WORK( J ) );
               } // 200
            }
         } // 210

         for (K = 1; K <= KB - 1; K++) { // 230
            J2 = I - K - 1 + max( 1, K-I0+2 )*KA1;

            // finish applying rotations in 1st set from the left

            DO 220 L = KB - K, 1, -1;
               NRT = ( N-J2+L ) / KA1;
               if (NRT > 0) dlartv( NRT, AB( L, J2+KA1-L ), INCA, AB( L+1, J2+KA1-L ), INCA, WORK( N+J2-M ), WORK( J2-M ), KA1 );
            } // 220
         } // 230

         if ( KB > 1 ) {
            DO 240 J = N - 1, I - KB + 2*KA + 1, -1;
               WORK( N+J-M ) = WORK( N+J-KA-M );
               WORK( J-M ) = WORK( J-KA-M );
            } // 240
         }

      } else {

         // Transform A, working with the lower triangle

         if ( UPDATE ) {

            // Form  inv(S(i))**T * A * inv(S(i))

            BII = BB( 1, I );
            for (J = I; J <= I1; J++) { // 250
               AB( J-I+1, I ) = AB( J-I+1, I ) / BII;
            } // 250
            DO 260 J = max( 1, I-KA ), I;
               AB( I-J+1, J ) = AB( I-J+1, J ) / BII;
            } // 260
            for (K = I - KBT; K <= I - 1; K++) { // 290
               for (J = I - KBT; J <= K; J++) { // 270
                  AB( K-J+1, J ) = AB( K-J+1, J ) - BB( I-J+1, J )*AB( I-K+1, K ) - BB( I-K+1, K )*AB( I-J+1, J ) + AB( 1, I )*BB( I-J+1, J )* BB( I-K+1, K );
               } // 270
               DO 280 J = max( 1, I-KA ), I - KBT - 1;
                  AB( K-J+1, J ) = AB( K-J+1, J ) - BB( I-K+1, K )*AB( I-J+1, J );
               } // 280
            } // 290
            for (J = I; J <= I1; J++) { // 310
               DO 300 K = max( J-KA, I-KBT ), I - 1;
                  AB( J-K+1, K ) = AB( J-K+1, K ) - BB( I-K+1, K )*AB( J-I+1, I );
               } // 300
            } // 310

            if ( WANTX ) {

               // post-multiply X by inv(S(i))

               dscal(N-M, ONE / BII, X( M+1, I ), 1 );
               if (KBT > 0) dger( N-M, KBT, -ONE, X( M+1, I ), 1, BB( KBT+1, I-KBT ), LDBB-1, X( M+1, I-KBT ), LDX );
            }

            // store a(i1,i) in RA1 for use in next loop over K

            RA1 = AB( I1-I+1, I );
         }

         // Generate and apply vectors of rotations to chase all the
         // existing bulges KA positions down toward the bottom of the
         // band

         for (K = 1; K <= KB - 1; K++) { // 360
            if ( UPDATE ) {

               // Determine the rotations which would annihilate the bulge
               // which has in theory just been created

               if ( I-K+KA < N && I-K > 1 ) {

                  // generate rotation to annihilate a(i-k+ka+1,i)

                  dlartg(AB( KA1-K, I ), RA1, WORK( N+I-K+KA-M ), WORK( I-K+KA-M ), RA );

                  // create nonzero element a(i-k+ka+1,i-k) outside the
                  // band and store it in WORK(i-k)

                  T = -BB( K+1, I-K )*RA1;
                  WORK( I-K ) = WORK( N+I-K+KA-M )*T - WORK( I-K+KA-M )*AB( KA1, I-K )                   AB( KA1, I-K ) = WORK( I-K+KA-M )*T + WORK( N+I-K+KA-M )*AB( KA1, I-K );
                  RA1 = RA;
               }
            }
            J2 = I - K - 1 + max( 1, K-I0+2 )*KA1;
            NR = ( N-J2+KA ) / KA1;
            J1 = J2 + ( NR-1 )*KA1;
            if ( UPDATE ) {
               J2T = max( J2, I+2*KA-K+1 );
            } else {
               J2T = J2;
            }
            NRT = ( N-J2T+KA ) / KA1;
            DO 320 J = J2T, J1, KA1;

               // create nonzero element a(j+1,j-ka) outside the band
               // and store it in WORK(j-m)

               WORK( J-M ) = WORK( J-M )*AB( KA1, J-KA+1 );
               AB( KA1, J-KA+1 ) = WORK( N+J-M )*AB( KA1, J-KA+1 );
            } // 320

            // generate rotations in 1st set to annihilate elements which
            // have been created outside the band

            if (NRT > 0) dlargv( NRT, AB( KA1, J2T-KA ), INCA, WORK( J2T-M ), KA1, WORK( N+J2T-M ), KA1 );
            if ( NR > 0 ) {

               // apply rotations in 1st set from the left

               for (L = 1; L <= KA - 1; L++) { // 330
                  dlartv(NR, AB( L+1, J2-L ), INCA, AB( L+2, J2-L ), INCA, WORK( N+J2-M ), WORK( J2-M ), KA1 );
               } // 330

               // apply rotations in 1st set from both sides to diagonal
               // blocks

               dlar2v(NR, AB( 1, J2 ), AB( 1, J2+1 ), AB( 2, J2 ), INCA, WORK( N+J2-M ), WORK( J2-M ), KA1 );

            }

            // start applying rotations in 1st set from the right

            DO 340 L = KA - 1, KB - K + 1, -1;
               NRT = ( N-J2+L ) / KA1;
               if (NRT > 0) dlartv( NRT, AB( KA1-L+1, J2 ), INCA, AB( KA1-L, J2+1 ), INCA, WORK( N+J2-M ), WORK( J2-M ), KA1 );
            } // 340

            if ( WANTX ) {

               // post-multiply X by product of rotations in 1st set

               DO 350 J = J2, J1, KA1;
                  drot(N-M, X( M+1, J ), 1, X( M+1, J+1 ), 1, WORK( N+J-M ), WORK( J-M ) );
               } // 350
            }
         } // 360

         if ( UPDATE ) {
            if ( I2 <= N && KBT > 0 ) {

               // create nonzero element a(i-kbt+ka+1,i-kbt) outside the
               // band and store it in WORK(i-kbt)

               WORK( I-KBT ) = -BB( KBT+1, I-KBT )*RA1;
            }
         }

         DO 400 K = KB, 1, -1;
            if ( UPDATE ) {
               J2 = I - K - 1 + max( 2, K-I0+1 )*KA1;
            } else {
               J2 = I - K - 1 + max( 1, K-I0+1 )*KA1;
            }

            // finish applying rotations in 2nd set from the right

            DO 370 L = KB - K, 1, -1;
               NRT = ( N-J2+KA+L ) / KA1;
               if (NRT > 0) dlartv( NRT, AB( KA1-L+1, J2-KA ), INCA, AB( KA1-L, J2-KA+1 ), INCA, WORK( N+J2-KA ), WORK( J2-KA ), KA1 );
            } // 370
            NR = ( N-J2+KA ) / KA1;
            J1 = J2 + ( NR-1 )*KA1;
            DO 380 J = J1, J2, -KA1;
               WORK( J ) = WORK( J-KA );
               WORK( N+J ) = WORK( N+J-KA );
            } // 380
            DO 390 J = J2, J1, KA1;

               // create nonzero element a(j+1,j-ka) outside the band
               // and store it in WORK(j)

               WORK( J ) = WORK( J )*AB( KA1, J-KA+1 );
               AB( KA1, J-KA+1 ) = WORK( N+J )*AB( KA1, J-KA+1 );
            } // 390
            if ( UPDATE ) {
               if (I-K < N-KA && K <= KBT) WORK( I-K+KA ) = WORK( I-K );
            }
         } // 400

         DO 440 K = KB, 1, -1;
            J2 = I - K - 1 + max( 1, K-I0+1 )*KA1;
            NR = ( N-J2+KA ) / KA1;
            J1 = J2 + ( NR-1 )*KA1;
            if ( NR > 0 ) {

               // generate rotations in 2nd set to annihilate elements
               // which have been created outside the band

               dlargv(NR, AB( KA1, J2-KA ), INCA, WORK( J2 ), KA1, WORK( N+J2 ), KA1 );

               // apply rotations in 2nd set from the left

               for (L = 1; L <= KA - 1; L++) { // 410
                  dlartv(NR, AB( L+1, J2-L ), INCA, AB( L+2, J2-L ), INCA, WORK( N+J2 ), WORK( J2 ), KA1 );
               } // 410

               // apply rotations in 2nd set from both sides to diagonal
               // blocks

               dlar2v(NR, AB( 1, J2 ), AB( 1, J2+1 ), AB( 2, J2 ), INCA, WORK( N+J2 ), WORK( J2 ), KA1 );

            }

            // start applying rotations in 2nd set from the right

            DO 420 L = KA - 1, KB - K + 1, -1;
               NRT = ( N-J2+L ) / KA1;
               if (NRT > 0) dlartv( NRT, AB( KA1-L+1, J2 ), INCA, AB( KA1-L, J2+1 ), INCA, WORK( N+J2 ), WORK( J2 ), KA1 );
            } // 420

            if ( WANTX ) {

               // post-multiply X by product of rotations in 2nd set

               DO 430 J = J2, J1, KA1;
                  drot(N-M, X( M+1, J ), 1, X( M+1, J+1 ), 1, WORK( N+J ), WORK( J ) );
               } // 430
            }
         } // 440

         for (K = 1; K <= KB - 1; K++) { // 460
            J2 = I - K - 1 + max( 1, K-I0+2 )*KA1;

            // finish applying rotations in 1st set from the right

            DO 450 L = KB - K, 1, -1;
               NRT = ( N-J2+L ) / KA1;
               if (NRT > 0) dlartv( NRT, AB( KA1-L+1, J2 ), INCA, AB( KA1-L, J2+1 ), INCA, WORK( N+J2-M ), WORK( J2-M ), KA1 );
            } // 450
         } // 460

         if ( KB > 1 ) {
            DO 470 J = N - 1, I - KB + 2*KA + 1, -1;
               WORK( N+J-M ) = WORK( N+J-KA-M );
               WORK( J-M ) = WORK( J-KA-M );
            } // 470
         }

      }

      GO TO 10;

      } // 480

      // **************************** Phase 2 *****************************

      // The logical structure of this phase is:

      // UPDATE = true;
      // DO I = 1, M
         // use S(i) to update A and create a new bulge
         // apply rotations to push all bulges KA positions upward
      // END DO
      // UPDATE = false;
      // DO I = M - KA - 1, 2, -1
         // apply rotations to push all bulges KA positions upward
      // END DO

      // To avoid duplicating code, the two loops are merged.

      UPDATE = true;
      I = 0;
      } // 490
      if ( UPDATE ) {
         I = I + 1;
         KBT = min( KB, M-I );
         I0 = I + 1;
         I1 = max( 1, I-KA );
         I2 = I + KBT - KA1;
         if ( I > M ) {
            UPDATE = false;
            I = I - 1;
            I0 = M + 1;
            if (KA == 0) return;
            GO TO 490;
         }
      } else {
         I = I - KA;
         if (I < 2) return;
      }

      if ( I < M-KBT ) {
         NX = M;
      } else {
         NX = N;
      }

      if ( UPPER ) {

         // Transform A, working with the upper triangle

         if ( UPDATE ) {

            // Form  inv(S(i))**T * A * inv(S(i))

            BII = BB( KB1, I );
            for (J = I1; J <= I; J++) { // 500
               AB( J-I+KA1, I ) = AB( J-I+KA1, I ) / BII;
            } // 500
            DO 510 J = I, min( N, I+KA );
               AB( I-J+KA1, J ) = AB( I-J+KA1, J ) / BII;
            } // 510
            for (K = I + 1; K <= I + KBT; K++) { // 540
               for (J = K; J <= I + KBT; J++) { // 520
                  AB( K-J+KA1, J ) = AB( K-J+KA1, J ) - BB( I-J+KB1, J )*AB( I-K+KA1, K ) - BB( I-K+KB1, K )*AB( I-J+KA1, J ) + AB( KA1, I )*BB( I-J+KB1, J )* BB( I-K+KB1, K );
               } // 520
               DO 530 J = I + KBT + 1, min( N, I+KA );
                  AB( K-J+KA1, J ) = AB( K-J+KA1, J ) - BB( I-K+KB1, K )*AB( I-J+KA1, J );
               } // 530
            } // 540
            for (J = I1; J <= I; J++) { // 560
               DO 550 K = I + 1, min( J+KA, I+KBT );
                  AB( J-K+KA1, K ) = AB( J-K+KA1, K ) - BB( I-K+KB1, K )*AB( J-I+KA1, I );
               } // 550
            } // 560

            if ( WANTX ) {

               // post-multiply X by inv(S(i))

               dscal(NX, ONE / BII, X( 1, I ), 1 );
               if (KBT > 0) dger( NX, KBT, -ONE, X( 1, I ), 1, BB( KB, I+1 ), LDBB-1, X( 1, I+1 ), LDX );
            }

            // store a(i1,i) in RA1 for use in next loop over K

            RA1 = AB( I1-I+KA1, I );
         }

         // Generate and apply vectors of rotations to chase all the
         // existing bulges KA positions up toward the top of the band

         for (K = 1; K <= KB - 1; K++) { // 610
            if ( UPDATE ) {

               // Determine the rotations which would annihilate the bulge
               // which has in theory just been created

               if ( I+K-KA1 > 0 && I+K < M ) {

                  // generate rotation to annihilate a(i+k-ka-1,i)

                  dlartg(AB( K+1, I ), RA1, WORK( N+I+K-KA ), WORK( I+K-KA ), RA );

                  // create nonzero element a(i+k-ka-1,i+k) outside the
                  // band and store it in WORK(m-kb+i+k)

                  T = -BB( KB1-K, I+K )*RA1;
                  WORK( M-KB+I+K ) = WORK( N+I+K-KA )*T - WORK( I+K-KA )*AB( 1, I+K )                   AB( 1, I+K ) = WORK( I+K-KA )*T + WORK( N+I+K-KA )*AB( 1, I+K );
                  RA1 = RA;
               }
            }
            J2 = I + K + 1 - max( 1, K+I0-M+1 )*KA1;
            NR = ( J2+KA-1 ) / KA1;
            J1 = J2 - ( NR-1 )*KA1;
            if ( UPDATE ) {
               J2T = min( J2, I-2*KA+K-1 );
            } else {
               J2T = J2;
            }
            NRT = ( J2T+KA-1 ) / KA1;
            DO 570 J = J1, J2T, KA1;

               // create nonzero element a(j-1,j+ka) outside the band
               // and store it in WORK(j)

               WORK( J ) = WORK( J )*AB( 1, J+KA-1 );
               AB( 1, J+KA-1 ) = WORK( N+J )*AB( 1, J+KA-1 );
            } // 570

            // generate rotations in 1st set to annihilate elements which
            // have been created outside the band

            if (NRT > 0) dlargv( NRT, AB( 1, J1+KA ), INCA, WORK( J1 ), KA1, WORK( N+J1 ), KA1 );
            if ( NR > 0 ) {

               // apply rotations in 1st set from the left

               for (L = 1; L <= KA - 1; L++) { // 580
                  dlartv(NR, AB( KA1-L, J1+L ), INCA, AB( KA-L, J1+L ), INCA, WORK( N+J1 ), WORK( J1 ), KA1 );
               } // 580

               // apply rotations in 1st set from both sides to diagonal
               // blocks

               dlar2v(NR, AB( KA1, J1 ), AB( KA1, J1-1 ), AB( KA, J1 ), INCA, WORK( N+J1 ), WORK( J1 ), KA1 );

            }

            // start applying rotations in 1st set from the right

            DO 590 L = KA - 1, KB - K + 1, -1;
               NRT = ( J2+L-1 ) / KA1;
               J1T = J2 - ( NRT-1 )*KA1;
               if (NRT > 0) dlartv( NRT, AB( L, J1T ), INCA, AB( L+1, J1T-1 ), INCA, WORK( N+J1T ), WORK( J1T ), KA1 );
            } // 590

            if ( WANTX ) {

               // post-multiply X by product of rotations in 1st set

               DO 600 J = J1, J2, KA1;
                  drot(NX, X( 1, J ), 1, X( 1, J-1 ), 1, WORK( N+J ), WORK( J ) );
               } // 600
            }
         } // 610

         if ( UPDATE ) {
            if ( I2 > 0 && KBT > 0 ) {

               // create nonzero element a(i+kbt-ka-1,i+kbt) outside the
               // band and store it in WORK(m-kb+i+kbt)

               WORK( M-KB+I+KBT ) = -BB( KB1-KBT, I+KBT )*RA1;
            }
         }

         DO 650 K = KB, 1, -1;
            if ( UPDATE ) {
               J2 = I + K + 1 - max( 2, K+I0-M )*KA1;
            } else {
               J2 = I + K + 1 - max( 1, K+I0-M )*KA1;
            }

            // finish applying rotations in 2nd set from the right

            DO 620 L = KB - K, 1, -1;
               NRT = ( J2+KA+L-1 ) / KA1;
               J1T = J2 - ( NRT-1 )*KA1;
               if (NRT > 0) dlartv( NRT, AB( L, J1T+KA ), INCA, AB( L+1, J1T+KA-1 ), INCA, WORK( N+M-KB+J1T+KA ), WORK( M-KB+J1T+KA ), KA1 );
            } // 620
            NR = ( J2+KA-1 ) / KA1;
            J1 = J2 - ( NR-1 )*KA1;
            DO 630 J = J1, J2, KA1;
               WORK( M-KB+J ) = WORK( M-KB+J+KA );
               WORK( N+M-KB+J ) = WORK( N+M-KB+J+KA );
            } // 630
            DO 640 J = J1, J2, KA1;

               // create nonzero element a(j-1,j+ka) outside the band
               // and store it in WORK(m-kb+j)

               WORK( M-KB+J ) = WORK( M-KB+J )*AB( 1, J+KA-1 );
               AB( 1, J+KA-1 ) = WORK( N+M-KB+J )*AB( 1, J+KA-1 );
            } // 640
            if ( UPDATE ) {
               if (I+K > KA1 && K <= KBT) WORK( M-KB+I+K-KA ) = WORK( M-KB+I+K );
            }
         } // 650

         DO 690 K = KB, 1, -1;
            J2 = I + K + 1 - max( 1, K+I0-M )*KA1;
            NR = ( J2+KA-1 ) / KA1;
            J1 = J2 - ( NR-1 )*KA1;
            if ( NR > 0 ) {

               // generate rotations in 2nd set to annihilate elements
               // which have been created outside the band

               dlargv(NR, AB( 1, J1+KA ), INCA, WORK( M-KB+J1 ), KA1, WORK( N+M-KB+J1 ), KA1 );

               // apply rotations in 2nd set from the left

               for (L = 1; L <= KA - 1; L++) { // 660
                  dlartv(NR, AB( KA1-L, J1+L ), INCA, AB( KA-L, J1+L ), INCA, WORK( N+M-KB+J1 ), WORK( M-KB+J1 ), KA1 );
               } // 660

               // apply rotations in 2nd set from both sides to diagonal
               // blocks

               dlar2v(NR, AB( KA1, J1 ), AB( KA1, J1-1 ), AB( KA, J1 ), INCA, WORK( N+M-KB+J1 ), WORK( M-KB+J1 ), KA1 );

            }

            // start applying rotations in 2nd set from the right

            DO 670 L = KA - 1, KB - K + 1, -1;
               NRT = ( J2+L-1 ) / KA1;
               J1T = J2 - ( NRT-1 )*KA1;
               if (NRT > 0) dlartv( NRT, AB( L, J1T ), INCA, AB( L+1, J1T-1 ), INCA, WORK( N+M-KB+J1T ), WORK( M-KB+J1T ), KA1 );
            } // 670

            if ( WANTX ) {

               // post-multiply X by product of rotations in 2nd set

               DO 680 J = J1, J2, KA1;
                  drot(NX, X( 1, J ), 1, X( 1, J-1 ), 1, WORK( N+M-KB+J ), WORK( M-KB+J ) );
               } // 680
            }
         } // 690

         for (K = 1; K <= KB - 1; K++) { // 710
            J2 = I + K + 1 - max( 1, K+I0-M+1 )*KA1;

            // finish applying rotations in 1st set from the right

            DO 700 L = KB - K, 1, -1;
               NRT = ( J2+L-1 ) / KA1;
               J1T = J2 - ( NRT-1 )*KA1;
               if (NRT > 0) dlartv( NRT, AB( L, J1T ), INCA, AB( L+1, J1T-1 ), INCA, WORK( N+J1T ), WORK( J1T ), KA1 );
            } // 700
         } // 710

         if ( KB > 1 ) {
            DO 720 J = 2, min( I+KB, M ) - 2*KA - 1;
               WORK( N+J ) = WORK( N+J+KA );
               WORK( J ) = WORK( J+KA );
            } // 720
         }

      } else {

         // Transform A, working with the lower triangle

         if ( UPDATE ) {

            // Form  inv(S(i))**T * A * inv(S(i))

            BII = BB( 1, I );
            for (J = I1; J <= I; J++) { // 730
               AB( I-J+1, J ) = AB( I-J+1, J ) / BII;
            } // 730
            DO 740 J = I, min( N, I+KA );
               AB( J-I+1, I ) = AB( J-I+1, I ) / BII;
            } // 740
            for (K = I + 1; K <= I + KBT; K++) { // 770
               for (J = K; J <= I + KBT; J++) { // 750
                  AB( J-K+1, K ) = AB( J-K+1, K ) - BB( J-I+1, I )*AB( K-I+1, I ) - BB( K-I+1, I )*AB( J-I+1, I ) + AB( 1, I )*BB( J-I+1, I )* BB( K-I+1, I );
               } // 750
               DO 760 J = I + KBT + 1, min( N, I+KA );
                  AB( J-K+1, K ) = AB( J-K+1, K ) - BB( K-I+1, I )*AB( J-I+1, I );
               } // 760
            } // 770
            for (J = I1; J <= I; J++) { // 790
               DO 780 K = I + 1, min( J+KA, I+KBT );
                  AB( K-J+1, J ) = AB( K-J+1, J ) - BB( K-I+1, I )*AB( I-J+1, J );
               } // 780
            } // 790

            if ( WANTX ) {

               // post-multiply X by inv(S(i))

               dscal(NX, ONE / BII, X( 1, I ), 1 );
               if (KBT > 0) dger( NX, KBT, -ONE, X( 1, I ), 1, BB( 2, I ), 1, X( 1, I+1 ), LDX );
            }

            // store a(i,i1) in RA1 for use in next loop over K

            RA1 = AB( I-I1+1, I1 );
         }

         // Generate and apply vectors of rotations to chase all the
         // existing bulges KA positions up toward the top of the band

         for (K = 1; K <= KB - 1; K++) { // 840
            if ( UPDATE ) {

               // Determine the rotations which would annihilate the bulge
               // which has in theory just been created

               if ( I+K-KA1 > 0 && I+K < M ) {

                  // generate rotation to annihilate a(i,i+k-ka-1)

                  dlartg(AB( KA1-K, I+K-KA ), RA1, WORK( N+I+K-KA ), WORK( I+K-KA ), RA );

                  // create nonzero element a(i+k,i+k-ka-1) outside the
                  // band and store it in WORK(m-kb+i+k)

                  T = -BB( K+1, I )*RA1;
                  WORK( M-KB+I+K ) = WORK( N+I+K-KA )*T - WORK( I+K-KA )*AB( KA1, I+K-KA )                   AB( KA1, I+K-KA ) = WORK( I+K-KA )*T + WORK( N+I+K-KA )*AB( KA1, I+K-KA );
                  RA1 = RA;
               }
            }
            J2 = I + K + 1 - max( 1, K+I0-M+1 )*KA1;
            NR = ( J2+KA-1 ) / KA1;
            J1 = J2 - ( NR-1 )*KA1;
            if ( UPDATE ) {
               J2T = min( J2, I-2*KA+K-1 );
            } else {
               J2T = J2;
            }
            NRT = ( J2T+KA-1 ) / KA1;
            DO 800 J = J1, J2T, KA1;

               // create nonzero element a(j+ka,j-1) outside the band
               // and store it in WORK(j)

               WORK( J ) = WORK( J )*AB( KA1, J-1 );
               AB( KA1, J-1 ) = WORK( N+J )*AB( KA1, J-1 );
            } // 800

            // generate rotations in 1st set to annihilate elements which
            // have been created outside the band

            if (NRT > 0) dlargv( NRT, AB( KA1, J1 ), INCA, WORK( J1 ), KA1, WORK( N+J1 ), KA1 );
            if ( NR > 0 ) {

               // apply rotations in 1st set from the right

               for (L = 1; L <= KA - 1; L++) { // 810
                  dlartv(NR, AB( L+1, J1 ), INCA, AB( L+2, J1-1 ), INCA, WORK( N+J1 ), WORK( J1 ), KA1 );
               } // 810

               // apply rotations in 1st set from both sides to diagonal
               // blocks

               dlar2v(NR, AB( 1, J1 ), AB( 1, J1-1 ), AB( 2, J1-1 ), INCA, WORK( N+J1 ), WORK( J1 ), KA1 );

            }

            // start applying rotations in 1st set from the left

            DO 820 L = KA - 1, KB - K + 1, -1;
               NRT = ( J2+L-1 ) / KA1;
               J1T = J2 - ( NRT-1 )*KA1;
               if (NRT > 0) dlartv( NRT, AB( KA1-L+1, J1T-KA1+L ), INCA, AB( KA1-L, J1T-KA1+L ), INCA, WORK( N+J1T ), WORK( J1T ), KA1 );
            } // 820

            if ( WANTX ) {

               // post-multiply X by product of rotations in 1st set

               DO 830 J = J1, J2, KA1;
                  drot(NX, X( 1, J ), 1, X( 1, J-1 ), 1, WORK( N+J ), WORK( J ) );
               } // 830
            }
         } // 840

         if ( UPDATE ) {
            if ( I2 > 0 && KBT > 0 ) {

               // create nonzero element a(i+kbt,i+kbt-ka-1) outside the
               // band and store it in WORK(m-kb+i+kbt)

               WORK( M-KB+I+KBT ) = -BB( KBT+1, I )*RA1;
            }
         }

         DO 880 K = KB, 1, -1;
            if ( UPDATE ) {
               J2 = I + K + 1 - max( 2, K+I0-M )*KA1;
            } else {
               J2 = I + K + 1 - max( 1, K+I0-M )*KA1;
            }

            // finish applying rotations in 2nd set from the left

            DO 850 L = KB - K, 1, -1;
               NRT = ( J2+KA+L-1 ) / KA1;
               J1T = J2 - ( NRT-1 )*KA1;
               if (NRT > 0) dlartv( NRT, AB( KA1-L+1, J1T+L-1 ), INCA, AB( KA1-L, J1T+L-1 ), INCA, WORK( N+M-KB+J1T+KA ), WORK( M-KB+J1T+KA ), KA1 );
            } // 850
            NR = ( J2+KA-1 ) / KA1;
            J1 = J2 - ( NR-1 )*KA1;
            DO 860 J = J1, J2, KA1;
               WORK( M-KB+J ) = WORK( M-KB+J+KA );
               WORK( N+M-KB+J ) = WORK( N+M-KB+J+KA );
            } // 860
            DO 870 J = J1, J2, KA1;

               // create nonzero element a(j+ka,j-1) outside the band
               // and store it in WORK(m-kb+j)

               WORK( M-KB+J ) = WORK( M-KB+J )*AB( KA1, J-1 );
               AB( KA1, J-1 ) = WORK( N+M-KB+J )*AB( KA1, J-1 );
            } // 870
            if ( UPDATE ) {
               if (I+K > KA1 && K <= KBT) WORK( M-KB+I+K-KA ) = WORK( M-KB+I+K );
            }
         } // 880

         DO 920 K = KB, 1, -1;
            J2 = I + K + 1 - max( 1, K+I0-M )*KA1;
            NR = ( J2+KA-1 ) / KA1;
            J1 = J2 - ( NR-1 )*KA1;
            if ( NR > 0 ) {

               // generate rotations in 2nd set to annihilate elements
               // which have been created outside the band

               dlargv(NR, AB( KA1, J1 ), INCA, WORK( M-KB+J1 ), KA1, WORK( N+M-KB+J1 ), KA1 );

               // apply rotations in 2nd set from the right

               for (L = 1; L <= KA - 1; L++) { // 890
                  dlartv(NR, AB( L+1, J1 ), INCA, AB( L+2, J1-1 ), INCA, WORK( N+M-KB+J1 ), WORK( M-KB+J1 ), KA1 );
               } // 890

               // apply rotations in 2nd set from both sides to diagonal
               // blocks

               dlar2v(NR, AB( 1, J1 ), AB( 1, J1-1 ), AB( 2, J1-1 ), INCA, WORK( N+M-KB+J1 ), WORK( M-KB+J1 ), KA1 );

            }

            // start applying rotations in 2nd set from the left

            DO 900 L = KA - 1, KB - K + 1, -1;
               NRT = ( J2+L-1 ) / KA1;
               J1T = J2 - ( NRT-1 )*KA1;
               if (NRT > 0) dlartv( NRT, AB( KA1-L+1, J1T-KA1+L ), INCA, AB( KA1-L, J1T-KA1+L ), INCA, WORK( N+M-KB+J1T ), WORK( M-KB+J1T ), KA1 );
            } // 900

            if ( WANTX ) {

               // post-multiply X by product of rotations in 2nd set

               DO 910 J = J1, J2, KA1;
                  drot(NX, X( 1, J ), 1, X( 1, J-1 ), 1, WORK( N+M-KB+J ), WORK( M-KB+J ) );
               } // 910
            }
         } // 920

         for (K = 1; K <= KB - 1; K++) { // 940
            J2 = I + K + 1 - max( 1, K+I0-M+1 )*KA1;

            // finish applying rotations in 1st set from the left

            DO 930 L = KB - K, 1, -1;
               NRT = ( J2+L-1 ) / KA1;
               J1T = J2 - ( NRT-1 )*KA1;
               if (NRT > 0) dlartv( NRT, AB( KA1-L+1, J1T-KA1+L ), INCA, AB( KA1-L, J1T-KA1+L ), INCA, WORK( N+J1T ), WORK( J1T ), KA1 );
            } // 930
         } // 940

         if ( KB > 1 ) {
            DO 950 J = 2, min( I+KB, M ) - 2*KA - 1;
               WORK( N+J ) = WORK( N+J+KA );
               WORK( J ) = WORK( J+KA );
            } // 950
         }

      }

      GO TO 490;
      }
