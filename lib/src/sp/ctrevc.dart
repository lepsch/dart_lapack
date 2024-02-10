      void ctrevc(SIDE, HOWMNY, SELECT, N, final Matrix<double> T, final int LDT, final Matrix<double> VL, final int LDVL, final Matrix<double> VR, final int LDVR, MM, M, WORK, RWORK, Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             HOWMNY, SIDE;
      int                INFO, LDT, LDVL, LDVR, M, MM, N;
      bool               SELECT( * );
      double               RWORK( * );
      Complex            T( LDT, * ), VL( LDVL, * ), VR( LDVR, * ), WORK( * );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      Complex            CMZERO, CMONE;
      const              CMZERO = ( 0.0, 0.0 ), CMONE = ( 1.0, 0.0 ) ;
      bool               ALLV, BOTHV, LEFTV, OVER, RIGHTV, SOMEV;
      int                I, II, IS, J, K, KI;
      double               OVFL, REMAX, SCALE, SMIN, SMLNUM, ULP, UNFL;
      Complex            CDUM;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                ICAMAX;
      //- REAL               SCASUM, SLAMCH;
      // EXTERNAL lsame, ICAMAX, SCASUM, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CCOPY, CGEMV, CLATRS, CSSCAL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, CMPLX, CONJG, MAX, REAL
      // ..
      // .. Statement Functions ..
      double               CABS1;
      // ..
      // .. Statement Function definitions ..
      CABS1[CDUM] = ( double( CDUM ) ).abs() + ( AIMAG( CDUM ) ).abs();

      // Decode and test the input parameters

      BOTHV = lsame( SIDE, 'B' );
      RIGHTV = lsame( SIDE, 'R' ) || BOTHV;
      LEFTV = lsame( SIDE, 'L' ) || BOTHV;

      ALLV = lsame( HOWMNY, 'A' );
      OVER = lsame( HOWMNY, 'B' );
      SOMEV = lsame( HOWMNY, 'S' );

      // Set M to the number of columns required to store the selected
      // eigenvectors.

      if ( SOMEV ) {
         M = 0;
         for (J = 1; J <= N; J++) { // 10
            if( SELECT( J ) ) M = M + 1;
         } // 10
      } else {
         M = N;
      }

      INFO = 0;
      if ( !RIGHTV && !LEFTV ) {
         INFO = -1;
      } else if ( !ALLV && !OVER && !SOMEV ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -4;
      } else if ( LDT < max( 1, N ) ) {
         INFO = -6;
      } else if ( LDVL < 1 || ( LEFTV && LDVL < N ) ) {
         INFO = -8;
      } else if ( LDVR < 1 || ( RIGHTV && LDVR < N ) ) {
         INFO = -10;
      } else if ( MM < M ) {
         INFO = -11;
      }
      if ( INFO != 0 ) {
         xerbla('CTREVC', -INFO );
         return;
      }

      // Quick return if possible.

      if (N == 0) return;

      // Set the constants to control overflow.

      UNFL = SLAMCH( 'Safe minimum' );
      OVFL = ONE / UNFL;
      ULP = SLAMCH( 'Precision' );
      SMLNUM = UNFL*( N / ULP );

      // Store the diagonal elements of T in working array WORK.

      for (I = 1; I <= N; I++) { // 20
         WORK[I+N] = T( I, I );
      } // 20

      // Compute 1-norm of each column of strictly upper triangular
      // part of T to control overflow in triangular solver.

      RWORK[1] = ZERO;
      for (J = 2; J <= N; J++) { // 30
         RWORK[J] = SCASUM( J-1, T( 1, J ), 1 );
      } // 30

      if ( RIGHTV ) {

         // Compute right eigenvectors.

         IS = M;
         for (KI = N; KI >= 1; KI--) { // 80

            if ( SOMEV ) {
               if( !SELECT( KI ) ) GO TO 80;
            }
            SMIN = max( ULP*( CABS1( T( KI, KI ) ) ), SMLNUM );

            WORK[1] = CMONE;

            // Form right-hand side.

            for (K = 1; K <= KI - 1; K++) { // 40
               WORK[K] = -T( K, KI );
            } // 40

            // Solve the triangular system:
            //    (T(1:KI-1,1:KI-1) - T(KI,KI))*X = SCALE*WORK.

            for (K = 1; K <= KI - 1; K++) { // 50
               T[K][K] = T( K, K ) - T( KI, KI );
               if[CABS1( T( K, K ) ) < SMIN ) T( K][K] = SMIN;
            } // 50

            if ( KI > 1 ) {
               clatrs('Upper', 'No transpose', 'Non-unit', 'Y', KI-1, T, LDT, WORK( 1 ), SCALE, RWORK, INFO );
               WORK[KI] = SCALE;
            }

            // Copy the vector x or Q*x to VR and normalize.

            if ( !OVER ) {
               ccopy(KI, WORK( 1 ), 1, VR( 1, IS ), 1 );

               II = ICAMAX( KI, VR( 1, IS ), 1 );
               REMAX = ONE / CABS1( VR( II, IS ) );
               csscal(KI, REMAX, VR( 1, IS ), 1 );

               for (K = KI + 1; K <= N; K++) { // 60
                  VR[K][IS] = CMZERO;
               } // 60
            } else {
               if (KI > 1) cgemv( 'N', N, KI-1, CMONE, VR, LDVR, WORK( 1 ), 1, CMPLX( SCALE ), VR( 1, KI ), 1 );

               II = ICAMAX( N, VR( 1, KI ), 1 );
               REMAX = ONE / CABS1( VR( II, KI ) );
               csscal(N, REMAX, VR( 1, KI ), 1 );
            }

            // Set back the original diagonal elements of T.

            for (K = 1; K <= KI - 1; K++) { // 70
               T[K][K] = WORK( K+N );
            } // 70

            IS = IS - 1;
         } // 80
      }

      if ( LEFTV ) {

         // Compute left eigenvectors.

         IS = 1;
         for (KI = 1; KI <= N; KI++) { // 130

            if ( SOMEV ) {
               if( !SELECT( KI ) ) GO TO 130;
            }
            SMIN = max( ULP*( CABS1( T( KI, KI ) ) ), SMLNUM );

            WORK[N] = CMONE;

            // Form right-hand side.

            for (K = KI + 1; K <= N; K++) { // 90
               WORK[K] = -CONJG( T( KI, K ) );
            } // 90

            // Solve the triangular system:
            //    (T(KI+1:N,KI+1:N) - T(KI,KI))**H*X = SCALE*WORK.

            for (K = KI + 1; K <= N; K++) { // 100
               T[K][K] = T( K, K ) - T( KI, KI );
               if[CABS1( T( K, K ) ) < SMIN ) T( K][K] = SMIN;
            } // 100

            if ( KI < N ) {
               clatrs('Upper', 'Conjugate transpose', 'Non-unit', 'Y', N-KI, T( KI+1, KI+1 ), LDT, WORK( KI+1 ), SCALE, RWORK, INFO );
               WORK[KI] = SCALE;
            }

            // Copy the vector x or Q*x to VL and normalize.

            if ( !OVER ) {
               ccopy(N-KI+1, WORK( KI ), 1, VL( KI, IS ), 1 );

               II = ICAMAX( N-KI+1, VL( KI, IS ), 1 ) + KI - 1;
               REMAX = ONE / CABS1( VL( II, IS ) );
               csscal(N-KI+1, REMAX, VL( KI, IS ), 1 );

               for (K = 1; K <= KI - 1; K++) { // 110
                  VL[K][IS] = CMZERO;
               } // 110
            } else {
               if (KI < N) cgemv( 'N', N, N-KI, CMONE, VL( 1, KI+1 ), LDVL, WORK( KI+1 ), 1, CMPLX( SCALE ), VL( 1, KI ), 1 );

               II = ICAMAX( N, VL( 1, KI ), 1 );
               REMAX = ONE / CABS1( VL( II, KI ) );
               csscal(N, REMAX, VL( 1, KI ), 1 );
            }

            // Set back the original diagonal elements of T.

            for (K = KI + 1; K <= N; K++) { // 120
               T[K][K] = WORK( K+N );
            } // 120

            IS = IS + 1;
         } // 130
      }

      }
