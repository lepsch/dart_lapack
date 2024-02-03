      SUBROUTINE ZTREVC3( SIDE, HOWMNY, SELECT, N, T, LDT, VL, LDVL, VR, LDVR, MM, M, WORK, LWORK, RWORK, LRWORK, INFO);
      IMPLICIT NONE;

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             HOWMNY, SIDE;
      int                INFO, LDT, LDVL, LDVR, LWORK, LRWORK, M, MM, N;
      // ..
      // .. Array Arguments ..
      bool               SELECT( * );
      double             RWORK( * );
      COMPLEX*16         T( LDT, * ), VL( LDVL, * ), VR( LDVR, * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      COMPLEX*16         CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE  = ( 1.0, 0.0 ) ;
      int                NBMIN, NBMAX;
      const              NBMIN = 8, NBMAX = 128 ;
      // ..
      // .. Local Scalars ..
      bool               ALLV, BOTHV, LEFTV, LQUERY, OVER, RIGHTV, SOMEV;
      int                I, II, IS, J, K, KI, IV, MAXWRK, NB;
      double             OVFL, REMAX, SCALE, SMIN, SMLNUM, ULP, UNFL;
      COMPLEX*16         CDUM;
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV, IZAMAX;
      double             DLAMCH, DZASUM;
      // EXTERNAL LSAME, ILAENV, IZAMAX, DLAMCH, DZASUM
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZCOPY, ZDSCAL, ZGEMV, ZLATRS, ZGEMM, ZLASET, ZLACPY
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DCMPLX, CONJG, DIMAG, MAX
      // ..
      // .. Statement Functions ..
      double             CABS1;
      // ..
      // .. Statement Function definitions ..
      CABS1( CDUM ) = ABS( DBLE( CDUM ) ) + ABS( DIMAG( CDUM ) );
      // ..
      // .. Executable Statements ..

      // Decode and test the input parameters

      BOTHV  = LSAME( SIDE, 'B' );
      RIGHTV = LSAME( SIDE, 'R' ) || BOTHV;
      LEFTV  = LSAME( SIDE, 'L' ) || BOTHV;

      ALLV  = LSAME( HOWMNY, 'A' );
      OVER  = LSAME( HOWMNY, 'B' );
      SOMEV = LSAME( HOWMNY, 'S' );

      // Set M to the number of columns required to store the selected
      // eigenvectors.

      if ( SOMEV ) {
         M = 0;
         for (J = 1; J <= N; J++) { // 10
            IF( SELECT( J ) ) M = M + 1;
         } // 10
      } else {
         M = N;
      }

      INFO = 0;
      NB = ILAENV( 1, 'ZTREVC', SIDE // HOWMNY, N, -1, -1, -1 );
      MAXWRK = MAX( 1, N + 2*N*NB );
      WORK(1) = MAXWRK;
      RWORK(1) = MAX( 1, N );
      LQUERY = ( LWORK == -1 || LRWORK == -1 );
      if ( !RIGHTV && !LEFTV ) {
         INFO = -1;
      } else if ( !ALLV && !OVER && !SOMEV ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -4;
      } else if ( LDT < MAX( 1, N ) ) {
         INFO = -6;
      } else if ( LDVL < 1 || ( LEFTV && LDVL < N ) ) {
         INFO = -8;
      } else if ( LDVR < 1 || ( RIGHTV && LDVR < N ) ) {
         INFO = -10;
      } else if ( MM < M ) {
         INFO = -11;
      } else if ( LWORK < MAX( 1, 2*N ) && !LQUERY ) {
         INFO = -14;
      } else if ( LRWORK < MAX( 1, N ) && !LQUERY ) {
         INFO = -16;
      }
      if ( INFO != 0 ) {
         xerbla('ZTREVC3', -INFO );
         RETURN;
      } else if ( LQUERY ) {
         RETURN;
      }

      // Quick return if possible.

      if (N == 0) RETURN;

      // Use blocked version of back-transformation if sufficient workspace.
      // Zero-out the workspace to avoid potential NaN propagation.

      if ( OVER && LWORK >= N + 2*N*NBMIN ) {
         NB = (LWORK - N) / (2*N);
         NB = MIN( NB, NBMAX );
         zlaset('F', N, 1+2*NB, CZERO, CZERO, WORK, N );
      } else {
         NB = 1;
      }

      // Set the constants to control overflow.

      UNFL = DLAMCH( 'Safe minimum' );
      OVFL = ONE / UNFL;
      ULP = DLAMCH( 'Precision' );
      SMLNUM = UNFL*( N / ULP );

      // Store the diagonal elements of T in working array WORK.

      for (I = 1; I <= N; I++) { // 20
         WORK( I ) = T( I, I );
      } // 20

      // Compute 1-norm of each column of strictly upper triangular
      // part of T to control overflow in triangular solver.

      RWORK( 1 ) = ZERO;
      for (J = 2; J <= N; J++) { // 30
         RWORK( J ) = DZASUM( J-1, T( 1, J ), 1 );
      } // 30

      if ( RIGHTV ) {

         // ============================================================
         // Compute right eigenvectors.

         // IV is index of column in current block.
         // Non-blocked version always uses IV=NB=1;
         // blocked     version starts with IV=NB, goes down to 1.
         // (Note the "0-th" column is used to store the original diagonal.)
         IV = NB;
         IS = M;
         DO 80 KI = N, 1, -1;
            if ( SOMEV ) {
               IF( !SELECT( KI ) ) GO TO 80;
            }
            SMIN = MAX( ULP*( CABS1( T( KI, KI ) ) ), SMLNUM );

            // --------------------------------------------------------
            // Complex right eigenvector

            WORK( KI + IV*N ) = CONE;

            // Form right-hand side.

            for (K = 1; K <= KI - 1; K++) { // 40
               WORK( K + IV*N ) = -T( K, KI );
            } // 40

            // Solve upper triangular system:
            // [ T(1:KI-1,1:KI-1) - T(KI,KI) ]*X = SCALE*WORK.

            for (K = 1; K <= KI - 1; K++) { // 50
               T( K, K ) = T( K, K ) - T( KI, KI );
               IF( CABS1( T( K, K ) ) < SMIN ) T( K, K ) = SMIN;
            } // 50

            if ( KI > 1 ) {
               zlatrs('Upper', 'No transpose', 'Non-unit', 'Y', KI-1, T, LDT, WORK( 1 + IV*N ), SCALE, RWORK, INFO );
               WORK( KI + IV*N ) = SCALE;
            }

            // Copy the vector x or Q*x to VR and normalize.

            if ( !OVER ) {
               // ------------------------------
               // no back-transform: copy x to VR and normalize.
               zcopy(KI, WORK( 1 + IV*N ), 1, VR( 1, IS ), 1 );

               II = IZAMAX( KI, VR( 1, IS ), 1 );
               REMAX = ONE / CABS1( VR( II, IS ) );
               zdscal(KI, REMAX, VR( 1, IS ), 1 );

               for (K = KI + 1; K <= N; K++) { // 60
                  VR( K, IS ) = CZERO;
               } // 60

            } else if ( NB == 1 ) {
               // ------------------------------
               // version 1: back-transform each vector with GEMV, Q*x.
               if (KI > 1) CALL ZGEMV( 'N', N, KI-1, CONE, VR, LDVR, WORK( 1 + IV*N ), 1, DCMPLX( SCALE ), VR( 1, KI ), 1 );

               II = IZAMAX( N, VR( 1, KI ), 1 );
               REMAX = ONE / CABS1( VR( II, KI ) );
               zdscal(N, REMAX, VR( 1, KI ), 1 );

            } else {
               // ------------------------------
               // version 2: back-transform block of vectors with GEMM
               // zero out below vector
               for (K = KI + 1; K <= N; K++) {
                  WORK( K + IV*N ) = CZERO;
               }

               // Columns IV:NB of work are valid vectors.
               // When the number of vectors stored reaches NB,
               // or if this was last vector, do the GEMM
               if ( (IV == 1) || (KI == 1) ) {
                  zgemm('N', 'N', N, NB-IV+1, KI+NB-IV, CONE, VR, LDVR, WORK( 1 + (IV)*N    ), N, CZERO, WORK( 1 + (NB+IV)*N ), N );
                  // normalize vectors
                  for (K = IV; K <= NB; K++) {
                     II = IZAMAX( N, WORK( 1 + (NB+K)*N ), 1 );
                     REMAX = ONE / CABS1( WORK( II + (NB+K)*N ) );
                     zdscal(N, REMAX, WORK( 1 + (NB+K)*N ), 1 );
                  }
                  zlacpy('F', N, NB-IV+1, WORK( 1 + (NB+IV)*N ), N, VR( 1, KI ), LDVR );
                  IV = NB;
               } else {
                  IV = IV - 1;
               }
            }

            // Restore the original diagonal elements of T.

            for (K = 1; K <= KI - 1; K++) { // 70
               T( K, K ) = WORK( K );
            } // 70

            IS = IS - 1;
         } // 80
      }

      if ( LEFTV ) {

         // ============================================================
         // Compute left eigenvectors.

         // IV is index of column in current block.
         // Non-blocked version always uses IV=1;
         // blocked     version starts with IV=1, goes up to NB.
         // (Note the "0-th" column is used to store the original diagonal.)
         IV = 1;
         IS = 1;
         for (KI = 1; KI <= N; KI++) { // 130

            if ( SOMEV ) {
               IF( !SELECT( KI ) ) GO TO 130;
            }
            SMIN = MAX( ULP*( CABS1( T( KI, KI ) ) ), SMLNUM );

            // --------------------------------------------------------
            // Complex left eigenvector

            WORK( KI + IV*N ) = CONE;

            // Form right-hand side.

            for (K = KI + 1; K <= N; K++) { // 90
               WORK( K + IV*N ) = -CONJG( T( KI, K ) );
            } // 90

            // Solve conjugate-transposed triangular system:
            // [ T(KI+1:N,KI+1:N) - T(KI,KI) ]**H * X = SCALE*WORK.

            for (K = KI + 1; K <= N; K++) { // 100
               T( K, K ) = T( K, K ) - T( KI, KI );
               IF( CABS1( T( K, K ) ) < SMIN ) T( K, K ) = SMIN;
            } // 100

            if ( KI < N ) {
               zlatrs('Upper', 'Conjugate transpose', 'Non-unit', 'Y', N-KI, T( KI+1, KI+1 ), LDT, WORK( KI+1 + IV*N ), SCALE, RWORK, INFO );
               WORK( KI + IV*N ) = SCALE;
            }

            // Copy the vector x or Q*x to VL and normalize.

            if ( !OVER ) {
               // ------------------------------
               // no back-transform: copy x to VL and normalize.
               zcopy(N-KI+1, WORK( KI + IV*N ), 1, VL(KI,IS), 1 );

               II = IZAMAX( N-KI+1, VL( KI, IS ), 1 ) + KI - 1;
               REMAX = ONE / CABS1( VL( II, IS ) );
               zdscal(N-KI+1, REMAX, VL( KI, IS ), 1 );

               for (K = 1; K <= KI - 1; K++) { // 110
                  VL( K, IS ) = CZERO;
               } // 110

            } else if ( NB == 1 ) {
               // ------------------------------
               // version 1: back-transform each vector with GEMV, Q*x.
               if (KI < N) CALL ZGEMV( 'N', N, N-KI, CONE, VL( 1, KI+1 ), LDVL, WORK( KI+1 + IV*N ), 1, DCMPLX( SCALE ), VL( 1, KI ), 1 );

               II = IZAMAX( N, VL( 1, KI ), 1 );
               REMAX = ONE / CABS1( VL( II, KI ) );
               zdscal(N, REMAX, VL( 1, KI ), 1 );

            } else {
               // ------------------------------
               // version 2: back-transform block of vectors with GEMM
               // zero out above vector
               // could go from KI-NV+1 to KI-1
               for (K = 1; K <= KI - 1; K++) {
                  WORK( K + IV*N ) = CZERO;
               }

               // Columns 1:IV of work are valid vectors.
               // When the number of vectors stored reaches NB,
               // or if this was last vector, do the GEMM
               if ( (IV == NB) || (KI == N) ) {
                  zgemm('N', 'N', N, IV, N-KI+IV, CONE, VL( 1, KI-IV+1 ), LDVL, WORK( KI-IV+1 + (1)*N ), N, CZERO, WORK( 1 + (NB+1)*N ), N );
                  // normalize vectors
                  for (K = 1; K <= IV; K++) {
                     II = IZAMAX( N, WORK( 1 + (NB+K)*N ), 1 );
                     REMAX = ONE / CABS1( WORK( II + (NB+K)*N ) );
                     zdscal(N, REMAX, WORK( 1 + (NB+K)*N ), 1 );
                  }
                  zlacpy('F', N, IV, WORK( 1 + (NB+1)*N ), N, VL( 1, KI-IV+1 ), LDVL );
                  IV = 1;
               } else {
                  IV = IV + 1;
               }
            }

            // Restore the original diagonal elements of T.

            for (K = KI + 1; K <= N; K++) { // 120
               T( K, K ) = WORK( K );
            } // 120

            IS = IS + 1;
         } // 130
      }

      RETURN;

      // End of ZTREVC3

      }
