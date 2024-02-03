      SUBROUTINE CTREVC( SIDE, HOWMNY, SELECT, N, T, LDT, VL, LDVL, VR, LDVR, MM, M, WORK, RWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             HOWMNY, SIDE;
      int                INFO, LDT, LDVL, LDVR, M, MM, N;
      // ..
      // .. Array Arguments ..
      bool               SELECT( * );
      REAL               RWORK( * )
      COMPLEX            T( LDT, * ), VL( LDVL, * ), VR( LDVR, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      COMPLEX            CMZERO, CMONE
      const              CMZERO = ( 0.0E+0, 0.0E+0 ), CMONE = ( 1.0E+0, 0.0E+0 ) ;
      // ..
      // .. Local Scalars ..
      bool               ALLV, BOTHV, LEFTV, OVER, RIGHTV, SOMEV;
      int                I, II, IS, J, K, KI;
      REAL               OVFL, REMAX, SCALE, SMIN, SMLNUM, ULP, UNFL
      COMPLEX            CDUM
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ICAMAX;
      REAL               SCASUM, SLAMCH
      // EXTERNAL LSAME, ICAMAX, SCASUM, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CCOPY, CGEMV, CLATRS, CSSCAL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, CMPLX, CONJG, MAX, REAL
      // ..
      // .. Statement Functions ..
      REAL               CABS1
      // ..
      // .. Statement Function definitions ..
      CABS1( CDUM ) = ABS( REAL( CDUM ) ) + ABS( AIMAG( CDUM ) )
      // ..
      // .. Executable Statements ..

      // Decode and test the input parameters

      BOTHV = LSAME( SIDE, 'B' )
      RIGHTV = LSAME( SIDE, 'R' ) .OR. BOTHV
      LEFTV = LSAME( SIDE, 'L' ) .OR. BOTHV

      ALLV = LSAME( HOWMNY, 'A' )
      OVER = LSAME( HOWMNY, 'B' )
      SOMEV = LSAME( HOWMNY, 'S' )

      // Set M to the number of columns required to store the selected
      // eigenvectors.

      if ( SOMEV ) {
         M = 0
         for (J = 1; J <= N; J++) { // 10
            IF( SELECT( J ) ) M = M + 1
         } // 10
      } else {
         M = N
      }

      INFO = 0
      if ( .NOT.RIGHTV .AND. .NOT.LEFTV ) {
         INFO = -1
      } else if ( .NOT.ALLV .AND. .NOT.OVER .AND. .NOT.SOMEV ) {
         INFO = -2
      } else if ( N.LT.0 ) {
         INFO = -4
      } else if ( LDT.LT.MAX( 1, N ) ) {
         INFO = -6
      } else if ( LDVL.LT.1 .OR. ( LEFTV .AND. LDVL.LT.N ) ) {
         INFO = -8
      } else if ( LDVR.LT.1 .OR. ( RIGHTV .AND. LDVR.LT.N ) ) {
         INFO = -10
      } else if ( MM.LT.M ) {
         INFO = -11
      }
      if ( INFO.NE.0 ) {
         xerbla('CTREVC', -INFO );
         RETURN
      }

      // Quick return if possible.

      IF( N.EQ.0 ) RETURN

      // Set the constants to control overflow.

      UNFL = SLAMCH( 'Safe minimum' )
      OVFL = ONE / UNFL
      ULP = SLAMCH( 'Precision' )
      SMLNUM = UNFL*( N / ULP )

      // Store the diagonal elements of T in working array WORK.

      for (I = 1; I <= N; I++) { // 20
         WORK( I+N ) = T( I, I )
      } // 20

      // Compute 1-norm of each column of strictly upper triangular
      // part of T to control overflow in triangular solver.

      RWORK( 1 ) = ZERO
      for (J = 2; J <= N; J++) { // 30
         RWORK( J ) = SCASUM( J-1, T( 1, J ), 1 )
      } // 30

      if ( RIGHTV ) {

         // Compute right eigenvectors.

         IS = M
         DO 80 KI = N, 1, -1

            if ( SOMEV ) {
               IF( .NOT.SELECT( KI ) ) GO TO 80
            }
            SMIN = MAX( ULP*( CABS1( T( KI, KI ) ) ), SMLNUM )

            WORK( 1 ) = CMONE

            // Form right-hand side.

            DO 40 K = 1, KI - 1
               WORK( K ) = -T( K, KI )
            } // 40

            // Solve the triangular system:
               // (T(1:KI-1,1:KI-1) - T(KI,KI))*X = SCALE*WORK.

            DO 50 K = 1, KI - 1
               T( K, K ) = T( K, K ) - T( KI, KI )
               IF( CABS1( T( K, K ) ).LT.SMIN ) T( K, K ) = SMIN
            } // 50

            if ( KI.GT.1 ) {
               clatrs('Upper', 'No transpose', 'Non-unit', 'Y', KI-1, T, LDT, WORK( 1 ), SCALE, RWORK, INFO );
               WORK( KI ) = SCALE
            }

            // Copy the vector x or Q*x to VR and normalize.

            if ( .NOT.OVER ) {
               ccopy(KI, WORK( 1 ), 1, VR( 1, IS ), 1 );

               II = ICAMAX( KI, VR( 1, IS ), 1 )
               REMAX = ONE / CABS1( VR( II, IS ) )
               csscal(KI, REMAX, VR( 1, IS ), 1 );

               DO 60 K = KI + 1, N
                  VR( K, IS ) = CMZERO
               } // 60
            } else {
               IF( KI.GT.1 ) CALL CGEMV( 'N', N, KI-1, CMONE, VR, LDVR, WORK( 1 ), 1, CMPLX( SCALE ), VR( 1, KI ), 1 )

               II = ICAMAX( N, VR( 1, KI ), 1 )
               REMAX = ONE / CABS1( VR( II, KI ) )
               csscal(N, REMAX, VR( 1, KI ), 1 );
            }

            // Set back the original diagonal elements of T.

            DO 70 K = 1, KI - 1
               T( K, K ) = WORK( K+N )
            } // 70

            IS = IS - 1
         } // 80
      }

      if ( LEFTV ) {

         // Compute left eigenvectors.

         IS = 1
         for (KI = 1; KI <= N; KI++) { // 130

            if ( SOMEV ) {
               IF( .NOT.SELECT( KI ) ) GO TO 130
            }
            SMIN = MAX( ULP*( CABS1( T( KI, KI ) ) ), SMLNUM )

            WORK( N ) = CMONE

            // Form right-hand side.

            DO 90 K = KI + 1, N
               WORK( K ) = -CONJG( T( KI, K ) )
            } // 90

            // Solve the triangular system:
               // (T(KI+1:N,KI+1:N) - T(KI,KI))**H*X = SCALE*WORK.

            DO 100 K = KI + 1, N
               T( K, K ) = T( K, K ) - T( KI, KI )
               IF( CABS1( T( K, K ) ).LT.SMIN ) T( K, K ) = SMIN
            } // 100

            if ( KI.LT.N ) {
               clatrs('Upper', 'Conjugate transpose', 'Non-unit', 'Y', N-KI, T( KI+1, KI+1 ), LDT, WORK( KI+1 ), SCALE, RWORK, INFO );
               WORK( KI ) = SCALE
            }

            // Copy the vector x or Q*x to VL and normalize.

            if ( .NOT.OVER ) {
               ccopy(N-KI+1, WORK( KI ), 1, VL( KI, IS ), 1 );

               II = ICAMAX( N-KI+1, VL( KI, IS ), 1 ) + KI - 1
               REMAX = ONE / CABS1( VL( II, IS ) )
               csscal(N-KI+1, REMAX, VL( KI, IS ), 1 );

               DO 110 K = 1, KI - 1
                  VL( K, IS ) = CMZERO
               } // 110
            } else {
               IF( KI.LT.N ) CALL CGEMV( 'N', N, N-KI, CMONE, VL( 1, KI+1 ), LDVL, WORK( KI+1 ), 1, CMPLX( SCALE ), VL( 1, KI ), 1 )

               II = ICAMAX( N, VL( 1, KI ), 1 )
               REMAX = ONE / CABS1( VL( II, KI ) )
               csscal(N, REMAX, VL( 1, KI ), 1 );
            }

            // Set back the original diagonal elements of T.

            DO 120 K = KI + 1, N
               T( K, K ) = WORK( K+N )
            } // 120

            IS = IS + 1
         } // 130
      }

      RETURN

      // End of CTREVC

      }
