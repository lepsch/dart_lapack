      SUBROUTINE ZTREVC( SIDE, HOWMNY, SELECT, N, T, LDT, VL, LDVL, VR, LDVR, MM, M, WORK, RWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             HOWMNY, SIDE;
      int                INFO, LDT, LDVL, LDVR, M, MM, N;
      // ..
      // .. Array Arguments ..
      bool               SELECT( * );
      double             RWORK( * );
      COMPLEX*16         T( LDT, * ), VL( LDVL, * ), VR( LDVR, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      COMPLEX*16         CMZERO, CMONE
      const              CMZERO = ( 0.0D+0, 0.0D+0 ), CMONE = ( 1.0D+0, 0.0D+0 ) ;
      // ..
      // .. Local Scalars ..
      bool               ALLV, BOTHV, LEFTV, OVER, RIGHTV, SOMEV;
      int                I, II, IS, J, K, KI;
      double             OVFL, REMAX, SCALE, SMIN, SMLNUM, ULP, UNFL;
      COMPLEX*16         CDUM
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                IZAMAX;
      double             DLAMCH, DZASUM;
      // EXTERNAL LSAME, IZAMAX, DLAMCH, DZASUM
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZCOPY, ZDSCAL, ZGEMV, ZLATRS
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DCMPLX, DCONJG, DIMAG, MAX
      // ..
      // .. Statement Functions ..
      double             CABS1;
      // ..
      // .. Statement Function definitions ..
      CABS1( CDUM ) = ABS( DBLE( CDUM ) ) + ABS( DIMAG( CDUM ) )
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
   10    CONTINUE
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
         xerbla('ZTREVC', -INFO );
         RETURN
      }

      // Quick return if possible.

      IF( N.EQ.0 ) RETURN

      // Set the constants to control overflow.

      UNFL = DLAMCH( 'Safe minimum' )
      OVFL = ONE / UNFL
      ULP = DLAMCH( 'Precision' )
      SMLNUM = UNFL*( N / ULP )

      // Store the diagonal elements of T in working array WORK.

      for (I = 1; I <= N; I++) { // 20
         WORK( I+N ) = T( I, I )
   20 CONTINUE

      // Compute 1-norm of each column of strictly upper triangular
      // part of T to control overflow in triangular solver.

      RWORK( 1 ) = ZERO
      for (J = 2; J <= N; J++) { // 30
         RWORK( J ) = DZASUM( J-1, T( 1, J ), 1 )
   30 CONTINUE

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
   40       CONTINUE

            // Solve the triangular system:
               // (T(1:KI-1,1:KI-1) - T(KI,KI))*X = SCALE*WORK.

            DO 50 K = 1, KI - 1
               T( K, K ) = T( K, K ) - T( KI, KI )
               IF( CABS1( T( K, K ) ).LT.SMIN ) T( K, K ) = SMIN
   50       CONTINUE

            if ( KI.GT.1 ) {
               zlatrs('Upper', 'No transpose', 'Non-unit', 'Y', KI-1, T, LDT, WORK( 1 ), SCALE, RWORK, INFO );
               WORK( KI ) = SCALE
            }

            // Copy the vector x or Q*x to VR and normalize.

            if ( .NOT.OVER ) {
               zcopy(KI, WORK( 1 ), 1, VR( 1, IS ), 1 );

               II = IZAMAX( KI, VR( 1, IS ), 1 )
               REMAX = ONE / CABS1( VR( II, IS ) )
               zdscal(KI, REMAX, VR( 1, IS ), 1 );

               DO 60 K = KI + 1, N
                  VR( K, IS ) = CMZERO
   60          CONTINUE
            } else {
               IF( KI.GT.1 ) CALL ZGEMV( 'N', N, KI-1, CMONE, VR, LDVR, WORK( 1 ), 1, DCMPLX( SCALE ), VR( 1, KI ), 1 )

               II = IZAMAX( N, VR( 1, KI ), 1 )
               REMAX = ONE / CABS1( VR( II, KI ) )
               zdscal(N, REMAX, VR( 1, KI ), 1 );
            }

            // Set back the original diagonal elements of T.

            DO 70 K = 1, KI - 1
               T( K, K ) = WORK( K+N )
   70       CONTINUE

            IS = IS - 1
   80    CONTINUE
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
               WORK( K ) = -DCONJG( T( KI, K ) )
   90       CONTINUE

            // Solve the triangular system:
               // (T(KI+1:N,KI+1:N) - T(KI,KI))**H * X = SCALE*WORK.

            DO 100 K = KI + 1, N
               T( K, K ) = T( K, K ) - T( KI, KI )
               IF( CABS1( T( K, K ) ).LT.SMIN ) T( K, K ) = SMIN
  100       CONTINUE

            if ( KI.LT.N ) {
               zlatrs('Upper', 'Conjugate transpose', 'Non-unit', 'Y', N-KI, T( KI+1, KI+1 ), LDT, WORK( KI+1 ), SCALE, RWORK, INFO );
               WORK( KI ) = SCALE
            }

            // Copy the vector x or Q*x to VL and normalize.

            if ( .NOT.OVER ) {
               zcopy(N-KI+1, WORK( KI ), 1, VL( KI, IS ), 1 );

               II = IZAMAX( N-KI+1, VL( KI, IS ), 1 ) + KI - 1
               REMAX = ONE / CABS1( VL( II, IS ) )
               zdscal(N-KI+1, REMAX, VL( KI, IS ), 1 );

               DO 110 K = 1, KI - 1
                  VL( K, IS ) = CMZERO
  110          CONTINUE
            } else {
               IF( KI.LT.N ) CALL ZGEMV( 'N', N, N-KI, CMONE, VL( 1, KI+1 ), LDVL, WORK( KI+1 ), 1, DCMPLX( SCALE ), VL( 1, KI ), 1 )

               II = IZAMAX( N, VL( 1, KI ), 1 )
               REMAX = ONE / CABS1( VL( II, KI ) )
               zdscal(N, REMAX, VL( 1, KI ), 1 );
            }

            // Set back the original diagonal elements of T.

            DO 120 K = KI + 1, N
               T( K, K ) = WORK( K+N )
  120       CONTINUE

            IS = IS + 1
  130    CONTINUE
      }

      RETURN

      // End of ZTREVC

      }
