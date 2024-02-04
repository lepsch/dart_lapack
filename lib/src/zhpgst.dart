      void zhpgst(ITYPE, UPLO, N, AP, BP, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, ITYPE, N;
      // ..
      // .. Array Arguments ..
      Complex         AP( * ), BP( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ONE, HALF;
      const              ONE = 1.0, HALF = 0.5 ;
      Complex         CONE;
      const              CONE = ( 1.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                J, J1, J1J1, JJ, K, K1, K1K1, KK;
      double             AJJ, AKK, BJJ, BKK;
      Complex         CT;
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZAXPY, ZDSCAL, ZHPMV, ZHPR2, ZTPMV, ZTPSV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE
      // ..
      // .. External Functions ..
      //- bool               LSAME;
      //- Complex         ZDOTC;
      // EXTERNAL LSAME, ZDOTC
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0;
      UPPER = LSAME( UPLO, 'U' );
      if ( ITYPE < 1 || ITYPE > 3 ) {
         INFO = -1;
      } else if ( !UPPER && !LSAME( UPLO, 'L' ) ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -3;
      }
      if ( INFO != 0 ) {
         xerbla('ZHPGST', -INFO );
         return;
      }

      if ( ITYPE == 1 ) {
         if ( UPPER ) {

            // Compute inv(U**H)*A*inv(U)

            // J1 and JJ are the indices of A(1,j) and A(j,j)

            JJ = 0;
            for (J = 1; J <= N; J++) { // 10
               J1 = JJ + 1;
               JJ = JJ + J;

               // Compute the j-th column of the upper triangle of A

               AP[JJ] = (AP( JJ )).toDouble();
               BJJ = (BP( JJ )).toDouble();
               ztpsv(UPLO, 'Conjugate transpose', 'Non-unit', J, BP, AP( J1 ), 1 );
               zhpmv(UPLO, J-1, -CONE, AP, BP( J1 ), 1, CONE, AP( J1 ), 1 );
               zdscal(J-1, ONE / BJJ, AP( J1 ), 1 );
               AP[JJ] = ( AP( JJ )-ZDOTC( J-1, AP( J1 ), 1, BP( J1 ), 1 ) ) / BJJ;
            } // 10
         } else {

            // Compute inv(L)*A*inv(L**H)

            // KK and K1K1 are the indices of A(k,k) and A(k+1,k+1)

            KK = 1;
            for (K = 1; K <= N; K++) { // 20
               K1K1 = KK + N - K + 1;

               // Update the lower triangle of A(k:n,k:n)

               AKK = (AP( KK )).toDouble();
               BKK = (BP( KK )).toDouble();
               AKK = AKK / BKK**2;
               AP[KK] = AKK;
               if ( K < N ) {
                  zdscal(N-K, ONE / BKK, AP( KK+1 ), 1 );
                  CT = -HALF*AKK;
                  zaxpy(N-K, CT, BP( KK+1 ), 1, AP( KK+1 ), 1 );
                  zhpr2(UPLO, N-K, -CONE, AP( KK+1 ), 1, BP( KK+1 ), 1, AP( K1K1 ) );
                  zaxpy(N-K, CT, BP( KK+1 ), 1, AP( KK+1 ), 1 );
                  ztpsv(UPLO, 'No transpose', 'Non-unit', N-K, BP( K1K1 ), AP( KK+1 ), 1 );
               }
               KK = K1K1;
            } // 20
         }
      } else {
         if ( UPPER ) {

            // Compute U*A*U**H

            // K1 and KK are the indices of A(1,k) and A(k,k)

            KK = 0;
            for (K = 1; K <= N; K++) { // 30
               K1 = KK + 1;
               KK = KK + K;

               // Update the upper triangle of A(1:k,1:k)

               AKK = (AP( KK )).toDouble();
               BKK = (BP( KK )).toDouble();
               ztpmv(UPLO, 'No transpose', 'Non-unit', K-1, BP, AP( K1 ), 1 );
               CT = HALF*AKK;
               zaxpy(K-1, CT, BP( K1 ), 1, AP( K1 ), 1 );
               zhpr2(UPLO, K-1, CONE, AP( K1 ), 1, BP( K1 ), 1, AP );
               zaxpy(K-1, CT, BP( K1 ), 1, AP( K1 ), 1 );
               zdscal(K-1, BKK, AP( K1 ), 1 );
               AP[KK] = AKK*BKK**2;
            } // 30
         } else {

            // Compute L**H *A*L

            // JJ and J1J1 are the indices of A(j,j) and A(j+1,j+1)

            JJ = 1;
            for (J = 1; J <= N; J++) { // 40
               J1J1 = JJ + N - J + 1;

               // Compute the j-th column of the lower triangle of A

               AJJ = (AP( JJ )).toDouble();
               BJJ = (BP( JJ )).toDouble();
               AP[JJ] = AJJ*BJJ + ZDOTC( N-J, AP( JJ+1 ), 1, BP( JJ+1 ), 1 );
               zdscal(N-J, BJJ, AP( JJ+1 ), 1 );
               zhpmv(UPLO, N-J, CONE, AP( J1J1 ), BP( JJ+1 ), 1, CONE, AP( JJ+1 ), 1 );
               ztpmv(UPLO, 'Conjugate transpose', 'Non-unit', N-J+1, BP( JJ ), AP( JJ ), 1 );
               JJ = J1J1;
            } // 40
         }
      }
      return;
      }
