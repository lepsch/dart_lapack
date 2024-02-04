      void chpgst(ITYPE, UPLO, N, AP, BP, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, ITYPE, N;
      // ..
      // .. Array Arguments ..
      COMPLEX            AP( * ), BP( * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ONE, HALF;
      const              ONE = 1.0, HALF = 0.5 ;
      COMPLEX            CONE;
      const              CONE = ( 1.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                J, J1, J1J1, JJ, K, K1, K1K1, KK;
      REAL               AJJ, AKK, BJJ, BKK;
      COMPLEX            CT;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CAXPY, CHPMV, CHPR2, CSSCAL, CTPMV, CTPSV, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC REAL
      // ..
      // .. External Functions ..
      //- bool               LSAME;
      //- COMPLEX            CDOTC;
      // EXTERNAL LSAME, CDOTC
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
         xerbla('CHPGST', -INFO );
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

               AP[JJ] = REAL( AP( JJ ) );
               BJJ = REAL( BP( JJ ) );
               ctpsv(UPLO, 'Conjugate transpose', 'Non-unit', J, BP, AP( J1 ), 1 );
               chpmv(UPLO, J-1, -CONE, AP, BP( J1 ), 1, CONE, AP( J1 ), 1 );
               csscal(J-1, ONE / BJJ, AP( J1 ), 1 );
               AP[JJ] = ( AP( JJ )-CDOTC( J-1, AP( J1 ), 1, BP( J1 ), 1 ) ) / BJJ;
            } // 10
         } else {

            // Compute inv(L)*A*inv(L**H)

            // KK and K1K1 are the indices of A(k,k) and A(k+1,k+1)

            KK = 1;
            for (K = 1; K <= N; K++) { // 20
               K1K1 = KK + N - K + 1;

               // Update the lower triangle of A(k:n,k:n)

               AKK = REAL( AP( KK ) );
               BKK = REAL( BP( KK ) );
               AKK = AKK / BKK**2;
               AP[KK] = AKK;
               if ( K < N ) {
                  csscal(N-K, ONE / BKK, AP( KK+1 ), 1 );
                  CT = -HALF*AKK;
                  caxpy(N-K, CT, BP( KK+1 ), 1, AP( KK+1 ), 1 );
                  chpr2(UPLO, N-K, -CONE, AP( KK+1 ), 1, BP( KK+1 ), 1, AP( K1K1 ) );
                  caxpy(N-K, CT, BP( KK+1 ), 1, AP( KK+1 ), 1 );
                  ctpsv(UPLO, 'No transpose', 'Non-unit', N-K, BP( K1K1 ), AP( KK+1 ), 1 );
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

               AKK = REAL( AP( KK ) );
               BKK = REAL( BP( KK ) );
               ctpmv(UPLO, 'No transpose', 'Non-unit', K-1, BP, AP( K1 ), 1 );
               CT = HALF*AKK;
               caxpy(K-1, CT, BP( K1 ), 1, AP( K1 ), 1 );
               chpr2(UPLO, K-1, CONE, AP( K1 ), 1, BP( K1 ), 1, AP );
               caxpy(K-1, CT, BP( K1 ), 1, AP( K1 ), 1 );
               csscal(K-1, BKK, AP( K1 ), 1 );
               AP[KK] = AKK*BKK**2;
            } // 30
         } else {

            // Compute L**H *A*L

            // JJ and J1J1 are the indices of A(j,j) and A(j+1,j+1)

            JJ = 1;
            for (J = 1; J <= N; J++) { // 40
               J1J1 = JJ + N - J + 1;

               // Compute the j-th column of the lower triangle of A

               AJJ = REAL( AP( JJ ) );
               BJJ = REAL( BP( JJ ) );
               AP[JJ] = AJJ*BJJ + CDOTC( N-J, AP( JJ+1 ), 1, BP( JJ+1 ), 1 );
               csscal(N-J, BJJ, AP( JJ+1 ), 1 );
               chpmv(UPLO, N-J, CONE, AP( J1J1 ), BP( JJ+1 ), 1, CONE, AP( JJ+1 ), 1 );
               ctpmv(UPLO, 'Conjugate transpose', 'Non-unit', N-J+1, BP( JJ ), AP( JJ ), 1 );
               JJ = J1J1;
            } // 40
         }
      }
      return;
      }
