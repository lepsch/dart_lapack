      SUBROUTINE SSPGST( ITYPE, UPLO, N, AP, BP, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, ITYPE, N;
      // ..
      // .. Array Arguments ..
      REAL               AP( * ), BP( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE, HALF
      const              ONE = 1.0, HALF = 0.5 ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                J, J1, J1J1, JJ, K, K1, K1K1, KK;
      REAL               AJJ, AKK, BJJ, BKK, CT
      // ..
      // .. External Subroutines ..
      // EXTERNAL SAXPY, SSCAL, SSPMV, SSPR2, STPMV, STPSV, XERBLA
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               SDOT
      // EXTERNAL LSAME, SDOT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      if ( ITYPE.LT.1 .OR. ITYPE.GT.3 ) {
         INFO = -1
      } else if ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -2
      } else if ( N.LT.0 ) {
         INFO = -3
      }
      if ( INFO.NE.0 ) {
         xerbla('SSPGST', -INFO );
         RETURN
      }

      if ( ITYPE.EQ.1 ) {
         if ( UPPER ) {

            // Compute inv(U**T)*A*inv(U)

            // J1 and JJ are the indices of A(1,j) and A(j,j)

            JJ = 0
            for (J = 1; J <= N; J++) { // 10
               J1 = JJ + 1
               JJ = JJ + J

               // Compute the j-th column of the upper triangle of A

               BJJ = BP( JJ )
               stpsv(UPLO, 'Transpose', 'Nonunit', J, BP, AP( J1 ), 1 )                CALL SSPMV( UPLO, J-1, -ONE, AP, BP( J1 ), 1, ONE, AP( J1 ), 1 );
               sscal(J-1, ONE / BJJ, AP( J1 ), 1 );
               AP( JJ ) = ( AP( JJ )-SDOT( J-1, AP( J1 ), 1, BP( J1 ), 1 ) ) / BJJ
            } // 10
         } else {

            // Compute inv(L)*A*inv(L**T)

            // KK and K1K1 are the indices of A(k,k) and A(k+1,k+1)

            KK = 1
            for (K = 1; K <= N; K++) { // 20
               K1K1 = KK + N - K + 1

               // Update the lower triangle of A(k:n,k:n)

               AKK = AP( KK )
               BKK = BP( KK )
               AKK = AKK / BKK**2
               AP( KK ) = AKK
               if ( K.LT.N ) {
                  sscal(N-K, ONE / BKK, AP( KK+1 ), 1 );
                  CT = -HALF*AKK
                  saxpy(N-K, CT, BP( KK+1 ), 1, AP( KK+1 ), 1 );
                  sspr2(UPLO, N-K, -ONE, AP( KK+1 ), 1, BP( KK+1 ), 1, AP( K1K1 ) );
                  saxpy(N-K, CT, BP( KK+1 ), 1, AP( KK+1 ), 1 );
                  stpsv(UPLO, 'No transpose', 'Non-unit', N-K, BP( K1K1 ), AP( KK+1 ), 1 );
               }
               KK = K1K1
            } // 20
         }
      } else {
         if ( UPPER ) {

            // Compute U*A*U**T

            // K1 and KK are the indices of A(1,k) and A(k,k)

            KK = 0
            for (K = 1; K <= N; K++) { // 30
               K1 = KK + 1
               KK = KK + K

               // Update the upper triangle of A(1:k,1:k)

               AKK = AP( KK )
               BKK = BP( KK )
               stpmv(UPLO, 'No transpose', 'Non-unit', K-1, BP, AP( K1 ), 1 );
               CT = HALF*AKK
               saxpy(K-1, CT, BP( K1 ), 1, AP( K1 ), 1 );
               sspr2(UPLO, K-1, ONE, AP( K1 ), 1, BP( K1 ), 1, AP );
               saxpy(K-1, CT, BP( K1 ), 1, AP( K1 ), 1 );
               sscal(K-1, BKK, AP( K1 ), 1 );
               AP( KK ) = AKK*BKK**2
            } // 30
         } else {

            // Compute L**T *A*L

            // JJ and J1J1 are the indices of A(j,j) and A(j+1,j+1)

            JJ = 1
            for (J = 1; J <= N; J++) { // 40
               J1J1 = JJ + N - J + 1

               // Compute the j-th column of the lower triangle of A

               AJJ = AP( JJ )
               BJJ = BP( JJ )
               AP( JJ ) = AJJ*BJJ + SDOT( N-J, AP( JJ+1 ), 1, BP( JJ+1 ), 1 )
               sscal(N-J, BJJ, AP( JJ+1 ), 1 );
               sspmv(UPLO, N-J, ONE, AP( J1J1 ), BP( JJ+1 ), 1, ONE, AP( JJ+1 ), 1 )                CALL STPMV( UPLO, 'Transpose', 'Non-unit', N-J+1, BP( JJ ), AP( JJ ), 1 );
               JJ = J1J1
            } // 40
         }
      }
      RETURN

      // End of SSPGST

      }
