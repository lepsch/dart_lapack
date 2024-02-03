      SUBROUTINE DSPGST( ITYPE, UPLO, N, AP, BP, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, ITYPE, N;
      // ..
      // .. Array Arguments ..
      double             AP( * ), BP( * );
      // ..
*
*  =====================================================================
*
      // .. Parameters ..
      double             ONE, HALF;
      PARAMETER          ( ONE = 1.0D0, HALF = 0.5D0 )
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                J, J1, J1J1, JJ, K, K1, K1K1, KK;
      double             AJJ, AKK, BJJ, BKK, CT;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DAXPY, DSCAL, DSPMV, DSPR2, DTPMV, DTPSV, XERBLA
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DDOT;
      // EXTERNAL LSAME, DDOT
      // ..
      // .. Executable Statements ..
*
      // Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( ITYPE.LT.1 .OR. ITYPE.GT.3 ) THEN
         INFO = -1
      ELSE IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSPGST', -INFO )
         RETURN
      END IF
*
      IF( ITYPE.EQ.1 ) THEN
         IF( UPPER ) THEN
*
            // Compute inv(U**T)*A*inv(U)
*
            // J1 and JJ are the indices of A(1,j) and A(j,j)
*
            JJ = 0
            DO 10 J = 1, N
               J1 = JJ + 1
               JJ = JJ + J
*
               // Compute the j-th column of the upper triangle of A
*
               BJJ = BP( JJ )
               CALL DTPSV( UPLO, 'Transpose', 'Nonunit', J, BP, AP( J1 ), 1 )                CALL DSPMV( UPLO, J-1, -ONE, AP, BP( J1 ), 1, ONE, AP( J1 ), 1 )
               CALL DSCAL( J-1, ONE / BJJ, AP( J1 ), 1 )
               AP( JJ ) = ( AP( JJ )-DDOT( J-1, AP( J1 ), 1, BP( J1 ), 1 ) ) / BJJ
   10       CONTINUE
         ELSE
*
            // Compute inv(L)*A*inv(L**T)
*
            // KK and K1K1 are the indices of A(k,k) and A(k+1,k+1)
*
            KK = 1
            DO 20 K = 1, N
               K1K1 = KK + N - K + 1
*
               // Update the lower triangle of A(k:n,k:n)
*
               AKK = AP( KK )
               BKK = BP( KK )
               AKK = AKK / BKK**2
               AP( KK ) = AKK
               IF( K.LT.N ) THEN
                  CALL DSCAL( N-K, ONE / BKK, AP( KK+1 ), 1 )
                  CT = -HALF*AKK
                  CALL DAXPY( N-K, CT, BP( KK+1 ), 1, AP( KK+1 ), 1 )
                  CALL DSPR2( UPLO, N-K, -ONE, AP( KK+1 ), 1, BP( KK+1 ), 1, AP( K1K1 ) )
                  CALL DAXPY( N-K, CT, BP( KK+1 ), 1, AP( KK+1 ), 1 )
                  CALL DTPSV( UPLO, 'No transpose', 'Non-unit', N-K, BP( K1K1 ), AP( KK+1 ), 1 )
               END IF
               KK = K1K1
   20       CONTINUE
         END IF
      ELSE
         IF( UPPER ) THEN
*
            // Compute U*A*U**T
*
            // K1 and KK are the indices of A(1,k) and A(k,k)
*
            KK = 0
            DO 30 K = 1, N
               K1 = KK + 1
               KK = KK + K
*
               // Update the upper triangle of A(1:k,1:k)
*
               AKK = AP( KK )
               BKK = BP( KK )
               CALL DTPMV( UPLO, 'No transpose', 'Non-unit', K-1, BP, AP( K1 ), 1 )
               CT = HALF*AKK
               CALL DAXPY( K-1, CT, BP( K1 ), 1, AP( K1 ), 1 )
               CALL DSPR2( UPLO, K-1, ONE, AP( K1 ), 1, BP( K1 ), 1, AP )
               CALL DAXPY( K-1, CT, BP( K1 ), 1, AP( K1 ), 1 )
               CALL DSCAL( K-1, BKK, AP( K1 ), 1 )
               AP( KK ) = AKK*BKK**2
   30       CONTINUE
         ELSE
*
            // Compute L**T *A*L
*
            // JJ and J1J1 are the indices of A(j,j) and A(j+1,j+1)
*
            JJ = 1
            DO 40 J = 1, N
               J1J1 = JJ + N - J + 1
*
               // Compute the j-th column of the lower triangle of A
*
               AJJ = AP( JJ )
               BJJ = BP( JJ )
               AP( JJ ) = AJJ*BJJ + DDOT( N-J, AP( JJ+1 ), 1, BP( JJ+1 ), 1 )
               CALL DSCAL( N-J, BJJ, AP( JJ+1 ), 1 )
               CALL DSPMV( UPLO, N-J, ONE, AP( J1J1 ), BP( JJ+1 ), 1, ONE, AP( JJ+1 ), 1 )                CALL DTPMV( UPLO, 'Transpose', 'Non-unit', N-J+1, BP( JJ ), AP( JJ ), 1 )
               JJ = J1J1
   40       CONTINUE
         END IF
      END IF
      RETURN
*
      // End of DSPGST
*
      END
