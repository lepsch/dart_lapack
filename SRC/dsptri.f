      SUBROUTINE DSPTRI( UPLO, N, AP, IPIV, WORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, N;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      double             AP( * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0D+0, ZERO = 0.0D+0 ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                J, K, KC, KCNEXT, KP, KPC, KSTEP, KX, NPP;
      double             AK, AKKP1, AKP1, D, T, TEMP;
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DDOT;
      // EXTERNAL LSAME, DDOT
      // ..
      // .. External Subroutines ..
      // EXTERNAL DCOPY, DSPMV, DSWAP, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      if ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      }
      if ( INFO.NE.0 ) {
         xerbla('DSPTRI', -INFO );
         RETURN
      }

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      // Check that the diagonal matrix D is nonsingular.

      if ( UPPER ) {

         // Upper triangular storage: examine D from bottom to top

         KP = N*( N+1 ) / 2
         DO 10 INFO = N, 1, -1
            IF( IPIV( INFO ).GT.0 .AND. AP( KP ).EQ.ZERO ) RETURN
            KP = KP - INFO
   10    CONTINUE
      } else {

         // Lower triangular storage: examine D from top to bottom.

         KP = 1
         DO 20 INFO = 1, N
            IF( IPIV( INFO ).GT.0 .AND. AP( KP ).EQ.ZERO ) RETURN
            KP = KP + N - INFO + 1
   20    CONTINUE
      }
      INFO = 0

      if ( UPPER ) {

         // Compute inv(A) from the factorization A = U*D*U**T.

         // K is the main loop index, increasing from 1 to N in steps of
         // 1 or 2, depending on the size of the diagonal blocks.

         K = 1
         KC = 1
   30    CONTINUE

         // If K > N, exit from loop.

         IF( K.GT.N ) GO TO 50

         KCNEXT = KC + K
         if ( IPIV( K ).GT.0 ) {

            // 1 x 1 diagonal block

            // Invert the diagonal block.

            AP( KC+K-1 ) = ONE / AP( KC+K-1 )

            // Compute column K of the inverse.

            if ( K.GT.1 ) {
               dcopy(K-1, AP( KC ), 1, WORK, 1 );
               dspmv(UPLO, K-1, -ONE, AP, WORK, 1, ZERO, AP( KC ), 1 )                AP( KC+K-1 ) = AP( KC+K-1 ) - DDOT( K-1, WORK, 1, AP( KC ), 1 );
            }
            KSTEP = 1
         } else {

            // 2 x 2 diagonal block

            // Invert the diagonal block.

            T = ABS( AP( KCNEXT+K-1 ) )
            AK = AP( KC+K-1 ) / T
            AKP1 = AP( KCNEXT+K ) / T
            AKKP1 = AP( KCNEXT+K-1 ) / T
            D = T*( AK*AKP1-ONE )
            AP( KC+K-1 ) = AKP1 / D
            AP( KCNEXT+K ) = AK / D
            AP( KCNEXT+K-1 ) = -AKKP1 / D

            // Compute columns K and K+1 of the inverse.

            if ( K.GT.1 ) {
               dcopy(K-1, AP( KC ), 1, WORK, 1 );
               dspmv(UPLO, K-1, -ONE, AP, WORK, 1, ZERO, AP( KC ), 1 )                AP( KC+K-1 ) = AP( KC+K-1 ) - DDOT( K-1, WORK, 1, AP( KC ), 1 )                AP( KCNEXT+K-1 ) = AP( KCNEXT+K-1 ) - DDOT( K-1, AP( KC ), 1, AP( KCNEXT ), 1 );
               dcopy(K-1, AP( KCNEXT ), 1, WORK, 1 );
               dspmv(UPLO, K-1, -ONE, AP, WORK, 1, ZERO, AP( KCNEXT ), 1 )                AP( KCNEXT+K ) = AP( KCNEXT+K ) - DDOT( K-1, WORK, 1, AP( KCNEXT ), 1 );
            }
            KSTEP = 2
            KCNEXT = KCNEXT + K + 1
         }

         KP = ABS( IPIV( K ) )
         if ( KP.NE.K ) {

            // Interchange rows and columns K and KP in the leading
            // submatrix A(1:k+1,1:k+1)

            KPC = ( KP-1 )*KP / 2 + 1
            dswap(KP-1, AP( KC ), 1, AP( KPC ), 1 );
            KX = KPC + KP - 1
            DO 40 J = KP + 1, K - 1
               KX = KX + J - 1
               TEMP = AP( KC+J-1 )
               AP( KC+J-1 ) = AP( KX )
               AP( KX ) = TEMP
   40       CONTINUE
            TEMP = AP( KC+K-1 )
            AP( KC+K-1 ) = AP( KPC+KP-1 )
            AP( KPC+KP-1 ) = TEMP
            if ( KSTEP.EQ.2 ) {
               TEMP = AP( KC+K+K-1 )
               AP( KC+K+K-1 ) = AP( KC+K+KP-1 )
               AP( KC+K+KP-1 ) = TEMP
            }
         }

         K = K + KSTEP
         KC = KCNEXT
         GO TO 30
   50    CONTINUE

      } else {

         // Compute inv(A) from the factorization A = L*D*L**T.

         // K is the main loop index, increasing from 1 to N in steps of
         // 1 or 2, depending on the size of the diagonal blocks.

         NPP = N*( N+1 ) / 2
         K = N
         KC = NPP
   60    CONTINUE

         // If K < 1, exit from loop.

         IF( K.LT.1 ) GO TO 80

         KCNEXT = KC - ( N-K+2 )
         if ( IPIV( K ).GT.0 ) {

            // 1 x 1 diagonal block

            // Invert the diagonal block.

            AP( KC ) = ONE / AP( KC )

            // Compute column K of the inverse.

            if ( K.LT.N ) {
               dcopy(N-K, AP( KC+1 ), 1, WORK, 1 );
               dspmv(UPLO, N-K, -ONE, AP( KC+N-K+1 ), WORK, 1, ZERO, AP( KC+1 ), 1 );
               AP( KC ) = AP( KC ) - DDOT( N-K, WORK, 1, AP( KC+1 ), 1 )
            }
            KSTEP = 1
         } else {

            // 2 x 2 diagonal block

            // Invert the diagonal block.

            T = ABS( AP( KCNEXT+1 ) )
            AK = AP( KCNEXT ) / T
            AKP1 = AP( KC ) / T
            AKKP1 = AP( KCNEXT+1 ) / T
            D = T*( AK*AKP1-ONE )
            AP( KCNEXT ) = AKP1 / D
            AP( KC ) = AK / D
            AP( KCNEXT+1 ) = -AKKP1 / D

            // Compute columns K-1 and K of the inverse.

            if ( K.LT.N ) {
               dcopy(N-K, AP( KC+1 ), 1, WORK, 1 );
               dspmv(UPLO, N-K, -ONE, AP( KC+( N-K+1 ) ), WORK, 1, ZERO, AP( KC+1 ), 1 );
               AP( KC ) = AP( KC ) - DDOT( N-K, WORK, 1, AP( KC+1 ), 1 )
               AP( KCNEXT+1 ) = AP( KCNEXT+1 ) - DDOT( N-K, AP( KC+1 ), 1, AP( KCNEXT+2 ), 1 )
               dcopy(N-K, AP( KCNEXT+2 ), 1, WORK, 1 );
               dspmv(UPLO, N-K, -ONE, AP( KC+( N-K+1 ) ), WORK, 1, ZERO, AP( KCNEXT+2 ), 1 )                AP( KCNEXT ) = AP( KCNEXT ) - DDOT( N-K, WORK, 1, AP( KCNEXT+2 ), 1 );
            }
            KSTEP = 2
            KCNEXT = KCNEXT - ( N-K+3 )
         }

         KP = ABS( IPIV( K ) )
         if ( KP.NE.K ) {

            // Interchange rows and columns K and KP in the trailing
            // submatrix A(k-1:n,k-1:n)

            KPC = NPP - ( N-KP+1 )*( N-KP+2 ) / 2 + 1
            IF( KP.LT.N ) CALL DSWAP( N-KP, AP( KC+KP-K+1 ), 1, AP( KPC+1 ), 1 )
            KX = KC + KP - K
            DO 70 J = K + 1, KP - 1
               KX = KX + N - J + 1
               TEMP = AP( KC+J-K )
               AP( KC+J-K ) = AP( KX )
               AP( KX ) = TEMP
   70       CONTINUE
            TEMP = AP( KC )
            AP( KC ) = AP( KPC )
            AP( KPC ) = TEMP
            if ( KSTEP.EQ.2 ) {
               TEMP = AP( KC-N+K-1 )
               AP( KC-N+K-1 ) = AP( KC-N+KP-1 )
               AP( KC-N+KP-1 ) = TEMP
            }
         }

         K = K - KSTEP
         KC = KCNEXT
         GO TO 60
   80    CONTINUE
      }

      RETURN

      // End of DSPTRI

      }
