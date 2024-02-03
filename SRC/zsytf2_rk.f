      SUBROUTINE ZSYTF2_RK( UPLO, N, A, LDA, E, IPIV, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, N;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      COMPLEX*16         A( LDA, * ), E( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      double             EIGHT, SEVTEN;
      const              EIGHT = 8.0D+0, SEVTEN = 17.0D+0 ;
      COMPLEX*16         CONE, CZERO
      const              CONE = ( 1.0D+0, 0.0D+0 ), CZERO = ( 0.0D+0, 0.0D+0 ) ;
      // ..
      // .. Local Scalars ..
      bool               UPPER, DONE;
      int                I, IMAX, J, JMAX, ITEMP, K, KK, KP, KSTEP, P, II;
      double             ABSAKK, ALPHA, COLMAX, ROWMAX, DTEMP, SFMIN;
      COMPLEX*16         D11, D12, D21, D22, T, WK, WKM1, WKP1, Z
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                IZAMAX;
      double             DLAMCH;
      // EXTERNAL LSAME, IZAMAX, DLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZSCAL, ZSWAP, ZSYR, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, SQRT, DIMAG, DBLE
      // ..
      // .. Statement Functions ..
      double             CABS1;
      // ..
      // .. Statement Function definitions ..
      CABS1( Z ) = ABS( DBLE( Z ) ) + ABS( DIMAG( Z ) )
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZSYTF2_RK', -INFO )
         RETURN
      END IF

      // Initialize ALPHA for use in choosing pivot block size.

      ALPHA = ( ONE+SQRT( SEVTEN ) ) / EIGHT

      // Compute machine safe minimum

      SFMIN = DLAMCH( 'S' )

      IF( UPPER ) THEN

         // Factorize A as U*D*U**T using the upper triangle of A

         // Initialize the first entry of array E, where superdiagonal
         // elements of D are stored

         E( 1 ) = CZERO

         // K is the main loop index, decreasing from N to 1 in steps of
         // 1 or 2

         K = N
   10    CONTINUE

         // If K < 1, exit from loop

         IF( K.LT.1 ) GO TO 34
         KSTEP = 1
         P = K

         // Determine rows and columns to be interchanged and whether
         // a 1-by-1 or 2-by-2 pivot block will be used

         ABSAKK = CABS1( A( K, K ) )

         // IMAX is the row-index of the largest off-diagonal element in
         // column K, and COLMAX is its absolute value.
         // Determine both COLMAX and IMAX.

         IF( K.GT.1 ) THEN
            IMAX = IZAMAX( K-1, A( 1, K ), 1 )
            COLMAX = CABS1( A( IMAX, K ) )
         } else {
            COLMAX = ZERO
         END IF

         IF( (MAX( ABSAKK, COLMAX ).EQ.ZERO) ) THEN

            // Column K is zero or underflow: set INFO and continue

            IF( INFO.EQ.0 ) INFO = K
            KP = K

            // Set E( K ) to zero

            IF( K.GT.1 ) E( K ) = CZERO

         } else {

            // Test for interchange

            // Equivalent to testing for (used to handle NaN and Inf)
            // ABSAKK.GE.ALPHA*COLMAX

            IF( .NOT.( ABSAKK.LT.ALPHA*COLMAX ) ) THEN

               // no interchange,
               // use 1-by-1 pivot block

               KP = K
            } else {

               DONE = .FALSE.

               // Loop until pivot found

   12          CONTINUE

                  // Begin pivot search loop body

                  // JMAX is the column-index of the largest off-diagonal
                  // element in row IMAX, and ROWMAX is its absolute value.
                  // Determine both ROWMAX and JMAX.

                  IF( IMAX.NE.K ) THEN
                     JMAX = IMAX + IZAMAX( K-IMAX, A( IMAX, IMAX+1 ), LDA )
                     ROWMAX = CABS1( A( IMAX, JMAX ) )
                  } else {
                     ROWMAX = ZERO
                  END IF

                  IF( IMAX.GT.1 ) THEN
                     ITEMP = IZAMAX( IMAX-1, A( 1, IMAX ), 1 )
                     DTEMP = CABS1( A( ITEMP, IMAX ) )
                     IF( DTEMP.GT.ROWMAX ) THEN
                        ROWMAX = DTEMP
                        JMAX = ITEMP
                     END IF
                  END IF

                  // Equivalent to testing for (used to handle NaN and Inf)
                  // ABS( A( IMAX, IMAX ) ).GE.ALPHA*ROWMAX

                  IF( .NOT.( CABS1( A( IMAX, IMAX ) ).LT.ALPHA*ROWMAX )) THEN

                     // interchange rows and columns K and IMAX,
                     // use 1-by-1 pivot block

                     KP = IMAX
                     DONE = .TRUE.

                  // Equivalent to testing for ROWMAX .EQ. COLMAX,
                  // used to handle NaN and Inf

                  ELSE IF( ( P.EQ.JMAX ).OR.( ROWMAX.LE.COLMAX ) ) THEN

                     // interchange rows and columns K+1 and IMAX,
                     // use 2-by-2 pivot block

                     KP = IMAX
                     KSTEP = 2
                     DONE = .TRUE.
                  } else {

                     // Pivot NOT found, set variables and repeat

                     P = IMAX
                     COLMAX = ROWMAX
                     IMAX = JMAX
                  END IF

                  // End pivot search loop body

               IF( .NOT. DONE ) GOTO 12

            END IF

            // Swap TWO rows and TWO columns

            // First swap

            IF( ( KSTEP.EQ.2 ) .AND. ( P.NE.K ) ) THEN

               // Interchange rows and column K and P in the leading
               // submatrix A(1:k,1:k) if we have a 2-by-2 pivot

               IF( P.GT.1 ) CALL ZSWAP( P-1, A( 1, K ), 1, A( 1, P ), 1 )                IF( P.LT.(K-1) ) CALL ZSWAP( K-P-1, A( P+1, K ), 1, A( P, P+1 ), LDA )
               T = A( K, K )
               A( K, K ) = A( P, P )
               A( P, P ) = T

               // Convert upper triangle of A into U form by applying
              t // he interchanges in columns k+1:N.

               IF( K.LT.N ) CALL ZSWAP( N-K, A( K, K+1 ), LDA, A( P, K+1 ), LDA )

            END IF

            // Second swap

            KK = K - KSTEP + 1
            IF( KP.NE.KK ) THEN

               // Interchange rows and columns KK and KP in the leading
               // submatrix A(1:k,1:k)

               IF( KP.GT.1 ) CALL ZSWAP( KP-1, A( 1, KK ), 1, A( 1, KP ), 1 )                IF( ( KK.GT.1 ) .AND. ( KP.LT.(KK-1) ) ) CALL ZSWAP( KK-KP-1, A( KP+1, KK ), 1, A( KP, KP+1 ), LDA )
               T = A( KK, KK )
               A( KK, KK ) = A( KP, KP )
               A( KP, KP ) = T
               IF( KSTEP.EQ.2 ) THEN
                  T = A( K-1, K )
                  A( K-1, K ) = A( KP, K )
                  A( KP, K ) = T
               END IF

               // Convert upper triangle of A into U form by applying
              t // he interchanges in columns k+1:N.

               IF( K.LT.N ) CALL ZSWAP( N-K, A( KK, K+1 ), LDA, A( KP, K+1 ), LDA )

            END IF

            // Update the leading submatrix

            IF( KSTEP.EQ.1 ) THEN

               // 1-by-1 pivot block D(k): column k now holds

               // W(k) = U(k)*D(k)

               // where U(k) is the k-th column of U

               IF( K.GT.1 ) THEN

                  // Perform a rank-1 update of A(1:k-1,1:k-1) and
                  // store U(k) in column k

                  IF( CABS1( A( K, K ) ).GE.SFMIN ) THEN

                     // Perform a rank-1 update of A(1:k-1,1:k-1) as
                     // A := A - U(k)*D(k)*U(k)**T
                        // = A - W(k)*1/D(k)*W(k)**T

                     D11 = CONE / A( K, K )
                     CALL ZSYR( UPLO, K-1, -D11, A( 1, K ), 1, A, LDA )

                     // Store U(k) in column k

                     CALL ZSCAL( K-1, D11, A( 1, K ), 1 )
                  } else {

                     // Store L(k) in column K

                     D11 = A( K, K )
                     DO 16 II = 1, K - 1
                        A( II, K ) = A( II, K ) / D11
   16                CONTINUE

                     // Perform a rank-1 update of A(k+1:n,k+1:n) as
                     // A := A - U(k)*D(k)*U(k)**T
                        // = A - W(k)*(1/D(k))*W(k)**T
                        // = A - (W(k)/D(k))*(D(k))*(W(k)/D(K))**T

                     CALL ZSYR( UPLO, K-1, -D11, A( 1, K ), 1, A, LDA )
                  END IF

                  // Store the superdiagonal element of D in array E

                  E( K ) = CZERO

               END IF

            } else {

               // 2-by-2 pivot block D(k): columns k and k-1 now hold

               // ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k)

               // where U(k) and U(k-1) are the k-th and (k-1)-th columns
               // of U

               // Perform a rank-2 update of A(1:k-2,1:k-2) as

               // A := A - ( U(k-1) U(k) )*D(k)*( U(k-1) U(k) )**T
                  // = A - ( ( A(k-1)A(k) )*inv(D(k)) ) * ( A(k-1)A(k) )**T

               // and store L(k) and L(k+1) in columns k and k+1

               IF( K.GT.2 ) THEN

                  D12 = A( K-1, K )
                  D22 = A( K-1, K-1 ) / D12
                  D11 = A( K, K ) / D12
                  T = CONE / ( D11*D22-CONE )

                  DO 30 J = K - 2, 1, -1

                     WKM1 = T*( D11*A( J, K-1 )-A( J, K ) )
                     WK = T*( D22*A( J, K )-A( J, K-1 ) )

                     DO 20 I = J, 1, -1
                        A( I, J ) = A( I, J ) - (A( I, K ) / D12 )*WK - ( A( I, K-1 ) / D12 )*WKM1
   20                CONTINUE

                     // Store U(k) and U(k-1) in cols k and k-1 for row J

                     A( J, K ) = WK / D12
                     A( J, K-1 ) = WKM1 / D12

   30             CONTINUE

               END IF

               // Copy superdiagonal elements of D(K) to E(K) and
               // ZERO out superdiagonal entry of A

               E( K ) = A( K-1, K )
               E( K-1 ) = CZERO
               A( K-1, K ) = CZERO

            END IF

            // End column K is nonsingular

         END IF

         // Store details of the interchanges in IPIV

         IF( KSTEP.EQ.1 ) THEN
            IPIV( K ) = KP
         } else {
            IPIV( K ) = -P
            IPIV( K-1 ) = -KP
         END IF

         // Decrease K and return to the start of the main loop

         K = K - KSTEP
         GO TO 10

   34    CONTINUE

      } else {

         // Factorize A as L*D*L**T using the lower triangle of A

         // Initialize the unused last entry of the subdiagonal array E.

         E( N ) = CZERO

         // K is the main loop index, increasing from 1 to N in steps of
         // 1 or 2

         K = 1
   40    CONTINUE

         // If K > N, exit from loop

         IF( K.GT.N ) GO TO 64
         KSTEP = 1
         P = K

         // Determine rows and columns to be interchanged and whether
         // a 1-by-1 or 2-by-2 pivot block will be used

         ABSAKK = CABS1( A( K, K ) )

         // IMAX is the row-index of the largest off-diagonal element in
         // column K, and COLMAX is its absolute value.
         // Determine both COLMAX and IMAX.

         IF( K.LT.N ) THEN
            IMAX = K + IZAMAX( N-K, A( K+1, K ), 1 )
            COLMAX = CABS1( A( IMAX, K ) )
         } else {
            COLMAX = ZERO
         END IF

         IF( ( MAX( ABSAKK, COLMAX ).EQ.ZERO ) ) THEN

            // Column K is zero or underflow: set INFO and continue

            IF( INFO.EQ.0 ) INFO = K
            KP = K

            // Set E( K ) to zero

            IF( K.LT.N ) E( K ) = CZERO

         } else {

            // Test for interchange

            // Equivalent to testing for (used to handle NaN and Inf)
            // ABSAKK.GE.ALPHA*COLMAX

            IF( .NOT.( ABSAKK.LT.ALPHA*COLMAX ) ) THEN

               // no interchange, use 1-by-1 pivot block

               KP = K

            } else {

               DONE = .FALSE.

               // Loop until pivot found

   42          CONTINUE

                  // Begin pivot search loop body

                  // JMAX is the column-index of the largest off-diagonal
                  // element in row IMAX, and ROWMAX is its absolute value.
                  // Determine both ROWMAX and JMAX.

                  IF( IMAX.NE.K ) THEN
                     JMAX = K - 1 + IZAMAX( IMAX-K, A( IMAX, K ), LDA )
                     ROWMAX = CABS1( A( IMAX, JMAX ) )
                  } else {
                     ROWMAX = ZERO
                  END IF

                  IF( IMAX.LT.N ) THEN
                     ITEMP = IMAX + IZAMAX( N-IMAX, A( IMAX+1, IMAX ), 1 )
                     DTEMP = CABS1( A( ITEMP, IMAX ) )
                     IF( DTEMP.GT.ROWMAX ) THEN
                        ROWMAX = DTEMP
                        JMAX = ITEMP
                     END IF
                  END IF

                  // Equivalent to testing for (used to handle NaN and Inf)
                  // ABS( A( IMAX, IMAX ) ).GE.ALPHA*ROWMAX

                  IF( .NOT.( CABS1( A( IMAX, IMAX ) ).LT.ALPHA*ROWMAX )) THEN

                     // interchange rows and columns K and IMAX,
                     // use 1-by-1 pivot block

                     KP = IMAX
                     DONE = .TRUE.

                  // Equivalent to testing for ROWMAX .EQ. COLMAX,
                  // used to handle NaN and Inf

                  ELSE IF( ( P.EQ.JMAX ).OR.( ROWMAX.LE.COLMAX ) ) THEN

                     // interchange rows and columns K+1 and IMAX,
                     // use 2-by-2 pivot block

                     KP = IMAX
                     KSTEP = 2
                     DONE = .TRUE.
                  } else {

                     // Pivot NOT found, set variables and repeat

                     P = IMAX
                     COLMAX = ROWMAX
                     IMAX = JMAX
                  END IF

                  // End pivot search loop body

               IF( .NOT. DONE ) GOTO 42

            END IF

            // Swap TWO rows and TWO columns

            // First swap

            IF( ( KSTEP.EQ.2 ) .AND. ( P.NE.K ) ) THEN

               // Interchange rows and column K and P in the trailing
               // submatrix A(k:n,k:n) if we have a 2-by-2 pivot

               IF( P.LT.N ) CALL ZSWAP( N-P, A( P+1, K ), 1, A( P+1, P ), 1 )                IF( P.GT.(K+1) ) CALL ZSWAP( P-K-1, A( K+1, K ), 1, A( P, K+1 ), LDA )
               T = A( K, K )
               A( K, K ) = A( P, P )
               A( P, P ) = T

               // Convert lower triangle of A into L form by applying
              t // he interchanges in columns 1:k-1.

               IF ( K.GT.1 ) CALL ZSWAP( K-1, A( K, 1 ), LDA, A( P, 1 ), LDA )

            END IF

            // Second swap

            KK = K + KSTEP - 1
            IF( KP.NE.KK ) THEN

               // Interchange rows and columns KK and KP in the trailing
               // submatrix A(k:n,k:n)

               IF( KP.LT.N ) CALL ZSWAP( N-KP, A( KP+1, KK ), 1, A( KP+1, KP ), 1 )                IF( ( KK.LT.N ) .AND. ( KP.GT.(KK+1) ) ) CALL ZSWAP( KP-KK-1, A( KK+1, KK ), 1, A( KP, KK+1 ), LDA )
               T = A( KK, KK )
               A( KK, KK ) = A( KP, KP )
               A( KP, KP ) = T
               IF( KSTEP.EQ.2 ) THEN
                  T = A( K+1, K )
                  A( K+1, K ) = A( KP, K )
                  A( KP, K ) = T
               END IF

               // Convert lower triangle of A into L form by applying
              t // he interchanges in columns 1:k-1.

               IF ( K.GT.1 ) CALL ZSWAP( K-1, A( KK, 1 ), LDA, A( KP, 1 ), LDA )

            END IF

            // Update the trailing submatrix

            IF( KSTEP.EQ.1 ) THEN

               // 1-by-1 pivot block D(k): column k now holds

               // W(k) = L(k)*D(k)

               // where L(k) is the k-th column of L

               IF( K.LT.N ) THEN

               // Perform a rank-1 update of A(k+1:n,k+1:n) and
               // store L(k) in column k

                  IF( CABS1( A( K, K ) ).GE.SFMIN ) THEN

                     // Perform a rank-1 update of A(k+1:n,k+1:n) as
                     // A := A - L(k)*D(k)*L(k)**T
                        // = A - W(k)*(1/D(k))*W(k)**T

                     D11 = CONE / A( K, K )
                     CALL ZSYR( UPLO, N-K, -D11, A( K+1, K ), 1, A( K+1, K+1 ), LDA )

                     // Store L(k) in column k

                     CALL ZSCAL( N-K, D11, A( K+1, K ), 1 )
                  } else {

                     // Store L(k) in column k

                     D11 = A( K, K )
                     DO 46 II = K + 1, N
                        A( II, K ) = A( II, K ) / D11
   46                CONTINUE

                     // Perform a rank-1 update of A(k+1:n,k+1:n) as
                     // A := A - L(k)*D(k)*L(k)**T
                        // = A - W(k)*(1/D(k))*W(k)**T
                        // = A - (W(k)/D(k))*(D(k))*(W(k)/D(K))**T

                     CALL ZSYR( UPLO, N-K, -D11, A( K+1, K ), 1, A( K+1, K+1 ), LDA )
                  END IF

                  // Store the subdiagonal element of D in array E

                  E( K ) = CZERO

               END IF

            } else {

               // 2-by-2 pivot block D(k): columns k and k+1 now hold

               // ( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k)

               // where L(k) and L(k+1) are the k-th and (k+1)-th columns
               // of L


               // Perform a rank-2 update of A(k+2:n,k+2:n) as

               // A := A - ( L(k) L(k+1) ) * D(k) * ( L(k) L(k+1) )**T
                  // = A - ( ( A(k)A(k+1) )*inv(D(k) ) * ( A(k)A(k+1) )**T

               // and store L(k) and L(k+1) in columns k and k+1

               IF( K.LT.N-1 ) THEN

                  D21 = A( K+1, K )
                  D11 = A( K+1, K+1 ) / D21
                  D22 = A( K, K ) / D21
                  T = CONE / ( D11*D22-CONE )

                  DO 60 J = K + 2, N

                     // Compute  D21 * ( W(k)W(k+1) ) * inv(D(k)) for row J

                     WK = T*( D11*A( J, K )-A( J, K+1 ) )
                     WKP1 = T*( D22*A( J, K+1 )-A( J, K ) )

                     // Perform a rank-2 update of A(k+2:n,k+2:n)

                     DO 50 I = J, N
                        A( I, J ) = A( I, J ) - ( A( I, K ) / D21 )*WK - ( A( I, K+1 ) / D21 )*WKP1
   50                CONTINUE

                     // Store L(k) and L(k+1) in cols k and k+1 for row J

                     A( J, K ) = WK / D21
                     A( J, K+1 ) = WKP1 / D21

   60             CONTINUE

               END IF

               // Copy subdiagonal elements of D(K) to E(K) and
               // ZERO out subdiagonal entry of A

               E( K ) = A( K+1, K )
               E( K+1 ) = CZERO
               A( K+1, K ) = CZERO

            END IF

            // End column K is nonsingular

         END IF

         // Store details of the interchanges in IPIV

         IF( KSTEP.EQ.1 ) THEN
            IPIV( K ) = KP
         } else {
            IPIV( K ) = -P
            IPIV( K+1 ) = -KP
         END IF

         // Increase K and return to the start of the main loop

         K = K + KSTEP
         GO TO 40

   64    CONTINUE

      END IF

      RETURN

      // End of ZSYTF2_RK

      }
