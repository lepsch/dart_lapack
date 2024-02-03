      SUBROUTINE ZHETF2( UPLO, N, A, LDA, IPIV, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, N;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      COMPLEX*16         A( LDA, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      double             EIGHT, SEVTEN;
      const              EIGHT = 8.0D+0, SEVTEN = 17.0D+0 ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                I, IMAX, J, JMAX, K, KK, KP, KSTEP;
      double             ABSAKK, ALPHA, COLMAX, D, D11, D22, R1, ROWMAX, TT;
      COMPLEX*16         D12, D21, T, WK, WKM1, WKP1, ZDUM
      // ..
      // .. External Functions ..
      bool               LSAME, DISNAN;
      int                IZAMAX;
      double             DLAPY2;
      // EXTERNAL LSAME, IZAMAX, DLAPY2, DISNAN
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZDSCAL, ZHER, ZSWAP
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DCMPLX, DCONJG, DIMAG, MAX, SQRT
      // ..
      // .. Statement Functions ..
      double             CABS1;
      // ..
      // .. Statement Function definitions ..
      CABS1( ZDUM ) = ABS( DBLE( ZDUM ) ) + ABS( DIMAG( ZDUM ) )
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      if ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -4
      }
      if ( INFO.NE.0 ) {
         xerbla('ZHETF2', -INFO );
         RETURN
      }

      // Initialize ALPHA for use in choosing pivot block size.

      ALPHA = ( ONE+SQRT( SEVTEN ) ) / EIGHT

      if ( UPPER ) {

         // Factorize A as U*D*U**H using the upper triangle of A

         // K is the main loop index, decreasing from N to 1 in steps of
         // 1 or 2

         K = N
         } // 10

         // If K < 1, exit from loop

         IF( K.LT.1 ) GO TO 90
         KSTEP = 1

         // Determine rows and columns to be interchanged and whether
         // a 1-by-1 or 2-by-2 pivot block will be used

         ABSAKK = ABS( DBLE( A( K, K ) ) )

         // IMAX is the row-index of the largest off-diagonal element in
         // column K, and COLMAX is its absolute value.
         // Determine both COLMAX and IMAX.

         if ( K.GT.1 ) {
            IMAX = IZAMAX( K-1, A( 1, K ), 1 )
            COLMAX = CABS1( A( IMAX, K ) )
         } else {
            COLMAX = ZERO
         }

         if ( (MAX( ABSAKK, COLMAX ).EQ.ZERO) .OR. DISNAN(ABSAKK) ) {

            // Column K is zero or underflow, or contains a NaN:
            // set INFO and continue

            IF( INFO.EQ.0 ) INFO = K
            KP = K
            A( K, K ) = DBLE( A( K, K ) )
         } else {

            // ============================================================

            // Test for interchange

            if ( ABSAKK.GE.ALPHA*COLMAX ) {

               // no interchange, use 1-by-1 pivot block

               KP = K
            } else {

               // JMAX is the column-index of the largest off-diagonal
               // element in row IMAX, and ROWMAX is its absolute value.
               // Determine only ROWMAX.

               JMAX = IMAX + IZAMAX( K-IMAX, A( IMAX, IMAX+1 ), LDA )
               ROWMAX = CABS1( A( IMAX, JMAX ) )
               if ( IMAX.GT.1 ) {
                  JMAX = IZAMAX( IMAX-1, A( 1, IMAX ), 1 )
                  ROWMAX = MAX( ROWMAX, CABS1( A( JMAX, IMAX ) ) )
               }

               if ( ABSAKK.GE.ALPHA*COLMAX*( COLMAX / ROWMAX ) ) {

                  // no interchange, use 1-by-1 pivot block

                  KP = K

               } else if ( ABS( DBLE( A( IMAX, IMAX ) ) ).GE.ALPHA*ROWMAX ) {

                  // interchange rows and columns K and IMAX, use 1-by-1
                  // pivot block

                  KP = IMAX
               } else {

                  // interchange rows and columns K-1 and IMAX, use 2-by-2
                  // pivot block

                  KP = IMAX
                  KSTEP = 2
               }

            }

            // ============================================================

            KK = K - KSTEP + 1
            if ( KP.NE.KK ) {

               // Interchange rows and columns KK and KP in the leading
               // submatrix A(1:k,1:k)

               zswap(KP-1, A( 1, KK ), 1, A( 1, KP ), 1 );
               for (J = KP + 1; J <= KK - 1; J++) { // 20
                  T = DCONJG( A( J, KK ) )
                  A( J, KK ) = DCONJG( A( KP, J ) )
                  A( KP, J ) = T
               } // 20
               A( KP, KK ) = DCONJG( A( KP, KK ) )
               R1 = DBLE( A( KK, KK ) )
               A( KK, KK ) = DBLE( A( KP, KP ) )
               A( KP, KP ) = R1
               if ( KSTEP.EQ.2 ) {
                  A( K, K ) = DBLE( A( K, K ) )
                  T = A( K-1, K )
                  A( K-1, K ) = A( KP, K )
                  A( KP, K ) = T
               }
            } else {
               A( K, K ) = DBLE( A( K, K ) )
               IF( KSTEP.EQ.2 ) A( K-1, K-1 ) = DBLE( A( K-1, K-1 ) )
            }

            // Update the leading submatrix

            if ( KSTEP.EQ.1 ) {

               // 1-by-1 pivot block D(k): column k now holds

               // W(k) = U(k)*D(k)

               // where U(k) is the k-th column of U

               // Perform a rank-1 update of A(1:k-1,1:k-1) as

               // A := A - U(k)*D(k)*U(k)**H = A - W(k)*1/D(k)*W(k)**H

               R1 = ONE / DBLE( A( K, K ) )
               zher(UPLO, K-1, -R1, A( 1, K ), 1, A, LDA );

               // Store U(k) in column k

               zdscal(K-1, R1, A( 1, K ), 1 );
            } else {

               // 2-by-2 pivot block D(k): columns k and k-1 now hold

               // ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k)

               // where U(k) and U(k-1) are the k-th and (k-1)-th columns
               // of U

               // Perform a rank-2 update of A(1:k-2,1:k-2) as

               // A := A - ( U(k-1) U(k) )*D(k)*( U(k-1) U(k) )**H
                  // = A - ( W(k-1) W(k) )*inv(D(k))*( W(k-1) W(k) )**H

               if ( K.GT.2 ) {

                  D = DLAPY2( DBLE( A( K-1, K ) ), DIMAG( A( K-1, K ) ) )
                  D22 = DBLE( A( K-1, K-1 ) ) / D
                  D11 = DBLE( A( K, K ) ) / D
                  TT = ONE / ( D11*D22-ONE )
                  D12 = A( K-1, K ) / D
                  D = TT / D

                  DO 40 J = K - 2, 1, -1
                     WKM1 = D*( D11*A( J, K-1 )-DCONJG( D12 )* A( J, K ) )
                     WK = D*( D22*A( J, K )-D12*A( J, K-1 ) )
                     DO 30 I = J, 1, -1
                        A( I, J ) = A( I, J ) - A( I, K )*DCONJG( WK ) - A( I, K-1 )*DCONJG( WKM1 )
                     } // 30
                     A( J, K ) = WK
                     A( J, K-1 ) = WKM1
                     A( J, J ) = DCMPLX( DBLE( A( J, J ) ), 0.0D+0 )
                  } // 40

               }

            }
         }

         // Store details of the interchanges in IPIV

         if ( KSTEP.EQ.1 ) {
            IPIV( K ) = KP
         } else {
            IPIV( K ) = -KP
            IPIV( K-1 ) = -KP
         }

         // Decrease K and return to the start of the main loop

         K = K - KSTEP
         GO TO 10

      } else {

         // Factorize A as L*D*L**H using the lower triangle of A

         // K is the main loop index, increasing from 1 to N in steps of
         // 1 or 2

         K = 1
         } // 50

         // If K > N, exit from loop

         IF( K.GT.N ) GO TO 90
         KSTEP = 1

         // Determine rows and columns to be interchanged and whether
         // a 1-by-1 or 2-by-2 pivot block will be used

         ABSAKK = ABS( DBLE( A( K, K ) ) )

         // IMAX is the row-index of the largest off-diagonal element in
         // column K, and COLMAX is its absolute value.
         // Determine both COLMAX and IMAX.

         if ( K.LT.N ) {
            IMAX = K + IZAMAX( N-K, A( K+1, K ), 1 )
            COLMAX = CABS1( A( IMAX, K ) )
         } else {
            COLMAX = ZERO
         }

         if ( (MAX( ABSAKK, COLMAX ).EQ.ZERO) .OR. DISNAN(ABSAKK) ) {

            // Column K is zero or underflow, or contains a NaN:
            // set INFO and continue

            IF( INFO.EQ.0 ) INFO = K
            KP = K
            A( K, K ) = DBLE( A( K, K ) )
         } else {

            // ============================================================

            // Test for interchange

            if ( ABSAKK.GE.ALPHA*COLMAX ) {

               // no interchange, use 1-by-1 pivot block

               KP = K
            } else {

               // JMAX is the column-index of the largest off-diagonal
               // element in row IMAX, and ROWMAX is its absolute value.
               // Determine only ROWMAX.

               JMAX = K - 1 + IZAMAX( IMAX-K, A( IMAX, K ), LDA )
               ROWMAX = CABS1( A( IMAX, JMAX ) )
               if ( IMAX.LT.N ) {
                  JMAX = IMAX + IZAMAX( N-IMAX, A( IMAX+1, IMAX ), 1 )
                  ROWMAX = MAX( ROWMAX, CABS1( A( JMAX, IMAX ) ) )
               }

               if ( ABSAKK.GE.ALPHA*COLMAX*( COLMAX / ROWMAX ) ) {

                  // no interchange, use 1-by-1 pivot block

                  KP = K

               } else if ( ABS( DBLE( A( IMAX, IMAX ) ) ).GE.ALPHA*ROWMAX ) {

                  // interchange rows and columns K and IMAX, use 1-by-1
                  // pivot block

                  KP = IMAX
               } else {

                  // interchange rows and columns K+1 and IMAX, use 2-by-2
                  // pivot block

                  KP = IMAX
                  KSTEP = 2
               }

            }

            // ============================================================

            KK = K + KSTEP - 1
            if ( KP.NE.KK ) {

               // Interchange rows and columns KK and KP in the trailing
               // submatrix A(k:n,k:n)

               IF( KP.LT.N ) CALL ZSWAP( N-KP, A( KP+1, KK ), 1, A( KP+1, KP ), 1 )
               for (J = KK + 1; J <= KP - 1; J++) { // 60
                  T = DCONJG( A( J, KK ) )
                  A( J, KK ) = DCONJG( A( KP, J ) )
                  A( KP, J ) = T
               } // 60
               A( KP, KK ) = DCONJG( A( KP, KK ) )
               R1 = DBLE( A( KK, KK ) )
               A( KK, KK ) = DBLE( A( KP, KP ) )
               A( KP, KP ) = R1
               if ( KSTEP.EQ.2 ) {
                  A( K, K ) = DBLE( A( K, K ) )
                  T = A( K+1, K )
                  A( K+1, K ) = A( KP, K )
                  A( KP, K ) = T
               }
            } else {
               A( K, K ) = DBLE( A( K, K ) )
               IF( KSTEP.EQ.2 ) A( K+1, K+1 ) = DBLE( A( K+1, K+1 ) )
            }

            // Update the trailing submatrix

            if ( KSTEP.EQ.1 ) {

               // 1-by-1 pivot block D(k): column k now holds

               // W(k) = L(k)*D(k)

               // where L(k) is the k-th column of L

               if ( K.LT.N ) {

                  // Perform a rank-1 update of A(k+1:n,k+1:n) as

                  // A := A - L(k)*D(k)*L(k)**H = A - W(k)*(1/D(k))*W(k)**H

                  R1 = ONE / DBLE( A( K, K ) )
                  zher(UPLO, N-K, -R1, A( K+1, K ), 1, A( K+1, K+1 ), LDA );

                  // Store L(k) in column K

                  zdscal(N-K, R1, A( K+1, K ), 1 );
               }
            } else {

               // 2-by-2 pivot block D(k)

               if ( K.LT.N-1 ) {

                  // Perform a rank-2 update of A(k+2:n,k+2:n) as

                  // A := A - ( L(k) L(k+1) )*D(k)*( L(k) L(k+1) )**H
                     // = A - ( W(k) W(k+1) )*inv(D(k))*( W(k) W(k+1) )**H

                  // where L(k) and L(k+1) are the k-th and (k+1)-th
                  // columns of L

                  D = DLAPY2( DBLE( A( K+1, K ) ), DIMAG( A( K+1, K ) ) )
                  D11 = DBLE( A( K+1, K+1 ) ) / D
                  D22 = DBLE( A( K, K ) ) / D
                  TT = ONE / ( D11*D22-ONE )
                  D21 = A( K+1, K ) / D
                  D = TT / D

                  for (J = K + 2; J <= N; J++) { // 80
                     WK = D*( D11*A( J, K )-D21*A( J, K+1 ) )
                     WKP1 = D*( D22*A( J, K+1 )-DCONJG( D21 )* A( J, K ) )
                     for (I = J; I <= N; I++) { // 70
                        A( I, J ) = A( I, J ) - A( I, K )*DCONJG( WK ) - A( I, K+1 )*DCONJG( WKP1 )
                     } // 70
                     A( J, K ) = WK
                     A( J, K+1 ) = WKP1
                     A( J, J ) = DCMPLX( DBLE( A( J, J ) ), 0.0D+0 )
                  } // 80
               }
            }
         }

         // Store details of the interchanges in IPIV

         if ( KSTEP.EQ.1 ) {
            IPIV( K ) = KP
         } else {
            IPIV( K ) = -KP
            IPIV( K+1 ) = -KP
         }

         // Increase K and return to the start of the main loop

         K = K + KSTEP
         GO TO 50

      }

      } // 90
      RETURN

      // End of ZHETF2

      }
