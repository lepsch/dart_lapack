      SUBROUTINE DLASYF_RK( UPLO, N, NB, KB, A, LDA, E, IPIV, W, LDW, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, KB, LDA, LDW, N, NB;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      double             A( LDA, * ), E( * ), W( LDW, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      double             EIGHT, SEVTEN;
      const              EIGHT = 8.0D+0, SEVTEN = 17.0D+0 ;
      // ..
      // .. Local Scalars ..
      bool               DONE;
      int                IMAX, ITEMP, J, JB, JJ, JMAX, K, KK, KW, KKW, KP, KSTEP, P, II;
      double             ABSAKK, ALPHA, COLMAX, D11, D12, D21, D22, DTEMP, R1, ROWMAX, T, SFMIN;
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                IDAMAX;
      double             DLAMCH;
      // EXTERNAL LSAME, IDAMAX, DLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL DCOPY, DGEMM, DGEMV, DSCAL, DSWAP
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN, SQRT
      // ..
      // .. Executable Statements ..

      INFO = 0

      // Initialize ALPHA for use in choosing pivot block size.

      ALPHA = ( ONE+SQRT( SEVTEN ) ) / EIGHT

      // Compute machine safe minimum

      SFMIN = DLAMCH( 'S' )

      if ( LSAME( UPLO, 'U' ) ) {

         // Factorize the trailing columns of A using the upper triangle
         // of A and working backwards, and compute the matrix W = U12*D
         // for use in updating A11

         // Initialize the first entry of array E, where superdiagonal
         // elements of D are stored

         E( 1 ) = ZERO

         // K is the main loop index, decreasing from N in steps of 1 or 2

         K = N
   10    CONTINUE

         // KW is the column of W which corresponds to column K of A

         KW = NB + K - N

         // Exit from loop

         IF( ( K.LE.N-NB+1 .AND. NB.LT.N ) .OR. K.LT.1 ) GO TO 30

         KSTEP = 1
         P = K

         // Copy column K of A to column KW of W and update it

         dcopy(K, A( 1, K ), 1, W( 1, KW ), 1 );
         IF( K.LT.N ) CALL DGEMV( 'No transpose', K, N-K, -ONE, A( 1, K+1 ), LDA, W( K, KW+1 ), LDW, ONE, W( 1, KW ), 1 )

         // Determine rows and columns to be interchanged and whether
         // a 1-by-1 or 2-by-2 pivot block will be used

         ABSAKK = ABS( W( K, KW ) )

         // IMAX is the row-index of the largest off-diagonal element in
         // column K, and COLMAX is its absolute value.
         // Determine both COLMAX and IMAX.

         if ( K.GT.1 ) {
            IMAX = IDAMAX( K-1, W( 1, KW ), 1 )
            COLMAX = ABS( W( IMAX, KW ) )
         } else {
            COLMAX = ZERO
         }

         if ( MAX( ABSAKK, COLMAX ).EQ.ZERO ) {

            // Column K is zero or underflow: set INFO and continue

            IF( INFO.EQ.0 ) INFO = K
            KP = K
            dcopy(K, W( 1, KW ), 1, A( 1, K ), 1 );

            // Set E( K ) to zero

            IF( K.GT.1 ) E( K ) = ZERO

         } else {

            // ============================================================

            // Test for interchange

            // Equivalent to testing for ABSAKK.GE.ALPHA*COLMAX
            // (used to handle NaN and Inf)

            if ( .NOT.( ABSAKK.LT.ALPHA*COLMAX ) ) {

               // no interchange, use 1-by-1 pivot block

               KP = K

            } else {

               DONE = .FALSE.

               // Loop until pivot found

   12          CONTINUE

                  // Begin pivot search loop body


                  // Copy column IMAX to column KW-1 of W and update it

                  dcopy(IMAX, A( 1, IMAX ), 1, W( 1, KW-1 ), 1 );
                  dcopy(K-IMAX, A( IMAX, IMAX+1 ), LDA, W( IMAX+1, KW-1 ), 1 );

                  IF( K.LT.N ) CALL DGEMV( 'No transpose', K, N-K, -ONE, A( 1, K+1 ), LDA, W( IMAX, KW+1 ), LDW, ONE, W( 1, KW-1 ), 1 )

                  // JMAX is the column-index of the largest off-diagonal
                  // element in row IMAX, and ROWMAX is its absolute value.
                  // Determine both ROWMAX and JMAX.

                  if ( IMAX.NE.K ) {
                     JMAX = IMAX + IDAMAX( K-IMAX, W( IMAX+1, KW-1 ), 1 )
                     ROWMAX = ABS( W( JMAX, KW-1 ) )
                  } else {
                     ROWMAX = ZERO
                  }

                  if ( IMAX.GT.1 ) {
                     ITEMP = IDAMAX( IMAX-1, W( 1, KW-1 ), 1 )
                     DTEMP = ABS( W( ITEMP, KW-1 ) )
                     if ( DTEMP.GT.ROWMAX ) {
                        ROWMAX = DTEMP
                        JMAX = ITEMP
                     }
                  }

                  // Equivalent to testing for
                  // ABS( W( IMAX, KW-1 ) ).GE.ALPHA*ROWMAX
                  // (used to handle NaN and Inf)

                  if ( .NOT.(ABS( W( IMAX, KW-1 ) ).LT.ALPHA*ROWMAX ) ) {

                     // interchange rows and columns K and IMAX,
                     // use 1-by-1 pivot block

                     KP = IMAX

                     // copy column KW-1 of W to column KW of W

                     dcopy(K, W( 1, KW-1 ), 1, W( 1, KW ), 1 );

                     DONE = .TRUE.

                  // Equivalent to testing for ROWMAX.EQ.COLMAX,
                  // (used to handle NaN and Inf)

                  } else if ( ( P.EQ.JMAX ) .OR. ( ROWMAX.LE.COLMAX ) ) {

                     // interchange rows and columns K-1 and IMAX,
                     // use 2-by-2 pivot block

                     KP = IMAX
                     KSTEP = 2
                     DONE = .TRUE.
                  } else {

                     // Pivot not found: set params and repeat

                     P = IMAX
                     COLMAX = ROWMAX
                     IMAX = JMAX

                     // Copy updated JMAXth (next IMAXth) column to Kth of W

                     dcopy(K, W( 1, KW-1 ), 1, W( 1, KW ), 1 );

                  }

                  // End pivot search loop body

               IF( .NOT. DONE ) GOTO 12

            }

            // ============================================================

            KK = K - KSTEP + 1

            // KKW is the column of W which corresponds to column KK of A

            KKW = NB + KK - N

            if ( ( KSTEP.EQ.2 ) .AND. ( P.NE.K ) ) {

               // Copy non-updated column K to column P

               dcopy(K-P, A( P+1, K ), 1, A( P, P+1 ), LDA );
               dcopy(P, A( 1, K ), 1, A( 1, P ), 1 );

               // Interchange rows K and P in last N-K+1 columns of A
               // and last N-K+2 columns of W

               dswap(N-K+1, A( K, K ), LDA, A( P, K ), LDA );
               dswap(N-KK+1, W( K, KKW ), LDW, W( P, KKW ), LDW );
            }

            // Updated column KP is already stored in column KKW of W

            if ( KP.NE.KK ) {

               // Copy non-updated column KK to column KP

               A( KP, K ) = A( KK, K )
               dcopy(K-1-KP, A( KP+1, KK ), 1, A( KP, KP+1 ), LDA );
               dcopy(KP, A( 1, KK ), 1, A( 1, KP ), 1 );

               // Interchange rows KK and KP in last N-KK+1 columns
               // of A and W

               dswap(N-KK+1, A( KK, KK ), LDA, A( KP, KK ), LDA );
               dswap(N-KK+1, W( KK, KKW ), LDW, W( KP, KKW ), LDW );
            }

            if ( KSTEP.EQ.1 ) {

               // 1-by-1 pivot block D(k): column KW of W now holds

               // W(k) = U(k)*D(k)

               // where U(k) is the k-th column of U

               // Store U(k) in column k of A

               dcopy(K, W( 1, KW ), 1, A( 1, K ), 1 );
               if ( K.GT.1 ) {
                  if ( ABS( A( K, K ) ).GE.SFMIN ) {
                     R1 = ONE / A( K, K )
                     dscal(K-1, R1, A( 1, K ), 1 );
                  } else if ( A( K, K ).NE.ZERO ) {
                     DO 14 II = 1, K - 1
                        A( II, K ) = A( II, K ) / A( K, K )
   14                CONTINUE
                  }

                  // Store the superdiagonal element of D in array E

                  E( K ) = ZERO

               }

            } else {

               // 2-by-2 pivot block D(k): columns KW and KW-1 of W now
               // hold

               // ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k)

               // where U(k) and U(k-1) are the k-th and (k-1)-th columns
               // of U

               if ( K.GT.2 ) {

                  // Store U(k) and U(k-1) in columns k and k-1 of A

                  D12 = W( K-1, KW )
                  D11 = W( K, KW ) / D12
                  D22 = W( K-1, KW-1 ) / D12
                  T = ONE / ( D11*D22-ONE )
                  DO 20 J = 1, K - 2
                     A( J, K-1 ) = T*( (D11*W( J, KW-1 )-W( J, KW ) ) / D12 )                      A( J, K ) = T*( ( D22*W( J, KW )-W( J, KW-1 ) ) / D12 )
   20             CONTINUE
               }

               // Copy diagonal elements of D(K) to A,
               // copy superdiagonal element of D(K) to E(K) and
               // ZERO out superdiagonal entry of A

               A( K-1, K-1 ) = W( K-1, KW-1 )
               A( K-1, K ) = ZERO
               A( K, K ) = W( K, KW )
               E( K ) = W( K-1, KW )
               E( K-1 ) = ZERO

            }

            // End column K is nonsingular

         }

         // Store details of the interchanges in IPIV

         if ( KSTEP.EQ.1 ) {
            IPIV( K ) = KP
         } else {
            IPIV( K ) = -P
            IPIV( K-1 ) = -KP
         }

         // Decrease K and return to the start of the main loop

         K = K - KSTEP
         GO TO 10

   30    CONTINUE

         // Update the upper triangle of A11 (= A(1:k,1:k)) as

         // A11 := A11 - U12*D*U12**T = A11 - U12*W**T

         // computing blocks of NB columns at a time

         DO 50 J = ( ( K-1 ) / NB )*NB + 1, 1, -NB
            JB = MIN( NB, K-J+1 )

            // Update the upper triangle of the diagonal block

            DO 40 JJ = J, J + JB - 1
               dgemv('No transpose', JJ-J+1, N-K, -ONE, A( J, K+1 ), LDA, W( JJ, KW+1 ), LDW, ONE, A( J, JJ ), 1 );
   40       CONTINUE

            // Update the rectangular superdiagonal block

            IF( J.GE.2 ) CALL DGEMM( 'No transpose', 'Transpose', J-1, JB, N-K, -ONE, A( 1, K+1 ), LDA, W( J, KW+1 ), LDW, ONE, A( 1, J ), LDA )
   50    CONTINUE

         // Set KB to the number of columns factorized

         KB = N - K

      } else {

         // Factorize the leading columns of A using the lower triangle
         // of A and working forwards, and compute the matrix W = L21*D
         // for use in updating A22

         // Initialize the unused last entry of the subdiagonal array E.

         E( N ) = ZERO

         // K is the main loop index, increasing from 1 in steps of 1 or 2

         K = 1
   70   CONTINUE

         // Exit from loop

         IF( ( K.GE.NB .AND. NB.LT.N ) .OR. K.GT.N ) GO TO 90

         KSTEP = 1
         P = K

         // Copy column K of A to column K of W and update it

         dcopy(N-K+1, A( K, K ), 1, W( K, K ), 1 );
         IF( K.GT.1 ) CALL DGEMV( 'No transpose', N-K+1, K-1, -ONE, A( K, 1 ), LDA, W( K, 1 ), LDW, ONE, W( K, K ), 1 )

         // Determine rows and columns to be interchanged and whether
         // a 1-by-1 or 2-by-2 pivot block will be used

         ABSAKK = ABS( W( K, K ) )

         // IMAX is the row-index of the largest off-diagonal element in
         // column K, and COLMAX is its absolute value.
         // Determine both COLMAX and IMAX.

         if ( K.LT.N ) {
            IMAX = K + IDAMAX( N-K, W( K+1, K ), 1 )
            COLMAX = ABS( W( IMAX, K ) )
         } else {
            COLMAX = ZERO
         }

         if ( MAX( ABSAKK, COLMAX ).EQ.ZERO ) {

            // Column K is zero or underflow: set INFO and continue

            IF( INFO.EQ.0 ) INFO = K
            KP = K
            dcopy(N-K+1, W( K, K ), 1, A( K, K ), 1 );

            // Set E( K ) to zero

            IF( K.LT.N ) E( K ) = ZERO

         } else {

            // ============================================================

            // Test for interchange

            // Equivalent to testing for ABSAKK.GE.ALPHA*COLMAX
            // (used to handle NaN and Inf)

            if ( .NOT.( ABSAKK.LT.ALPHA*COLMAX ) ) {

               // no interchange, use 1-by-1 pivot block

               KP = K

            } else {

               DONE = .FALSE.

               // Loop until pivot found

   72          CONTINUE

                  // Begin pivot search loop body


                  // Copy column IMAX to column K+1 of W and update it

                  dcopy(IMAX-K, A( IMAX, K ), LDA, W( K, K+1 ), 1);
                  dcopy(N-IMAX+1, A( IMAX, IMAX ), 1, W( IMAX, K+1 ), 1 )                   IF( K.GT.1 ) CALL DGEMV( 'No transpose', N-K+1, K-1, -ONE, A( K, 1 ), LDA, W( IMAX, 1 ), LDW, ONE, W( K, K+1 ), 1 );

                  // JMAX is the column-index of the largest off-diagonal
                  // element in row IMAX, and ROWMAX is its absolute value.
                  // Determine both ROWMAX and JMAX.

                  if ( IMAX.NE.K ) {
                     JMAX = K - 1 + IDAMAX( IMAX-K, W( K, K+1 ), 1 )
                     ROWMAX = ABS( W( JMAX, K+1 ) )
                  } else {
                     ROWMAX = ZERO
                  }

                  if ( IMAX.LT.N ) {
                     ITEMP = IMAX + IDAMAX( N-IMAX, W( IMAX+1, K+1 ), 1)
                     DTEMP = ABS( W( ITEMP, K+1 ) )
                     if ( DTEMP.GT.ROWMAX ) {
                        ROWMAX = DTEMP
                        JMAX = ITEMP
                     }
                  }

                  // Equivalent to testing for
                  // ABS( W( IMAX, K+1 ) ).GE.ALPHA*ROWMAX
                  // (used to handle NaN and Inf)

                  if ( .NOT.( ABS( W( IMAX, K+1 ) ).LT.ALPHA*ROWMAX ) ) {

                     // interchange rows and columns K and IMAX,
                     // use 1-by-1 pivot block

                     KP = IMAX

                     // copy column K+1 of W to column K of W

                     dcopy(N-K+1, W( K, K+1 ), 1, W( K, K ), 1 );

                     DONE = .TRUE.

                  // Equivalent to testing for ROWMAX.EQ.COLMAX,
                  // (used to handle NaN and Inf)

                  } else if ( ( P.EQ.JMAX ) .OR. ( ROWMAX.LE.COLMAX ) ) {

                     // interchange rows and columns K+1 and IMAX,
                     // use 2-by-2 pivot block

                     KP = IMAX
                     KSTEP = 2
                     DONE = .TRUE.
                  } else {

                     // Pivot not found: set params and repeat

                     P = IMAX
                     COLMAX = ROWMAX
                     IMAX = JMAX

                     // Copy updated JMAXth (next IMAXth) column to Kth of W

                     dcopy(N-K+1, W( K, K+1 ), 1, W( K, K ), 1 );

                  }

                  // End pivot search loop body

               IF( .NOT. DONE ) GOTO 72

            }

            // ============================================================

            KK = K + KSTEP - 1

            if ( ( KSTEP.EQ.2 ) .AND. ( P.NE.K ) ) {

               // Copy non-updated column K to column P

               dcopy(P-K, A( K, K ), 1, A( P, K ), LDA );
               dcopy(N-P+1, A( P, K ), 1, A( P, P ), 1 );

               // Interchange rows K and P in first K columns of A
               // and first K+1 columns of W

               dswap(K, A( K, 1 ), LDA, A( P, 1 ), LDA );
               dswap(KK, W( K, 1 ), LDW, W( P, 1 ), LDW );
            }

            // Updated column KP is already stored in column KK of W

            if ( KP.NE.KK ) {

               // Copy non-updated column KK to column KP

               A( KP, K ) = A( KK, K )
               dcopy(KP-K-1, A( K+1, KK ), 1, A( KP, K+1 ), LDA );
               dcopy(N-KP+1, A( KP, KK ), 1, A( KP, KP ), 1 );

               // Interchange rows KK and KP in first KK columns of A and W

               dswap(KK, A( KK, 1 ), LDA, A( KP, 1 ), LDA );
               dswap(KK, W( KK, 1 ), LDW, W( KP, 1 ), LDW );
            }

            if ( KSTEP.EQ.1 ) {

               // 1-by-1 pivot block D(k): column k of W now holds

               // W(k) = L(k)*D(k)

               // where L(k) is the k-th column of L

               // Store L(k) in column k of A

               dcopy(N-K+1, W( K, K ), 1, A( K, K ), 1 );
               if ( K.LT.N ) {
                  if ( ABS( A( K, K ) ).GE.SFMIN ) {
                     R1 = ONE / A( K, K )
                     dscal(N-K, R1, A( K+1, K ), 1 );
                  } else if ( A( K, K ).NE.ZERO ) {
                     DO 74 II = K + 1, N
                        A( II, K ) = A( II, K ) / A( K, K )
   74                CONTINUE
                  }

                  // Store the subdiagonal element of D in array E

                  E( K ) = ZERO

               }

            } else {

               // 2-by-2 pivot block D(k): columns k and k+1 of W now hold

               // ( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k)

               // where L(k) and L(k+1) are the k-th and (k+1)-th columns
               // of L

               if ( K.LT.N-1 ) {

                  // Store L(k) and L(k+1) in columns k and k+1 of A

                  D21 = W( K+1, K )
                  D11 = W( K+1, K+1 ) / D21
                  D22 = W( K, K ) / D21
                  T = ONE / ( D11*D22-ONE )
                  DO 80 J = K + 2, N
                     A( J, K ) = T*( ( D11*W( J, K )-W( J, K+1 ) ) / D21 )                      A( J, K+1 ) = T*( ( D22*W( J, K+1 )-W( J, K ) ) / D21 )
   80             CONTINUE
               }

               // Copy diagonal elements of D(K) to A,
               // copy subdiagonal element of D(K) to E(K) and
               // ZERO out subdiagonal entry of A

               A( K, K ) = W( K, K )
               A( K+1, K ) = ZERO
               A( K+1, K+1 ) = W( K+1, K+1 )
               E( K ) = W( K+1, K )
               E( K+1 ) = ZERO

            }

            // End column K is nonsingular

         }

         // Store details of the interchanges in IPIV

         if ( KSTEP.EQ.1 ) {
            IPIV( K ) = KP
         } else {
            IPIV( K ) = -P
            IPIV( K+1 ) = -KP
         }

         // Increase K and return to the start of the main loop

         K = K + KSTEP
         GO TO 70

   90    CONTINUE

         // Update the lower triangle of A22 (= A(k:n,k:n)) as

         // A22 := A22 - L21*D*L21**T = A22 - L21*W**T

         // computing blocks of NB columns at a time

         DO 110 J = K, N, NB
            JB = MIN( NB, N-J+1 )

            // Update the lower triangle of the diagonal block

            DO 100 JJ = J, J + JB - 1
               dgemv('No transpose', J+JB-JJ, K-1, -ONE, A( JJ, 1 ), LDA, W( JJ, 1 ), LDW, ONE, A( JJ, JJ ), 1 );
  100       CONTINUE

            // Update the rectangular subdiagonal block

            IF( J+JB.LE.N ) CALL DGEMM( 'No transpose', 'Transpose', N-J-JB+1, JB, K-1, -ONE, A( J+JB, 1 ), LDA, W( J, 1 ), LDW, ONE, A( J+JB, J ), LDA )
  110    CONTINUE

         // Set KB to the number of columns factorized

         KB = K - 1

      }

      RETURN

      // End of DLASYF_RK

      }
