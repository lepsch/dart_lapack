import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

      void dlasyf_rook(UPLO, N, NB, KB, final Matrix<double> A, final int LDA, IPIV, final Matrix<double> W, final int LDW, final Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                INFO, KB, LDA, LDW, N, NB;
      int                IPIV( * );
      double             A( LDA, * ), W( LDW, * );
      // ..

      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      double             EIGHT, SEVTEN;
      const              EIGHT = 8.0, SEVTEN = 17.0 ;
      bool               DONE;
      int                IMAX, ITEMP, J, JB, JJ, JMAX, JP1, JP2, K, KK, KW, KKW, KP, KSTEP, P, II;
       double             ABSAKK, ALPHA, COLMAX, D11, D12, D21, D22, DTEMP, R1, ROWMAX, T, SFMIN;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                idamax;
      //- double             DLAMCH;
      // EXTERNAL lsame, idamax, DLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL DCOPY, DGEMM, DGEMV, DSCAL, DSWAP
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN, SQRT

      INFO = 0;

      // Initialize ALPHA for use in choosing pivot block size.

      ALPHA = ( ONE+sqrt( SEVTEN ) ) / EIGHT;

      // Compute machine safe minimum

      SFMIN = dlamch( 'S' );

      if ( lsame( UPLO, 'U' ) ) {

         // Factorize the trailing columns of A using the upper triangle
         // of A and working backwards, and compute the matrix W = U12*D
         // for use in updating A11

         // K is the main loop index, decreasing from N in steps of 1 or 2

         K = N;
        //  } // 10

         // KW is the column of W which corresponds to column K of A

         KW = NB + K - N;

         // Exit from loop

         if( ( K <= N-NB+1 && NB < N ) || K < 1 ) GO TO 30;

         KSTEP = 1;
         P = K;

         // Copy column K of A to column KW of W and update it

         dcopy(K, A( 1, K ), 1, W( 1, KW ), 1 );
         if (K < N) dgemv( 'No transpose', K, N-K, -ONE, A( 1, K+1 ), LDA, W( K, KW+1 ), LDW, ONE, W( 1, KW ), 1 );

         // Determine rows and columns to be interchanged and whether
         // a 1-by-1 or 2-by-2 pivot block will be used

         ABSAKK = ( W( K, KW ) ).abs();

         // IMAX is the row-index of the largest off-diagonal element in
         // column K, and COLMAX is its absolute value.
         // Determine both COLMAX and IMAX.

         if ( K > 1 ) {
            IMAX = idamax( K-1, W( 1, KW ), 1 );
            COLMAX = ( W( IMAX, KW ) ).abs();
         } else {
            COLMAX = ZERO;
         }

         if ( max( ABSAKK, COLMAX ) == ZERO ) {

            // Column K is zero or underflow: set INFO and continue

            if (INFO == 0) INFO = K;
            KP = K;
            dcopy(K, W( 1, KW ), 1, A( 1, K ), 1 );
         } else {

            // ============================================================

            // Test for interchange

            // Equivalent to testing for ABSAKK >= ALPHA*COLMAX
            // (used to handle NaN and Inf)

            if ( !( ABSAKK < ALPHA*COLMAX ) ) {

               // no interchange, use 1-by-1 pivot block

               KP = K;

            } else {

               DONE = false;

               // Loop until pivot found

              //  } // 12

                  // Begin pivot search loop body


                  // Copy column IMAX to column KW-1 of W and update it

                  dcopy(IMAX, A( 1, IMAX ), 1, W( 1, KW-1 ), 1 );
                  dcopy(K-IMAX, A( IMAX, IMAX+1 ), LDA, W( IMAX+1, KW-1 ), 1 );

                  if (K < N) dgemv( 'No transpose', K, N-K, -ONE, A( 1, K+1 ), LDA, W( IMAX, KW+1 ), LDW, ONE, W( 1, KW-1 ), 1 );

                  // JMAX is the column-index of the largest off-diagonal
                  // element in row IMAX, and ROWMAX is its absolute value.
                  // Determine both ROWMAX and JMAX.

                  if ( IMAX != K ) {
                     JMAX = IMAX + idamax( K-IMAX, W( IMAX+1, KW-1 ), 1 );
                     ROWMAX = ( W( JMAX, KW-1 ) ).abs();
                  } else {
                     ROWMAX = ZERO;
                  }

                  if ( IMAX > 1 ) {
                     ITEMP = idamax( IMAX-1, W( 1, KW-1 ), 1 );
                     DTEMP = ( W( ITEMP, KW-1 ) ).abs();
                     if ( DTEMP > ROWMAX ) {
                        ROWMAX = DTEMP;
                        JMAX = ITEMP;
                     }
                  }

                  // Equivalent to testing for
                  // ABS( W( IMAX, KW-1 ) ) >= ALPHA*ROWMAX
                  // (used to handle NaN and Inf)

                  if ( !(( W( IMAX, KW-1 ) ).abs() < ALPHA*ROWMAX ) ) {

                     // interchange rows and columns K and IMAX,
                     // use 1-by-1 pivot block

                     KP = IMAX;

                     // copy column KW-1 of W to column KW of W

                     dcopy(K, W( 1, KW-1 ), 1, W( 1, KW ), 1 );

                     DONE = true;

                  // Equivalent to testing for ROWMAX == COLMAX,
                  // (used to handle NaN and Inf)

                  } else if ( ( P == JMAX ) || ( ROWMAX <= COLMAX ) ) {

                     // interchange rows and columns K-1 and IMAX,
                     // use 2-by-2 pivot block

                     KP = IMAX;
                     KSTEP = 2;
                     DONE = true;
                  } else {

                     // Pivot not found: set params and repeat

                     P = IMAX;
                     COLMAX = ROWMAX;
                     IMAX = JMAX;

                     // Copy updated JMAXth (next IMAXth) column to Kth of W

                     dcopy(K, W( 1, KW-1 ), 1, W( 1, KW ), 1 );

                  }

                  // End pivot search loop body

               if ( !DONE) GOTO 12;

            }

            // ============================================================

            KK = K - KSTEP + 1;

            // KKW is the column of W which corresponds to column KK of A

            KKW = NB + KK - N;

            if ( ( KSTEP == 2 ) && ( P != K ) ) {

               // Copy non-updated column K to column P

               dcopy(K-P, A( P+1, K ), 1, A( P, P+1 ), LDA );
               dcopy(P, A( 1, K ), 1, A( 1, P ), 1 );

               // Interchange rows K and P in last N-K+1 columns of A
               // and last N-K+2 columns of W

               dswap(N-K+1, A( K, K ), LDA, A( P, K ), LDA );
               dswap(N-KK+1, W( K, KKW ), LDW, W( P, KKW ), LDW );
            }

            // Updated column KP is already stored in column KKW of W

            if ( KP != KK ) {

               // Copy non-updated column KK to column KP

               A[KP][K] = A( KK, K );
               dcopy(K-1-KP, A( KP+1, KK ), 1, A( KP, KP+1 ), LDA );
               dcopy(KP, A( 1, KK ), 1, A( 1, KP ), 1 );

               // Interchange rows KK and KP in last N-KK+1 columns
               // of A and W

               dswap(N-KK+1, A( KK, KK ), LDA, A( KP, KK ), LDA );
               dswap(N-KK+1, W( KK, KKW ), LDW, W( KP, KKW ), LDW );
            }

            if ( KSTEP == 1 ) {

               // 1-by-1 pivot block D(k): column KW of W now holds

               // W(k) = U(k)*D(k)

               // where U(k) is the k-th column of U

               // Store U(k) in column k of A

               dcopy(K, W( 1, KW ), 1, A( 1, K ), 1 );
               if ( K > 1 ) {
                  if ( ( A( K, K ) ).abs() >= SFMIN ) {
                     R1 = ONE / A( K, K );
                     dscal(K-1, R1, A( 1, K ), 1 );
                  } else if ( A( K, K ) != ZERO ) {
                     for (II = 1; II <= K - 1; II++) { // 14
                        A[II][K] = A( II, K ) / A( K, K );
                     } // 14
                  }
               }

            } else {

               // 2-by-2 pivot block D(k): columns KW and KW-1 of W now
               // hold

               // ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k)

               // where U(k) and U(k-1) are the k-th and (k-1)-th columns
               // of U

               if ( K > 2 ) {

                  // Store U(k) and U(k-1) in columns k and k-1 of A

                  D12 = W( K-1, KW );
                  D11 = W( K, KW ) / D12;
                  D22 = W( K-1, KW-1 ) / D12;
                  T = ONE / ( D11*D22-ONE );
                  for (J = 1; J <= K - 2; J++) { // 20
                     A[J][K-1] = T*( (D11*W( J, KW-1 )-W( J, KW ) ) / D12 )                      A( J, K ) = T*( ( D22*W( J, KW )-W( J, KW-1 ) ) / D12 );
                  } // 20
               }

               // Copy D(k) to A

               A[K-1][K-1] = W( K-1, KW-1 );
               A[K-1][K] = W( K-1, KW );
               A[K][K] = W( K, KW );
            }
         }

         // Store details of the interchanges in IPIV

         if ( KSTEP == 1 ) {
            IPIV[K] = KP;
         } else {
            IPIV[K] = -P;
            IPIV[K-1] = -KP;
         }

         // Decrease K and return to the start of the main loop

         K = K - KSTEP;
         GO TO 10;

        //  } // 30

         // Update the upper triangle of A11 (= A(1:k,1:k)) as

         // A11 := A11 - U12*D*U12**T = A11 - U12*W**T

         // computing blocks of NB columns at a time

         for (J = ( ( K-1 ) / NB )*NB + 1; -NB < 0 ? J >= 1 : J <= 1; J += -NB) { // 50
            JB = min( NB, K-J+1 );

            // Update the upper triangle of the diagonal block

            for (JJ = J; JJ <= J + JB - 1; JJ++) { // 40
               dgemv('No transpose', JJ-J+1, N-K, -ONE, A( J, K+1 ), LDA, W( JJ, KW+1 ), LDW, ONE, A( J, JJ ), 1 );
            } // 40

            // Update the rectangular superdiagonal block

            if (J >= 2) dgemm( 'No transpose', 'Transpose', J-1, JB, N-K, -ONE, A( 1, K+1 ), LDA, W( J, KW+1 ), LDW, ONE, A( 1, J ), LDA );
         } // 50

         // Put U12 in standard form by partially undoing the interchanges
         // in columns k+1:n

         J = K + 1;
        //  } // 60

            KSTEP = 1;
            JP1 = 1;
            JJ = J;
            JP2 = IPIV( J );
            if ( JP2 < 0 ) {
               JP2 = -JP2;
               J = J + 1;
               JP1 = -IPIV( J );
               KSTEP = 2;
            }

            J = J + 1;
            if (JP2 != JJ && J <= N) dswap( N-J+1, A( JP2, J ), LDA, A( JJ, J ), LDA );
            JJ = J - 1;
            if (JP1 != JJ && KSTEP == 2) dswap( N-J+1, A( JP1, J ), LDA, A( JJ, J ), LDA );
            if( J <= N ) GO TO 60;

         // Set KB to the number of columns factorized

         KB = N - K;

      } else {

         // Factorize the leading columns of A using the lower triangle
         // of A and working forwards, and compute the matrix W = L21*D
         // for use in updating A22

         // K is the main loop index, increasing from 1 in steps of 1 or 2

         K = 1;
        // } // 70

         // Exit from loop

         if( ( K >= NB && NB < N ) || K > N ) GO TO 90;

         KSTEP = 1;
         P = K;

         // Copy column K of A to column K of W and update it

         dcopy(N-K+1, A( K, K ), 1, W( K, K ), 1 );
         if (K > 1) dgemv( 'No transpose', N-K+1, K-1, -ONE, A( K, 1 ), LDA, W( K, 1 ), LDW, ONE, W( K, K ), 1 );

         // Determine rows and columns to be interchanged and whether
         // a 1-by-1 or 2-by-2 pivot block will be used

         ABSAKK = ( W( K, K ) ).abs();

         // IMAX is the row-index of the largest off-diagonal element in
         // column K, and COLMAX is its absolute value.
         // Determine both COLMAX and IMAX.

         if ( K < N ) {
            IMAX = K + idamax( N-K, W( K+1, K ), 1 );
            COLMAX = ( W( IMAX, K ) ).abs();
         } else {
            COLMAX = ZERO;
         }

         if ( max( ABSAKK, COLMAX ) == ZERO ) {

            // Column K is zero or underflow: set INFO and continue

            if (INFO == 0) INFO = K;
            KP = K;
            dcopy(N-K+1, W( K, K ), 1, A( K, K ), 1 );
         } else {

            // ============================================================

            // Test for interchange

            // Equivalent to testing for ABSAKK >= ALPHA*COLMAX
            // (used to handle NaN and Inf)

            if ( !( ABSAKK < ALPHA*COLMAX ) ) {

               // no interchange, use 1-by-1 pivot block

               KP = K;

            } else {

               DONE = false;

               // Loop until pivot found

              //  } // 72

                  // Begin pivot search loop body


                  // Copy column IMAX to column K+1 of W and update it

                  dcopy(IMAX-K, A( IMAX, K ), LDA, W( K, K+1 ), 1);
                  dcopy(N-IMAX+1, A( IMAX, IMAX ), 1, W( IMAX, K+1 ), 1 )                   IF( K > 1 ) CALL DGEMV( 'No transpose', N-K+1, K-1, -ONE, A( K, 1 ), LDA, W( IMAX, 1 ), LDW, ONE, W( K, K+1 ), 1 );

                  // JMAX is the column-index of the largest off-diagonal
                  // element in row IMAX, and ROWMAX is its absolute value.
                  // Determine both ROWMAX and JMAX.

                  if ( IMAX != K ) {
                     JMAX = K - 1 + idamax( IMAX-K, W( K, K+1 ), 1 );
                     ROWMAX = ( W( JMAX, K+1 ) ).abs();
                  } else {
                     ROWMAX = ZERO;
                  }

                  if ( IMAX < N ) {
                     ITEMP = IMAX + idamax( N-IMAX, W( IMAX+1, K+1 ), 1);
                     DTEMP = ( W( ITEMP, K+1 ) ).abs();
                     if ( DTEMP > ROWMAX ) {
                        ROWMAX = DTEMP;
                        JMAX = ITEMP;
                     }
                  }

                  // Equivalent to testing for
                  // ABS( W( IMAX, K+1 ) ) >= ALPHA*ROWMAX
                  // (used to handle NaN and Inf)

                  if ( !( ( W( IMAX, K+1 ) ).abs() < ALPHA*ROWMAX ) ) {

                     // interchange rows and columns K and IMAX,
                     // use 1-by-1 pivot block

                     KP = IMAX;

                     // copy column K+1 of W to column K of W

                     dcopy(N-K+1, W( K, K+1 ), 1, W( K, K ), 1 );

                     DONE = true;

                  // Equivalent to testing for ROWMAX == COLMAX,
                  // (used to handle NaN and Inf)

                  } else if ( ( P == JMAX ) || ( ROWMAX <= COLMAX ) ) {

                     // interchange rows and columns K+1 and IMAX,
                     // use 2-by-2 pivot block

                     KP = IMAX;
                     KSTEP = 2;
                     DONE = true;
                  } else {

                     // Pivot not found: set params and repeat

                     P = IMAX;
                     COLMAX = ROWMAX;
                     IMAX = JMAX;

                     // Copy updated JMAXth (next IMAXth) column to Kth of W

                     dcopy(N-K+1, W( K, K+1 ), 1, W( K, K ), 1 );

                  }

                  // End pivot search loop body

               if ( !DONE) GOTO 72;

            }

            // ============================================================

            KK = K + KSTEP - 1;

            if ( ( KSTEP == 2 ) && ( P != K ) ) {

               // Copy non-updated column K to column P

               dcopy(P-K, A( K, K ), 1, A( P, K ), LDA );
               dcopy(N-P+1, A( P, K ), 1, A( P, P ), 1 );

               // Interchange rows K and P in first K columns of A
               // and first K+1 columns of W

               dswap(K, A( K, 1 ), LDA, A( P, 1 ), LDA );
               dswap(KK, W( K, 1 ), LDW, W( P, 1 ), LDW );
            }

            // Updated column KP is already stored in column KK of W

            if ( KP != KK ) {

               // Copy non-updated column KK to column KP

               A[KP][K] = A( KK, K );
               dcopy(KP-K-1, A( K+1, KK ), 1, A( KP, K+1 ), LDA );
               dcopy(N-KP+1, A( KP, KK ), 1, A( KP, KP ), 1 );

               // Interchange rows KK and KP in first KK columns of A and W

               dswap(KK, A( KK, 1 ), LDA, A( KP, 1 ), LDA );
               dswap(KK, W( KK, 1 ), LDW, W( KP, 1 ), LDW );
            }

            if ( KSTEP == 1 ) {

               // 1-by-1 pivot block D(k): column k of W now holds

               // W(k) = L(k)*D(k)

               // where L(k) is the k-th column of L

               // Store L(k) in column k of A

               dcopy(N-K+1, W( K, K ), 1, A( K, K ), 1 );
               if ( K < N ) {
                  if ( ( A( K, K ) ).abs() >= SFMIN ) {
                     R1 = ONE / A( K, K );
                     dscal(N-K, R1, A( K+1, K ), 1 );
                  } else if ( A( K, K ) != ZERO ) {
                     for (II = K + 1; II <= N; II++) { // 74
                        A[II][K] = A( II, K ) / A( K, K );
                     } // 74
                  }
               }

            } else {

               // 2-by-2 pivot block D(k): columns k and k+1 of W now hold

               // ( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k)

               // where L(k) and L(k+1) are the k-th and (k+1)-th columns
               // of L

               if ( K < N-1 ) {

                  // Store L(k) and L(k+1) in columns k and k+1 of A

                  D21 = W( K+1, K );
                  D11 = W( K+1, K+1 ) / D21;
                  D22 = W( K, K ) / D21;
                  T = ONE / ( D11*D22-ONE );
                  for (J = K + 2; J <= N; J++) { // 80
                     A[J][K] = T*( ( D11*W( J, K )-W( J, K+1 ) ) / D21 )                      A( J, K+1 ) = T*( ( D22*W( J, K+1 )-W( J, K ) ) / D21 );
                  } // 80
               }

               // Copy D(k) to A

               A[K][K] = W( K, K );
               A[K+1][K] = W( K+1, K );
               A[K+1][K+1] = W( K+1, K+1 );
            }
         }

         // Store details of the interchanges in IPIV

         if ( KSTEP == 1 ) {
            IPIV[K] = KP;
         } else {
            IPIV[K] = -P;
            IPIV[K+1] = -KP;
         }

         // Increase K and return to the start of the main loop

         K = K + KSTEP;
         GO TO 70;

        //  } // 90

         // Update the lower triangle of A22 (= A(k:n,k:n)) as

         // A22 := A22 - L21*D*L21**T = A22 - L21*W**T

         // computing blocks of NB columns at a time

         for (J = K; NB < 0 ? J >= N : J <= N; J += NB) { // 110
            JB = min( NB, N-J+1 );

            // Update the lower triangle of the diagonal block

            for (JJ = J; JJ <= J + JB - 1; JJ++) { // 100
               dgemv('No transpose', J+JB-JJ, K-1, -ONE, A( JJ, 1 ), LDA, W( JJ, 1 ), LDW, ONE, A( JJ, JJ ), 1 );
            } // 100

            // Update the rectangular subdiagonal block

            if (J+JB <= N) dgemm( 'No transpose', 'Transpose', N-J-JB+1, JB, K-1, -ONE, A( J+JB, 1 ), LDA, W( J, 1 ), LDW, ONE, A( J+JB, J ), LDA );
         } // 110

         // Put L21 in standard form by partially undoing the interchanges
         // in columns 1:k-1

         J = K - 1;
        //  } // 120

            KSTEP = 1;
            JP1 = 1;
            JJ = J;
            JP2 = IPIV( J );
            if ( JP2 < 0 ) {
               JP2 = -JP2;
               J = J - 1;
               JP1 = -IPIV( J );
               KSTEP = 2;
            }

            J = J - 1;
            if (JP2 != JJ && J >= 1) dswap( J, A( JP2, 1 ), LDA, A( JJ, 1 ), LDA );
            JJ = J + 1;
            if (JP1 != JJ && KSTEP == 2) dswap( J, A( JP1, 1 ), LDA, A( JJ, 1 ), LDA );
         IF( J >= 1 ) GO TO 120;

         // Set KB to the number of columns factorized

         KB = K - 1;

      }
      }
