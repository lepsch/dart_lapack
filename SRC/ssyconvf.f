      SUBROUTINE SSYCONVF( UPLO, WAY, N, A, LDA, E, IPIV, INFO );

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO, WAY;
      int                INFO, LDA, N;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      REAL               A( LDA, * ), E( * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ZERO;
      const              ZERO = 0.0 ;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME

      // .. External Subroutines ..
      // EXTERNAL SSWAP, XERBLA
      // .. Local Scalars ..
      bool               UPPER, CONVERT;
      int                I, IP;
      // ..
      // .. Executable Statements ..

      INFO = 0;
      UPPER = LSAME( UPLO, 'U' );
      CONVERT = LSAME( WAY, 'C' );
      if ( !UPPER && !LSAME( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( !CONVERT && !LSAME( WAY, 'R' ) ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -3;
      } else if ( LDA < MAX( 1, N ) ) {
         INFO = -5;

      }
      if ( INFO != 0 ) {
         xerbla('SSYCONVF', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      if ( UPPER ) {

         // Begin A is UPPER

         if ( CONVERT ) {

            // Convert A (A is upper)


            // Convert VALUE

            // Assign superdiagonal entries of D to array E and zero out
            // corresponding entries in input storage A

            I = N;
            E( 1 ) = ZERO;
            DO WHILE ( I > 1 );
               if ( IPIV( I ) < 0 ) {
                  E( I ) = A( I-1, I );
                  E( I-1 ) = ZERO;
                  A( I-1, I ) = ZERO;
                  I = I - 1;
               } else {
                  E( I ) = ZERO;
               }
               I = I - 1;
            }

            // Convert PERMUTATIONS and IPIV

            // Apply permutations to submatrices of upper part of A
            // in factorization order where i decreases from N to 1

            I = N;
            DO WHILE ( I >= 1 );
               if ( IPIV( I ) > 0 ) {

                  // 1-by-1 pivot interchange

                  // Swap rows i and IPIV(i) in A(1:i,N-i:N)

                  IP = IPIV( I );
                  if ( I < N ) {
                     if ( IP != I ) {
                        sswap(N-I, A( I, I+1 ), LDA, A( IP, I+1 ), LDA );
                     }
                  }

               } else {

                  // 2-by-2 pivot interchange

                  // Swap rows i-1 and IPIV(i) in A(1:i,N-i:N)

                  IP = -IPIV( I );
                  if ( I < N ) {
                     if ( IP != (I-1) ) {
                        sswap(N-I, A( I-1, I+1 ), LDA, A( IP, I+1 ), LDA );
                     }
                  }

                  // Convert IPIV
                  // There is no interchange of rows i and and IPIV(i),
                  // so this should be reflected in IPIV format for
                  // *SYTRF_RK ( or *SYTRF_BK)

                  IPIV( I ) = I;

                  I = I - 1;

               }
               I = I - 1;
            }

         } else {

            // Revert A (A is upper)


            // Revert PERMUTATIONS and IPIV

            // Apply permutations to submatrices of upper part of A
            // in reverse factorization order where i increases from 1 to N

            I = 1;
            DO WHILE ( I <= N );
               if ( IPIV( I ) > 0 ) {

                  // 1-by-1 pivot interchange

                  // Swap rows i and IPIV(i) in A(1:i,N-i:N)

                  IP = IPIV( I );
                  if ( I < N ) {
                     if ( IP != I ) {
                        sswap(N-I, A( IP, I+1 ), LDA, A( I, I+1 ), LDA );
                     }
                  }

               } else {

                  // 2-by-2 pivot interchange

                  // Swap rows i-1 and IPIV(i) in A(1:i,N-i:N)

                  I = I + 1;
                  IP = -IPIV( I );
                  if ( I < N ) {
                     if ( IP != (I-1) ) {
                        sswap(N-I, A( IP, I+1 ), LDA, A( I-1, I+1 ), LDA );
                     }
                  }

                  // Convert IPIV
                  // There is one interchange of rows i-1 and IPIV(i-1),
                  // so this should be recorded in two consecutive entries
                  // in IPIV format for *SYTRF

                  IPIV( I ) = IPIV( I-1 );

               }
               I = I + 1;
            }

            // Revert VALUE
            // Assign superdiagonal entries of D from array E to
            // superdiagonal entries of A.

            I = N;
            DO WHILE ( I > 1 );
               if ( IPIV( I ) < 0 ) {
                  A( I-1, I ) = E( I );
                  I = I - 1;
               }
               I = I - 1;
            }

         // End A is UPPER

         }

      } else {

         // Begin A is LOWER

         if ( CONVERT ) {

            // Convert A (A is lower)


            // Convert VALUE
            // Assign subdiagonal entries of D to array E and zero out
            // corresponding entries in input storage A

            I = 1;
            E( N ) = ZERO;
            DO WHILE ( I <= N );
               if ( I < N && IPIV(I) < 0 ) {
                  E( I ) = A( I+1, I );
                  E( I+1 ) = ZERO;
                  A( I+1, I ) = ZERO;
                  I = I + 1;
               } else {
                  E( I ) = ZERO;
               }
               I = I + 1;
            }

            // Convert PERMUTATIONS and IPIV

            // Apply permutations to submatrices of lower part of A
            // in factorization order where k increases from 1 to N

            I = 1;
            DO WHILE ( I <= N );
               if ( IPIV( I ) > 0 ) {

                  // 1-by-1 pivot interchange

                  // Swap rows i and IPIV(i) in A(i:N,1:i-1)

                  IP = IPIV( I );
                  if ( I > 1 ) {
                     if ( IP != I ) {
                        sswap(I-1, A( I, 1 ), LDA, A( IP, 1 ), LDA );
                     }
                  }

               } else {

                  // 2-by-2 pivot interchange

                  // Swap rows i+1 and IPIV(i) in A(i:N,1:i-1)

                  IP = -IPIV( I );
                  if ( I > 1 ) {
                     if ( IP != (I+1) ) {
                        sswap(I-1, A( I+1, 1 ), LDA, A( IP, 1 ), LDA );
                     }
                  }

                  // Convert IPIV
                  // There is no interchange of rows i and and IPIV(i),
                  // so this should be reflected in IPIV format for
                  // *SYTRF_RK ( or *SYTRF_BK)

                  IPIV( I ) = I;

                  I = I + 1;

               }
               I = I + 1;
            }

         } else {

            // Revert A (A is lower)


            // Revert PERMUTATIONS and IPIV

            // Apply permutations to submatrices of lower part of A
            // in reverse factorization order where i decreases from N to 1

            I = N;
            DO WHILE ( I >= 1 );
               if ( IPIV( I ) > 0 ) {

                  // 1-by-1 pivot interchange

                  // Swap rows i and IPIV(i) in A(i:N,1:i-1)

                  IP = IPIV( I );
                  if ( I > 1 ) {
                     if ( IP != I ) {
                        sswap(I-1, A( IP, 1 ), LDA, A( I, 1 ), LDA );
                     }
                  }

               } else {

                  // 2-by-2 pivot interchange

                  // Swap rows i+1 and IPIV(i) in A(i:N,1:i-1)

                  I = I - 1;
                  IP = -IPIV( I );
                  if ( I > 1 ) {
                     if ( IP != (I+1) ) {
                        sswap(I-1, A( IP, 1 ), LDA, A( I+1, 1 ), LDA );
                     }
                  }

                  // Convert IPIV
                  // There is one interchange of rows i+1 and IPIV(i+1),
                  // so this should be recorded in consecutive entries
                  // in IPIV format for *SYTRF

                  IPIV( I ) = IPIV( I+1 );

               }
               I = I - 1;
            }

            // Revert VALUE
            // Assign subdiagonal entries of D from array E to
            // subdiagonal entries of A.

            I = 1;
            DO WHILE ( I <= N-1 );
               if ( IPIV( I ) < 0 ) {
                  A( I + 1, I ) = E( I );
                  I = I + 1;
               }
               I = I + 1;
            }

         }

         // End A is LOWER

      }

      return;

      // End of SSYCONVF

      }
