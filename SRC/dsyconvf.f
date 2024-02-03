      SUBROUTINE DSYCONVF( UPLO, WAY, N, A, LDA, E, IPIV, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO, WAY;
      int                INFO, LDA, N;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      double             A( LDA, * ), E( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO;
      const              ZERO = 0.0D+0 ;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME

      // .. External Subroutines ..
      // EXTERNAL DSWAP, XERBLA
      // .. Local Scalars ..
      bool               UPPER, CONVERT;
      int                I, IP;
      // ..
      // .. Executable Statements ..

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      CONVERT = LSAME( WAY, 'C' )
      if ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( .NOT.CONVERT .AND. .NOT.LSAME( WAY, 'R' ) ) {
         INFO = -2
      } else if ( N.LT.0 ) {
         INFO = -3
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -5

      }
      if ( INFO.NE.0 ) {
         CALL XERBLA( 'DSYCONVF', -INFO )
         RETURN
      }

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      if ( UPPER ) {

         // Begin A is UPPER

         if ( CONVERT ) {

            // Convert A (A is upper)


            // Convert VALUE

            // Assign superdiagonal entries of D to array E and zero out
            // corresponding entries in input storage A

            I = N
            E( 1 ) = ZERO
            DO WHILE ( I.GT.1 )
               if ( IPIV( I ).LT.0 ) {
                  E( I ) = A( I-1, I )
                  E( I-1 ) = ZERO
                  A( I-1, I ) = ZERO
                  I = I - 1
               } else {
                  E( I ) = ZERO
               }
               I = I - 1
            END DO

            // Convert PERMUTATIONS and IPIV

            // Apply permutations to submatrices of upper part of A
            // in factorization order where i decreases from N to 1

            I = N
            DO WHILE ( I.GE.1 )
               if ( IPIV( I ).GT.0 ) {

                  // 1-by-1 pivot interchange

                  // Swap rows i and IPIV(i) in A(1:i,N-i:N)

                  IP = IPIV( I )
                  if ( I.LT.N ) {
                     if ( IP.NE.I ) {
                        CALL DSWAP( N-I, A( I, I+1 ), LDA, A( IP, I+1 ), LDA )
                     }
                  }

               } else {

                  // 2-by-2 pivot interchange

                  // Swap rows i-1 and IPIV(i) in A(1:i,N-i:N)

                  IP = -IPIV( I )
                  if ( I.LT.N ) {
                     if ( IP.NE.(I-1) ) {
                        CALL DSWAP( N-I, A( I-1, I+1 ), LDA, A( IP, I+1 ), LDA )
                     }
                  }

                  // Convert IPIV
                  // There is no interchange of rows i and and IPIV(i),
                  // so this should be reflected in IPIV format for
                  // *SYTRF_RK ( or *SYTRF_BK)

                  IPIV( I ) = I

                  I = I - 1

               }
               I = I - 1
            END DO

         } else {

            // Revert A (A is upper)


            // Revert PERMUTATIONS and IPIV

            // Apply permutations to submatrices of upper part of A
            // in reverse factorization order where i increases from 1 to N

            I = 1
            DO WHILE ( I.LE.N )
               if ( IPIV( I ).GT.0 ) {

                  // 1-by-1 pivot interchange

                  // Swap rows i and IPIV(i) in A(1:i,N-i:N)

                  IP = IPIV( I )
                  if ( I.LT.N ) {
                     if ( IP.NE.I ) {
                        CALL DSWAP( N-I, A( IP, I+1 ), LDA, A( I, I+1 ), LDA )
                     }
                  }

               } else {

                  // 2-by-2 pivot interchange

                  // Swap rows i-1 and IPIV(i) in A(1:i,N-i:N)

                  I = I + 1
                  IP = -IPIV( I )
                  if ( I.LT.N ) {
                     if ( IP.NE.(I-1) ) {
                        CALL DSWAP( N-I, A( IP, I+1 ), LDA, A( I-1, I+1 ), LDA )
                     }
                  }

                  // Convert IPIV
                  // There is one interchange of rows i-1 and IPIV(i-1),
                  // so this should be recorded in two consecutive entries
                  // in IPIV format for *SYTRF

                  IPIV( I ) = IPIV( I-1 )

               }
               I = I + 1
            END DO

            // Revert VALUE
            // Assign superdiagonal entries of D from array E to
            // superdiagonal entries of A.

            I = N
            DO WHILE ( I.GT.1 )
               if ( IPIV( I ).LT.0 ) {
                  A( I-1, I ) = E( I )
                  I = I - 1
               }
               I = I - 1
            END DO

         // End A is UPPER

         }

      } else {

         // Begin A is LOWER

         if ( CONVERT ) {

            // Convert A (A is lower)


            // Convert VALUE
            // Assign subdiagonal entries of D to array E and zero out
            // corresponding entries in input storage A

            I = 1
            E( N ) = ZERO
            DO WHILE ( I.LE.N )
               if ( I.LT.N .AND. IPIV(I).LT.0 ) {
                  E( I ) = A( I+1, I )
                  E( I+1 ) = ZERO
                  A( I+1, I ) = ZERO
                  I = I + 1
               } else {
                  E( I ) = ZERO
               }
               I = I + 1
            END DO

            // Convert PERMUTATIONS and IPIV

            // Apply permutations to submatrices of lower part of A
            // in factorization order where k increases from 1 to N

            I = 1
            DO WHILE ( I.LE.N )
               if ( IPIV( I ).GT.0 ) {

                  // 1-by-1 pivot interchange

                  // Swap rows i and IPIV(i) in A(i:N,1:i-1)

                  IP = IPIV( I )
                  if ( I.GT.1 ) {
                     if ( IP.NE.I ) {
                        CALL DSWAP( I-1, A( I, 1 ), LDA, A( IP, 1 ), LDA )
                     }
                  }

               } else {

                  // 2-by-2 pivot interchange

                  // Swap rows i+1 and IPIV(i) in A(i:N,1:i-1)

                  IP = -IPIV( I )
                  if ( I.GT.1 ) {
                     if ( IP.NE.(I+1) ) {
                        CALL DSWAP( I-1, A( I+1, 1 ), LDA, A( IP, 1 ), LDA )
                     }
                  }

                  // Convert IPIV
                  // There is no interchange of rows i and and IPIV(i),
                  // so this should be reflected in IPIV format for
                  // *SYTRF_RK ( or *SYTRF_BK)

                  IPIV( I ) = I

                  I = I + 1

               }
               I = I + 1
            END DO

         } else {

            // Revert A (A is lower)


            // Revert PERMUTATIONS and IPIV

            // Apply permutations to submatrices of lower part of A
            // in reverse factorization order where i decreases from N to 1

            I = N
            DO WHILE ( I.GE.1 )
               if ( IPIV( I ).GT.0 ) {

                  // 1-by-1 pivot interchange

                  // Swap rows i and IPIV(i) in A(i:N,1:i-1)

                  IP = IPIV( I )
                  if ( I.GT.1 ) {
                     if ( IP.NE.I ) {
                        CALL DSWAP( I-1, A( IP, 1 ), LDA, A( I, 1 ), LDA )
                     }
                  }

               } else {

                  // 2-by-2 pivot interchange

                  // Swap rows i+1 and IPIV(i) in A(i:N,1:i-1)

                  I = I - 1
                  IP = -IPIV( I )
                  if ( I.GT.1 ) {
                     if ( IP.NE.(I+1) ) {
                        CALL DSWAP( I-1, A( IP, 1 ), LDA, A( I+1, 1 ), LDA )
                     }
                  }

                  // Convert IPIV
                  // There is one interchange of rows i+1 and IPIV(i+1),
                  // so this should be recorded in consecutive entries
                  // in IPIV format for *SYTRF

                  IPIV( I ) = IPIV( I+1 )

               }
               I = I - 1
            END DO

            // Revert VALUE
            // Assign subdiagonal entries of D from array E to
            // subdiagonal entries of A.

            I = 1
            DO WHILE ( I.LE.N-1 )
               if ( IPIV( I ).LT.0 ) {
                  A( I + 1, I ) = E( I )
                  I = I + 1
               }
               I = I + 1
            END DO

         }

         // End A is LOWER

      }

      RETURN

      // End of DSYCONVF

      }
