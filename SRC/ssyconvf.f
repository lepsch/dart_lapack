      SUBROUTINE SSYCONVF( UPLO, WAY, N, A, LDA, E, IPIV, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO, WAY;
      int                INFO, LDA, N;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      REAL               A( LDA, * ), E( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E+0 )
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

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      CONVERT = LSAME( WAY, 'C' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.CONVERT .AND. .NOT.LSAME( WAY, 'R' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5

      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SSYCONVF', -INFO )
         RETURN
      END IF

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      IF( UPPER ) THEN

         // Begin A is UPPER

         IF ( CONVERT ) THEN

            // Convert A (A is upper)


            // Convert VALUE

            // Assign superdiagonal entries of D to array E and zero out
            // corresponding entries in input storage A

            I = N
            E( 1 ) = ZERO
            DO WHILE ( I.GT.1 )
               IF( IPIV( I ).LT.0 ) THEN
                  E( I ) = A( I-1, I )
                  E( I-1 ) = ZERO
                  A( I-1, I ) = ZERO
                  I = I - 1
               ELSE
                  E( I ) = ZERO
               END IF
               I = I - 1
            END DO

            // Convert PERMUTATIONS and IPIV

            // Apply permutations to submatrices of upper part of A
            // in factorization order where i decreases from N to 1

            I = N
            DO WHILE ( I.GE.1 )
               IF( IPIV( I ).GT.0 ) THEN

                  // 1-by-1 pivot interchange

                  // Swap rows i and IPIV(i) in A(1:i,N-i:N)

                  IP = IPIV( I )
                  IF( I.LT.N ) THEN
                     IF( IP.NE.I ) THEN
                        CALL SSWAP( N-I, A( I, I+1 ), LDA, A( IP, I+1 ), LDA )
                     END IF
                  END IF

               ELSE

                  // 2-by-2 pivot interchange

                  // Swap rows i-1 and IPIV(i) in A(1:i,N-i:N)

                  IP = -IPIV( I )
                  IF( I.LT.N ) THEN
                     IF( IP.NE.(I-1) ) THEN
                        CALL SSWAP( N-I, A( I-1, I+1 ), LDA, A( IP, I+1 ), LDA )
                     END IF
                  END IF

                  // Convert IPIV
                  // There is no interchange of rows i and and IPIV(i),
                  // so this should be reflected in IPIV format for
                  // *SYTRF_RK ( or *SYTRF_BK)

                  IPIV( I ) = I

                  I = I - 1

               END IF
               I = I - 1
            END DO

         ELSE

            // Revert A (A is upper)


            // Revert PERMUTATIONS and IPIV

            // Apply permutations to submatrices of upper part of A
            // in reverse factorization order where i increases from 1 to N

            I = 1
            DO WHILE ( I.LE.N )
               IF( IPIV( I ).GT.0 ) THEN

                  // 1-by-1 pivot interchange

                  // Swap rows i and IPIV(i) in A(1:i,N-i:N)

                  IP = IPIV( I )
                  IF( I.LT.N ) THEN
                     IF( IP.NE.I ) THEN
                        CALL SSWAP( N-I, A( IP, I+1 ), LDA, A( I, I+1 ), LDA )
                     END IF
                  END IF

               ELSE

                  // 2-by-2 pivot interchange

                  // Swap rows i-1 and IPIV(i) in A(1:i,N-i:N)

                  I = I + 1
                  IP = -IPIV( I )
                  IF( I.LT.N ) THEN
                     IF( IP.NE.(I-1) ) THEN
                        CALL SSWAP( N-I, A( IP, I+1 ), LDA, A( I-1, I+1 ), LDA )
                     END IF
                  END IF

                  // Convert IPIV
                  // There is one interchange of rows i-1 and IPIV(i-1),
                  // so this should be recorded in two consecutive entries
                  // in IPIV format for *SYTRF

                  IPIV( I ) = IPIV( I-1 )

               END IF
               I = I + 1
            END DO

            // Revert VALUE
            // Assign superdiagonal entries of D from array E to
            // superdiagonal entries of A.

            I = N
            DO WHILE ( I.GT.1 )
               IF( IPIV( I ).LT.0 ) THEN
                  A( I-1, I ) = E( I )
                  I = I - 1
               END IF
               I = I - 1
            END DO

         // End A is UPPER

         END IF

      ELSE

         // Begin A is LOWER

         IF ( CONVERT ) THEN

            // Convert A (A is lower)


            // Convert VALUE
            // Assign subdiagonal entries of D to array E and zero out
            // corresponding entries in input storage A

            I = 1
            E( N ) = ZERO
            DO WHILE ( I.LE.N )
               IF( I.LT.N .AND. IPIV(I).LT.0 ) THEN
                  E( I ) = A( I+1, I )
                  E( I+1 ) = ZERO
                  A( I+1, I ) = ZERO
                  I = I + 1
               ELSE
                  E( I ) = ZERO
               END IF
               I = I + 1
            END DO

            // Convert PERMUTATIONS and IPIV

            // Apply permutations to submatrices of lower part of A
            // in factorization order where k increases from 1 to N

            I = 1
            DO WHILE ( I.LE.N )
               IF( IPIV( I ).GT.0 ) THEN

                  // 1-by-1 pivot interchange

                  // Swap rows i and IPIV(i) in A(i:N,1:i-1)

                  IP = IPIV( I )
                  IF ( I.GT.1 ) THEN
                     IF( IP.NE.I ) THEN
                        CALL SSWAP( I-1, A( I, 1 ), LDA, A( IP, 1 ), LDA )
                     END IF
                  END IF

               ELSE

                  // 2-by-2 pivot interchange

                  // Swap rows i+1 and IPIV(i) in A(i:N,1:i-1)

                  IP = -IPIV( I )
                  IF ( I.GT.1 ) THEN
                     IF( IP.NE.(I+1) ) THEN
                        CALL SSWAP( I-1, A( I+1, 1 ), LDA, A( IP, 1 ), LDA )
                     END IF
                  END IF

                  // Convert IPIV
                  // There is no interchange of rows i and and IPIV(i),
                  // so this should be reflected in IPIV format for
                  // *SYTRF_RK ( or *SYTRF_BK)

                  IPIV( I ) = I

                  I = I + 1

               END IF
               I = I + 1
            END DO

         ELSE

            // Revert A (A is lower)


            // Revert PERMUTATIONS and IPIV

            // Apply permutations to submatrices of lower part of A
            // in reverse factorization order where i decreases from N to 1

            I = N
            DO WHILE ( I.GE.1 )
               IF( IPIV( I ).GT.0 ) THEN

                  // 1-by-1 pivot interchange

                  // Swap rows i and IPIV(i) in A(i:N,1:i-1)

                  IP = IPIV( I )
                  IF ( I.GT.1 ) THEN
                     IF( IP.NE.I ) THEN
                        CALL SSWAP( I-1, A( IP, 1 ), LDA, A( I, 1 ), LDA )
                     END IF
                  END IF

               ELSE

                  // 2-by-2 pivot interchange

                  // Swap rows i+1 and IPIV(i) in A(i:N,1:i-1)

                  I = I - 1
                  IP = -IPIV( I )
                  IF ( I.GT.1 ) THEN
                     IF( IP.NE.(I+1) ) THEN
                        CALL SSWAP( I-1, A( IP, 1 ), LDA, A( I+1, 1 ), LDA )
                     END IF
                  END IF

                  // Convert IPIV
                  // There is one interchange of rows i+1 and IPIV(i+1),
                  // so this should be recorded in consecutive entries
                  // in IPIV format for *SYTRF

                  IPIV( I ) = IPIV( I+1 )

               END IF
               I = I - 1
            END DO

            // Revert VALUE
            // Assign subdiagonal entries of D from array E to
            // subdiagonal entries of A.

            I = 1
            DO WHILE ( I.LE.N-1 )
               IF( IPIV( I ).LT.0 ) THEN
                  A( I + 1, I ) = E( I )
                  I = I + 1
               END IF
               I = I + 1
            END DO

         END IF

         // End A is LOWER

      END IF

      RETURN

      // End of SSYCONVF

      END
