      SUBROUTINE CHETRI_3X( UPLO, N, A, LDA, E, IPIV, WORK, NB, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, N, NB;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      COMPLEX            A( LDA, * ), E( * ), WORK( N+NB+1, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE
      const              ONE = 1.0E+0 ;
      COMPLEX            CONE, CZERO
      const              CONE = ( 1.0E+0, 0.0E+0 ), CZERO = ( 0.0E+0, 0.0E+0 ) ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                CUT, I, ICOUNT, INVD, IP, K, NNB, J, U11;
      REAL               AK, AKP1, T
      COMPLEX            AKKP1, D, U01_I_J, U01_IP1_J, U11_I_J, U11_IP1_J
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMM, CHESWAPR, CTRTRI, CTRMM, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, CONJG, MAX, REAL
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

      // Quick return if possible

      if ( INFO.NE.0 ) {
         CALL XERBLA( 'CHETRI_3X', -INFO )
         RETURN
      }
      IF( N.EQ.0 ) RETURN

      // Workspace got Non-diag elements of D

      DO K = 1, N
         WORK( K, 1 ) = E( K )
      END DO

      // Check that the diagonal matrix D is nonsingular.

      if ( UPPER ) {

         // Upper triangular storage: examine D from bottom to top

         DO INFO = N, 1, -1
            IF( IPIV( INFO ).GT.0 .AND. A( INFO, INFO ).EQ.CZERO ) RETURN
         END DO
      } else {

         // Lower triangular storage: examine D from top to bottom.

         DO INFO = 1, N
            IF( IPIV( INFO ).GT.0 .AND. A( INFO, INFO ).EQ.CZERO ) RETURN
         END DO
      }

      INFO = 0

      // Splitting Workspace
      // U01 is a block ( N, NB+1 )
      // The first element of U01 is in WORK( 1, 1 )
      // U11 is a block ( NB+1, NB+1 )
      // The first element of U11 is in WORK( N+1, 1 )

      U11 = N

      // INVD is a block ( N, 2 )
      // The first element of INVD is in WORK( 1, INVD )

      INVD = NB + 2

      if ( UPPER ) {

         // Begin Upper

         // invA = P * inv(U**H) * inv(D) * inv(U) * P**T.

         CALL CTRTRI( UPLO, 'U', N, A, LDA, INFO )

         // inv(D) and inv(D) * inv(U)

         K = 1
         DO WHILE( K.LE.N )
            if ( IPIV( K ).GT.0 ) {
               // 1 x 1 diagonal NNB
               WORK( K, INVD ) = ONE / REAL( A( K, K ) )
               WORK( K, INVD+1 ) = CZERO
            } else {
               // 2 x 2 diagonal NNB
               T = ABS( WORK( K+1, 1 ) )
               AK = REAL( A( K, K ) ) / T
               AKP1 = REAL( A( K+1, K+1 ) ) / T
               AKKP1 = WORK( K+1, 1 )  / T
               D = T*( AK*AKP1-CONE )
               WORK( K, INVD ) = AKP1 / D
               WORK( K+1, INVD+1 ) = AK / D
               WORK( K, INVD+1 ) = -AKKP1 / D
               WORK( K+1, INVD ) = CONJG( WORK( K, INVD+1 ) )
               K = K + 1
            }
            K = K + 1
         END DO

         // inv(U**H) = (inv(U))**H

         // inv(U**H) * inv(D) * inv(U)

         CUT = N
         DO WHILE( CUT.GT.0 )
            NNB = NB
            if ( CUT.LE.NNB ) {
               NNB = CUT
            } else {
               ICOUNT = 0
               // count negative elements,
               DO I = CUT+1-NNB, CUT
                  IF( IPIV( I ).LT.0 ) ICOUNT = ICOUNT + 1
               END DO
               // need a even number for a clear cut
               IF( MOD( ICOUNT, 2 ).EQ.1 ) NNB = NNB + 1
            }

            CUT = CUT - NNB

            // U01 Block

            DO I = 1, CUT
               DO J = 1, NNB
                  WORK( I, J ) = A( I, CUT+J )
               END DO
            END DO

            // U11 Block

            DO I = 1, NNB
               WORK( U11+I, I ) = CONE
               DO J = 1, I-1
                  WORK( U11+I, J ) = CZERO
                END DO
                DO J = I+1, NNB
                   WORK( U11+I, J ) = A( CUT+I, CUT+J )
                END DO
            END DO

            // invD * U01

            I = 1
            DO WHILE( I.LE.CUT )
               if ( IPIV( I ).GT.0 ) {
                  DO J = 1, NNB
                     WORK( I, J ) = WORK( I, INVD ) * WORK( I, J )
                  END DO
               } else {
                  DO J = 1, NNB
                     U01_I_J = WORK( I, J )
                     U01_IP1_J = WORK( I+1, J )
                     WORK( I, J ) = WORK( I, INVD ) * U01_I_J + WORK( I, INVD+1 ) * U01_IP1_J                      WORK( I+1, J ) = WORK( I+1, INVD ) * U01_I_J + WORK( I+1, INVD+1 ) * U01_IP1_J
                  END DO
                  I = I + 1
               }
               I = I + 1
            END DO

            // invD1 * U11

            I = 1
            DO WHILE ( I.LE.NNB )
               if ( IPIV( CUT+I ).GT.0 ) {
                  DO J = I, NNB
                     WORK( U11+I, J ) = WORK(CUT+I,INVD) * WORK(U11+I,J)
                  END DO
               } else {
                  DO J = I, NNB
                     U11_I_J = WORK(U11+I,J)
                     U11_IP1_J = WORK(U11+I+1,J)
                     WORK( U11+I, J ) = WORK(CUT+I,INVD) * WORK(U11+I,J) + WORK(CUT+I,INVD+1) * WORK(U11+I+1,J)                      WORK( U11+I+1, J ) = WORK(CUT+I+1,INVD) * U11_I_J + WORK(CUT+I+1,INVD+1) * U11_IP1_J
                  END DO
                  I = I + 1
               }
               I = I + 1
            END DO

            // U11**H * invD1 * U11 -> U11

            CALL CTRMM( 'L', 'U', 'C', 'U', NNB, NNB, CONE, A( CUT+1, CUT+1 ), LDA, WORK( U11+1, 1 ), N+NB+1 )

            DO I = 1, NNB
               DO J = I, NNB
                  A( CUT+I, CUT+J ) = WORK( U11+I, J )
               END DO
            END DO

            // U01**H * invD * U01 -> A( CUT+I, CUT+J )

            CALL CGEMM( 'C', 'N', NNB, NNB, CUT, CONE, A( 1, CUT+1 ), LDA, WORK, N+NB+1, CZERO, WORK(U11+1,1), N+NB+1 )


            // U11 =  U11**H * invD1 * U11 + U01**H * invD * U01

            DO I = 1, NNB
               DO J = I, NNB
                  A( CUT+I, CUT+J ) = A( CUT+I, CUT+J ) + WORK(U11+I,J)
               END DO
            END DO

            // U01 =  U00**H * invD0 * U01

            CALL CTRMM( 'L', UPLO, 'C', 'U', CUT, NNB, CONE, A, LDA, WORK, N+NB+1 )


            // Update U01

            DO I = 1, CUT
               DO J = 1, NNB
                  A( I, CUT+J ) = WORK( I, J )
               END DO
            END DO

            // Next Block

         END DO

         // Apply PERMUTATIONS P and P**T:
         // P * inv(U**H) * inv(D) * inv(U) * P**T.
         // Interchange rows and columns I and IPIV(I) in reverse order
         // from the formation order of IPIV vector for Upper case.

         // ( We can use a loop over IPIV with increment 1,
         // since the ABS value of IPIV(I) represents the row (column)
         // index of the interchange with row (column) i in both 1x1
         // and 2x2 pivot cases, i.e. we don't need separate code branches
         // for 1x1 and 2x2 pivot cases )

         DO I = 1, N
             IP = ABS( IPIV( I ) )
             if ( IP.NE.I ) {
                IF (I .LT. IP) CALL CHESWAPR( UPLO, N, A, LDA, I ,IP )
                IF (I .GT. IP) CALL CHESWAPR( UPLO, N, A, LDA, IP ,I )
             }
         END DO

      } else {

         // Begin Lower

         // inv A = P * inv(L**H) * inv(D) * inv(L) * P**T.

         CALL CTRTRI( UPLO, 'U', N, A, LDA, INFO )

         // inv(D) and inv(D) * inv(L)

         K = N
         DO WHILE ( K .GE. 1 )
            if ( IPIV( K ).GT.0 ) {
               // 1 x 1 diagonal NNB
               WORK( K, INVD ) = ONE / REAL( A( K, K ) )
               WORK( K, INVD+1 ) = CZERO
            } else {
               // 2 x 2 diagonal NNB
               T = ABS( WORK( K-1, 1 ) )
               AK = REAL( A( K-1, K-1 ) ) / T
               AKP1 = REAL( A( K, K ) ) / T
               AKKP1 = WORK( K-1, 1 ) / T
               D = T*( AK*AKP1-CONE )
               WORK( K-1, INVD ) = AKP1 / D
               WORK( K, INVD ) = AK / D
               WORK( K, INVD+1 ) = -AKKP1 / D
               WORK( K-1, INVD+1 ) = CONJG( WORK( K, INVD+1 ) )
               K = K - 1
            }
            K = K - 1
         END DO

         // inv(L**H) = (inv(L))**H

         // inv(L**H) * inv(D) * inv(L)

         CUT = 0
         DO WHILE( CUT.LT.N )
            NNB = NB
            if ( (CUT + NNB).GT.N ) {
               NNB = N - CUT
            } else {
               ICOUNT = 0
               // count negative elements,
               DO I = CUT + 1, CUT+NNB
                  IF ( IPIV( I ).LT.0 ) ICOUNT = ICOUNT + 1
               END DO
               // need a even number for a clear cut
               IF( MOD( ICOUNT, 2 ).EQ.1 ) NNB = NNB + 1
            }

            // L21 Block

            DO I = 1, N-CUT-NNB
               DO J = 1, NNB
                 WORK( I, J ) = A( CUT+NNB+I, CUT+J )
               END DO
            END DO

            // L11 Block

            DO I = 1, NNB
               WORK( U11+I, I) = CONE
               DO J = I+1, NNB
                  WORK( U11+I, J ) = CZERO
               END DO
               DO J = 1, I-1
                  WORK( U11+I, J ) = A( CUT+I, CUT+J )
               END DO
            END DO

            // invD*L21

            I = N-CUT-NNB
            DO WHILE( I.GE.1 )
               if ( IPIV( CUT+NNB+I ).GT.0 ) {
                  DO J = 1, NNB
                     WORK( I, J ) = WORK( CUT+NNB+I, INVD) * WORK( I, J)
                  END DO
               } else {
                  DO J = 1, NNB
                     U01_I_J = WORK(I,J)
                     U01_IP1_J = WORK(I-1,J)
                     WORK(I,J)=WORK(CUT+NNB+I,INVD)*U01_I_J+ WORK(CUT+NNB+I,INVD+1)*U01_IP1_J                      WORK(I-1,J)=WORK(CUT+NNB+I-1,INVD+1)*U01_I_J+ WORK(CUT+NNB+I-1,INVD)*U01_IP1_J
                  END DO
                  I = I - 1
               }
               I = I - 1
            END DO

            // invD1*L11

            I = NNB
            DO WHILE( I.GE.1 )
               if ( IPIV( CUT+I ).GT.0 ) {
                  DO J = 1, NNB
                     WORK( U11+I, J ) = WORK( CUT+I, INVD)*WORK(U11+I,J)
                  END DO

               } else {
                  DO J = 1, NNB
                     U11_I_J = WORK( U11+I, J )
                     U11_IP1_J = WORK( U11+I-1, J )
                     WORK( U11+I, J ) = WORK(CUT+I,INVD) * WORK(U11+I,J) + WORK(CUT+I,INVD+1) * U11_IP1_J                      WORK( U11+I-1, J ) = WORK(CUT+I-1,INVD+1) * U11_I_J + WORK(CUT+I-1,INVD) * U11_IP1_J
                  END DO
                  I = I - 1
               }
               I = I - 1
            END DO

            // L11**H * invD1 * L11 -> L11

            CALL CTRMM( 'L', UPLO, 'C', 'U', NNB, NNB, CONE, A( CUT+1, CUT+1 ), LDA, WORK( U11+1, 1 ), N+NB+1 )


            DO I = 1, NNB
               DO J = 1, I
                  A( CUT+I, CUT+J ) = WORK( U11+I, J )
               END DO
            END DO

            if ( (CUT+NNB).LT.N ) {

               // L21**H * invD2*L21 -> A( CUT+I, CUT+J )

               CALL CGEMM( 'C', 'N', NNB, NNB, N-NNB-CUT, CONE, A( CUT+NNB+1, CUT+1 ), LDA, WORK, N+NB+1, CZERO, WORK( U11+1, 1 ), N+NB+1 )


               // L11 =  L11**H * invD1 * L11 + U01**H * invD * U01

               DO I = 1, NNB
                  DO J = 1, I
                     A( CUT+I, CUT+J ) = A( CUT+I, CUT+J )+WORK(U11+I,J)
                  END DO
               END DO

               // L01 =  L22**H * invD2 * L21

               CALL CTRMM( 'L', UPLO, 'C', 'U', N-NNB-CUT, NNB, CONE, A( CUT+NNB+1, CUT+NNB+1 ), LDA, WORK, N+NB+1 )

               // Update L21

               DO I = 1, N-CUT-NNB
                  DO J = 1, NNB
                     A( CUT+NNB+I, CUT+J ) = WORK( I, J )
                  END DO
               END DO

            } else {

               // L11 =  L11**H * invD1 * L11

               DO I = 1, NNB
                  DO J = 1, I
                     A( CUT+I, CUT+J ) = WORK( U11+I, J )
                  END DO
               END DO
            }

            // Next Block

            CUT = CUT + NNB

         END DO

         // Apply PERMUTATIONS P and P**T:
         // P * inv(L**H) * inv(D) * inv(L) * P**T.
         // Interchange rows and columns I and IPIV(I) in reverse order
         // from the formation order of IPIV vector for Lower case.

         // ( We can use a loop over IPIV with increment -1,
         // since the ABS value of IPIV(I) represents the row (column)
         // index of the interchange with row (column) i in both 1x1
         // and 2x2 pivot cases, i.e. we don't need separate code branches
         // for 1x1 and 2x2 pivot cases )

         DO I = N, 1, -1
             IP = ABS( IPIV( I ) )
             if ( IP.NE.I ) {
                IF (I .LT. IP) CALL CHESWAPR( UPLO, N, A, LDA, I ,IP )
                IF (I .GT. IP) CALL CHESWAPR( UPLO, N, A, LDA, IP ,I )
             }
         END DO

      }

      RETURN

      // End of CHETRI_3X

      }
