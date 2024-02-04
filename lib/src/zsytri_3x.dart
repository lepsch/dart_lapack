      void zsytri_3x(UPLO, N, A, LDA, E, IPIV, WORK, NB, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, N, NB;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      Complex         A( LDA, * ), E( * ), WORK( N+NB+1, * );
      // ..

// =====================================================================

      // .. Parameters ..
      Complex         CONE, CZERO;
      const              CONE = ( 1.0, 0.0 ), CZERO = ( 0.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                CUT, I, ICOUNT, INVD, IP, K, NNB, J, U11;
      Complex         AK, AKKP1, AKP1, D, T, U01_I_J, U01_IP1_J, U11_I_J, U11_IP1_J;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZGEMM, ZSYSWAPR, ZTRTRI, ZTRMM, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MOD
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0;
      UPPER = lsame( UPLO, 'U' );
      if ( !UPPER && !lsame( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -4;
      }

      // Quick return if possible

      if ( INFO != 0 ) {
         xerbla('ZSYTRI_3X', -INFO );
         return;
      }
      if (N == 0) return;

      // Workspace got Non-diag elements of D

      for (K = 1; K <= N; K++) {
         WORK[K, 1] = E( K );
      }

      // Check that the diagonal matrix D is nonsingular.

      if ( UPPER ) {

         // Upper triangular storage: examine D from bottom to top

         for (INFO = N; INFO >= 1; INFO--) {
            if( IPIV( INFO ) > 0 && A( INFO, INFO ) == CZERO ) return;
         }
      } else {

         // Lower triangular storage: examine D from top to bottom.

         for (INFO = 1; INFO <= N; INFO++) {
            if( IPIV( INFO ) > 0 && A( INFO, INFO ) == CZERO ) return;
         }
      }

      INFO = 0;

      // Splitting Workspace
      // U01 is a block ( N, NB+1 )
      // The first element of U01 is in WORK( 1, 1 )
      // U11 is a block ( NB+1, NB+1 )
      // The first element of U11 is in WORK( N+1, 1 )

      U11 = N;

      // INVD is a block ( N, 2 )
      // The first element of INVD is in WORK( 1, INVD )

      INVD = NB + 2;

      if ( UPPER ) {

         // Begin Upper

         // invA = P * inv(U**T) * inv(D) * inv(U) * P**T.

         ztrtri(UPLO, 'U', N, A, LDA, INFO );

         // inv(D) and inv(D) * inv(U)

         K = 1;
         while (K <= N) {
            if ( IPIV( K ) > 0 ) {
               // 1 x 1 diagonal NNB
               WORK[K, INVD] = CONE /  A( K, K );
               WORK[K, INVD+1] = CZERO;
            } else {
               // 2 x 2 diagonal NNB
               T = WORK( K+1, 1 );
               AK = A( K, K ) / T;
               AKP1 = A( K+1, K+1 ) / T;
               AKKP1 = WORK( K+1, 1 )  / T;
               D = T*( AK*AKP1-CONE );
               WORK[K, INVD] = AKP1 / D;
               WORK[K+1, INVD+1] = AK / D;
               WORK[K, INVD+1] = -AKKP1 / D;
               WORK[K+1, INVD] = WORK( K, INVD+1 );
               K = K + 1;
            }
            K = K + 1;
         }

         // inv(U**T) = (inv(U))**T

         // inv(U**T) * inv(D) * inv(U)

         CUT = N;
         while (CUT > 0) {
            NNB = NB;
            if ( CUT <= NNB ) {
               NNB = CUT;
            } else {
               ICOUNT = 0;
               // count negative elements,
               for (I = CUT+1-NNB; I <= CUT; I++) {
                  if( IPIV( I ) < 0 ) ICOUNT = ICOUNT + 1;
               }
               // need a even number for a clear cut
               if( (ICOUNT % 2) == 1 ) NNB = NNB + 1;
            }

            CUT = CUT - NNB;

            // U01 Block

            for (I = 1; I <= CUT; I++) {
               for (J = 1; J <= NNB; J++) {
                  WORK[I, J] = A( I, CUT+J );
               }
            }

            // U11 Block

            for (I = 1; I <= NNB; I++) {
               WORK[U11+I, I] = CONE;
               for (J = 1; J <= I-1; J++) {
                  WORK[U11+I, J] = CZERO;
                }
                for (J = I+1; J <= NNB; J++) {
                   WORK[U11+I, J] = A( CUT+I, CUT+J );
                }
            }

            // invD * U01

            I = 1;
            while (I <= CUT) {
               if ( IPIV( I ) > 0 ) {
                  for (J = 1; J <= NNB; J++) {
                     WORK[I, J] = WORK( I, INVD ) * WORK( I, J );
                  }
               } else {
                  for (J = 1; J <= NNB; J++) {
                     U01_I_J = WORK( I, J );
                     U01_IP1_J = WORK( I+1, J );
                     WORK[I, J] = WORK( I, INVD ) * U01_I_J + WORK( I, INVD+1 ) * U01_IP1_J                      WORK( I+1, J ) = WORK( I+1, INVD ) * U01_I_J + WORK( I+1, INVD+1 ) * U01_IP1_J;
                  }
                  I = I + 1;
               }
               I = I + 1;
            }

            // invD1 * U11

            I = 1;
            while (I <= NNB) {
               if ( IPIV( CUT+I ) > 0 ) {
                  for (J = I; J <= NNB; J++) {
                     WORK[U11+I, J] = WORK(CUT+I,INVD) * WORK(U11+I,J);
                  }
               } else {
                  for (J = I; J <= NNB; J++) {
                     U11_I_J = WORK(U11+I,J);
                     U11_IP1_J = WORK(U11+I+1,J);
                     WORK[U11+I, J] = WORK(CUT+I,INVD) * WORK(U11+I,J) + WORK(CUT+I,INVD+1) * WORK(U11+I+1,J)                      WORK( U11+I+1, J ) = WORK(CUT+I+1,INVD) * U11_I_J + WORK(CUT+I+1,INVD+1) * U11_IP1_J;
                  }
                  I = I + 1;
               }
               I = I + 1;
            }

            // U11**T * invD1 * U11 -> U11

            ztrmm('L', 'U', 'T', 'U', NNB, NNB, CONE, A( CUT+1, CUT+1 ), LDA, WORK( U11+1, 1 ), N+NB+1 );

            for (I = 1; I <= NNB; I++) {
               for (J = I; J <= NNB; J++) {
                  A[CUT+I, CUT+J] = WORK( U11+I, J );
               }
            }

            // U01**T * invD * U01 -> A( CUT+I, CUT+J )

            zgemm('T', 'N', NNB, NNB, CUT, CONE, A( 1, CUT+1 ), LDA, WORK, N+NB+1, CZERO, WORK(U11+1,1), N+NB+1 );


            // U11 =  U11**T * invD1 * U11 + U01**T * invD * U01

            for (I = 1; I <= NNB; I++) {
               for (J = I; J <= NNB; J++) {
                  A[CUT+I, CUT+J] = A( CUT+I, CUT+J ) + WORK(U11+I,J);
               }
            }

            // U01 =  U00**T * invD0 * U01

            ztrmm('L', UPLO, 'T', 'U', CUT, NNB, CONE, A, LDA, WORK, N+NB+1 );


            // Update U01

            for (I = 1; I <= CUT; I++) {
               for (J = 1; J <= NNB; J++) {
                  A[I, CUT+J] = WORK( I, J );
               }
            }

            // Next Block

         }

         // Apply PERMUTATIONS P and P**T:
         // P * inv(U**T) * inv(D) * inv(U) * P**T.
         // Interchange rows and columns I and IPIV(I) in reverse order
         // from the formation order of IPIV vector for Upper case.

         // ( We can use a loop over IPIV with increment 1,
         // since the ABS value of IPIV(I) represents the row (column)
         // index of the interchange with row (column) i in both 1x1
         // and 2x2 pivot cases, i.e. we don't need separate code branches
         // for 1x1 and 2x2 pivot cases )

         for (I = 1; I <= N; I++) {
             IP = ( IPIV( I ) ).abs();
             if ( IP != I ) {
                if (I < IP) zsyswapr( UPLO, N, A, LDA, I ,IP );
                if (I > IP) zsyswapr( UPLO, N, A, LDA, IP ,I );
             }
         }

      } else {

         // Begin Lower

         // inv A = P * inv(L**T) * inv(D) * inv(L) * P**T.

         ztrtri(UPLO, 'U', N, A, LDA, INFO );

         // inv(D) and inv(D) * inv(L)

         K = N;
         while (K >= 1) {
            if ( IPIV( K ) > 0 ) {
               // 1 x 1 diagonal NNB
               WORK[K, INVD] = CONE /  A( K, K );
               WORK[K, INVD+1] = CZERO;
            } else {
               // 2 x 2 diagonal NNB
               T = WORK( K-1, 1 );
               AK = A( K-1, K-1 ) / T;
               AKP1 = A( K, K ) / T;
               AKKP1 = WORK( K-1, 1 ) / T;
               D = T*( AK*AKP1-CONE );
               WORK[K-1, INVD] = AKP1 / D;
               WORK[K, INVD] = AK / D;
               WORK[K, INVD+1] = -AKKP1 / D;
               WORK[K-1, INVD+1] = WORK( K, INVD+1 );
               K = K - 1;
            }
            K = K - 1;
         }

         // inv(L**T) = (inv(L))**T

         // inv(L**T) * inv(D) * inv(L)

         CUT = 0;
         while (CUT < N) {
            NNB = NB;
            if ( (CUT + NNB) > N ) {
               NNB = N - CUT;
            } else {
               ICOUNT = 0;
               // count negative elements,
               for (I = CUT + 1; I <= CUT+NNB; I++) {
                  if ( IPIV( I ) < 0 ) ICOUNT = ICOUNT + 1;
               }
               // need a even number for a clear cut
               if( (ICOUNT % 2) == 1 ) NNB = NNB + 1;
            }

            // L21 Block

            for (I = 1; I <= N-CUT-NNB; I++) {
               for (J = 1; J <= NNB; J++) {
                 WORK[I, J] = A( CUT+NNB+I, CUT+J );
               }
            }

            // L11 Block

            for (I = 1; I <= NNB; I++) {
               WORK[U11+I, I] = CONE;
               for (J = I+1; J <= NNB; J++) {
                  WORK[U11+I, J] = CZERO;
               }
               for (J = 1; J <= I-1; J++) {
                  WORK[U11+I, J] = A( CUT+I, CUT+J );
               }
            }

            // invD*L21

            I = N-CUT-NNB;
            while (I >= 1) {
               if ( IPIV( CUT+NNB+I ) > 0 ) {
                  for (J = 1; J <= NNB; J++) {
                     WORK[I, J] = WORK( CUT+NNB+I, INVD) * WORK( I, J);
                  }
               } else {
                  for (J = 1; J <= NNB; J++) {
                     U01_I_J = WORK(I,J);
                     U01_IP1_J = WORK(I-1,J);
                     WORK(I,J)=WORK(CUT+NNB+I,INVD)*U01_I_J+ WORK(CUT+NNB+I,INVD+1)*U01_IP1_J                      WORK(I-1,J)=WORK(CUT+NNB+I-1,INVD+1)*U01_I_J+ WORK(CUT+NNB+I-1,INVD)*U01_IP1_J;
                  }
                  I = I - 1;
               }
               I = I - 1;
            }

            // invD1*L11

            I = NNB;
            while (I >= 1) {
               if ( IPIV( CUT+I ) > 0 ) {
                  for (J = 1; J <= NNB; J++) {
                     WORK[U11+I, J] = WORK( CUT+I, INVD)*WORK(U11+I,J);
                  }

               } else {
                  for (J = 1; J <= NNB; J++) {
                     U11_I_J = WORK( U11+I, J );
                     U11_IP1_J = WORK( U11+I-1, J );
                     WORK[U11+I, J] = WORK(CUT+I,INVD) * WORK(U11+I,J) + WORK(CUT+I,INVD+1) * U11_IP1_J                      WORK( U11+I-1, J ) = WORK(CUT+I-1,INVD+1) * U11_I_J + WORK(CUT+I-1,INVD) * U11_IP1_J;
                  }
                  I = I - 1;
               }
               I = I - 1;
            }

            // L11**T * invD1 * L11 -> L11

            ztrmm('L', UPLO, 'T', 'U', NNB, NNB, CONE, A( CUT+1, CUT+1 ), LDA, WORK( U11+1, 1 ), N+NB+1 );


            for (I = 1; I <= NNB; I++) {
               for (J = 1; J <= I; J++) {
                  A[CUT+I, CUT+J] = WORK( U11+I, J );
               }
            }

            if ( (CUT+NNB) < N ) {

               // L21**T * invD2*L21 -> A( CUT+I, CUT+J )

               zgemm('T', 'N', NNB, NNB, N-NNB-CUT, CONE, A( CUT+NNB+1, CUT+1 ), LDA, WORK, N+NB+1, CZERO, WORK( U11+1, 1 ), N+NB+1 );


               // L11 =  L11**T * invD1 * L11 + U01**T * invD * U01

               for (I = 1; I <= NNB; I++) {
                  for (J = 1; J <= I; J++) {
                     A[CUT+I, CUT+J] = A( CUT+I, CUT+J )+WORK(U11+I,J);
                  }
               }

               // L01 =  L22**T * invD2 * L21

               ztrmm('L', UPLO, 'T', 'U', N-NNB-CUT, NNB, CONE, A( CUT+NNB+1, CUT+NNB+1 ), LDA, WORK, N+NB+1 );

               // Update L21

               for (I = 1; I <= N-CUT-NNB; I++) {
                  for (J = 1; J <= NNB; J++) {
                     A[CUT+NNB+I, CUT+J] = WORK( I, J );
                  }
               }

            } else {

               // L11 =  L11**T * invD1 * L11

               for (I = 1; I <= NNB; I++) {
                  for (J = 1; J <= I; J++) {
                     A[CUT+I, CUT+J] = WORK( U11+I, J );
                  }
               }
            }

            // Next Block

            CUT = CUT + NNB;

         }

         // Apply PERMUTATIONS P and P**T:
         // P * inv(L**T) * inv(D) * inv(L) * P**T.
         // Interchange rows and columns I and IPIV(I) in reverse order
         // from the formation order of IPIV vector for Lower case.

         // ( We can use a loop over IPIV with increment -1,
         // since the ABS value of IPIV(I) represents the row (column)
         // index of the interchange with row (column) i in both 1x1
         // and 2x2 pivot cases, i.e. we don't need separate code branches
         // for 1x1 and 2x2 pivot cases )

         for (I = N; I >= 1; I--) {
             IP = ( IPIV( I ) ).abs();
             if ( IP != I ) {
                if (I < IP) zsyswapr( UPLO, N, A, LDA, I ,IP );
                if (I > IP) zsyswapr( UPLO, N, A, LDA, IP ,I );
             }
         }

      }

      return;
      }