      void ssytri2x(final int UPLO, final int N, final Matrix<double> A_, final int LDA, final Array<int> IPIV_, final Array<double> _WORK_, final int NB, final Box<int> INFO,) {
  final A = A_.dim();
  final IPIV = IPIV_.dim();
  final _WORK = _WORK_.dim();

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                INFO, LDA, N, NB;
      int                IPIV( * );
      double               A( LDA, * ), WORK( N+NB+1,* );
      // ..

      double               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      bool               UPPER;
      int                I, IINFO, IP, K, CUT, NNB;
      int                COUNT;
      int                J, U11, INVD;

      double               AK, AKKP1, AKP1, D, T;
      double               U01_I_J, U01_IP1_J;
      double               U11_I_J, U11_IP1_J;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL SSYCONV, XERBLA, STRTRI
      // EXTERNAL SGEMM, STRMM, SSYSWAPR
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX

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
         xerbla('SSYTRI2X', -INFO );
         return;
      }
      if (N == 0) return;

      // Convert A
      // Workspace got Non-diag elements of D

      ssyconv(UPLO, 'C', N, A, LDA, IPIV, WORK, IINFO );

      // Check that the diagonal matrix D is nonsingular.

      if ( UPPER ) {

         // Upper triangular storage: examine D from bottom to top

         for (INFO = N; INFO >= 1; INFO--) {
            if( IPIV( INFO ) > 0 && A( INFO, INFO ) == ZERO ) return;
         }
      } else {

         // Lower triangular storage: examine D from top to bottom.

         for (INFO = 1; INFO <= N; INFO++) {
            if( IPIV( INFO ) > 0 && A( INFO, INFO ) == ZERO ) return;
         }
      }
      INFO = 0;

// Splitting Workspace
      // U01 is a block (N,NB+1)
      // The first element of U01 is in WORK(1,1)
      // U11 is a block (NB+1,NB+1)
      // The first element of U11 is in WORK(N+1,1)
      U11 = N;
      // INVD is a block (N,2)
      // The first element of INVD is in WORK(1,INVD)
      INVD = NB+2;

      if ( UPPER ) {

         // invA = P * inv(U**T)*inv(D)*inv(U)*P**T.

        strtri(UPLO, 'U', N, A, LDA, INFO );

        // inv(D) and inv(D)*inv(U)

        K=1;
        while (K <= N) {
         if ( IPIV( K ) > 0 ) {
            // 1 x 1 diagonal NNB
             WORK[K][INVD] = ONE /  A( K, K );
             WORK[K][INVD+1] = 0;
            K=K+1;
         } else {
            // 2 x 2 diagonal NNB
             T = WORK(K+1,1);
             AK = A( K, K ) / T;
             AKP1 = A( K+1, K+1 ) / T;
             AKKP1 = WORK(K+1,1)  / T;
             D = T*( AK*AKP1-ONE );
             WORK[K][INVD] = AKP1 / D;
             WORK[K+1][INVD+1] = AK / D;
             WORK[K][INVD+1] = -AKKP1 / D;
             WORK[K+1][INVD] = -AKKP1 / D;
            K=K+2;
         }
        }

        // inv(U**T) = (inv(U))**T

        // inv(U**T)*inv(D)*inv(U)

        CUT=N;
        while (CUT > 0) {
           NNB=NB;
           if (CUT <= NNB) {
              NNB=CUT;
           } else {
              COUNT = 0;
              // count negative elements,
              for (I = CUT+1-NNB; I <= CUT; I++) {
                  if (IPIV(I) < 0) COUNT=COUNT+1;
              }
              // need a even number for a clear cut
              if ((COUNT % 2) == 1) NNB=NNB+1;
           }

           CUT=CUT-NNB;

           // U01 Block

           for (I = 1; I <= CUT; I++) {
             for (J = 1; J <= NNB; J++) {
              WORK(I,J)=A(I,CUT+J);
             }
           }

           // U11 Block

           for (I = 1; I <= NNB; I++) {
             WORK(U11+I,I)=ONE;
             for (J = 1; J <= I-1; J++) {
                WORK(U11+I,J)=ZERO;
             }
             for (J = I+1; J <= NNB; J++) {
                WORK(U11+I,J)=A(CUT+I,CUT+J);
             }
           }

           // invD*U01

           I=1;
           while (I <= CUT) {
             if (IPIV(I) > 0) {
                for (J = 1; J <= NNB; J++) {
                    WORK(I,J)=WORK(I,INVD)*WORK(I,J);
                }
                I=I+1;
             } else {
                for (J = 1; J <= NNB; J++) {
                   U01_I_J = WORK(I,J);
                   U01_IP1_J = WORK(I+1,J);
                   WORK(I,J)=WORK(I,INVD)*U01_I_J+ WORK(I,INVD+1)*U01_IP1_J                    WORK(I+1,J)=WORK(I+1,INVD)*U01_I_J+ WORK(I+1,INVD+1)*U01_IP1_J;
                }
                I=I+2;
             }
           }

         // invD1*U11

           I=1;
           while (I <= NNB) {
             if (IPIV(CUT+I) > 0) {
                for (J = I; J <= NNB; J++) {
                    WORK(U11+I,J)=WORK(CUT+I,INVD)*WORK(U11+I,J);
                }
                I=I+1;
             } else {
                for (J = I; J <= NNB; J++) {
                   U11_I_J = WORK(U11+I,J);
                   U11_IP1_J = WORK(U11+I+1,J);
                WORK(U11+I,J)=WORK(CUT+I,INVD)*WORK(U11+I,J) + WORK(CUT+I,INVD+1)*WORK(U11+I+1,J)                 WORK(U11+I+1,J)=WORK(CUT+I+1,INVD)*U11_I_J+ WORK(CUT+I+1,INVD+1)*U11_IP1_J;
                }
                I=I+2;
             }
           }

        // U11**T*invD1*U11->U11

        strmm('L','U','T','U',NNB, NNB, ONE,A(CUT+1,CUT+1),LDA,WORK(U11+1,1),N+NB+1);

         for (I = 1; I <= NNB; I++) {
            for (J = I; J <= NNB; J++) {
              A(CUT+I,CUT+J)=WORK(U11+I,J);
            }
         }

           // U01**T*invD*U01->A(CUT+I,CUT+J)

         sgemm('T','N',NNB,NNB,CUT,ONE,A(1,CUT+1),LDA, WORK,N+NB+1, ZERO, WORK(U11+1,1), N+NB+1);

         // U11 =  U11**T*invD1*U11 + U01**T*invD*U01

         for (I = 1; I <= NNB; I++) {
            for (J = I; J <= NNB; J++) {
              A(CUT+I,CUT+J)=A(CUT+I,CUT+J)+WORK(U11+I,J);
            }
         }

         // U01 =  U00**T*invD0*U01

         strmm('L',UPLO,'T','U',CUT, NNB, ONE,A,LDA,WORK,N+NB+1);


         // Update U01

         for (I = 1; I <= CUT; I++) {
           for (J = 1; J <= NNB; J++) {
            A(I,CUT+J)=WORK(I,J);
           }
         }

       // Next Block

       }

         // Apply PERMUTATIONS P and P**T: P * inv(U**T)*inv(D)*inv(U) *P**T

            I=1;
            while (I <= N) {
               if ( IPIV(I) > 0 ) {
                  IP=IPIV(I);
                 if (I < IP) ssyswapr( UPLO, N, A, LDA, I ,IP );
                 if (I > IP) ssyswapr( UPLO, N, A, LDA, IP ,I );
               } else {
                 IP=-IPIV(I);
                 I=I+1;
                 if ( (I-1) < IP) ssyswapr( UPLO, N, A, LDA, I-1 ,IP );
                 IF ( (I-1) > IP) ssyswapr( UPLO, N, A, LDA, IP ,I-1 );
              }
               I=I+1;
            }
      } else {

         // LOWER...

         // invA = P * inv(U**T)*inv(D)*inv(U)*P**T.

         strtri(UPLO, 'U', N, A, LDA, INFO );

        // inv(D) and inv(D)*inv(U)

        K=N;
        while (K >= 1) {
         if ( IPIV( K ) > 0 ) {
            // 1 x 1 diagonal NNB
             WORK[K][INVD] = ONE /  A( K, K );
             WORK[K][INVD+1] = 0;
            K=K-1;
         } else {
            // 2 x 2 diagonal NNB
             T = WORK(K-1,1);
             AK = A( K-1, K-1 ) / T;
             AKP1 = A( K, K ) / T;
             AKKP1 = WORK(K-1,1) / T;
             D = T*( AK*AKP1-ONE );
             WORK[K-1][INVD] = AKP1 / D;
             WORK[K][INVD] = AK / D;
             WORK[K][INVD+1] = -AKKP1 / D;
             WORK[K-1][INVD+1] = -AKKP1 / D;
            K=K-2;
         }
        }

        // inv(U**T) = (inv(U))**T

        // inv(U**T)*inv(D)*inv(U)

        CUT=0;
        while (CUT < N) {
           NNB=NB;
           if (CUT + NNB > N) {
              NNB=N-CUT;
           } else {
              COUNT = 0;
              // count negative elements,
              for (I = CUT+1; I <= CUT+NNB; I++) {
                  if (IPIV(I) < 0) COUNT=COUNT+1;
              }
              // need a even number for a clear cut
              if ((COUNT % 2) == 1) NNB=NNB+1;
           }
      // L21 Block
           for (I = 1; I <= N-CUT-NNB; I++) {
             for (J = 1; J <= NNB; J++) {
              WORK(I,J)=A(CUT+NNB+I,CUT+J);
             }
           }
      // L11 Block
           for (I = 1; I <= NNB; I++) {
             WORK(U11+I,I)=ONE;
             for (J = I+1; J <= NNB; J++) {
                WORK(U11+I,J)=ZERO;
             }
             for (J = 1; J <= I-1; J++) {
                WORK(U11+I,J)=A(CUT+I,CUT+J);
             }
           }

           // invD*L21

           I=N-CUT-NNB;
           while (I >= 1) {
             if (IPIV(CUT+NNB+I) > 0) {
                for (J = 1; J <= NNB; J++) {
                    WORK(I,J)=WORK(CUT+NNB+I,INVD)*WORK(I,J);
                }
                I=I-1;
             } else {
                for (J = 1; J <= NNB; J++) {
                   U01_I_J = WORK(I,J);
                   U01_IP1_J = WORK(I-1,J);
                   WORK(I,J)=WORK(CUT+NNB+I,INVD)*U01_I_J+ WORK(CUT+NNB+I,INVD+1)*U01_IP1_J                    WORK(I-1,J)=WORK(CUT+NNB+I-1,INVD+1)*U01_I_J+ WORK(CUT+NNB+I-1,INVD)*U01_IP1_J;
                }
                I=I-2;
             }
           }

         // invD1*L11

           I=NNB;
           while (I >= 1) {
             if (IPIV(CUT+I) > 0) {
                for (J = 1; J <= NNB; J++) {
                    WORK(U11+I,J)=WORK(CUT+I,INVD)*WORK(U11+I,J);
                }
                I=I-1;
             } else {
                for (J = 1; J <= NNB; J++) {
                   U11_I_J = WORK(U11+I,J);
                   U11_IP1_J = WORK(U11+I-1,J);
                WORK(U11+I,J)=WORK(CUT+I,INVD)*WORK(U11+I,J) + WORK(CUT+I,INVD+1)*U11_IP1_J                 WORK(U11+I-1,J)=WORK(CUT+I-1,INVD+1)*U11_I_J+ WORK(CUT+I-1,INVD)*U11_IP1_J;
                }
                I=I-2;
             }
           }

        // L11**T*invD1*L11->L11

        strmm('L',UPLO,'T','U',NNB, NNB, ONE,A(CUT+1,CUT+1),LDA,WORK(U11+1,1),N+NB+1);


         for (I = 1; I <= NNB; I++) {
            for (J = 1; J <= I; J++) {
              A(CUT+I,CUT+J)=WORK(U11+I,J);
            }
         }

        if ( (CUT+NNB) < N ) {

           // L21**T*invD2*L21->A(CUT+I,CUT+J)

         sgemm('T','N',NNB,NNB,N-NNB-CUT,ONE,A(CUT+NNB+1,CUT+1) ,LDA,WORK,N+NB+1, ZERO, WORK(U11+1,1), N+NB+1);


         // L11 =  L11**T*invD1*L11 + U01**T*invD*U01

         for (I = 1; I <= NNB; I++) {
            for (J = 1; J <= I; J++) {
              A(CUT+I,CUT+J)=A(CUT+I,CUT+J)+WORK(U11+I,J);
            }
         }

         // L01 =  L22**T*invD2*L21

         strmm('L',UPLO,'T','U', N-NNB-CUT, NNB, ONE,A(CUT+NNB+1,CUT+NNB+1),LDA,WORK,N+NB+1);

       // Update L21

         for (I = 1; I <= N-CUT-NNB; I++) {
           for (J = 1; J <= NNB; J++) {
              A(CUT+NNB+I,CUT+J)=WORK(I,J);
           }
         }

       } else {

         // L11 =  L11**T*invD1*L11

         for (I = 1; I <= NNB; I++) {
            for (J = 1; J <= I; J++) {
              A(CUT+I,CUT+J)=WORK(U11+I,J);
            }
         }
       }

       // Next Block

           CUT=CUT+NNB;
       }

         // Apply PERMUTATIONS P and P**T: P * inv(U**T)*inv(D)*inv(U) *P**T

            I=N;
            while (I >= 1) {
               if ( IPIV(I) > 0 ) {
                  IP=IPIV(I);
                 if (I < IP) ssyswapr( UPLO, N, A, LDA, I ,IP  );
                 if (I > IP) ssyswapr( UPLO, N, A, LDA, IP ,I );
               } else {
                 IP=-IPIV(I);
                 if (I < IP) ssyswapr( UPLO, N, A, LDA, I ,IP );
                 if (I > IP) ssyswapr( UPLO, N, A, LDA, IP ,I );
                 I=I-1;
               }
               I=I-1;
            }
      }

      }
