      SUBROUTINE CHETRI2X( UPLO, N, A, LDA, IPIV, WORK, NB, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, N, NB;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      COMPLEX            A( LDA, * ), WORK( N+NB+1,* )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE
      COMPLEX            CONE, ZERO
      const              ONE = 1.0E+0, CONE = ( 1.0E+0, 0.0E+0 ), ZERO = ( 0.0E+0, 0.0E+0 ) ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                I, IINFO, IP, K, CUT, NNB;
      int                COUNT;
      int                J, U11, INVD;

      COMPLEX   AK, AKKP1, AKP1, D, T
      COMPLEX   U01_I_J, U01_IP1_J
      COMPLEX   U11_I_J, U11_IP1_J
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL CSYCONV, XERBLA, CTRTRI
      // EXTERNAL CGEMM, CTRMM, CHESWAPR
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
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
         xerbla('CHETRI2X', -INFO );
         RETURN
      }
      IF( N.EQ.0 ) RETURN

      // Convert A
      // Workspace got Non-diag elements of D

      csyconv(UPLO, 'C', N, A, LDA, IPIV, WORK, IINFO );

      // Check that the diagonal matrix D is nonsingular.

      if ( UPPER ) {

         // Upper triangular storage: examine D from bottom to top

         DO INFO = N, 1, -1
            IF( IPIV( INFO ).GT.0 .AND. A( INFO, INFO ).EQ.ZERO ) RETURN
         }
      } else {

         // Lower triangular storage: examine D from top to bottom.

         for (INFO = 1; INFO <= N; INFO++) {
            IF( IPIV( INFO ).GT.0 .AND. A( INFO, INFO ).EQ.ZERO ) RETURN
         }
      }
      INFO = 0

*  Splitting Workspace
      // U01 is a block (N,NB+1)
      // The first element of U01 is in WORK(1,1)
      // U11 is a block (NB+1,NB+1)
      // The first element of U11 is in WORK(N+1,1)
      U11 = N
      // INVD is a block (N,2)
      // The first element of INVD is in WORK(1,INVD)
      INVD = NB+2

      if ( UPPER ) {

         // invA = P * inv(U**H)*inv(D)*inv(U)*P**H.

        ctrtri(UPLO, 'U', N, A, LDA, INFO );

        // inv(D) and inv(D)*inv(U)

        K=1
        DO WHILE ( K .LE. N )
         if ( IPIV( K ).GT.0 ) {
            // 1 x 1 diagonal NNB
             WORK(K,INVD) = ONE / REAL ( A( K, K ) )
             WORK(K,INVD+1) = 0
            K=K+1
         } else {
            // 2 x 2 diagonal NNB
             T = ABS ( WORK(K+1,1) )
             AK = REAL ( A( K, K ) ) / T
             AKP1 = REAL ( A( K+1, K+1 ) ) / T
             AKKP1 = WORK(K+1,1)  / T
             D = T*( AK*AKP1-ONE )
             WORK(K,INVD) = AKP1 / D
             WORK(K+1,INVD+1) = AK / D
             WORK(K,INVD+1) = -AKKP1 / D
             WORK(K+1,INVD) = CONJG (WORK(K,INVD+1) )
            K=K+2
         }
        }

        // inv(U**H) = (inv(U))**H

        // inv(U**H)*inv(D)*inv(U)

        CUT=N
        DO WHILE (CUT .GT. 0)
           NNB=NB
           if (CUT .LE. NNB) {
              NNB=CUT
           } else {
              COUNT = 0
              // count negative elements,
              for (I = CUT+1-NNB; I <= CUT; I++) {
                  IF (IPIV(I) .LT. 0) COUNT=COUNT+1
              }
              // need a even number for a clear cut
              IF (MOD(COUNT,2) .EQ. 1) NNB=NNB+1
           }

           CUT=CUT-NNB

           // U01 Block

           for (I = 1; I <= CUT; I++) {
             for (J = 1; J <= NNB; J++) {
              WORK(I,J)=A(I,CUT+J)
             }
           }

           // U11 Block

           for (I = 1; I <= NNB; I++) {
             WORK(U11+I,I)=CONE
             for (J = 1; J <= I-1; J++) {
                WORK(U11+I,J)=ZERO
             }
             for (J = I+1; J <= NNB; J++) {
                WORK(U11+I,J)=A(CUT+I,CUT+J)
             }
           }

           // invD*U01

           I=1
           DO WHILE (I .LE. CUT)
             if (IPIV(I) > 0) {
                for (J = 1; J <= NNB; J++) {
                    WORK(I,J)=WORK(I,INVD)*WORK(I,J)
                }
                I=I+1
             } else {
                for (J = 1; J <= NNB; J++) {
                   U01_I_J = WORK(I,J)
                   U01_IP1_J = WORK(I+1,J)
                   WORK(I,J)=WORK(I,INVD)*U01_I_J+ WORK(I,INVD+1)*U01_IP1_J                    WORK(I+1,J)=WORK(I+1,INVD)*U01_I_J+ WORK(I+1,INVD+1)*U01_IP1_J
                }
                I=I+2
             }
           }

         // invD1*U11

           I=1
           DO WHILE (I .LE. NNB)
             if (IPIV(CUT+I) > 0) {
                for (J = I; J <= NNB; J++) {
                    WORK(U11+I,J)=WORK(CUT+I,INVD)*WORK(U11+I,J)
                }
                I=I+1
             } else {
                for (J = I; J <= NNB; J++) {
                   U11_I_J = WORK(U11+I,J)
                   U11_IP1_J = WORK(U11+I+1,J)
                WORK(U11+I,J)=WORK(CUT+I,INVD)*WORK(U11+I,J) + WORK(CUT+I,INVD+1)*WORK(U11+I+1,J)                 WORK(U11+I+1,J)=WORK(CUT+I+1,INVD)*U11_I_J+ WORK(CUT+I+1,INVD+1)*U11_IP1_J
                }
                I=I+2
             }
           }

        // U11**H*invD1*U11->U11

        ctrmm('L','U','C','U',NNB, NNB, CONE,A(CUT+1,CUT+1),LDA,WORK(U11+1,1),N+NB+1);

         for (I = 1; I <= NNB; I++) {
            for (J = I; J <= NNB; J++) {
              A(CUT+I,CUT+J)=WORK(U11+I,J)
            }
         }

           // U01**H*invD*U01->A(CUT+I,CUT+J)

         cgemm('C','N',NNB,NNB,CUT,CONE,A(1,CUT+1),LDA, WORK,N+NB+1, ZERO, WORK(U11+1,1), N+NB+1);

         // U11 =  U11**H*invD1*U11 + U01**H*invD*U01

         for (I = 1; I <= NNB; I++) {
            for (J = I; J <= NNB; J++) {
              A(CUT+I,CUT+J)=A(CUT+I,CUT+J)+WORK(U11+I,J)
            }
         }

         // U01 =  U00**H*invD0*U01

         ctrmm('L',UPLO,'C','U',CUT, NNB, CONE,A,LDA,WORK,N+NB+1);


         // Update U01

         for (I = 1; I <= CUT; I++) {
           for (J = 1; J <= NNB; J++) {
            A(I,CUT+J)=WORK(I,J)
           }
         }

       // Next Block

       }

         // Apply PERMUTATIONS P and P**H: P * inv(U**H)*inv(D)*inv(U) *P**H

            I=1
            DO WHILE ( I .LE. N )
               if ( IPIV(I) .GT. 0 ) {
                  IP=IPIV(I)
                 IF (I .LT. IP) CALL CHESWAPR( UPLO, N, A, LDA, I ,IP )
                 IF (I .GT. IP) CALL CHESWAPR( UPLO, N, A, LDA, IP ,I )
               } else {
                 IP=-IPIV(I)
                 I=I+1
                 IF ( (I-1) .LT. IP) CALL CHESWAPR( UPLO, N, A, LDA, I-1 ,IP )                  IF ( (I-1) .GT. IP) CALL CHESWAPR( UPLO, N, A, LDA, IP ,I-1 )
              }
               I=I+1
            }
      } else {

         // LOWER...

         // invA = P * inv(U**H)*inv(D)*inv(U)*P**H.

         ctrtri(UPLO, 'U', N, A, LDA, INFO );

        // inv(D) and inv(D)*inv(U)

        K=N
        DO WHILE ( K .GE. 1 )
         if ( IPIV( K ).GT.0 ) {
            // 1 x 1 diagonal NNB
             WORK(K,INVD) = ONE / REAL ( A( K, K ) )
             WORK(K,INVD+1) = 0
            K=K-1
         } else {
            // 2 x 2 diagonal NNB
             T = ABS ( WORK(K-1,1) )
             AK = REAL ( A( K-1, K-1 ) ) / T
             AKP1 = REAL ( A( K, K ) ) / T
             AKKP1 = WORK(K-1,1) / T
             D = T*( AK*AKP1-ONE )
             WORK(K-1,INVD) = AKP1 / D
             WORK(K,INVD) = AK / D
             WORK(K,INVD+1) = -AKKP1 / D
             WORK(K-1,INVD+1) = CONJG (WORK(K,INVD+1) )
            K=K-2
         }
        }

        // inv(U**H) = (inv(U))**H

        // inv(U**H)*inv(D)*inv(U)

        CUT=0
        DO WHILE (CUT .LT. N)
           NNB=NB
           if (CUT + NNB .GE. N) {
              NNB=N-CUT
           } else {
              COUNT = 0
              // count negative elements,
              for (I = CUT+1; I <= CUT+NNB; I++) {
                  IF (IPIV(I) .LT. 0) COUNT=COUNT+1
              }
              // need a even number for a clear cut
              IF (MOD(COUNT,2) .EQ. 1) NNB=NNB+1
           }
       // L21 Block
           for (I = 1; I <= N-CUT-NNB; I++) {
             for (J = 1; J <= NNB; J++) {
              WORK(I,J)=A(CUT+NNB+I,CUT+J)
             }
           }
      // L11 Block
           for (I = 1; I <= NNB; I++) {
             WORK(U11+I,I)=CONE
             for (J = I+1; J <= NNB; J++) {
                WORK(U11+I,J)=ZERO
             }
             for (J = 1; J <= I-1; J++) {
                WORK(U11+I,J)=A(CUT+I,CUT+J)
             }
           }

           // invD*L21

           I=N-CUT-NNB
           DO WHILE (I .GE. 1)
             if (IPIV(CUT+NNB+I) > 0) {
                for (J = 1; J <= NNB; J++) {
                    WORK(I,J)=WORK(CUT+NNB+I,INVD)*WORK(I,J)
                }
                I=I-1
             } else {
                for (J = 1; J <= NNB; J++) {
                   U01_I_J = WORK(I,J)
                   U01_IP1_J = WORK(I-1,J)
                   WORK(I,J)=WORK(CUT+NNB+I,INVD)*U01_I_J+ WORK(CUT+NNB+I,INVD+1)*U01_IP1_J                    WORK(I-1,J)=WORK(CUT+NNB+I-1,INVD+1)*U01_I_J+ WORK(CUT+NNB+I-1,INVD)*U01_IP1_J
                }
                I=I-2
             }
           }

         // invD1*L11

           I=NNB
           DO WHILE (I .GE. 1)
             if (IPIV(CUT+I) > 0) {
                for (J = 1; J <= NNB; J++) {
                    WORK(U11+I,J)=WORK(CUT+I,INVD)*WORK(U11+I,J)
                }
                I=I-1
             } else {
                for (J = 1; J <= NNB; J++) {
                   U11_I_J = WORK(U11+I,J)
                   U11_IP1_J = WORK(U11+I-1,J)
                WORK(U11+I,J)=WORK(CUT+I,INVD)*WORK(U11+I,J) + WORK(CUT+I,INVD+1)*U11_IP1_J                 WORK(U11+I-1,J)=WORK(CUT+I-1,INVD+1)*U11_I_J+ WORK(CUT+I-1,INVD)*U11_IP1_J
                }
                I=I-2
             }
           }

        // L11**H*invD1*L11->L11

        ctrmm('L',UPLO,'C','U',NNB, NNB, CONE,A(CUT+1,CUT+1),LDA,WORK(U11+1,1),N+NB+1);

         for (I = 1; I <= NNB; I++) {
            for (J = 1; J <= I; J++) {
              A(CUT+I,CUT+J)=WORK(U11+I,J)
            }
         }

        if ( (CUT+NNB) .LT. N ) {

           // L21**H*invD2*L21->A(CUT+I,CUT+J)

         cgemm('C','N',NNB,NNB,N-NNB-CUT,CONE,A(CUT+NNB+1,CUT+1) ,LDA,WORK,N+NB+1, ZERO, WORK(U11+1,1), N+NB+1);


         // L11 =  L11**H*invD1*L11 + U01**H*invD*U01

         for (I = 1; I <= NNB; I++) {
            for (J = 1; J <= I; J++) {
              A(CUT+I,CUT+J)=A(CUT+I,CUT+J)+WORK(U11+I,J)
            }
         }

         // L01 =  L22**H*invD2*L21

         ctrmm('L',UPLO,'C','U', N-NNB-CUT, NNB, CONE,A(CUT+NNB+1,CUT+NNB+1),LDA,WORK,N+NB+1);

       // Update L21
         for (I = 1; I <= N-CUT-NNB; I++) {
           for (J = 1; J <= NNB; J++) {
              A(CUT+NNB+I,CUT+J)=WORK(I,J)
           }
         }
       } else {

         // L11 =  L11**H*invD1*L11

         for (I = 1; I <= NNB; I++) {
            for (J = 1; J <= I; J++) {
              A(CUT+I,CUT+J)=WORK(U11+I,J)
            }
         }
       }

       // Next Block

           CUT=CUT+NNB
       }

         // Apply PERMUTATIONS P and P**H: P * inv(U**H)*inv(D)*inv(U) *P**H

            I=N
            DO WHILE ( I .GE. 1 )
               if ( IPIV(I) .GT. 0 ) {
                  IP=IPIV(I)
                 IF (I .LT. IP) CALL CHESWAPR( UPLO, N, A, LDA, I ,IP  )
                 IF (I .GT. IP) CALL CHESWAPR( UPLO, N, A, LDA, IP ,I )
               } else {
                 IP=-IPIV(I)
                 IF ( I .LT. IP) CALL CHESWAPR( UPLO, N, A, LDA, I ,IP )
                 IF ( I .GT. IP) CALL CHESWAPR( UPLO, N, A, LDA, IP ,I )
                 I=I-1
               }
               I=I-1
            }
      }

      RETURN

      // End of CHETRI2X

      }
