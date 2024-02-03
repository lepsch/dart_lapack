      SUBROUTINE DSYTRI2X( UPLO, N, A, LDA, IPIV, WORK, NB, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, N, NB;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      double             A( LDA, * ), WORK( N+NB+1,* );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0D+0, ZERO = 0.0D+0 ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                I, IINFO, IP, K, CUT, NNB;
      int                COUNT;
      int                J, U11, INVD;

      double             AK, AKKP1, AKP1, D, T;
      double             U01_I_J, U01_IP1_J;
      double             U11_I_J, U11_IP1_J;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL DSYCONV, XERBLA, DTRTRI
      // EXTERNAL DGEMM, DTRMM, DSYSWAPR
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
         xerbla('DSYTRI2X', -INFO );
         RETURN
      }
      IF( N.EQ.0 ) RETURN

      // Convert A
      // Workspace got Non-diag elements of D

      dsyconv(UPLO, 'C', N, A, LDA, IPIV, WORK, IINFO );

      // Check that the diagonal matrix D is nonsingular.

      if ( UPPER ) {

         // Upper triangular storage: examine D from bottom to top

         DO INFO = N, 1, -1
            IF( IPIV( INFO ).GT.0 .AND. A( INFO, INFO ).EQ.ZERO ) RETURN
         END DO
      } else {

         // Lower triangular storage: examine D from top to bottom.

         DO INFO = 1, N
            IF( IPIV( INFO ).GT.0 .AND. A( INFO, INFO ).EQ.ZERO ) RETURN
         END DO
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

         // invA = P * inv(U**T)*inv(D)*inv(U)*P**T.

        dtrtri(UPLO, 'U', N, A, LDA, INFO );

        // inv(D) and inv(D)*inv(U)

        K=1
        DO WHILE ( K .LE. N )
         if ( IPIV( K ).GT.0 ) {
            // 1 x 1 diagonal NNB
             WORK(K,INVD) = ONE /  A( K, K )
             WORK(K,INVD+1) = 0
            K=K+1
         } else {
            // 2 x 2 diagonal NNB
             T = WORK(K+1,1)
             AK = A( K, K ) / T
             AKP1 = A( K+1, K+1 ) / T
             AKKP1 = WORK(K+1,1)  / T
             D = T*( AK*AKP1-ONE )
             WORK(K,INVD) = AKP1 / D
             WORK(K+1,INVD+1) = AK / D
             WORK(K,INVD+1) = -AKKP1 / D
             WORK(K+1,INVD) = -AKKP1 / D
            K=K+2
         }
        END DO

        // inv(U**T) = (inv(U))**T

        // inv(U**T)*inv(D)*inv(U)

        CUT=N
        DO WHILE (CUT .GT. 0)
           NNB=NB
           if (CUT .LE. NNB) {
              NNB=CUT
           } else {
              COUNT = 0
              // count negative elements,
              DO I=CUT+1-NNB,CUT
                  IF (IPIV(I) .LT. 0) COUNT=COUNT+1
              END DO
              // need a even number for a clear cut
              IF (MOD(COUNT,2) .EQ. 1) NNB=NNB+1
           }

           CUT=CUT-NNB

           // U01 Block

           DO I=1,CUT
             DO J=1,NNB
              WORK(I,J)=A(I,CUT+J)
             END DO
           END DO

           // U11 Block

           DO I=1,NNB
             WORK(U11+I,I)=ONE
             DO J=1,I-1
                WORK(U11+I,J)=ZERO
             END DO
             DO J=I+1,NNB
                WORK(U11+I,J)=A(CUT+I,CUT+J)
             END DO
           END DO

           // invD*U01

           I=1
           DO WHILE (I .LE. CUT)
             if (IPIV(I) > 0) {
                DO J=1,NNB
                    WORK(I,J)=WORK(I,INVD)*WORK(I,J)
                END DO
                I=I+1
             } else {
                DO J=1,NNB
                   U01_I_J = WORK(I,J)
                   U01_IP1_J = WORK(I+1,J)
                   WORK(I,J)=WORK(I,INVD)*U01_I_J+ WORK(I,INVD+1)*U01_IP1_J                    WORK(I+1,J)=WORK(I+1,INVD)*U01_I_J+ WORK(I+1,INVD+1)*U01_IP1_J
                END DO
                I=I+2
             }
           END DO

         // invD1*U11

           I=1
           DO WHILE (I .LE. NNB)
             if (IPIV(CUT+I) > 0) {
                DO J=I,NNB
                    WORK(U11+I,J)=WORK(CUT+I,INVD)*WORK(U11+I,J)
                END DO
                I=I+1
             } else {
                DO J=I,NNB
                   U11_I_J = WORK(U11+I,J)
                   U11_IP1_J = WORK(U11+I+1,J)
                WORK(U11+I,J)=WORK(CUT+I,INVD)*WORK(U11+I,J) + WORK(CUT+I,INVD+1)*WORK(U11+I+1,J)                 WORK(U11+I+1,J)=WORK(CUT+I+1,INVD)*U11_I_J+ WORK(CUT+I+1,INVD+1)*U11_IP1_J
                END DO
                I=I+2
             }
           END DO

        // U11**T*invD1*U11->U11

        dtrmm('L','U','T','U',NNB, NNB, ONE,A(CUT+1,CUT+1),LDA,WORK(U11+1,1),N+NB+1);

         DO I=1,NNB
            DO J=I,NNB
              A(CUT+I,CUT+J)=WORK(U11+I,J)
            END DO
         END DO

           // U01**T*invD*U01->A(CUT+I,CUT+J)

         dgemm('T','N',NNB,NNB,CUT,ONE,A(1,CUT+1),LDA, WORK,N+NB+1, ZERO, WORK(U11+1,1), N+NB+1);


         // U11 =  U11**T*invD1*U11 + U01**T*invD*U01

         DO I=1,NNB
            DO J=I,NNB
              A(CUT+I,CUT+J)=A(CUT+I,CUT+J)+WORK(U11+I,J)
            END DO
         END DO

         // U01 =  U00**T*invD0*U01

         dtrmm('L',UPLO,'T','U',CUT, NNB, ONE,A,LDA,WORK,N+NB+1);


         // Update U01

         DO I=1,CUT
           DO J=1,NNB
            A(I,CUT+J)=WORK(I,J)
           END DO
         END DO

       // Next Block

       END DO

         // Apply PERMUTATIONS P and P**T: P * inv(U**T)*inv(D)*inv(U) *P**T

            I=1
            DO WHILE ( I .LE. N )
               if ( IPIV(I) .GT. 0 ) {
                  IP=IPIV(I)
                 IF (I .LT. IP) CALL DSYSWAPR( UPLO, N, A, LDA, I ,IP )
                 IF (I .GT. IP) CALL DSYSWAPR( UPLO, N, A, LDA, IP ,I )
               } else {
                 IP=-IPIV(I)
                 I=I+1
                 IF ( (I-1) .LT. IP) CALL DSYSWAPR( UPLO, N, A, LDA, I-1 ,IP )                  IF ( (I-1) .GT. IP) CALL DSYSWAPR( UPLO, N, A, LDA, IP ,I-1 )
              ENDIF
               I=I+1
            END DO
      } else {

         // LOWER...

         // invA = P * inv(U**T)*inv(D)*inv(U)*P**T.

         dtrtri(UPLO, 'U', N, A, LDA, INFO );

        // inv(D) and inv(D)*inv(U)

        K=N
        DO WHILE ( K .GE. 1 )
         if ( IPIV( K ).GT.0 ) {
            // 1 x 1 diagonal NNB
             WORK(K,INVD) = ONE /  A( K, K )
             WORK(K,INVD+1) = 0
            K=K-1
         } else {
            // 2 x 2 diagonal NNB
             T = WORK(K-1,1)
             AK = A( K-1, K-1 ) / T
             AKP1 = A( K, K ) / T
             AKKP1 = WORK(K-1,1) / T
             D = T*( AK*AKP1-ONE )
             WORK(K-1,INVD) = AKP1 / D
             WORK(K,INVD) = AK / D
             WORK(K,INVD+1) = -AKKP1 / D
             WORK(K-1,INVD+1) = -AKKP1 / D
            K=K-2
         }
        END DO

        // inv(U**T) = (inv(U))**T

        // inv(U**T)*inv(D)*inv(U)

        CUT=0
        DO WHILE (CUT .LT. N)
           NNB=NB
           if (CUT + NNB .GT. N) {
              NNB=N-CUT
           } else {
              COUNT = 0
              // count negative elements,
              DO I=CUT+1,CUT+NNB
                  IF (IPIV(I) .LT. 0) COUNT=COUNT+1
              END DO
              // need a even number for a clear cut
              IF (MOD(COUNT,2) .EQ. 1) NNB=NNB+1
           }
      // L21 Block
           DO I=1,N-CUT-NNB
             DO J=1,NNB
              WORK(I,J)=A(CUT+NNB+I,CUT+J)
             END DO
           END DO
      // L11 Block
           DO I=1,NNB
             WORK(U11+I,I)=ONE
             DO J=I+1,NNB
                WORK(U11+I,J)=ZERO
             END DO
             DO J=1,I-1
                WORK(U11+I,J)=A(CUT+I,CUT+J)
             END DO
           END DO

           // invD*L21

           I=N-CUT-NNB
           DO WHILE (I .GE. 1)
             if (IPIV(CUT+NNB+I) > 0) {
                DO J=1,NNB
                    WORK(I,J)=WORK(CUT+NNB+I,INVD)*WORK(I,J)
                END DO
                I=I-1
             } else {
                DO J=1,NNB
                   U01_I_J = WORK(I,J)
                   U01_IP1_J = WORK(I-1,J)
                   WORK(I,J)=WORK(CUT+NNB+I,INVD)*U01_I_J+ WORK(CUT+NNB+I,INVD+1)*U01_IP1_J                    WORK(I-1,J)=WORK(CUT+NNB+I-1,INVD+1)*U01_I_J+ WORK(CUT+NNB+I-1,INVD)*U01_IP1_J
                END DO
                I=I-2
             }
           END DO

         // invD1*L11

           I=NNB
           DO WHILE (I .GE. 1)
             if (IPIV(CUT+I) > 0) {
                DO J=1,NNB
                    WORK(U11+I,J)=WORK(CUT+I,INVD)*WORK(U11+I,J)
                END DO
                I=I-1
             } else {
                DO J=1,NNB
                   U11_I_J = WORK(U11+I,J)
                   U11_IP1_J = WORK(U11+I-1,J)
                WORK(U11+I,J)=WORK(CUT+I,INVD)*WORK(U11+I,J) + WORK(CUT+I,INVD+1)*U11_IP1_J                 WORK(U11+I-1,J)=WORK(CUT+I-1,INVD+1)*U11_I_J+ WORK(CUT+I-1,INVD)*U11_IP1_J
                END DO
                I=I-2
             }
           END DO

        // L11**T*invD1*L11->L11

        dtrmm('L',UPLO,'T','U',NNB, NNB, ONE,A(CUT+1,CUT+1),LDA,WORK(U11+1,1),N+NB+1);


         DO I=1,NNB
            DO J=1,I
              A(CUT+I,CUT+J)=WORK(U11+I,J)
            END DO
         END DO

        if ( (CUT+NNB) .LT. N ) {

           // L21**T*invD2*L21->A(CUT+I,CUT+J)

         dgemm('T','N',NNB,NNB,N-NNB-CUT,ONE,A(CUT+NNB+1,CUT+1) ,LDA,WORK,N+NB+1, ZERO, WORK(U11+1,1), N+NB+1);


         // L11 =  L11**T*invD1*L11 + U01**T*invD*U01

         DO I=1,NNB
            DO J=1,I
              A(CUT+I,CUT+J)=A(CUT+I,CUT+J)+WORK(U11+I,J)
            END DO
         END DO

         // L01 =  L22**T*invD2*L21

         dtrmm('L',UPLO,'T','U', N-NNB-CUT, NNB, ONE,A(CUT+NNB+1,CUT+NNB+1),LDA,WORK,N+NB+1);

       // Update L21

         DO I=1,N-CUT-NNB
           DO J=1,NNB
              A(CUT+NNB+I,CUT+J)=WORK(I,J)
           END DO
         END DO

       } else {

         // L11 =  L11**T*invD1*L11

         DO I=1,NNB
            DO J=1,I
              A(CUT+I,CUT+J)=WORK(U11+I,J)
            END DO
         END DO
       }

       // Next Block

           CUT=CUT+NNB
       END DO

         // Apply PERMUTATIONS P and P**T: P * inv(U**T)*inv(D)*inv(U) *P**T

            I=N
            DO WHILE ( I .GE. 1 )
               if ( IPIV(I) .GT. 0 ) {
                  IP=IPIV(I)
                 IF (I .LT. IP) CALL DSYSWAPR( UPLO, N, A, LDA, I ,IP  )
                 IF (I .GT. IP) CALL DSYSWAPR( UPLO, N, A, LDA, IP ,I )
               } else {
                 IP=-IPIV(I)
                 IF ( I .LT. IP) CALL DSYSWAPR( UPLO, N, A, LDA, I ,IP )
                 IF ( I .GT. IP) CALL DSYSWAPR( UPLO, N, A, LDA, IP, I )
                 I=I-1
               ENDIF
               I=I-1
            END DO
      }

      RETURN

      // End of DSYTRI2X

      }
