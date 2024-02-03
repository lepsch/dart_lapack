      SUBROUTINE ZSYTRF_AA_2STAGE( UPLO, N, A, LDA, TB, LTB, IPIV, IPIV2, WORK, LWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      IMPLICIT NONE

      // .. Scalar Arguments ..
      String             UPLO;
      int                N, LDA, LTB, LWORK, INFO;
      // ..
      // .. Array Arguments ..
      int                IPIV( * ), IPIV2( * );
      COMPLEX*16         A( LDA, * ), TB( * ), WORK( * )
      // ..

*  =====================================================================
      // .. Parameters ..
      COMPLEX*16         CZERO, CONE
      const              CZERO = ( 0.0D+0, 0.0D+0 ), CONE  = ( 1.0D+0, 0.0D+0 ) ;

      // .. Local Scalars ..
      bool               UPPER, TQUERY, WQUERY;
      int                I, J, K, I1, I2, TD;
      int                LDTB, NB, KB, JB, NT, IINFO;
      COMPLEX*16         PIV
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      // EXTERNAL LSAME, ILAENV
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZCOPY, ZGBTRF, ZGEMM, ZGETRF,   ZLACPY, ZLASET, ZLASWP, ZTRSM, ZSWAP
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MIN, MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      WQUERY = ( LWORK.EQ.-1 )
      TQUERY = ( LTB.EQ.-1 )
      if ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -4
      } else if ( LTB .LT. 4*N .AND. .NOT.TQUERY ) {
         INFO = -6
      } else if ( LWORK .LT. N .AND. .NOT.WQUERY ) {
         INFO = -10
      }

      if ( INFO.NE.0 ) {
         xerbla('ZSYTRF_AA_2STAGE', -INFO );
         RETURN
      }

      // Answer the query

      NB = ILAENV( 1, 'ZSYTRF_AA_2STAGE', UPLO, N, -1, -1, -1 )
      if ( INFO.EQ.0 ) {
         if ( TQUERY ) {
            TB( 1 ) = (3*NB+1)*N
         }
         if ( WQUERY ) {
            WORK( 1 ) = N*NB
         }
      }
      if ( TQUERY .OR. WQUERY ) {
         RETURN
      }

      // Quick return

      if ( N.EQ.0 ) {
         RETURN
      ENDIF

      // Determine the number of the block size

      LDTB = LTB/N
      if ( LDTB .LT. 3*NB+1 ) {
         NB = (LDTB-1)/3
      }
      if ( LWORK .LT. NB*N ) {
         NB = LWORK/N
      }

      // Determine the number of the block columns

      NT = (N+NB-1)/NB
      TD = 2*NB
      KB = MIN(NB, N)

      // Initialize vectors/matrices

      for (J = 1; J <= KB; J++) {
         IPIV( J ) = J
      END DO

      // Save NB

      TB( 1 ) = NB

      if ( UPPER ) {

         // .....................................................
         // Factorize A as U**T*D*U using the upper triangle of A
         // .....................................................

         DO J = 0, NT-1

            // Generate Jth column of W and H

            KB = MIN(NB, N-J*NB)
            DO I = 1, J-1
               if ( I.EQ.1 ) {
                   // H(I,J) = T(I,I)*U(I,J) + T(I+1,I)*U(I+1,J)
                  if ( I .EQ. (J-1) ) {
                     JB = NB+KB
                  } else {
                     JB = 2*NB
                  }
                  zgemm('NoTranspose', 'NoTranspose', NB, KB, JB, CONE,  TB( TD+1 + (I*NB)*LDTB ), LDTB-1, A( (I-1)*NB+1, J*NB+1 ), LDA, CZERO, WORK( I*NB+1 ), N );
               } else {
                  // H(I,J) = T(I,I-1)*U(I-1,J) + T(I,I)*U(I,J) + T(I,I+1)*U(I+1,J)
                  if ( I .EQ. (J-1) ) {
                     JB = 2*NB+KB
                  } else {
                     JB = 3*NB
                  }
                  zgemm('NoTranspose', 'NoTranspose', NB, KB, JB, CONE,  TB( TD+NB+1 + ((I-1)*NB)*LDTB ), LDTB-1, A( (I-2)*NB+1, J*NB+1 ), LDA, CZERO, WORK( I*NB+1 ), N );
               }
            END DO

            // Compute T(J,J)

            zlacpy('Upper', KB, KB, A( J*NB+1, J*NB+1 ), LDA, TB( TD+1 + (J*NB)*LDTB ), LDTB-1 );
            if ( J.GT.1 ) {
               // T(J,J) = U(1:J,J)'*H(1:J)
               zgemm('Transpose', 'NoTranspose', KB, KB, (J-1)*NB, -CONE, A( 1, J*NB+1 ), LDA, WORK( NB+1 ), N, CONE, TB( TD+1 + (J*NB)*LDTB ), LDTB-1 );
               // T(J,J) += U(J,J)'*T(J,J-1)*U(J-1,J)
               zgemm('Transpose', 'NoTranspose', KB, NB, KB, CONE,  A( (J-1)*NB+1, J*NB+1 ), LDA, TB( TD+NB+1 + ((J-1)*NB)*LDTB ), LDTB-1, CZERO, WORK( 1 ), N )                CALL ZGEMM( 'NoTranspose', 'NoTranspose', KB, KB, NB, -CONE, WORK( 1 ), N, A( (J-2)*NB+1, J*NB+1 ), LDA, CONE, TB( TD+1 + (J*NB)*LDTB ), LDTB-1 );
            }

            // Expand T(J,J) into full format

            for (I = 1; I <= KB; I++) {
               DO K = I+1, KB
                  TB( TD+(K-I)+1 + (J*NB+I-1)*LDTB ) = TB( TD-(K-(I+1)) + (J*NB+K-1)*LDTB )
               END DO
            END DO
            if ( J.GT.0 ) {
                // CALL CHEGST( 1, 'Upper', KB,
      // $                      TB( TD+1 + (J*NB)*LDTB ), LDTB-1,
      // $                      A( (J-1)*NB+1, J*NB+1 ), LDA, IINFO )
               ztrsm('L', 'U', 'T', 'N', KB, KB, CONE, A( (J-1)*NB+1, J*NB+1 ), LDA, TB( TD+1 + (J*NB)*LDTB ), LDTB-1 )                CALL ZTRSM( 'R', 'U', 'N', 'N', KB, KB, CONE, A( (J-1)*NB+1, J*NB+1 ), LDA, TB( TD+1 + (J*NB)*LDTB ), LDTB-1 );
            }

            if ( J.LT.NT-1 ) {
               if ( J.GT.0 ) {

                  // Compute H(J,J)

                  if ( J.EQ.1 ) {
                     zgemm('NoTranspose', 'NoTranspose', KB, KB, KB, CONE,  TB( TD+1 + (J*NB)*LDTB ), LDTB-1, A( (J-1)*NB+1, J*NB+1 ), LDA, CZERO, WORK( J*NB+1 ), N );
                  } else {
                     zgemm('NoTranspose', 'NoTranspose', KB, KB, NB+KB, CONE, TB( TD+NB+1 + ((J-1)*NB)*LDTB ), LDTB-1, A( (J-2)*NB+1, J*NB+1 ), LDA, CZERO, WORK( J*NB+1 ), N );
                  }

                  // Update with the previous column

                  zgemm('Transpose', 'NoTranspose', NB, N-(J+1)*NB, J*NB, -CONE, WORK( NB+1 ), N, A( 1, (J+1)*NB+1 ), LDA, CONE, A( J*NB+1, (J+1)*NB+1 ), LDA );
               }

               // Copy panel to workspace to call ZGETRF

               for (K = 1; K <= NB; K++) {
                   zcopy(N-(J+1)*NB, A( J*NB+K, (J+1)*NB+1 ), LDA, WORK( 1+(K-1)*N ), 1 );
               END DO

               // Factorize panel

               zgetrf(N-(J+1)*NB, NB,  WORK, N, IPIV( (J+1)*NB+1 ), IINFO );
                // IF (IINFO.NE.0 .AND. INFO.EQ.0) THEN
                   // INFO = IINFO+(J+1)*NB
                // END IF

               // Copy panel back

               for (K = 1; K <= NB; K++) {
                   zcopy(N-(J+1)*NB, WORK( 1+(K-1)*N ), 1, A( J*NB+K, (J+1)*NB+1 ), LDA );
               END DO

               // Compute T(J+1, J), zero out for GEMM update

               KB = MIN(NB, N-(J+1)*NB)
               zlaset('Full', KB, NB, CZERO, CZERO,  TB( TD+NB+1 + (J*NB)*LDTB), LDTB-1 )                CALL ZLACPY( 'Upper', KB, NB, WORK, N, TB( TD+NB+1 + (J*NB)*LDTB ), LDTB-1 );
               if ( J.GT.0 ) {
                  ztrsm('R', 'U', 'N', 'U', KB, NB, CONE, A( (J-1)*NB+1, J*NB+1 ), LDA, TB( TD+NB+1 + (J*NB)*LDTB ), LDTB-1 );
               }

               // Copy T(J,J+1) into T(J+1, J), both upper/lower for GEMM
               // updates

               for (K = 1; K <= NB; K++) {
                  for (I = 1; I <= KB; I++) {
                     TB( TD-NB+K-I+1 + (J*NB+NB+I-1)*LDTB ) = TB( TD+NB+I-K+1 + (J*NB+K-1)*LDTB )
                  END DO
               END DO
               zlaset('Lower', KB, NB, CZERO, CONE,  A( J*NB+1, (J+1)*NB+1), LDA );

               // Apply pivots to trailing submatrix of A

               for (K = 1; K <= KB; K++) {
                  // > Adjust ipiv
                  IPIV( (J+1)*NB+K ) = IPIV( (J+1)*NB+K ) + (J+1)*NB

                  I1 = (J+1)*NB+K
                  I2 = IPIV( (J+1)*NB+K )
                  if ( I1.NE.I2 ) {
                     // > Apply pivots to previous columns of L
                     zswap(K-1, A( (J+1)*NB+1, I1 ), 1,  A( (J+1)*NB+1, I2 ), 1 );
                     // > Swap A(I1+1:M, I1) with A(I2, I1+1:M)
                     IF( I2.GT.(I1+1) ) CALL ZSWAP( I2-I1-1, A( I1, I1+1 ), LDA, A( I1+1, I2 ), 1 )
                     // > Swap A(I2+1:M, I1) with A(I2+1:M, I2)
                     IF( I2.LT.N ) CALL ZSWAP( N-I2, A( I1, I2+1 ), LDA, A( I2, I2+1 ), LDA )
                     // > Swap A(I1, I1) with A(I2, I2)
                     PIV = A( I1, I1 )
                     A( I1, I1 ) = A( I2, I2 )
                     A( I2, I2 ) = PIV
                     // > Apply pivots to previous columns of L
                     if ( J.GT.0 ) {
                        zswap(J*NB, A( 1, I1 ), 1, A( 1, I2 ), 1 );
                     }
                  ENDIF
               END DO
            }
         END DO
      } else {

         // .....................................................
         // Factorize A as L*D*L**T using the lower triangle of A
         // .....................................................

         DO J = 0, NT-1

            // Generate Jth column of W and H

            KB = MIN(NB, N-J*NB)
            DO I = 1, J-1
               if ( I.EQ.1 ) {
                   // H(I,J) = T(I,I)*L(J,I)' + T(I+1,I)'*L(J,I+1)'
                  if ( I .EQ. (J-1) ) {
                     JB = NB+KB
                  } else {
                     JB = 2*NB
                  }
                  zgemm('NoTranspose', 'Transpose', NB, KB, JB, CONE, TB( TD+1 + (I*NB)*LDTB ), LDTB-1, A( J*NB+1, (I-1)*NB+1 ), LDA, CZERO, WORK( I*NB+1 ), N );
               } else {
                  // H(I,J) = T(I,I-1)*L(J,I-1)' + T(I,I)*L(J,I)' + T(I,I+1)*L(J,I+1)'
                  if ( I .EQ. (J-1) ) {
                     JB = 2*NB+KB
                  } else {
                     JB = 3*NB
                  }
                  zgemm('NoTranspose', 'Transpose', NB, KB, JB, CONE,  TB( TD+NB+1 + ((I-1)*NB)*LDTB ), LDTB-1, A( J*NB+1, (I-2)*NB+1 ), LDA, CZERO, WORK( I*NB+1 ), N );
               }
            END DO

            // Compute T(J,J)

            zlacpy('Lower', KB, KB, A( J*NB+1, J*NB+1 ), LDA, TB( TD+1 + (J*NB)*LDTB ), LDTB-1 );
            if ( J.GT.1 ) {
               // T(J,J) = L(J,1:J)*H(1:J)
               zgemm('NoTranspose', 'NoTranspose', KB, KB, (J-1)*NB, -CONE, A( J*NB+1, 1 ), LDA, WORK( NB+1 ), N, CONE, TB( TD+1 + (J*NB)*LDTB ), LDTB-1 );
               // T(J,J) += L(J,J)*T(J,J-1)*L(J,J-1)'
               zgemm('NoTranspose', 'NoTranspose', KB, NB, KB, CONE,  A( J*NB+1, (J-1)*NB+1 ), LDA, TB( TD+NB+1 + ((J-1)*NB)*LDTB ), LDTB-1, CZERO, WORK( 1 ), N )                CALL ZGEMM( 'NoTranspose', 'Transpose', KB, KB, NB, -CONE, WORK( 1 ), N, A( J*NB+1, (J-2)*NB+1 ), LDA, CONE, TB( TD+1 + (J*NB)*LDTB ), LDTB-1 );
            }

            // Expand T(J,J) into full format

            for (I = 1; I <= KB; I++) {
               DO K = I+1, KB
                  TB( TD-(K-(I+1)) + (J*NB+K-1)*LDTB ) = TB( TD+(K-I)+1 + (J*NB+I-1)*LDTB )
               END DO
            END DO
            if ( J.GT.0 ) {
                // CALL CHEGST( 1, 'Lower', KB,
      // $                      TB( TD+1 + (J*NB)*LDTB ), LDTB-1,
      // $                      A( J*NB+1, (J-1)*NB+1 ), LDA, IINFO )
               ztrsm('L', 'L', 'N', 'N', KB, KB, CONE, A( J*NB+1, (J-1)*NB+1 ), LDA, TB( TD+1 + (J*NB)*LDTB ), LDTB-1 )                CALL ZTRSM( 'R', 'L', 'T', 'N', KB, KB, CONE, A( J*NB+1, (J-1)*NB+1 ), LDA, TB( TD+1 + (J*NB)*LDTB ), LDTB-1 );
            }

            // Symmetrize T(J,J)

            for (I = 1; I <= KB; I++) {
               DO K = I+1, KB
                  TB( TD-(K-(I+1)) + (J*NB+K-1)*LDTB ) = TB( TD+(K-I)+1 + (J*NB+I-1)*LDTB )
               END DO
            END DO

            if ( J.LT.NT-1 ) {
               if ( J.GT.0 ) {

                  // Compute H(J,J)

                  if ( J.EQ.1 ) {
                     zgemm('NoTranspose', 'Transpose', KB, KB, KB, CONE,  TB( TD+1 + (J*NB)*LDTB ), LDTB-1, A( J*NB+1, (J-1)*NB+1 ), LDA, CZERO, WORK( J*NB+1 ), N );
                  } else {
                     zgemm('NoTranspose', 'Transpose', KB, KB, NB+KB, CONE, TB( TD+NB+1 + ((J-1)*NB)*LDTB ), LDTB-1, A( J*NB+1, (J-2)*NB+1 ), LDA, CZERO, WORK( J*NB+1 ), N );
                  }

                  // Update with the previous column

                  zgemm('NoTranspose', 'NoTranspose', N-(J+1)*NB, NB, J*NB, -CONE, A( (J+1)*NB+1, 1 ), LDA, WORK( NB+1 ), N, CONE, A( (J+1)*NB+1, J*NB+1 ), LDA );
               }

               // Factorize panel

               zgetrf(N-(J+1)*NB, NB,  A( (J+1)*NB+1, J*NB+1 ), LDA, IPIV( (J+1)*NB+1 ), IINFO );
                // IF (IINFO.NE.0 .AND. INFO.EQ.0) THEN
                   // INFO = IINFO+(J+1)*NB
                // END IF

               // Compute T(J+1, J), zero out for GEMM update

               KB = MIN(NB, N-(J+1)*NB)
               zlaset('Full', KB, NB, CZERO, CZERO,  TB( TD+NB+1 + (J*NB)*LDTB), LDTB-1 )                CALL ZLACPY( 'Upper', KB, NB, A( (J+1)*NB+1, J*NB+1 ), LDA, TB( TD+NB+1 + (J*NB)*LDTB ), LDTB-1 );
               if ( J.GT.0 ) {
                  ztrsm('R', 'L', 'T', 'U', KB, NB, CONE, A( J*NB+1, (J-1)*NB+1 ), LDA, TB( TD+NB+1 + (J*NB)*LDTB ), LDTB-1 );
               }

               // Copy T(J+1,J) into T(J, J+1), both upper/lower for GEMM
               // updates

               for (K = 1; K <= NB; K++) {
                  for (I = 1; I <= KB; I++) {
                     TB( TD-NB+K-I+1 + (J*NB+NB+I-1)*LDTB ) = TB( TD+NB+I-K+1 + (J*NB+K-1)*LDTB )
                  END DO
               END DO
               zlaset('Upper', KB, NB, CZERO, CONE,  A( (J+1)*NB+1, J*NB+1 ), LDA );

               // Apply pivots to trailing submatrix of A

               for (K = 1; K <= KB; K++) {
                  // > Adjust ipiv
                  IPIV( (J+1)*NB+K ) = IPIV( (J+1)*NB+K ) + (J+1)*NB

                  I1 = (J+1)*NB+K
                  I2 = IPIV( (J+1)*NB+K )
                  if ( I1.NE.I2 ) {
                     // > Apply pivots to previous columns of L
                     zswap(K-1, A( I1, (J+1)*NB+1 ), LDA,  A( I2, (J+1)*NB+1 ), LDA );
                     // > Swap A(I1+1:M, I1) with A(I2, I1+1:M)
                     IF( I2.GT.(I1+1) ) CALL ZSWAP( I2-I1-1, A( I1+1, I1 ), 1, A( I2, I1+1 ), LDA )
                     // > Swap A(I2+1:M, I1) with A(I2+1:M, I2)
                     IF( I2.LT.N ) CALL ZSWAP( N-I2, A( I2+1, I1 ), 1, A( I2+1, I2 ), 1 )
                     // > Swap A(I1, I1) with A(I2, I2)
                     PIV = A( I1, I1 )
                     A( I1, I1 ) = A( I2, I2 )
                     A( I2, I2 ) = PIV
                     // > Apply pivots to previous columns of L
                     if ( J.GT.0 ) {
                        zswap(J*NB, A( I1, 1 ), LDA, A( I2, 1 ), LDA );
                     }
                  ENDIF
               END DO

               // Apply pivots to previous columns of L

                // CALL ZLASWP( J*NB, A( 1, 1 ), LDA,
      // $                     (J+1)*NB+1, (J+1)*NB+KB, IPIV, 1 )
            }
         END DO
      }

      // Factor the band matrix
      zgbtrf(N, N, NB, NB, TB, LDTB, IPIV2, INFO );

      RETURN

      // End of ZSYTRF_AA_2STAGE

      }
